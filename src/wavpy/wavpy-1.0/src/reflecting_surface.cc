 /*
    Copyright (C) 2017  Fran Fabra

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details at:

    <http://www.gnu.org/licenses/>
*/

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include "wavpy_global_variables.h"
#include "reflecting_surface.h"
#include "ancillary_functions.h"

extern "C" {
void compute_epsilon_sea_water_(double *, double *, double *, double *, double *);
void compute_epsilon_sea_ice_(double *, double *, double *);
void compute_epsilon_dry_snow_(double *, double *, double *);
void compute_epsilon_wet_snow_(double *, double *, double *, double *, double *);
void compute_fresnel_linear_refl_(double *, double *, double *, double *, double *, double *, double *, double *, double *);
void compute_fresnel_linear_trans_(double *, double *, double *, double *, double *, double *, double *, double *, double *);
void sea_spectrum_elfouhaily_(double psi[], double kx[], double ky[], double *, int *, double *, double *, double *);
}

using namespace std;

Reflecting_surface::Reflecting_surface( void ){
	//Set default values
	mss_x = 0.0075;
	mss_y = 0.0075;
	sigma_z = 0.00669*exp(156.1*(mss_x + mss_y));
	wind_U10_azimuth = 0.; //Degrees
	wind_U10_speed = 5.; //m/sec
	//Gram-Charlier C21 and C03 coefficients for 5 m/s wind
	c21_coeff = -0.033;
	c03_coeff = -0.125;
	epsilon_real = 73.423;
	epsilon_imag = 56.067;
	medium = "Sea water with T=15C and sal=35psu";
	freq_GHz = FREQ_GPS_L1/1000000000.0;
	n_spectrum = 0;
	omnidirectional_spec = false;
	k_threshold = 2.0*PI_NUM*FREQ_GPS_L1/(3.0*C_LIGHT); //Brown 1978 for theta=0
	size_lon_wgrid = 0;
	size_lat_wgrid = 0;
	use_wind_grid = false;
	return;
}

void Reflecting_surface::dump_parameters( void ){
	printf("\n"); 
	printf("====== Reflecting surface Parameters ======\n");
	printf("| MSS_x = %f \n", mss_x);
	printf("| MSS_y = %f \n", mss_y);
	if(!use_wind_grid){
		printf("| Wind Azimuth = %f degrees\n", wind_U10_azimuth);
		printf("| Wind Speed = %f m/s\n", wind_U10_speed);
	}else{
		printf("| Using a wind grid\n");
	}
	printf("| C21 coeff = %f \n", c21_coeff);
	printf("| C03 coeff = %f \n", c03_coeff);
	printf("| DielConst (real) = %f \n", epsilon_real);
	printf("| DielConst (imag) = %f \n", epsilon_imag);
	printf("| Medium = %s \n", medium.c_str());
	printf("| Frequency = %f GHz\n", freq_GHz);
	printf("| k threshold = %f m^(-1)\n", k_threshold);
	printf("===========================================\n");
	return;
}

void Reflecting_surface::set_frequency( double frequency_GHz ){
	freq_GHz = frequency_GHz;
	return;
}

void Reflecting_surface::set_k_threshold( double k_threshold_in ){
	if(k_threshold_in < 0.0){
		printf("ERROR! k_threshold not valid (< 0)\n");
		return;
	}
	k_threshold = k_threshold_in;
	return;
}

void Reflecting_surface::set_k_threshold_Brown( double incidence_deg ){
	k_threshold = 2.0*PI_NUM*freq_GHz*1000000000.0*cos(incidence_deg*PI_NUM/180.0)/(3.0*C_LIGHT);
	return;
}

double Reflecting_surface::get_k_threshold( void ){
	return k_threshold;
}

void Reflecting_surface::epsilon_sea_water( double salinity, double temperature ){
	double sal, temp;
	string sal_str, temp_str;
	stringstream number;
	sal = salinity;
	temp = temperature;
	number << sal;
	number >> sal_str;
	number.clear();
	number << temp;
	number >> temp_str;
	medium = "Sea water with temperature = " + temp_str + "C and salinity = " + sal_str + "psu";
	compute_epsilon_sea_water_(&sal, &temp, &freq_GHz, &epsilon_real, &epsilon_imag);
	return;
}

void Reflecting_surface::epsilon_sea_ice( double brine_volume ){
	double v_b;
	string v_b_str;
	stringstream number;
	v_b = brine_volume;
	number << v_b;
	number >> v_b_str;
	medium = "Sea ice with brine volume = " + v_b_str + " 1/1000";
	compute_epsilon_sea_ice_(&v_b, &epsilon_real, &epsilon_imag);
	return;
}

void Reflecting_surface::epsilon_dry_snow( double snow_density ){
	double rho_s;
	string rho_s_str;
	stringstream number;
	rho_s = snow_density;
	number << rho_s;
	number >> rho_s_str;
	medium = "Dry snow with snow density = " + rho_s_str + "gr/cm^3";
	compute_epsilon_dry_snow_(&rho_s, &epsilon_real, &epsilon_imag);
	return;
}

void Reflecting_surface::epsilon_wet_snow( double snow_density, double water_volume ){
	double rho_s, m_v;
	string rho_s_str, m_v_str;
	stringstream number;
	rho_s = snow_density;
	m_v = water_volume;
	number << rho_s;
	number >> rho_s_str;
	number.clear();
	number << m_v;
	number >> m_v_str;
	medium = "Wet snow with snow density = " + rho_s_str + "gr/cm^3 and water volume = " + m_v_str + "%";
	compute_epsilon_wet_snow_(&rho_s, &m_v, &freq_GHz, &epsilon_real, &epsilon_imag);
	return;
}

void Reflecting_surface::compute_Rfresnel_linear( double incidence_deg, double epsilon_upLayer[2], double rvv[2], double rhh[2] ){
	double theta_rad, rvv_real, rvv_imag, rhh_real, rhh_imag, epsilon_real_UP, epsilon_imag_UP;
	theta_rad = incidence_deg*PI_NUM/180.0;
	epsilon_real_UP = epsilon_upLayer[0];
	epsilon_imag_UP = epsilon_upLayer[1];
	compute_fresnel_linear_refl_(&theta_rad, &epsilon_real_UP, &epsilon_imag_UP, &epsilon_real, &epsilon_imag, &rvv_real, &rvv_imag, &rhh_real, &rhh_imag);
	rvv[0] = rvv_real;
	rvv[1] = rvv_imag;
	rhh[0] = rhh_real;
	rhh[1] = rhh_imag;
	return;
}

void Reflecting_surface::compute_Rfresnel_circular( double incidence_deg, double epsilon_upLayer[2], double rco[2], double rcross[2] ){
	double theta_rad, rvv_real, rvv_imag, rhh_real, rhh_imag, epsilon_real_UP, epsilon_imag_UP;
	theta_rad = incidence_deg*PI_NUM/180.0;
	epsilon_real_UP = epsilon_upLayer[0];
	epsilon_imag_UP = epsilon_upLayer[1];
	compute_fresnel_linear_refl_(&theta_rad, &epsilon_real_UP, &epsilon_imag_UP, &epsilon_real, &epsilon_imag, &rvv_real, &rvv_imag, &rhh_real, &rhh_imag);
	rco[0] = 0.5*(rvv_real + rhh_real);
	rco[1] = 0.5*(rvv_imag + rhh_imag);
	rcross[0] = 0.5*(rvv_real - rhh_real);
	rcross[1] = 0.5*(rvv_imag - rhh_imag);
	return;
}

void Reflecting_surface::compute_Tfresnel_linear( double incidence_deg, double epsilon_upLayer[2], double tvv[2], double thh[2] ){
	double theta_rad, tvv_real, tvv_imag, thh_real, thh_imag, epsilon_real_UP, epsilon_imag_UP;
	theta_rad = incidence_deg*PI_NUM/180.0;
	epsilon_real_UP = epsilon_upLayer[0];
	epsilon_imag_UP = epsilon_upLayer[1];
	compute_fresnel_linear_trans_(&theta_rad, &epsilon_real_UP, &epsilon_imag_UP, &epsilon_real, &epsilon_imag, &tvv_real, &tvv_imag, &thh_real, &thh_imag);
	tvv[0] = tvv_real;
	tvv[1] = tvv_imag;
	thh[0] = thh_real;
	thh[1] = thh_imag;
	return;
}

void Reflecting_surface::compute_Tfresnel_circular( double incidence_deg, double epsilon_upLayer[2], double tco[2], double tcross[2] ){
	double theta_rad, tvv_real, tvv_imag, thh_real, thh_imag, epsilon_real_UP, epsilon_imag_UP;
	theta_rad = incidence_deg*PI_NUM/180.0;
	epsilon_real_UP = epsilon_upLayer[0];
	epsilon_imag_UP = epsilon_upLayer[1];
	compute_fresnel_linear_trans_(&theta_rad, &epsilon_real_UP, &epsilon_imag_UP, &epsilon_real, &epsilon_imag, &tvv_real, &tvv_imag, &thh_real, &thh_imag);
	tco[0] = 0.5*(tvv_real + thh_real);
	tco[1] = 0.5*(tvv_imag + thh_imag);
	tcross[0] = 0.5*(tvv_real - thh_real);
	tcross[1] = 0.5*(tvv_imag - thh_imag);
	return;
}

void Reflecting_surface::compute_sea_spectrum( int num_samples, double delta_k_in, double theta_deg, double omega ){
	int i, j, n;
	double thetaD, ome;
	double *psi;
	if(num_samples <= 0){
		printf("ERROR! num_samples not valid (<= 0)\n");
		return;
	}
	//Free memory
	if(n_spectrum > 0){
		for(i=0; i<n_spectrum; i++){
			free(surface_spectrum[i]);
		}
		free(surface_spectrum);
		free(kx_spec);
		free(ky_spec);
	}
	omnidirectional_spec = false;
	n = num_samples/2;
	n_spectrum = num_samples;
	delta_k = delta_k_in;
	thetaD = theta_deg;
	ome = omega;
	surface_spectrum = (double **) malloc(n_spectrum*sizeof(double *));
	for(i=0;i<n_spectrum;i++)
		surface_spectrum[i] = (double *) malloc(n_spectrum*sizeof(double));
	psi = (double*) malloc (n_spectrum*n_spectrum*sizeof(double));
	kx_spec = (double*) malloc (n_spectrum*sizeof(double));
	ky_spec = (double*) malloc (n_spectrum*sizeof(double));
	sea_spectrum_elfouhaily_((double*)psi, (double*)kx_spec, (double*)ky_spec, &delta_k, &n, &wind_U10_speed, &thetaD, &ome);
	for(i=0; i<n_spectrum; i++){
		for(j=0; j<n_spectrum; j++){
			surface_spectrum[i][j] = psi[j*n_spectrum + i];
		}
	}
	free(psi);
	return;
}

void Reflecting_surface::set_surf_spectrum( double* spectrum_in, int dim_x, int dim_y, double* kx_spec_in, int dim_kx, double* ky_spec_in, int dim_ky, double delta_k_in ){
	int i, j;
	if(dim_x <= 0){
		printf("ERROR! Dimension not valid (<= 0)\n");
		return;
	}
	if(!((dim_x == dim_y)&&(dim_x == dim_kx)&&(dim_x == dim_ky))){
		printf("ERROR! Inputs must have the same dim (spectrum[dim, dim], kx[dim] and ky[dim])\n");
		return;
	}
	//Free memory
	if(n_spectrum > 0){
		for(i=0; i<n_spectrum; i++){
			free(surface_spectrum[i]);
		}
		free(surface_spectrum);
		free(kx_spec);
		free(ky_spec);
	}
	n_spectrum = dim_x;
	delta_k = delta_k_in;
	omnidirectional_spec = false;
	surface_spectrum = (double **) malloc(n_spectrum*sizeof(double *));
	for(i=0; i<n_spectrum; i++)
		surface_spectrum[i] = (double *) malloc(n_spectrum*sizeof(double));
	kx_spec = (double*) malloc (n_spectrum*sizeof(double));
	ky_spec = (double*) malloc (n_spectrum*sizeof(double));
	for(i=0; i<n_spectrum; i++){
		for(j=0; j<n_spectrum; j++){
			surface_spectrum[i][j] = spectrum_in[i*n_spectrum + j];
		}
		kx_spec[i] = kx_spec_in[i];
		ky_spec[i] = ky_spec_in[i];
	}
	return;
}

void Reflecting_surface::set_surf_spectrum_omnidir( double* spectrum_in_omnidir, int dim_x_omnidir, double* k_spec_in_omnidir, int dim_k_omnidir, double delta_k_in ){
	int i;
	if(dim_x_omnidir <= 0){
		printf("ERROR! Dimension not valid (<= 0)\n");
		return;
	}
	if(dim_x_omnidir != dim_k_omnidir){
		printf("ERROR! Inputs must have the same dim (spectrum[dim] and k[dim])\n");
		return;
	}
	//Free memory
	if(n_spectrum > 0){
		for(i=0; i<n_spectrum; i++){
			free(surface_spectrum[i]);
		}
		free(surface_spectrum);
		free(kx_spec);
		free(ky_spec);
	}
	n_spectrum = dim_x_omnidir;
	delta_k = delta_k_in;
	omnidirectional_spec = true;
	surface_spectrum = (double **) malloc(n_spectrum*sizeof(double *));
	for(i=0; i<n_spectrum; i++)
		surface_spectrum[i] = (double *) malloc(1*sizeof(double));
	kx_spec = (double*) malloc (n_spectrum*sizeof(double));
	for(i=0; i<n_spectrum; i++){
		surface_spectrum[i][0] = spectrum_in_omnidir[i];
		kx_spec[i] = k_spec_in_omnidir[i];
	}
	return;
}

void Reflecting_surface::get_surf_spectrum( int x, int y, double kx_ky_spectrum[3] ){
	if(n_spectrum == 0){
		printf("ERROR! Spectrum not in memory\n");
		return;
	}
	if(omnidirectional_spec){
		printf("ERROR! Omnidirectional spectrum. Use get_surf_spectrum_omnidir(int pos)\n");
		return;
	}
	if((x>=n_spectrum)||(y>=n_spectrum)||(x<0)||(y<0)){
		printf("ERROR! Index out of spectrum limits (max=%d)\n", n_spectrum);
		return;
	}
	kx_ky_spectrum[0] = kx_spec[x];
	kx_ky_spectrum[1] = ky_spec[y];
	kx_ky_spectrum[2] = surface_spectrum[x][y];
	return;
}

void Reflecting_surface::get_surf_spectrum_omnidir( int x, double k_spectrum[2] ){
	if(n_spectrum == 0){
		printf("ERROR! Spectrum not in memory\n");
		return;
	}
	if(!omnidirectional_spec){
		printf("ERROR! Not omnidirectional spectrum. Use get_surf_spectrum(int pos_x, int pos_y)\n");
		return;
	}
	if((x>=n_spectrum)||(x<0)){
		printf("ERROR! Index out of spectrum limits (max=%d)\n", n_spectrum);
		return;
	}
	k_spectrum[0] = kx_spec[x];
	k_spectrum[1] = surface_spectrum[x][0];
	return;
}

void Reflecting_surface::compute_mss_from_spectrum( void ){
	int i, j;
	double accum_mss, accum_rmsz2;
	if(n_spectrum == 0){
		printf("ERROR! Spectrum not in memory\n");
		return;
	}
	accum_rmsz2 = 0.0;
	if(omnidirectional_spec){
		accum_mss = 0.0;
		for(i=0; i<n_spectrum; i++){
			if(kx_spec[i] < k_threshold){
				accum_rmsz2 = accum_rmsz2 + surface_spectrum[i][0]*delta_k;
				accum_mss = accum_mss + surface_spectrum[i][0]*kx_spec[i]*kx_spec[i]*delta_k;
			}
		}
		mss_x = accum_mss/2.0;
		mss_y = accum_mss/2.0;
	}else{
		mss_x = 0.0;
		mss_y = 0.0;
		for(i=0; i<n_spectrum; i++){
			for(j=0; j<n_spectrum; j++){
				if(sqrt(kx_spec[i]*kx_spec[i] + ky_spec[j]*ky_spec[j]) < k_threshold){
					accum_rmsz2 = accum_rmsz2 + surface_spectrum[i][j]*delta_k*delta_k;
					mss_x = mss_x + surface_spectrum[i][j]*kx_spec[i]*kx_spec[i]*delta_k*delta_k;
					mss_y = mss_y + surface_spectrum[i][j]*ky_spec[j]*ky_spec[j]*delta_k*delta_k;
				}
			}
		}
	}
	sigma_z = sqrt(accum_rmsz2);
	return;
}

void Reflecting_surface::compute_mss_from_wind( void ){
	double fU;
	if(wind_U10_speed <= 0.0){
		printf("ERROR! Wind speed not valid (<= 0.0)\n");
		return;
	}
	//Based on Katzberg et al 2006 for sea surface
	if(wind_U10_speed < 3.49){
		fU = wind_U10_speed;
	}else{
		if(wind_U10_speed < 46.0){
			fU = 6.0*log(wind_U10_speed) - 4.0;
		}else{
			fU = 0.411*wind_U10_speed;
		}
	}
	mss_x = 0.45*0.00316*fU;
	mss_y = 0.45*(0.003 + 0.00192*fU);
	//============================================
	sigma_z = 0.00669*exp(156.1*(mss_x + mss_y));
	return;
}

void Reflecting_surface::set_wind_grid( double* wspeed_grid_in, int dim_wsx, int dim_wsy, double* wazim_grid_in, int dim_wax, int dim_way, double* lon_in, int dim_lon, double* lat_in, int dim_lat ){
	int i, j, i_lon, j_lat;
	bool reversed_lon, reversed_lat;
	if(dim_lon < 2){
		printf("ERROR! Dimension not valid (< 2)\n");
		return;
	}
	if(!((dim_lon == dim_wsx)&&(dim_lon == dim_wax))){
		printf("ERROR! Inputs must have the same dim (wind_speed[dim_x, dim_y], wind_azimuth[dim_x, dim_y], lon[dim_x] and lat[dim_y])\n");
		return;
	}
	if(dim_lat < 2){
		printf("ERROR! Dimension not valid (< 2)\n");
		return;
	}
	if(!((dim_lat == dim_wsy)&&(dim_lat == dim_way))){
		printf("ERROR! Inputs must have the same dim (wind_speed[dim_x, dim_y], wind_azimuth[dim_x, dim_y], lon[dim_x] and lat[dim_y])\n");
		return;
	}
	//Check longitudes and latitudes
	if(range180(lon_in[1]) > range180(lon_in[0])){
		reversed_lon = false;
		for(i=1; i<dim_lon; i++){
			if(range180(lon_in[i]) <= range180(lon_in[i - 1])){
				printf("ERROR! Longitudes must be strictly increasing or decreasing.\n");
				return;
			}
		}
	}else{
		reversed_lon = true;
		for(i=1; i<dim_lon; i++){
			if(range180(lon_in[i]) >= range180(lon_in[i - 1])){
				printf("ERROR! Longitudes must be strictly increasing or decreasing.\n");
				return;
			}
		}
	}
	if(range180(lat_in[1]) > range180(lat_in[0])){
		reversed_lat = false;
		for(i=1; i<dim_lat; i++){
			if(range180(lat_in[i]) <= range180(lat_in[i - 1])){
				printf("ERROR! Latitudes must be strictly increasing or decreasing.\n");
				return;
			}
		}
	}else{
		reversed_lat = true;
		for(i=1; i<dim_lat; i++){
			if(range180(lat_in[i]) >= range180(lat_in[i - 1])){
				printf("ERROR! Latitudes must be strictly increasing or decreasing.\n");
				return;
			}
		}
	}
	//Free memory
	if(size_lon_wgrid > 0){
		for(i=0; i<size_lon_wgrid; i++){
			free(wind_U10_speed_grid[i]);
			free(wind_U10_azimuth_grid[i]);
		}
		free(wind_U10_speed_grid);
		free(wind_U10_azimuth_grid);
		free(lon_wgrid);
		free(lat_wgrid);
	}
	//Assign new values
	size_lon_wgrid = dim_lon;
	size_lat_wgrid = dim_lat;
	use_wind_grid = true;
	wind_U10_speed_grid = (double **) malloc(size_lon_wgrid*sizeof(double *));
	for(i=0; i<size_lon_wgrid; i++)
		wind_U10_speed_grid[i] = (double *) malloc(size_lat_wgrid*sizeof(double));
	wind_U10_azimuth_grid = (double **) malloc(size_lon_wgrid*sizeof(double *));
	for(i=0; i<size_lon_wgrid; i++)
		wind_U10_azimuth_grid[i] = (double *) malloc(size_lat_wgrid*sizeof(double));
	lon_wgrid = (double*) malloc (size_lon_wgrid*sizeof(double));
	lat_wgrid = (double*) malloc (size_lat_wgrid*sizeof(double));
	for(i=0; i<size_lon_wgrid; i++){
		i_lon = i;
		if(reversed_lon){
			i_lon = size_lon_wgrid - 1 - i;
		}
		for(j=0; j<size_lat_wgrid; j++){
			j_lat = j;
			if(reversed_lat){
				j_lat = size_lat_wgrid - 1 - j;
			}
			wind_U10_speed_grid[i_lon][j_lat] = wspeed_grid_in[i*size_lat_wgrid + j];
			wind_U10_azimuth_grid[i_lon][j_lat] = wazim_grid_in[i*size_lat_wgrid + j];
			if(i==0){
				lat_wgrid[j_lat] = range180(lat_in[j]);
			}
		}
		lon_wgrid[i_lon] = range180(lon_in[i]);
	}
	return;
}

void Reflecting_surface::interp_wind_grid( double lon_in, double lat_in )
{
	int i, j;
	double longitude, latitude;
	double *ws_grid, *wa_grid_I, *wa_grid_Q;
	double interp_wa_I, interp_wa_Q, coeff_up, coeff_dw;
	const gsl_interp2d_type *T = gsl_interp2d_bilinear;
	longitude = range180(lon_in);
	latitude = range180(lat_in);
	if((longitude >= lon_wgrid[0])&&(longitude <= lon_wgrid[size_lon_wgrid - 1])&&(latitude >= lat_wgrid[0])&&(latitude <= lat_wgrid[size_lat_wgrid - 1])){
		gsl_spline2d *spline = gsl_spline2d_alloc(T, size_lon_wgrid, size_lat_wgrid);
		gsl_interp_accel *xacc = gsl_interp_accel_alloc();
		gsl_interp_accel *yacc = gsl_interp_accel_alloc();
		ws_grid = (double*) malloc (size_lon_wgrid*size_lat_wgrid*sizeof(double));
		wa_grid_I = (double*) malloc (size_lon_wgrid*size_lat_wgrid*sizeof(double));
		wa_grid_Q = (double*) malloc (size_lon_wgrid*size_lat_wgrid*sizeof(double));
		for(i=0; i<size_lon_wgrid; i++){
			for(j=0; j<size_lat_wgrid; j++){
				ws_grid[j*size_lon_wgrid + i] = wind_U10_speed_grid[i][j];
				wa_grid_I[j*size_lon_wgrid + i] = cos(wind_U10_azimuth_grid[i][j]*PI_NUM/180.0);
				wa_grid_Q[j*size_lon_wgrid + i] = sin(wind_U10_azimuth_grid[i][j]*PI_NUM/180.0);
			}
		}
		gsl_spline2d_init(spline, lon_wgrid, lat_wgrid, ws_grid, size_lon_wgrid, size_lat_wgrid);
		wind_U10_speed = gsl_spline2d_eval(spline, longitude, latitude, xacc, yacc);
		free(ws_grid);
		gsl_spline2d_init(spline, lon_wgrid, lat_wgrid, wa_grid_I, size_lon_wgrid, size_lat_wgrid);
		interp_wa_I = gsl_spline2d_eval(spline, longitude, latitude, xacc, yacc);
		free(wa_grid_I);
		gsl_spline2d_init(spline, lon_wgrid, lat_wgrid, wa_grid_Q, size_lon_wgrid, size_lat_wgrid);
		interp_wa_Q = gsl_spline2d_eval(spline, longitude, latitude, xacc, yacc);
		free(wa_grid_Q);
		wind_U10_azimuth = atan2(interp_wa_Q, interp_wa_I)*180.0/PI_NUM;
		gsl_spline2d_free(spline);
		gsl_interp_accel_free(xacc);
		gsl_interp_accel_free(yacc);
	}else{ //If input coordinates are outside grid, then the closest value is taken (applying linear interpolation at the limits)
		if((longitude < lon_wgrid[0])&&(latitude < lat_wgrid[0])){
			wind_U10_speed = wind_U10_speed_grid[0][0];
			wind_U10_azimuth = wind_U10_azimuth_grid[0][0];
		}
		if((longitude < lon_wgrid[0])&&(latitude > lat_wgrid[size_lat_wgrid - 1])){
			wind_U10_speed = wind_U10_speed_grid[0][size_lat_wgrid - 1];
			wind_U10_azimuth = wind_U10_azimuth_grid[0][size_lat_wgrid - 1];
		}
		if((longitude > lon_wgrid[size_lon_wgrid - 1])&&(latitude < lat_wgrid[0])){
			wind_U10_speed = wind_U10_speed_grid[size_lon_wgrid - 1][0];
			wind_U10_azimuth = wind_U10_azimuth_grid[size_lon_wgrid - 1][0];
		}
		if((longitude > lon_wgrid[size_lon_wgrid - 1])&&(latitude > lat_wgrid[size_lat_wgrid - 1])){
			wind_U10_speed = wind_U10_speed_grid[size_lon_wgrid - 1][size_lat_wgrid - 1];
			wind_U10_azimuth = wind_U10_azimuth_grid[size_lon_wgrid - 1][size_lat_wgrid - 1];
		}
		if((longitude < lon_wgrid[0])&&(latitude >= lat_wgrid[0])&&(latitude <= lat_wgrid[size_lat_wgrid - 1])){
			i = 1;
			while(latitude > lat_wgrid[i]) i++;
			coeff_up = (latitude - lat_wgrid[i - 1])/(lat_wgrid[i] - lat_wgrid[i - 1]);
			coeff_dw = (lat_wgrid[i] - latitude)/(lat_wgrid[i] - lat_wgrid[i - 1]);
			wind_U10_speed = coeff_up*wind_U10_speed_grid[0][i] + coeff_dw*wind_U10_speed_grid[0][i - 1];
			interp_wa_I = coeff_up*cos(wind_U10_azimuth_grid[0][i]*PI_NUM/180.0) + coeff_dw*cos(wind_U10_azimuth_grid[0][i - 1]*PI_NUM/180.0);
			interp_wa_Q = coeff_up*sin(wind_U10_azimuth_grid[0][i]*PI_NUM/180.0) + coeff_dw*sin(wind_U10_azimuth_grid[0][i - 1]*PI_NUM/180.0);
			wind_U10_azimuth = atan2(interp_wa_Q, interp_wa_I)*180.0/PI_NUM;
		}
		if((longitude > lon_wgrid[size_lon_wgrid - 1])&&(latitude >= lat_wgrid[0])&&(latitude <= lat_wgrid[size_lat_wgrid - 1])){
			i = 1;
			while(latitude > lat_wgrid[i]) i++;
			coeff_up = (latitude - lat_wgrid[i - 1])/(lat_wgrid[i] - lat_wgrid[i - 1]);
			coeff_dw = (lat_wgrid[i] - latitude)/(lat_wgrid[i] - lat_wgrid[i - 1]);
			wind_U10_speed = coeff_up*wind_U10_speed_grid[size_lon_wgrid - 1][i] + coeff_dw*wind_U10_speed_grid[size_lon_wgrid - 1][i - 1];
			interp_wa_I = coeff_up*cos(wind_U10_azimuth_grid[size_lon_wgrid - 1][i]*PI_NUM/180.0) + coeff_dw*cos(wind_U10_azimuth_grid[size_lon_wgrid - 1][i - 1]*PI_NUM/180.0);
			interp_wa_Q = coeff_up*sin(wind_U10_azimuth_grid[size_lon_wgrid - 1][i]*PI_NUM/180.0) + coeff_dw*sin(wind_U10_azimuth_grid[size_lon_wgrid - 1][i - 1]*PI_NUM/180.0);
			wind_U10_azimuth = atan2(interp_wa_Q, interp_wa_I)*180.0/PI_NUM;
		}
		if((longitude >= lon_wgrid[0])&&(longitude <= lon_wgrid[size_lon_wgrid - 1])&&(latitude < lat_wgrid[0])){
			i = 1;
			while(longitude > lon_wgrid[i]) i++;
			coeff_up = (longitude - lon_wgrid[i - 1])/(lon_wgrid[i] - lon_wgrid[i - 1]);
			coeff_dw = (lon_wgrid[i] - longitude)/(lon_wgrid[i] - lon_wgrid[i - 1]);
			wind_U10_speed = coeff_up*wind_U10_speed_grid[i][0] + coeff_dw*wind_U10_speed_grid[i - 1][0];
			interp_wa_I = coeff_up*cos(wind_U10_azimuth_grid[i][0]*PI_NUM/180.0) + coeff_dw*cos(wind_U10_azimuth_grid[i - 1][0]*PI_NUM/180.0);
			interp_wa_Q = coeff_up*sin(wind_U10_azimuth_grid[i][0]*PI_NUM/180.0) + coeff_dw*sin(wind_U10_azimuth_grid[i - 1][0]*PI_NUM/180.0);
			wind_U10_azimuth = atan2(interp_wa_Q, interp_wa_I)*180.0/PI_NUM;
		}
		if((longitude >= lon_wgrid[0])&&(longitude <= lon_wgrid[size_lon_wgrid - 1])&&(latitude > lat_wgrid[size_lat_wgrid - 1])){
			i = 1;
			while(longitude > lon_wgrid[i]) i++;
			coeff_up = (longitude - lon_wgrid[i - 1])/(lon_wgrid[i] - lon_wgrid[i - 1]);
			coeff_dw = (lon_wgrid[i] - longitude)/(lon_wgrid[i] - lon_wgrid[i - 1]);
			wind_U10_speed = coeff_up*wind_U10_speed_grid[i][size_lat_wgrid - 1] + coeff_dw*wind_U10_speed_grid[i - 1][size_lat_wgrid - 1];
			interp_wa_I = coeff_up*cos(wind_U10_azimuth_grid[i][size_lat_wgrid - 1]*PI_NUM/180.0) + coeff_dw*cos(wind_U10_azimuth_grid[i - 1][size_lat_wgrid - 1]*PI_NUM/180.0);
			interp_wa_Q = coeff_up*sin(wind_U10_azimuth_grid[i][size_lat_wgrid - 1]*PI_NUM/180.0) + coeff_dw*sin(wind_U10_azimuth_grid[i - 1][size_lat_wgrid - 1]*PI_NUM/180.0);
			wind_U10_azimuth = atan2(interp_wa_Q, interp_wa_I)*180.0/PI_NUM;
		}
	}
	return;
}

void Reflecting_surface::disable_wind_grid( void ){
	int i;
	//Free memory
	if(size_lon_wgrid > 0){
		for(i=0; i<size_lon_wgrid; i++){
			free(wind_U10_speed_grid[i]);
			free(wind_U10_azimuth_grid[i]);
		}
		free(wind_U10_speed_grid);
		free(wind_U10_azimuth_grid);
		free(lon_wgrid);
		free(lat_wgrid);
	}
	size_lon_wgrid = 0;
	size_lat_wgrid = 0;
	use_wind_grid = false;
	return;
}

bool Reflecting_surface::get_wind_grid_status( void ){
	return use_wind_grid;
}