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

#include "wavpy_global_variables.h"
#include "ancillary_functions.h"
#include "specular_geometry.h"
#include "reflecting_surface.h"
#include "rf_front_end.h"
#include "gnss_composite.h"
#include "waveform_power.h"
#include "ZavorotnyVoronovichModel.h"
#include "gsl/gsl_fft_complex.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

ZaVoModel_GNSSR::ZaVoModel_GNSSR( void ){
	double antenna_vector_E_up[3], antenna_vector_H_up[3];
	dump_isolines_data = false;
	size_ddm_stored[0] = 0;
	size_ddm_stored[1] = 0;
	len_cov_stored = 0;
	polarization = 'L';
	exponent_wav_model_length = 10;
	num_angles = 120;
	wav_length = 256;
	ddm_half_dopplers = 0;
	sampling_rate = 80000000.0;
	delta_doppler = 0.0;
	delta_freq = 0.0;
	coherent_integration = 0.001; // 1msec
	receiver_Down.set_receiver_params(3.0, 200.0, 3.0, 12000000.0, 1);
	receiver_Up.set_receiver_params(3.0, 10.0, 3.0, 12000000.0, 1);
	antenna_vector_E_up[0] = 1.0;
	antenna_vector_E_up[1] = 0.0;
	antenna_vector_E_up[2] = 0.0;
	antenna_vector_H_up[0] = 0.0;
	antenna_vector_H_up[1] = -1.0;
	antenna_vector_H_up[2] = 0.0;
	receiver_Up.set_antenna_orientation_BF_EH(antenna_vector_E_up, antenna_vector_H_up); //Up-looking antenna
	gnss_signal.set_instrumental_params(sampling_rate, receiver_Down.get_filter_BB_BW(), 1);
	waveform_POW.set_sampling_rate(sampling_rate);
	return;
}

void ZaVoModel_GNSSR::enable_isolines_data_dump( const char* namefile ){
	dump_isolines_data = true;
	isolines_data_namefile = namefile;
	return;
}

void ZaVoModel_GNSSR::disable_isolines_data_dump( void ){
	dump_isolines_data = false;
	return;
}

void ZaVoModel_GNSSR::compute_waveform( int interferometric, int apply_curvature, int add_coherent_pow, int compute_wav_cov, int compute_ddm_cov ){
	int i, j, k, n, doppler_lines, central_doppler;
	double tau, tau_specular, tau_absolute, doppler_0, doppler_x, incidence, cosIncidence, tanIncidence;
	double a, b, c, theta, dtheta, diff_angle, prev_x1, prev_x0, first_x0, first_x1;
	double *power_tau, *ft, *prev_ft, *first_ft, *ftm, *sinc_Doppler;
	double dr, dt, d_eta_dx, d_eta_dy, d_theta_dx, d_theta_dy, jacobian, max_model, power_wav, dist_TR;
	double sinc_arg, sigma_0, directivity, isotropic_factor, noise_floor, eta_factor, power_trans, pow_coherent, reflect_pow, power_direct;
	double x[3], xR[3], xT[3], n_incident[3], n_scattered[3], vec_TR[3];
	double posR_ECEF_km[3], posT_ECEF_km[3], velR_ECEF_kms[3], velT_ECEF_kms[3], posS_ECEF_km[3];
	double posR_local[3], velR_local[3], posT_local[3], velT_local[3];
	double ant_E_BF[3], ant_H_BF[3], ant_k_BF[3], ant_E_local[3], ant_H_local[3], ant_k_local[3], angles_PhiTheta[2], rco[2], rcross[2];
	int wav_model_length = int(pow(2.0, double(exponent_wav_model_length)));
	double waveform[wav_length];
	double waveform_model[wav_model_length];
	double **power;
	double power_single_doppler[wav_model_length];
	double *lambda_func, *range_lambda, *lambda_conv;
	//Curvature approximation
	double delay1, theta_origin, rad1, rad2, delta;
	//Local constants
	double tau_approx_zero = 0.000001; // It might cause problems when computing sigma_0
	double atmospheric_loss = pow(10.0, (0.5/10.0));
	double epsilon_air[2];
	epsilon_air[0] = 1.0;
	epsilon_air[1] = 0.0;
	//Wind grid
	double ECEF2ENU_S[3][3];
	double x_ECEF_km[3], x_ECEF_rot[3], PhiLambdaH_S[3];
	double x_lon, x_lat, x_height;
	//Binary surface file
	struct toDump{
		float x, y, lon, lat, tau, sigma_0, directivity, doppler_x;
	} isolineRecord;
	FILE *fpo_data;
	if(ddm_half_dopplers < 0){
		printf("ERROR! Doppler lines not valid (< 0).");
		return;
	}
	bool cov_mode = false;
	bool cov_ddm_comp = false;
	if((compute_wav_cov == 1)||(compute_ddm_cov == 1)){
		cov_mode = true;
		if(compute_ddm_cov == 1){
			cov_ddm_comp = true;
		}
	}
	int cov_ddm_factor = 1;
	if(cov_ddm_comp){
		cov_ddm_factor = 5;
		doppler_lines = 1 + ddm_half_dopplers*4*cov_ddm_factor; //= 1 + ddm_half_dopplers*2 + 4*int(1.0/(coherent_integration*delta_doppler));
	}else{
		doppler_lines = 1 + ddm_half_dopplers*2;
	}
	central_doppler = doppler_lines/2;
	power = (double **) calloc(doppler_lines, sizeof(double *));
	for(i=0; i<doppler_lines; i++){
		power[i] = (double *) calloc(wav_model_length, sizeof(double));
	}
	if(size_ddm_stored[0] > 0){
		for(i=0; i<size_ddm_stored[0]; i++){
			free(ddm[i]);
		}
		free(ddm);
		size_ddm_stored[0] = 0;
		size_ddm_stored[1] = 0;
	}
	if(ddm_half_dopplers > 0){
		ddm = (double **) calloc(wav_length, sizeof(double *));
		for(i=0; i<wav_length; i++){
			ddm[i] = (double *) calloc(2*ddm_half_dopplers, sizeof(double));
		}
		size_ddm_stored[0] = wav_length;
		size_ddm_stored[1] = 2*ddm_half_dopplers;
	}
	power_tau = (double *) malloc (doppler_lines*sizeof(double));
	ft = (double *) calloc(doppler_lines, sizeof(double));
	prev_ft = (double *) calloc(doppler_lines, sizeof(double));
	first_ft = (double *) calloc(doppler_lines, sizeof(double));
	ftm = (double *) calloc(doppler_lines, sizeof(double));
	sinc_Doppler = (double *) calloc(doppler_lines, sizeof(double));
	if(((!cov_mode) && (len_cov_stored > 0)) || ((cov_mode) && (!cov_ddm_comp) && (len_cov_stored != wav_length) && (len_cov_stored > 0)) || ((cov_ddm_comp) && (len_cov_stored != wav_length*(1 + ddm_half_dopplers*2)) && (len_cov_stored > 0))){
		for(i=0; i<len_cov_stored; i++){
			free(cov[i]);
			free(chol[i]);
		}
		free(cov);
		free(chol);
		len_cov_stored = 0;
	}
	if(cov_mode){
		if(ddm_half_dopplers < 5){
			printf("ERROR! Not enough Doppler lines for covariance computation (< 5).\n");
			return;
		}
		if(cov_ddm_comp){
			n = 1 + ddm_half_dopplers*2;
		}else{
			n = 1;
		}
		if(len_cov_stored != wav_length*n){
			cov = (double **) calloc(wav_length*n, sizeof(double *));
			chol = (double **) calloc(wav_length*n, sizeof(double *));
			for(i=0; i<wav_length*n; i++){
				cov[i] = (double *) calloc(wav_length*n, sizeof(double));
				chol[i] = (double *) calloc(wav_length*n, sizeof(double));
			}
			len_cov_stored = wav_length*n;
		}
	}
	//==================================================================================================
	if(dump_isolines_data){
		if((fpo_data = fopen(isolines_data_namefile.c_str(), "wb")) == NULL){
			printf("ERROR! isolines_data_file can't be opened.\n");
			return;
		}
	}
	if((polarization != 'L')&&(polarization != 'R')){
		printf("ERROR! Polarization not implemented. Only L (for LHCP) or R (for RHCP).");
		return;
	}
	surface.set_frequency(gnss_signal.frequency/1000000000.0);
	gnss_signal.set_instrumental_params(sampling_rate, receiver_Down.get_filter_BB_BW(), 1);
	if((wav_model_length < 3*gnss_signal.lambda_size)||(wav_model_length < wav_length)){
		printf("ERROR! wav_model_length is too short. Increase exponent_wav_model_length.\n");
		return;
	}
	lambda_func = (double *) malloc (gnss_signal.lambda_size*sizeof(double));
	range_lambda = (double *) malloc (gnss_signal.lambda_size*sizeof(double));
	lambda_conv = (double *) malloc (gnss_signal.lambda_size*sizeof(double));
	waveform_POW.set_sampling_rate(sampling_rate);
	waveform_POW.set_init_range(geometry.geometric_delay*1000.0 + double(-(gnss_signal.lambda_size/2) - 1)*C_LIGHT/sampling_rate);
	if(geometry.elevation < 0.0){
		printf("ERROR! Elevation must be > 0.0\n");
		return;
	}
	isotropic_factor = pow((C_LIGHT/gnss_signal.frequency), 2.0)/pow((4.0*PI_NUM), 3.0);
	incidence = (90.0 - geometry.elevation);
	cosIncidence = cos(incidence*PI_NUM/180.0);
	tanIncidence = tan(incidence*PI_NUM/180.0);
	geometry.get_ECEFpos_Receiver(posR_ECEF_km);
	geometry.get_ECEFvel_Receiver(velR_ECEF_kms);
	geometry.get_ECEFpos_Transmitter(posT_ECEF_km);
	geometry.get_ECEFvel_Transmitter(velT_ECEF_kms);
	geometry.get_ECEFpos_Specular(posS_ECEF_km);
	get_local_geometry_vectors(geometry.azimuthT, posS_ECEF_km, posR_ECEF_km, velR_ECEF_kms, posT_ECEF_km, velT_ECEF_kms, posR_local, velR_local, posT_local, velT_local);
	for(k=0; k<3; k++){
		vec_TR[k] = posR_local[k] - posT_local[k];
	}
	dist_TR = norm3vec(vec_TR);
	if(surface.get_wind_grid_status()){ //There is a wind grid to interpolate
		XYZ2GEOID(posS_ECEF_km, PhiLambdaH_S, ECEF2ENU_S);
	}
	//Computation of transmitted power PtGt
	power_trans = Compute_power_trans(gnss_signal, geometry.elevation, atmospheric_loss, norm3vec(posT_local));
	if(power_trans == -1.0){
		return;
	}
	receiver_Down.get_antenna_orientation_BF(ant_E_BF, ant_H_BF, ant_k_BF);
	geometry.rotate_vector_BF_to_local(ant_E_BF, ant_E_local);
	geometry.rotate_vector_BF_to_local(ant_H_BF, ant_H_local);
	geometry.rotate_vector_BF_to_local(ant_k_BF, ant_k_local);
	for(k=0; k<3; k++){
		n_incident[k] = -posT_local[k]/norm3vec(posT_local);
		n_scattered[k] = posR_local[k]/norm3vec(posR_local);
	}
	doppler_0 = Doppler_func(n_scattered, velR_local, n_incident, velT_local, gnss_signal.frequency);
	tau_specular = geometry.geometric_delay*1000.0;
	//Coherent reflection power component
	pow_coherent = 0.0;
	if(add_coherent_pow == 1){
		AntGain_PhiTheta_angles(ant_k_local, ant_E_local, ant_H_local, n_scattered, angles_PhiTheta);
		directivity = pow(10.0, (receiver_Down.get_PhiTheta_gain_dB(angles_PhiTheta[0], angles_PhiTheta[1])/10.0));
		surface.compute_Rfresnel_circular(incidence, epsilon_air, rco, rcross);
		reflect_pow = 0.0;
		if(polarization == 'L'){
			reflect_pow = rcross[0]*rcross[0] + rcross[1]*rcross[1];
		}
		if(polarization == 'R'){
			reflect_pow = rco[0]*rco[0] + rco[1]*rco[1];
		}
		sinc_arg = PI_NUM*delta_freq*coherent_integration;
		if(sinc_arg == 0.0){
			sinc_Doppler[central_doppler] = 1.0;
		}else{
			sinc_Doppler[central_doppler] = pow((sin(sinc_arg)/sinc_arg), 2.0);
		}
		pow_coherent = power_trans*reflect_pow*exp(-1.0*pow(4.0*PI_NUM*surface.sigma_z*cosIncidence*gnss_signal.frequency/C_LIGHT, 2.0))*directivity*sinc_Doppler[central_doppler]*pow(C_LIGHT*coherent_integration/(4.0*PI_NUM*gnss_signal.frequency*(norm3vec(posT_local) + norm3vec(posR_local))), 2.0)/atmospheric_loss;
	}
	//Computation power distribution over reflecting surface
	//origin located at sample = gnss_signal.lambda_size
	for(i=gnss_signal.lambda_size; i<(wav_model_length - gnss_signal.lambda_size); i++){
		tau = double(i - gnss_signal.lambda_size)*C_LIGHT/sampling_rate;
		if(i == gnss_signal.lambda_size){
			tau = tau_approx_zero;
		}
		memset(power_tau, 0, sizeof(double)*doppler_lines);
		tau_absolute = tau + tau_specular;
		a = (tau_absolute/(cosIncidence*cosIncidence))*sqrt(1.0 - tau_specular/tau_absolute);
		b = a*cosIncidence;
		c = (tanIncidence/cosIncidence)*tau;
		x[2] = 0.0;
		for(j=0; j<num_angles; j++){
			theta = double(j)*2.0*PI_NUM/num_angles;
			x[0] = b*cos(theta);
			x[1] = a*sin(theta) + c;
			//Curvature approximation
			if(apply_curvature == 1){
				rad1 = sqrt(x[0]*x[0] + x[1]*x[1]);
				theta_origin = atan2(x[1], x[0]);
				x[2] = sqrt((A_EARTH_SEMIAXIS_KM*1000.0)*(A_EARTH_SEMIAXIS_KM*1000.0) - rad1*rad1) - A_EARTH_SEMIAXIS_KM*1000.0;
				for(k=0; k<3; k++){
					xR[k] = posR_local[k] - x[k];
					xT[k] = posT_local[k] - x[k];
				}
				dr = norm3vec(xR);
				dt = norm3vec(xT);
				delay1 = dr + dt - dist_TR - tau_specular;
				rad2 = sqrt(tau/delay1)*rad1;
				x[0] = rad2*cos(theta_origin);
				x[1] = rad2*sin(theta_origin);
				x[2] = sqrt((A_EARTH_SEMIAXIS_KM*1000.0)*(A_EARTH_SEMIAXIS_KM*1000.0) - rad2*rad2) - A_EARTH_SEMIAXIS_KM*1000.0;
			}
			//
			for(k=0; k<3; k++){
				xR[k] = posR_local[k] - x[k];
				xT[k] = posT_local[k] - x[k];
			}
			dr = norm3vec(xR);
			dt = norm3vec(xT);
			for(k=0; k<3; k++){
				n_incident[k] = -xT[k]/dt;
				n_scattered[k] = xR[k]/dr;
			}
			d_eta_dx = x[0]*(1.0/dt + 1.0/dr);
			d_eta_dy = (x[1] - posR_local[1])/dr + (x[1] - posT_local[1])/dt;
			d_theta_dx = x[1]/(x[0]*x[0] + x[1]*x[1]);
			d_theta_dy = -x[0]/(x[0]*x[0] + x[1]*x[1]);
			jacobian = fabs(1./(d_eta_dx*d_theta_dy - d_eta_dy*d_theta_dx));
			//Compute S (doppler sinc)
			doppler_x = Doppler_func(n_scattered, velR_local, n_incident, velT_local, gnss_signal.frequency) - doppler_0 - delta_freq;
			for(n=0; n<doppler_lines; n++){
				if(cov_mode){
					if(fabs(doppler_x - double(n - central_doppler)*delta_doppler/(double(cov_ddm_factor))) < (delta_doppler/(double(cov_ddm_factor)))/2){
						sinc_Doppler[n] = 1.0;
					}else{
						sinc_Doppler[n] = 0.0;
					}
				}else{
					sinc_arg = PI_NUM*(doppler_x - double(n - central_doppler)*delta_doppler)*coherent_integration;
					if(sinc_arg == 0.0){
						sinc_Doppler[n] = 1.0;
					}else{
						sinc_Doppler[n] = pow((sin(sinc_arg)/sinc_arg), 2.0);
					}
				}
			}
			//Compute sigma0
			if(surface.get_wind_grid_status() || dump_isolines_data){
				Local2ECEF_Rot(geometry.azimuthT, x, x_ECEF_rot, ECEF2ENU_S);
				for(k=0; k<3; k++){
					x_ECEF_km[k] = posS_ECEF_km[k] + x_ECEF_rot[k]/1000.0;
				}
				Compute_LatLonH_from_ECEF(x_ECEF_km, &x_lon, &x_lat, &x_height);
			}
			if(surface.get_wind_grid_status()){ //There is a wind grid to interpolate
				surface.interp_wind_grid(x_lon, x_lat);
				surface.compute_mss_from_wind();
			}
			sigma_0 = Sigma0_func(surface, n_scattered, n_incident, polarization, geometry.azimuthT);
			//Compute antenna gain
			//AntGain_EH_angles(ant_k_local, ant_E_local, ant_H_local, n_scattered, angles_PhiTheta);
			//directivity = pow(10.0, (receiver_Down.get_anglesEH_gain_dB(angles_PhiTheta[0], angles_PhiTheta[1])/10.0));
			AntGain_PhiTheta_angles(ant_k_local, ant_E_local, ant_H_local, n_scattered, angles_PhiTheta);
			directivity = pow(10.0, (receiver_Down.get_PhiTheta_gain_dB(angles_PhiTheta[0], angles_PhiTheta[1])/10.0));
			//Last step: contribution from differential annuli angular section (dtheta x dtau)
			if(dump_isolines_data){
				isolineRecord.x = float(x[0]);
				isolineRecord.y = float(x[1]);
				isolineRecord.lon = float(x_lon);
				isolineRecord.lat = float(x_lat);
				isolineRecord.tau = float(tau);
				isolineRecord.sigma_0 = float(sigma_0);
				isolineRecord.directivity = float(directivity);
				isolineRecord.doppler_x = float(doppler_x); 
				fwrite(&isolineRecord, sizeof(toDump), 1, fpo_data);
			}
			if(j>0){
				diff_angle = fabs(atan2(x[1], x[0]) - atan2(prev_x1, prev_x0));
				dtheta = std::min(diff_angle, (2.0*PI_NUM - diff_angle));
				prev_x0 = x[0];
				prev_x1 = x[1];
				for(n=0; n<doppler_lines; n++){
					ft[n] = directivity*isotropic_factor*(1.0/(dt*dt))*(1.0/(dr*dr))*sigma_0*sinc_Doppler[n]*jacobian;
					ftm[n] = 0.5*(ft[n] + prev_ft[n]);
					prev_ft[n] = ft[n];
					power_tau[n] = power_tau[n] + dtheta*ftm[n];
					if(j==(num_angles - 1)){
						diff_angle = fabs(atan2(first_x1, first_x0) - atan2(x[1], x[0]));
						dtheta = std::min(diff_angle, (2.0*PI_NUM - diff_angle));
						ftm[n] = 0.5*(first_ft[n] + ft[n]);
						power_tau[n] = power_tau[n] + dtheta*ftm[n];
					}
				}
			}else{
				prev_x0 = x[0];
				prev_x1 = x[1];
				first_x0 = x[0];
				first_x1 = x[1];
				for(n=0; n<doppler_lines; n++){
					first_ft[n] = directivity*isotropic_factor*(1.0/(dt*dt))*(1.0/(dr*dr))*sigma_0*sinc_Doppler[n]*jacobian;
					prev_ft[n] = first_ft[n];
				}
			}
		}
		for(n=0; n<doppler_lines; n++){
			power[n][i] = power_tau[n];
		}
	}
	if(dump_isolines_data){
		fclose(fpo_data);
	}
	//Convolution between power and lambda function
	gnss_signal.get_lambda_func(range_lambda, gnss_signal.lambda_size, lambda_func, gnss_signal.lambda_size);
	if((interferometric == 1)||(cov_mode)){ //Compute direct signal power at receiver
		receiver_Up.get_antenna_orientation_BF(ant_E_BF, ant_H_BF, ant_k_BF);
		geometry.rotate_vector_BF_to_local(ant_E_BF, ant_E_local);
		geometry.rotate_vector_BF_to_local(ant_H_BF, ant_H_local);
		geometry.rotate_vector_BF_to_local(ant_k_BF, ant_k_local);
		for(k=0; k<3; k++){
			n_incident[k] = vec_TR[k]/dist_TR;
		}
		//AntGain_EH_angles(ant_k_local, ant_E_local, ant_H_local, n_incident, angles_PhiTheta);
		//directivity = pow(10.0, (receiver_Up.get_anglesEH_gain_dB(angles_PhiTheta[0], angles_PhiTheta[1])/10.0));
		AntGain_PhiTheta_angles(ant_k_local, ant_E_local, ant_H_local, n_incident, angles_PhiTheta);
		directivity = pow(10.0, (receiver_Up.get_PhiTheta_gain_dB(angles_PhiTheta[0], angles_PhiTheta[1])/10.0));
		power_direct = power_trans*directivity*pow(C_LIGHT/(gnss_signal.frequency*4.0*PI_NUM*dist_TR), 2.0)/atmospheric_loss;
	}
	if(cov_mode){
		Compute_Covariance_DDM(power_trans*C_LIGHT/(sampling_rate*atmospheric_loss), power_direct, sampling_rate, coherent_integration, receiver_Down.get_filter_BB_BW(), pow(10.0, (receiver_Up.get_noise_pow_dBW()/10.0)), pow(10.0, (receiver_Down.get_noise_pow_dBW()/10.0)), delta_doppler, power, doppler_lines, wav_model_length, lambda_func, gnss_signal.lambda_size, waveform, ddm, cov, chol, wav_length, interferometric, cov_ddm_factor, cov_ddm_comp);
		waveform_POW.set_waveform(waveform, wav_length);
		waveform_POW.compute_delays();
	}else{
		for(i=0; i<gnss_signal.lambda_size; i++){
			if(i<(gnss_signal.lambda_size/2)){
				lambda_conv[i + (gnss_signal.lambda_size/2) + 1] = lambda_func[i];
			}else{
				lambda_conv[i - (gnss_signal.lambda_size/2)] = lambda_func[i];
			}
		}
		for(i=0; i<wav_model_length; i++){
			power_single_doppler[i] = power[central_doppler][i];
		}
		convolution_gsl(power_single_doppler, wav_model_length, lambda_conv, gnss_signal.lambda_size, waveform_model);
		max_model = 0.0;
		for(i=0; i<wav_length; i++){
			power_wav = waveform_model[i + gnss_signal.lambda_size/2]*power_trans*C_LIGHT/(sampling_rate*atmospheric_loss);
			if(max_model < power_wav){
				max_model = power_wav;
			}
		}
		//Compute floor noise power
		if(interferometric == 1){ //Interferometric case
			eta_factor = 1.0 + (1.0 + (max_model/pow(10.0, (receiver_Down.get_noise_pow_dBW()/10.0))))/(power_direct/pow(10.0, (receiver_Up.get_noise_pow_dBW()/10.0))); // eta = 1 + (1 + SNR_refl)/SNR_direct
		}else{ //Clean-replica case
			eta_factor = 1.0;
		}
		noise_floor = eta_factor*pow(10.0, (receiver_Down.get_noise_pow_dBW()/10.0))/(receiver_Down.get_filter_BB_BW()*coherent_integration);
		//Compute waveform
		for(i=0; i<wav_length; i++){
			if(i < gnss_signal.lambda_size){
				waveform[i] = fabs(waveform_model[i + gnss_signal.lambda_size/2]*power_trans*C_LIGHT/(sampling_rate*atmospheric_loss) + noise_floor + pow_coherent*lambda_func[i]);
			}else{
				waveform[i] = fabs(waveform_model[i + gnss_signal.lambda_size/2]*power_trans*C_LIGHT/(sampling_rate*atmospheric_loss) + noise_floor);
			}
		}
		waveform_POW.set_waveform(waveform, wav_length);
		waveform_POW.compute_delays();
		//Compute DDM
		for(n=0; n<doppler_lines; n++){
			if(n != central_doppler){
				for(i=0; i<wav_model_length; i++){
					power_single_doppler[i] = power[n][i];
				}
				convolution_gsl(power_single_doppler, wav_model_length, lambda_conv, gnss_signal.lambda_size, waveform_model);
				for(i=0; i<wav_length; i++){
					if(n < central_doppler){
						if(i < gnss_signal.lambda_size){
							ddm[i][n] = fabs(waveform_model[i + gnss_signal.lambda_size/2]*power_trans*C_LIGHT/(sampling_rate*atmospheric_loss) + noise_floor + pow_coherent*lambda_func[i]);
						}else{
							ddm[i][n] = fabs(waveform_model[i + gnss_signal.lambda_size/2]*power_trans*C_LIGHT/(sampling_rate*atmospheric_loss) + noise_floor);
						}
					}else{
						if(i < gnss_signal.lambda_size){
							ddm[i][n - 1] = fabs(waveform_model[i + gnss_signal.lambda_size/2]*power_trans*C_LIGHT/(sampling_rate*atmospheric_loss) + noise_floor + pow_coherent*lambda_func[i]);
						}else{
							ddm[i][n - 1] = fabs(waveform_model[i + gnss_signal.lambda_size/2]*power_trans*C_LIGHT/(sampling_rate*atmospheric_loss) + noise_floor);
						}
					}
				}
			}
		}
	}
	for(i=0; i<doppler_lines; i++){
		free(power[i]);
	}
	free(power);
	free(power_tau);
	free(ft);
	free(prev_ft); 
	free(first_ft); 
	free(ftm);
	free(sinc_Doppler);
	free(lambda_func);
	free(range_lambda);
	free(lambda_conv);
	return;
}

void ZaVoModel_GNSSR::get_DDM_doppler_slice( int doppler_index, double* dm_out, int dm_out_length ){
	int i, wav_length_stored;
	wav_length_stored = waveform_POW.get_wav_length();
	if((doppler_index == 0) && (dm_out_length != wav_length_stored)){
		printf("ERROR! DM length not valid: %d (it has to be %d)\n", dm_out_length, wav_length_stored);
		return;
	}
	if((doppler_index != 0) && (dm_out_length != size_ddm_stored[0])){
		printf("ERROR! DM length not valid: %d (it has to be %d)\n", dm_out_length, size_ddm_stored[0]);
		return;
	}
	if(abs(doppler_index) > size_ddm_stored[1]/2){
		printf("ERROR! Doppler slice not stored: %d (max absolute value is %d)\n", doppler_index, size_ddm_stored[1]/2);
		return;
	}
	if(doppler_index == 0){
		waveform_POW.get_waveform(dm_out, dm_out_length);
		return;
	}else{
		if(doppler_index < 0){
			for(i=0; i<size_ddm_stored[0]; i++){
				dm_out[i] = ddm[i][doppler_index + size_ddm_stored[1]/2];
			}
		}else{
			for(i=0; i<size_ddm_stored[0]; i++){
				dm_out[i] = ddm[i][doppler_index - 1 + size_ddm_stored[1]/2];
			}
		}
		return;
	}
}

void ZaVoModel_GNSSR::get_cov_slice( int cov_index, double* cov_out, int cov_out_length ){
	int i;
	if(cov_out_length != len_cov_stored){
		printf("ERROR! COV length not valid: %d (it has to be %d)\n", cov_out_length, len_cov_stored);
		return;
	}
	if((cov_index < 0)||(cov_index >= cov_out_length)){
		printf("ERROR! COV index out of range.\n");
		return;
	}
	for(i=0; i<len_cov_stored; i++){
		cov_out[i] = cov[cov_index][i];
	}
	return;
}

void ZaVoModel_GNSSR::get_noisy_waveform( double* wav_out, int wav_out_length, unsigned long int seed_in ){
	int i, j;
	double *mean_ddm;
	double *noisy_ddm;
	double *mean_wav;
	if(waveform_POW.get_wav_length() == 0){
		printf("ERROR! Mean waveform not stored.\n");
		return;
	}
	if(len_cov_stored == 0){
		printf("ERROR! Waveform covariance not stored.\n");
		return;
	}
	if((len_cov_stored != waveform_POW.get_wav_length())&&(len_cov_stored != waveform_POW.get_wav_length()*(size_ddm_stored[1] + 1))){
		printf("ERROR! Mean waveform and covariance stored have different length.\n");
		return;
	}
	if(wav_out_length != waveform_POW.get_wav_length()){
		printf("ERROR! Waveform length not valid: %d (it has to be %d).\n", wav_out_length, waveform_POW.get_wav_length());
		return;
	}
	mean_wav = (double *) malloc (wav_out_length*sizeof(double));
	waveform_POW.get_waveform(mean_wav, wav_out_length);
	if(len_cov_stored == waveform_POW.get_wav_length()){
		Get_noise_string(wav_out, mean_wav, wav_out_length, chol, seed_in);
	}else{
		mean_ddm = (double *) malloc (len_cov_stored*sizeof(double));
		noisy_ddm = (double *) malloc (len_cov_stored*sizeof(double));
		//Negative Doppler offset
		for(j=0; j<size_ddm_stored[1]/2; j++){
			for(i=0; i<size_ddm_stored[0]; i++){
				mean_ddm[i + j*size_ddm_stored[0]] = ddm[i][j];
			}
		}
		//Central Doppler
		for(i=0; i<size_ddm_stored[0]; i++){
			mean_ddm[i + size_ddm_stored[0]*size_ddm_stored[1]/2] = mean_wav[i];
		}
		//Positive Doppler offset
		for(j=(size_ddm_stored[1]/2); j<size_ddm_stored[1]; j++){
			for(i=0; i<size_ddm_stored[0]; i++){
				mean_ddm[i + (j + 1)*size_ddm_stored[0]] = ddm[i][j];
			}
		}
		Get_noise_string(noisy_ddm, mean_ddm, len_cov_stored, chol, seed_in);
		free(mean_ddm);
		for(i=0; i<wav_out_length; i++){
			wav_out[i] = noisy_ddm[i + wav_out_length*size_ddm_stored[1]/2];
		}
		free(noisy_ddm);
	}
	free(mean_wav);
	return;
}

void ZaVoModel_GNSSR::get_noisy_DDM( double* ddm_row_out, int ddm_row_out_length, unsigned long int seed_in ){
	int i, j;
	double *mean_ddm;
	double *mean_wav;
	if(waveform_POW.get_wav_length() == 0){
		printf("ERROR! Mean waveform not stored.\n");
		return;
	}
	if((size_ddm_stored[0] == 0) || (size_ddm_stored[1] == 0)){
		printf("ERROR! Mean DDM not stored.\n");
		return;
	}
	if(size_ddm_stored[0] != waveform_POW.get_wav_length()){
		printf("ERROR! Mean waveform and mean DDM stored have different length.\n");
		return;
	}
	if(len_cov_stored == 0){
		printf("ERROR! DDM covariance not stored.\n");
		return;
	}
	if(len_cov_stored != size_ddm_stored[0]*(size_ddm_stored[1] + 1)){
		printf("ERROR! Mean DDM and covariance stored have different length.\n");
		return;
	}
	if(ddm_row_out_length != len_cov_stored){
		printf("ERROR! DDM length not valid: %d (it has to be %d)\n", ddm_row_out_length, len_cov_stored);
		return;
	}
	mean_ddm = (double *) malloc (len_cov_stored*sizeof(double));
	//Negative Doppler offset
	for(j=0; j<size_ddm_stored[1]/2; j++){
		for(i=0; i<size_ddm_stored[0]; i++){
			mean_ddm[i + j*size_ddm_stored[0]] = ddm[i][j];
		}
	}
	//Central Doppler
	mean_wav = (double *) malloc (waveform_POW.get_wav_length()*sizeof(double));
	waveform_POW.get_waveform(mean_wav, waveform_POW.get_wav_length());
	for(i=0; i<size_ddm_stored[0]; i++){
		mean_ddm[i + size_ddm_stored[0]*size_ddm_stored[1]/2] = mean_wav[i];
	}
	free(mean_wav);
	//Positive Doppler offset
	for(j=(size_ddm_stored[1]/2); j<size_ddm_stored[1]; j++){
		for(i=0; i<size_ddm_stored[0]; i++){
			mean_ddm[i + (j + 1)*size_ddm_stored[0]] = ddm[i][j];
		}
	}
	Get_noise_string(ddm_row_out, mean_ddm, len_cov_stored, chol, seed_in);
	free(mean_ddm);
	return;
}


double Compute_power_trans( GNSS_composite signal, double elev, double atm_loss, double posTlocal_norm )
{
	double powT, ant_gain_GPS_dB;
	double coeff_elev_to_gain[3];
	if(((signal.weight_CA > 0.0)||(signal.weight_PY > 0.0)||(signal.weight_M > 0.0)||(signal.weight_IM > 0.0))&&(signal.weight_L1C == 0.0)){ //GPS
		//Coefficients to obtain GPS antenna gain [dB] as a function of the elevation angle [deg]
		coeff_elev_to_gain[0] = 11.805;
		coeff_elev_to_gain[1] = 0.063546;
		coeff_elev_to_gain[2] = -0.00090420;
		ant_gain_GPS_dB = elev*elev*coeff_elev_to_gain[2] + elev*coeff_elev_to_gain[1] + coeff_elev_to_gain[0];
		powT = pow((double(signal.weight_CA)*pow(10.0, ((POW_CA_Trans_dBW + ant_gain_GPS_dB)/20.0)) + double(signal.weight_PY)*pow(10.0, ((POW_PY_Trans_dBW + ant_gain_GPS_dB)/20.0)) + double(signal.weight_M)*pow(10.0, ((POW_M_Trans_dBW + ant_gain_GPS_dB)/20.0)) + double(signal.weight_IM)*pow(10.0, ((POW_IM_Trans_dBW + ant_gain_GPS_dB)/20.0))), 2.0); // power_trans = SUM_codes(Pt*Gt*weight_code^2)
		if(signal.frequency == FREQ_GPS_L5){
			powT = pow(10.0, ((POW_CA_Trans_dBW + ant_gain_GPS_dB + 3.0)/10.0)); //Power is 2 times L1 CA signal
		}
	}else{
		if((signal.weight_E1A > 0.0)||(signal.weight_E1B > 0.0)||(signal.weight_E1C > 0.0)){ //Galileo
			powT = pow((double(signal.weight_E1A)*pow(10.0, (POW_E1A_RecEarth_dBW/20.0)) + double(signal.weight_E1B)*pow(10.0, (POW_E1B_RecEarth_dBW/20.0)) + double(signal.weight_E1C)*pow(10.0, (POW_E1C_RecEarth_dBW/20.0))), 2.0);
			powT = powT*atm_loss*pow((4.0*PI_NUM*posTlocal_norm*signal.frequency/C_LIGHT), 2.0); // power_trans = Pt*Gt*weight_code^2, Pt*Gt = Pr*La*(4*pi*Rt/lambda)^2
		}else{
			if(signal.weight_B1I > 0.0){ //BeiDou
				powT = double(signal.weight_B1I*signal.weight_B1I)*atm_loss*pow(10.0, (POW_B1I_RecEarth_dBW/10.0))*pow((4.0*PI_NUM*posTlocal_norm*signal.frequency/C_LIGHT), 2.0); // power_trans = Pt*Gt*weight_code^2, Pt*Gt = Pr*La*(4*pi*Rt/lambda)^2
			}else{
				if(signal.weight_L1C > 0.0){ //QZSS
					powT = pow((double(signal.weight_CA)*pow(10.0, (POW_CA_RecEarth_dBW/20.0)) + double(signal.weight_L1C)*pow(10.0, (POW_L1C_RecEarth_dBW/20.0))), 2.0);
					powT = powT*atm_loss*pow((4.0*PI_NUM*posTlocal_norm*signal.frequency/C_LIGHT), 2.0); // power_trans = Pt*Gt*weight_code^2, Pt*Gt = Pr*La*(4*pi*Rt/lambda)^2
				}else{
					printf("ERROR! All code weigths are <= 0.0! Waveform is not computed.\n");
					return -1.0;
				}
			}
		}
	}
	return powT;
}

double Doppler_func( double n_scat[3], double velR[3], double n_inc[3], double velT[3], double frequency )
{
	int i;
	double doppler = 0.0;
	for(i=0; i<3; i++){
		doppler = doppler + (velT[i]*n_inc[i] - velR[i]*n_scat[i])*frequency/C_LIGHT;
	}
	return doppler;
}

double Sigma0_func( Reflecting_surface surf, double n_scat[3], double n_inc[3], char pol, double azimuthT )
{
	int i;
	double sigma, q_mod, theta_cos, reflect_pow, pdf_q;
	double q[3], rco[2], rcross[2], epsilon_air[2];
	epsilon_air[0] = 1.0;
	epsilon_air[1] = 0.0;
	for(i=0; i<3; i++){
		q[i] = n_scat[i] - n_inc[i];
	}
	q_mod = norm3vec(q);
	theta_cos = 0.0;
	if(q_mod > 0.0){
		for(i=0; i<3; i++){
			theta_cos = theta_cos + n_scat[i]*q[i]/q_mod;
		}
	}
	if(theta_cos > 1.0){
		printf("ERROR! Problems with approximation at tau = 0.0 when computing sigma_0\n");
		return 0.0;
	}
	surf.compute_Rfresnel_circular(acos(theta_cos)*180.0/PI_NUM, epsilon_air, rco, rcross);
	reflect_pow = 0.0;
	if(pol == 'L'){
		reflect_pow = rcross[0]*rcross[0] + rcross[1]*rcross[1];
	}
	if(pol == 'R'){
		reflect_pow = rco[0]*rco[0] + rco[1]*rco[1];
	}
	pdf_q = PDFunction(surf, q, azimuthT);
	sigma = PI_NUM*reflect_pow*pdf_q*pow((q_mod/q[2]), 4.0);
	return sigma;
}

double PDFunction( Reflecting_surface surf, double q[3], double azimuthT )
{
	double pdf, phi;
	double qx[2], qxiw[2];
	//Anisotropy as a function of wind azimuth. If mss_x = mss_y (isotropy), phi has no effect over pdf.
	qx[0] = -q[0]/q[2];
	qx[1] = -q[1]/q[2];
	phi = (90.0 + azimuthT - surf.wind_U10_azimuth)*PI_NUM/180.0;
	qxiw[0] = (qx[0]*cos(phi) + qx[1]*sin(phi))/sqrt(surf.mss_x);
	qxiw[1] = (qx[1]*cos(phi) - qx[0]*sin(phi))/sqrt(surf.mss_y);
	//Two-dimensional Gaussian distribution
	pdf = (1.0/(2.0*PI_NUM*sqrt(surf.mss_x*surf.mss_y)))*exp(-1.0*(qxiw[0]*qxiw[0] + qxiw[1]*qxiw[1])/2.0);
	if(strncmp(surf.medium.c_str(), "Sea water", 9) == 0){ //Gram-Charlier distribution (only for sea water)
		pdf = pdf*(1.0 + (surf.c21_coeff/2.0)*(qxiw[1]*qxiw[1] - 1.0)*qxiw[0] + (surf.c03_coeff/6.0)*(qxiw[0]*qxiw[0]*qxiw[0] - 3.0*qxiw[0]));
	}
	return pdf;
}

void AntGain_EH_angles( double ant_k[3], double ant_E[3], double ant_H[3], double n_target_to_ant[3], double out_angles[2] )
{
	int i;
	double angle_E_ant, angle_H_ant, sign_E, sign_H;
	double vec_x_Ek_plane[3], vec_x_Hk_plane[3], n_ant_to_target[3];
	double x_over_k, x_over_E, x_over_H, cos_angle_E_ant, cos_angle_H_ant;
	for(i=0; i<3; i++){
		n_ant_to_target[i] = -n_target_to_ant[i];
	}
	x_over_k = scalar3Prod(ant_k, n_ant_to_target);
	x_over_E = scalar3Prod(ant_E, n_ant_to_target);
	if(x_over_E >= 0.0){
		sign_E = 1.0;
	}else{
		sign_E = -1.0;
	}
	x_over_H = scalar3Prod(ant_H, n_ant_to_target);
	if(x_over_H >= 0.0){
		sign_H = 1.0;
	}else{
		sign_H = -1.0;
	}
	for(i=0; i<3; i++){
		vec_x_Ek_plane[i] = ant_E[i]*x_over_E + ant_k[i]*x_over_k;
		vec_x_Hk_plane[i] = ant_H[i]*x_over_H + ant_k[i]*x_over_k;
	}
	cos_angle_E_ant = scalar3Prod(ant_k, vec_x_Ek_plane)/(norm3vec(vec_x_Ek_plane)*norm3vec(ant_k));
	if(fabs(cos_angle_E_ant) > 1.0){
		angle_E_ant = sign_E*180.0;
	}else{
		angle_E_ant = sign_E*acos(cos_angle_E_ant)*180.0/PI_NUM;
	}
	cos_angle_H_ant = scalar3Prod(ant_k, vec_x_Hk_plane)/(norm3vec(vec_x_Hk_plane)*norm3vec(ant_k));
	if(fabs(cos_angle_H_ant) > 1.0){
		angle_H_ant = sign_H*180.0;
	}else{
		angle_H_ant = sign_H*acos(cos_angle_H_ant)*180.0/PI_NUM;
	}
	out_angles[0] = angle_E_ant;
	out_angles[1] = angle_H_ant;
	return;
}

void AntGain_PhiTheta_angles( double ant_k[3], double ant_E[3], double ant_H[3], double n_target_to_ant[3], double out_angles[2] )
{
	int i;
	double angle_phi_ant, angle_theta_ant;
	double n_ant_to_target[3];
	double x_over_k, x_over_E, x_over_H;
	for(i=0; i<3; i++){
		n_ant_to_target[i] = -n_target_to_ant[i];
	}
	x_over_k = scalar3Prod(ant_k, n_ant_to_target);
	x_over_E = scalar3Prod(ant_E, n_ant_to_target);
	x_over_H = scalar3Prod(ant_H, n_ant_to_target);
	angle_theta_ant = acos(x_over_k)*180.0/PI_NUM;
	angle_phi_ant = atan2(x_over_H, x_over_E)*180.0/PI_NUM;
	out_angles[0] = angle_phi_ant;
	out_angles[1] = angle_theta_ant;
	return;
}

void Compute_Covariance_DDM( double power_direct_surf, double power_direct_rec, double sampling_rate, double coh_time, double BW, double pow_nd, double pow_nr, double delta_doppler, double** power_surf, int doppler_len, int delay_len, double* lambda_func, int lambda_len, double* wav_out, double** ddm_out, double** cov_out, double** chol_out, int wav_len, int interferometric, int ddm_freq_factor, bool compute_ddm )
{
	int i, j, k, dtau, dk1, dk2, doppler_len_ddm, cov_len, chol_status;
	double sinc_dk1, sinc_dk2, sinc_arg, max_val, accum_cov_val;
	double* power_single_doppler;
	double* cov_sd_sr;
	double* cov_sd_nr;
	double* cov_nd_nr;
	double* lambda_conv;
	double* sinc_rec;
	double* fft_cov;
	double*** accum_cov_1; //accum_cov_1[doppler_len][wav_len][wav_len]
	double*** accum_cov_2;
	sinc_rec = (double *) malloc (delay_len*sizeof(double));
	power_single_doppler = (double *) malloc (delay_len*sizeof(double));
	cov_sd_sr = (double *) malloc (delay_len*sizeof(double));
	cov_sd_nr = (double *) malloc (delay_len*sizeof(double));
	cov_nd_nr = (double *) malloc (delay_len*sizeof(double));
	lambda_conv = (double *) malloc (lambda_len*sizeof(double));
	fft_cov = (double *) malloc (2*doppler_len*sizeof(double));
	accum_cov_1 = (double ***) calloc(doppler_len, sizeof(double **));
	//doppler_len_ddm = doppler_len - 4*int(1.0/(coh_time*delta_doppler));
	doppler_len_ddm = 1 + (doppler_len - doppler_len/2)/ddm_freq_factor;
	accum_cov_2 = (double ***) calloc(doppler_len, sizeof(double **));
	for(k=0; k<doppler_len; k++){
		accum_cov_1[k] = (double **) calloc(wav_len, sizeof(double *));
		accum_cov_2[k] = (double **) calloc(wav_len, sizeof(double *));
		for(i=0; i<wav_len; i++){
			accum_cov_1[k][i] = (double *) calloc(wav_len, sizeof(double));
			accum_cov_2[k][i] = (double *) calloc(wav_len, sizeof(double));
		}
	}
	//Correlation between direct noise and reflected noise
	for(i=0; i<delay_len; i++){
		sinc_arg = PI_NUM*double(i - delay_len/2)*2*BW/sampling_rate;
		if(sinc_arg == 0.0){
			sinc_rec[i] = 1.0;
			cov_nd_nr[i] = 1.0;
		}else{
			sinc_rec[i] = sin(sinc_arg)/sinc_arg;
			cov_nd_nr[i] = fabs(sinc_rec[i]);
		}
	}
	//Correlation between direct signal and reflected noise
	for(i=0; i<lambda_len; i++){
		if(i<(lambda_len/2)){
			lambda_conv[i + (lambda_len/2) + 1] = sqrt(lambda_func[i]);
		}else{
			lambda_conv[i - (lambda_len/2)] = sqrt(lambda_func[i]);
		}
	}
	convolution_gsl(sinc_rec, delay_len, lambda_conv, lambda_len, cov_sd_nr);
	free(sinc_rec);
	max_val = 0.0;
	for(i=0; i<delay_len; i++){
		if(fabs(cov_sd_nr[i]) > max_val){
			max_val = fabs(cov_sd_nr[i]);
		}
	}
	if(max_val > 0.0){
		for(i=0; i<delay_len; i++){
			cov_sd_nr[i] = fabs(cov_sd_nr[i])/max_val;
		}
	}
	//Correlation between direct and reflected signal
	for(k=0; k<doppler_len; k++){
		for(i=0; i<delay_len; i++){
			power_single_doppler[i] = power_surf[k][i];
		}
		for(dtau=0; dtau<(2*lambda_len - 1); dtau++){
			memset(lambda_conv, 0, sizeof(double)*lambda_len);
			for(i=0; i<lambda_len; i++){
				j = i - dtau + (lambda_len - 1);
				if(i<(lambda_len/2)){
					if((j>=0)&&(j<lambda_len)){
						lambda_conv[i + (lambda_len/2) + 1] = sqrt(lambda_func[i]*lambda_func[j]);
					}
				}else{
					if((j>=0)&&(j<lambda_len)){
						lambda_conv[i - (lambda_len/2)] = sqrt(lambda_func[i]*lambda_func[j]);
					}
				}
			}
			convolution_gsl(power_single_doppler, delay_len, lambda_conv, lambda_len, cov_sd_sr);
			for(i=0; i<wav_len; i++){
				j = i - dtau + (lambda_len - 1);
				if((j>=0)&&(j<wav_len)){
					accum_cov_1[k][i][j] = fabs(cov_sd_sr[i + lambda_len/2])*power_direct_surf*power_direct_rec;
				}
			}
		}
	}
	free(power_single_doppler);
	free(lambda_conv);
	free(cov_sd_sr);
	//Add previous noise contributions
	for(k=0; k<doppler_len; k++){
		for(i=0; i<wav_len; i++){
			for(j=0; j<wav_len; j++){
				if(((j - i + delay_len/2)>=0)&&((j - i + delay_len/2)<delay_len)){
					accum_cov_1[k][i][j] = accum_cov_1[k][i][j] + double(interferometric)*2.0*(delta_doppler/double(ddm_freq_factor))*pow_nd*pow_nr*cov_nd_nr[j - i + delay_len/2]/BW + (delta_doppler/double(ddm_freq_factor))*power_direct_rec*pow_nr*cov_sd_nr[j - i + delay_len/2]/BW;
				}
			}
		}
	}
	free(cov_nd_nr);
	free(cov_sd_nr);
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc(doppler_len);
	workspace = gsl_fft_complex_workspace_alloc(doppler_len);
	if(compute_ddm){
		cov_len = wav_len*doppler_len_ddm;
	}else{
		cov_len = wav_len;
	}
	gsl_matrix *cov_init = gsl_matrix_alloc (cov_len, cov_len);
	//Doppler intersections and FFT-1 to convert to time domain
	for(dk1=0; dk1<doppler_len_ddm; dk1++){
		for(dk2=0; dk2<doppler_len_ddm; dk2++){
			if(((dk1==(doppler_len_ddm/2))&&(dk2==(doppler_len_ddm/2)))||(compute_ddm)){
				for(i=0; i<wav_len; i++){
					for(j=0; j<wav_len; j++){
						for(k=0; k<doppler_len; k++){
							sinc_arg = PI_NUM*(delta_doppler*double(doppler_len_ddm/2 - dk1) + delta_doppler*double(k - doppler_len/2)/double(ddm_freq_factor))*coh_time;
							if(sinc_arg == 0.0){
								sinc_dk1 = 1.0;
							}else{
								sinc_dk1 = sin(sinc_arg)/sinc_arg;
							}
							sinc_arg = PI_NUM*(delta_doppler*double(doppler_len_ddm/2 - dk2) + delta_doppler*double(k - doppler_len/2)/double(ddm_freq_factor))*coh_time;
							if(sinc_arg == 0.0){
								sinc_dk2 = 1.0;
							}else{
								sinc_dk2 = sin(sinc_arg)/sinc_arg;
							}
							if(k<doppler_len/2){
								fft_cov[2*(k + 1 + doppler_len/2)] = sinc_dk1*sinc_dk2*accum_cov_1[k][i][j];
							}else{
								fft_cov[2*(k - doppler_len/2)] = sinc_dk1*sinc_dk2*accum_cov_1[k][i][j];
							}
							fft_cov[2*k + 1] = 0.0;
						}
						gsl_fft_complex_inverse(fft_cov, 1, doppler_len, wavetable, workspace);
						accum_cov_val = 0.0;
						for(k=0; k<doppler_len/2; k++){
							if(k%(int(coh_time*(double(doppler_len - 1)*delta_doppler/double(ddm_freq_factor)))) == 0){
								accum_cov_2[k][i][j] = fft_cov[2*k]*fft_cov[2*k] + fft_cov[2*k + 1]*fft_cov[2*k + 1];
								accum_cov_val = accum_cov_val + accum_cov_2[k][i][j]*coh_time;
							}
						}
						for(k=doppler_len/2; k<doppler_len; k++){
							if((doppler_len - k)%(int(coh_time*(double(doppler_len - 1)*delta_doppler/double(ddm_freq_factor)))) == 0){
								accum_cov_2[k][i][j] = fft_cov[2*k]*fft_cov[2*k] + fft_cov[2*k + 1]*fft_cov[2*k + 1];
								accum_cov_val = accum_cov_val + accum_cov_2[k][i][j]*coh_time;
							}
						}
						if(compute_ddm){
							gsl_matrix_set(cov_init, i + dk1*wav_len, j + dk2*wav_len, accum_cov_val*0.5);
						}else{
							gsl_matrix_set(cov_init, i, j, accum_cov_val*0.5);
						}
					}
				}
				if((dk1==dk2)&&(dk1!=(doppler_len_ddm/2))){
					if(dk1 < (doppler_len_ddm/2)){
						for(i=0; i<wav_len; i++){
							ddm_out[i][dk1] = sqrt(accum_cov_2[0][i][i]);
						}
					}else{
						for(i=0; i<wav_len; i++){
							ddm_out[i][dk1 - 1] = sqrt(accum_cov_2[0][i][i]);
						}
					}
				}
				if((dk1==(doppler_len_ddm/2))&&(dk2==(doppler_len_ddm/2))){
					for(i=0; i<wav_len; i++){
						wav_out[i] = sqrt(accum_cov_2[0][i][i]);
					}
				}
			}
		}
	}
	for(k=0; k<doppler_len; k++){
		for(i=0; i<wav_len; i++){
			free(accum_cov_1[k][i]);
			free(accum_cov_2[k][i]);
		}
		free(accum_cov_1[k]);
		free(accum_cov_2[k]);
	}
	free(accum_cov_1);
	free(accum_cov_2);
	free(fft_cov);
	gsl_fft_complex_wavetable_free(wavetable);
	gsl_fft_complex_workspace_free(workspace);
	//Cov_sym = (Cov + Cov')/2
	gsl_matrix *cov_init_trans = gsl_matrix_alloc (cov_len, cov_len);
	gsl_matrix_transpose_memcpy(cov_init_trans, cov_init);
	gsl_matrix_add(cov_init, cov_init_trans);
	gsl_matrix_free(cov_init_trans);
	gsl_matrix *cov_end = gsl_matrix_alloc(cov_len, cov_len);
	gsl_matrix_memcpy(cov_end, cov_init);
	//Try Cholesky decomposition
	gsl_error_handler_t *prev_handler;
	prev_handler = gsl_set_error_handler_off();
	chol_status = gsl_linalg_cholesky_decomp(cov_end);
	gsl_set_error_handler(prev_handler);
	if(chol_status != GSL_SUCCESS){
		//[V,D] = eig(Cov_sym)
		gsl_vector *eigen_vals = gsl_vector_alloc(cov_len);
		gsl_matrix *V_evecs = gsl_matrix_alloc(cov_len, cov_len);
		gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(cov_len);
		gsl_eigen_symmv(cov_init, eigen_vals, V_evecs, w);
		gsl_eigen_symmv_free(w);
		gsl_eigen_symmv_sort(eigen_vals, V_evecs, GSL_EIGEN_SORT_ABS_ASC);
		//D = abs(D)
		gsl_matrix *D_evals = gsl_matrix_alloc(cov_len, cov_len);
		gsl_matrix_set_zero(D_evals);
		for(i=0; i<cov_len; i++){
			gsl_matrix_set(D_evals, i, i, fabs(gsl_vector_get(eigen_vals, i)));
		}
		gsl_vector_free(eigen_vals);
		//Cov = V*D*V'
		gsl_matrix *DVtrans = gsl_matrix_alloc(cov_len, cov_len);
		gsl_matrix_set_zero(DVtrans);
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, D_evals, V_evecs, 0.0, DVtrans);
		gsl_matrix_free(D_evals);
		gsl_matrix_set_zero(cov_end);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V_evecs, DVtrans, 0.0, cov_end);
		gsl_matrix_free(V_evecs);
		gsl_matrix_free(DVtrans);
		gsl_matrix_memcpy(cov_init, cov_end);
		//Cholesky decomposition
		gsl_linalg_cholesky_decomp(cov_end);
	}
	for(i=0; i<cov_len; i++){
		for(j=0; j<cov_len; j++){
			cov_out[i][j] = gsl_matrix_get(cov_init, i, j);
			if(j <= i){
				chol_out[i][j] = gsl_matrix_get(cov_end, i, j);
			}else{
				chol_out[i][j] = 0.0;
			}
		}
	}
	gsl_matrix_free(cov_init);
	gsl_matrix_free(cov_end);
	return;
}

void Get_noise_string( double* out_string, double* mean_string, int len_string, double** chol_mat, int seed_rn )
{
	int i, j;
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set (r, seed_rn);
	gsl_vector *noise_string = gsl_vector_alloc(len_string);
	gsl_vector *rand_vec = gsl_vector_alloc(len_string);
	gsl_matrix *chol_lowerT = gsl_matrix_alloc(len_string, len_string);
	for(i=0; i<len_string; i++){
		gsl_vector_set(noise_string, i, mean_string[i]);
		gsl_vector_set(rand_vec, i, gsl_ran_gaussian(r, 1.0));
		for(j=0; j<len_string; j++){
			gsl_matrix_set(chol_lowerT, i, j, chol_mat[i][j]);
		}
	}
	gsl_rng_free(r);
	gsl_blas_dgemv(CblasNoTrans, 1.0, chol_lowerT, rand_vec, 1.0, noise_string);
	for(i=0; i<len_string; i++){
		out_string[i] = gsl_vector_get(noise_string, i);
	}
	gsl_vector_free(noise_string);
	gsl_vector_free(rand_vec);
	gsl_matrix_free(chol_lowerT);
	return;
}
