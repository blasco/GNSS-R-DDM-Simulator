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

class Reflecting_surface
{
	double freq_GHz;
	double **surface_spectrum;
	double *kx_spec;
	double *ky_spec;
	bool omnidirectional_spec;
	int n_spectrum;
	double delta_k;
	double k_threshold;
	double **wind_U10_speed_grid;
	double **wind_U10_azimuth_grid;
	double *lon_wgrid;
	double *lat_wgrid;
	int size_lon_wgrid;
	int size_lat_wgrid;
	bool use_wind_grid;
  public:
	//Relative permittivity
	double epsilon_real;
	double epsilon_imag;
	//Roughness
	double mss_x;
	double mss_y;
	double sigma_z;
	double c21_coeff;
	double c03_coeff;
	double wind_U10_speed;
	double wind_U10_azimuth;
	std::string medium;
	Reflecting_surface( void );
	void dump_parameters( void );
	void set_frequency( double frequency_GHz );
	void set_k_threshold( double k_threshold_in );
	void set_k_threshold_Brown( double incidence_deg );
	double get_k_threshold( void );
	void epsilon_sea_water( double salinity, double temperature );
	void epsilon_sea_ice( double brine_volume );
	void epsilon_dry_snow( double snow_density );
	void epsilon_wet_snow( double snow_density, double water_volume );
	void compute_Rfresnel_linear( double incidence_deg, double epsilon_upLayer[2], double rvv[2], double rhh[2] );
	void compute_Rfresnel_circular( double incidence_deg, double epsilon_upLayer[2], double rco[2], double rcross[2] );
	void compute_Tfresnel_linear( double incidence_deg, double epsilon_upLayer[2], double tvv[2], double thh[2] );
	void compute_Tfresnel_circular( double incidence_deg, double epsilon_upLayer[2], double tco[2], double tcross[2] );
	void compute_sea_spectrum( int num_samples, double delta_k_in, double theta_deg, double omega );
	void set_surf_spectrum( double* spectrum_in, int dim_x, int dim_y, double* kx_spec_in, int dim_kx, double* ky_spec_in, int dim_ky, double delta_k_in );
	void set_surf_spectrum_omnidir( double* spectrum_in_omnidir, int dim_x_omnidir, double* k_spec_in_omnidir, int dim_k_omnidir, double delta_k_in );
	void get_surf_spectrum( int x, int y, double kx_ky_spectrum[3] );
	void get_surf_spectrum_omnidir( int x, double k_spectrum[2] );
	void compute_mss_from_spectrum( void );
	void compute_mss_from_wind( void );
	void set_wind_grid( double* wspeed_grid_in, int dim_wsx, int dim_wsy, double* wazim_grid_in, int dim_wax, int dim_way, double* lon_in, int dim_lon, double* lat_in, int dim_lat );
	void interp_wind_grid( double lon_in, double lat_in );
	void disable_wind_grid();
	bool get_wind_grid_status();
};
