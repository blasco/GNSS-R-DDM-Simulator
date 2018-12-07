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

class MRSR_Model
{
	int num_layers;
	//XYZ system centered at the surface level (plane XY) under the receiver and Z positive towards it
	double height_z;
	double *depth_layer;
	double *alpha_x; //pitch in degrees
	double *alpha_y; //roll in degrees
	double *epsilon_r;
	double *epsilon_i;
  public:
	RF_FrontEnd receiver;
	GNSS_composite gnss_signal;
	Waveform_complex_cluster waveforms;
	MRSR_Model( void );
	void set_general_scenario( double height_in, double* depths_in, int num_depths, double* alpha_x_in, int num_alpha_x, double* alpha_y_in, int num_alpha_y, double* epsilon_r_in, int num_epsilon_r, double* epsilon_i_in, int num_epsilon_i );
	void set_planar_layers_scenario( double height_in, double* depths_in, int num_depths, double* epsilon_r_in, int num_epsilon_r, double* epsilon_i_in, int num_epsilon_i );
	void set_dry_snow_planar_layers_scenario( double height_in, double* depths_in, int num_depths, double* snow_dens_in, int num_snow_dens );
	void mod_height_depths( double height_in, double* depths_in, int num_depths );
	void mod_alphas( double* alpha_x_in, int num_alpha_x, double* alpha_y_in, int num_alpha_y );
	void mod_epsilon( double* epsilon_r_in, int num_epsilon_r, double* epsilon_i_in, int num_epsilon_i );
	void compute_GNSS_wavcluster( int wav_lags, int lag_direct_pos, double sampling_rate, double* elevations, int size_elevs, double* yaws, int size_yaws );
	void compute_LH_freqs_and_depths( double elev_range[2], double azim_range[2], double time_range[2], double* freq_LH, int samples_freq_LH, double* depth_LH, int samples_depth_LH );
	void compute_pow_linearPol( double elevation, double yaw, double freq, double pow_HV[2] );
};

int compute_elevs_out( double elev_in, double *elevs_out, bool *valid_elevs_out, double *depth_layer, double *alpha_plane, double *epsilon_r, double *epsilon_i, int num_layers );
double theta_layer2_Snell( double theta_layer1, double eps_r_layer1, double eps_i_layer1, double eps_r_layer2, double eps_i_layer2 );
int get_layered_powrange( double elev_in, double height_rec, double lambda, double *depth_layer, double *alpha_plane, double *epsilon_r, double *epsilon_i, double *elevs_out, bool *valid_elevs_out, int num_layers, double *range_out, double *amp_pol1, double *amp_pol2, bool lin_pol );
bool traverse_single_layer( double theta_normal, double orig_thickness, double alpha_layer, double alpha_beneath_layer, double epsilon_r, double epsilon_i, double lambda, bool up2down, double &horizontal_pos, double &delta_range, double &amp_atten );
