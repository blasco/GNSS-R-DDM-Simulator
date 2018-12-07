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
#include "reflecting_surface.h"
#include "rf_front_end.h"
#include "gnss_composite.h"
#include "waveform_complex.h"
#include "waveform_power.h"
#include "MultiRaySingleRefl_model.h"
#include <gsl/gsl_multifit.h>

MRSR_Model::MRSR_Model( void ){
	num_layers = 0;
	height_z = 0.0;
	receiver.set_receiver_params(3.0, 200.0, 3.0, 12000000.0, 1);
	gnss_signal.set_instrumental_params(80000000.0, receiver.get_filter_BB_BW(), 1);
	return;
}

void MRSR_Model::set_general_scenario( double height_in, double* depths_in, int num_depths, double* alpha_x_in, int num_alpha_x, double* alpha_y_in, int num_alpha_y, double* epsilon_r_in, int num_epsilon_r, double* epsilon_i_in, int num_epsilon_i ){
	int i;
	if((num_depths <= 0)||(num_depths != num_alpha_x)||(num_depths != num_alpha_y)||(num_depths != num_epsilon_r)||(num_depths != num_epsilon_i)){
		printf("ERROR! Incorrect sizes.\n");
		return;
	}
	if(height_in <= 0.0){
		printf("ERROR! Height must be > 0.\n");
		return;
	}
	for(i=0; i<num_depths; i++){
		if(i==0){
			if(depths_in[i] < 0.0){
				printf("ERROR! Incorrect depth at layer %d.\n", i);
				return;
			}
		}else{
			if(depths_in[i] <= depths_in[i - 1]){
				printf("ERROR! Incorrect depth at layer %d.\n", i);
				return;
			}
		}
		if((fabs(alpha_x_in[i]) >= 90.0)||(fabs(alpha_y_in[i]) >= 90.0)){
			printf("ERROR! Incorrect alpha at layer %d.\n", i);
			return;
		}
		if((epsilon_r_in[i] < 0.0)||(epsilon_i_in[i] < 0.0)){
			printf("ERROR! Incorrect epsilon at layer %d.\n", i);
			return;
		}
	}
	if(num_layers != num_depths){
		if(num_layers != 0){
			free(depth_layer);
			free(alpha_x);
			free(alpha_y);
			free(epsilon_r);
			free(epsilon_i);
		}
		num_layers = num_depths;
		depth_layer = (double*) malloc (num_layers*sizeof(double));
		alpha_x = (double*) malloc (num_layers*sizeof(double));
		alpha_y = (double*) malloc (num_layers*sizeof(double));
		epsilon_r = (double*) malloc (num_layers*sizeof(double));
		epsilon_i = (double*) malloc (num_layers*sizeof(double));
	}
	height_z = height_in + depths_in[0];
	for(i=0; i<num_layers; i++){
		depth_layer[i] = depths_in[i] - depths_in[0];
		alpha_x[i] = alpha_x_in[i];
		alpha_y[i] = alpha_y_in[i];
		epsilon_r[i] = epsilon_r_in[i];
		epsilon_i[i] = epsilon_i_in[i];
	}
	return;
}

void MRSR_Model::set_planar_layers_scenario( double height_in, double* depths_in, int num_depths, double* epsilon_r_in, int num_epsilon_r, double* epsilon_i_in, int num_epsilon_i ){
	int i;
	if((num_depths <= 0)||(num_depths != num_epsilon_r)||(num_depths != num_epsilon_i)){
		printf("ERROR! Incorrect sizes.\n");
		return;
	}
	if(height_in <= 0.0){
		printf("ERROR! Height must be > 0.\n");
		return;
	}
	for(i=0; i<num_depths; i++){
		if(i==0){
			if(depths_in[i] < 0.0){
				printf("ERROR! Incorrect depth at layer %d.\n", i);
				return;
			}
		}else{
			if(depths_in[i] <= depths_in[i - 1]){
				printf("ERROR! Incorrect depth at layer %d.\n", i);
				return;
			}
		}
		if((epsilon_r_in[i] < 0.0)||(epsilon_i_in[i] < 0.0)){
			printf("ERROR! Incorrect epsilon at layer %d.\n", i);
			return;
		}
	}
	if(num_layers != num_depths){
		if(num_layers != 0){
			free(depth_layer);
			free(alpha_x);
			free(alpha_y);
			free(epsilon_r);
			free(epsilon_i);
		}
		num_layers = num_depths;
		depth_layer = (double*) malloc (num_layers*sizeof(double));
		alpha_x = (double*) malloc (num_layers*sizeof(double));
		alpha_y = (double*) malloc (num_layers*sizeof(double));
		epsilon_r = (double*) malloc (num_layers*sizeof(double));
		epsilon_i = (double*) malloc (num_layers*sizeof(double));
	}
	height_z = height_in + depths_in[0];
	for(i=0; i<num_layers; i++){
		depth_layer[i] = depths_in[i] - depths_in[0];
		alpha_x[i] = 0.0;
		alpha_y[i] = 0.0;
		epsilon_r[i] = epsilon_r_in[i];
		epsilon_i[i] = epsilon_i_in[i];
	}
	return;
}

void MRSR_Model::set_dry_snow_planar_layers_scenario( double height_in, double* depths_in, int num_depths, double* snow_dens_in, int num_snow_dens ){
	int i;
	Reflecting_surface surface;
	surface.medium = "Dry snow";
	surface.set_frequency(gnss_signal.frequency);
	if((num_depths <= 0)||(num_depths != num_snow_dens)){
		printf("ERROR! Incorrect sizes.\n");
		return;
	}
	if(height_in <= 0.0){
		printf("ERROR! Height must be > 0.\n");
		return;
	}
	for(i=0; i<num_depths; i++){
		if(i==0){
			if(depths_in[i] < 0.0){
				printf("ERROR! Incorrect depth at layer %d.\n", i);
				return;
			}
		}else{
			if(depths_in[i] <= depths_in[i - 1]){
				printf("ERROR! Incorrect depth at layer %d.\n", i);
				return;
			}
		}
		if(snow_dens_in[i] < 0.0){
			printf("ERROR! Incorrect snow density at layer %d.\n", i);
			return;
		}
	}
	if(num_layers != num_depths){
		if(num_layers != 0){
			free(depth_layer);
			free(alpha_x);
			free(alpha_y);
			free(epsilon_r);
			free(epsilon_i);
		}
		num_layers = num_depths;
		depth_layer = (double*) malloc (num_layers*sizeof(double));
		alpha_x = (double*) malloc (num_layers*sizeof(double));
		alpha_y = (double*) malloc (num_layers*sizeof(double));
		epsilon_r = (double*) malloc (num_layers*sizeof(double));
		epsilon_i = (double*) malloc (num_layers*sizeof(double));
	}
	height_z = height_in + depths_in[0];
	for(i=0; i<num_layers; i++){
		depth_layer[i] = depths_in[i] - depths_in[0];
		alpha_x[i] = 0.0;
		alpha_y[i] = 0.0;
		surface.epsilon_dry_snow(snow_dens_in[i]);
		epsilon_r[i] = surface.epsilon_real;
		epsilon_i[i] = surface.epsilon_imag;
	}
	return;
}

void MRSR_Model::mod_height_depths( double height_in, double* depths_in, int num_depths ){
	int i;
	if(height_in <= 0.0){
		printf("ERROR! Height must be > 0.\n");
		return;
	}
	if(num_depths != num_layers){
		printf("ERROR! Incorrect size.\n");
		return;
	}
	for(i=0; i<num_layers; i++){
		if(i==0){
			if(depths_in[i] < 0.0){
				printf("ERROR! Incorrect depth at layer %d.\n", i);
				return;
			}
		}else{
			if(depths_in[i] <= depths_in[i - 1]){
				printf("ERROR! Incorrect depth at layer %d.\n", i);
				return;
			}
		}
	}
	height_z = height_in + depths_in[0];
	for(i=0; i<num_layers; i++){
		depth_layer[i] = depths_in[i] - depths_in[0];
	}
	return;
}

void MRSR_Model::mod_alphas( double* alpha_x_in, int num_alpha_x, double* alpha_y_in, int num_alpha_y ){
	int i;
	if((num_alpha_x != num_layers)||(num_alpha_y != num_layers)){
		printf("ERROR! Incorrect size.\n");
		return;
	}
	for(i=0; i<num_layers; i++){
		if((fabs(alpha_x_in[i]) >= 90.0)||(fabs(alpha_y_in[i]) >= 90.0)){
			printf("ERROR! Incorrect alpha at layer %d.\n", i);
			return;
		}
	}
	for(i=0; i<num_layers; i++){
		alpha_x[i] = alpha_x_in[i];
		alpha_y[i] = alpha_y_in[i];
	}
	return;
}

void MRSR_Model::mod_epsilon( double* epsilon_r_in, int num_epsilon_r, double* epsilon_i_in, int num_epsilon_i ){
	int i;
	if((num_epsilon_r != num_layers)||(num_epsilon_i != num_layers)){
		printf("ERROR! Incorrect size.\n");
		return;
	}
	for(i=0; i<num_layers; i++){
		if((epsilon_r_in[i] < 0.0)||(epsilon_i_in[i] < 0.0)){
			printf("ERROR! Incorrect epsilon at layer %d.\n", i);
			return;
		}
	}
	for(i=0; i<num_layers; i++){
		epsilon_r[i] = epsilon_r_in[i];
		epsilon_i[i] = epsilon_i_in[i];
	}
	return;
}

void MRSR_Model::compute_GNSS_wavcluster( int wav_lags, int lag_direct_pos, double sampling_rate, double* elevations, int size_elevs, double* yaws, int size_yaws ){
	int i, j, k, num_valid_layers, num_crossed_layers, accum_wavs;
	double lambda, phase, gain_ant, range_direct;
	double rec_inc_vec[3];
	double *alpha_plane, *elevs_out, *range_out, *amp_co, *amp_cross;
	double *wav_I, *wav_Q, *lambda_func, *range_lambda;
	bool *valid_elevs_out;
	Waveform_power accum_wav_I;
	Waveform_power accum_wav_Q;
	if(wav_lags <= 0){
		printf("ERROR! Num lags not valid.\n");
		return;
	}
	if(num_layers == 0){
		printf("ERROR! There are no layers stored.\n");
		return;
	}
	if(size_elevs != size_yaws){
		printf("ERROR! Elevations and yaws must have the same length.\n");
		return;
	}
	if(size_elevs <= 0){
		printf("ERROR! Elevations and yaws must have length > 0.\n");
		return;
	}
	if((lag_direct_pos < 0)||(lag_direct_pos >= wav_lags)){
		printf("ERROR! Lag of the direct signal outside of waveform's range.\n");
		return;
	}
	lambda = C_LIGHT/gnss_signal.frequency;
	gnss_signal.set_instrumental_params(sampling_rate, receiver.get_filter_BB_BW(), 1);
	range_direct = double(lag_direct_pos - gnss_signal.lambda_size/2)*C_LIGHT/sampling_rate;
	if(gnss_signal.lambda_size > wav_lags){
		printf("ERROR! Waveform lenght (%d) too short (min=%d).\n", wav_lags, gnss_signal.lambda_size);
		return;
	}
	waveforms.initialize(size_elevs, wav_lags);
	accum_wav_I.set_sampling_rate(sampling_rate);
	accum_wav_Q.set_sampling_rate(sampling_rate);
	lambda_func = (double *) malloc (gnss_signal.lambda_size*sizeof(double));
	range_lambda = (double *) malloc (gnss_signal.lambda_size*sizeof(double));
	gnss_signal.get_lambda_func(range_lambda, gnss_signal.lambda_size, lambda_func, gnss_signal.lambda_size);
	free(range_lambda);
	alpha_plane = (double *) malloc(num_layers*sizeof(double));
	elevs_out = (double *) malloc(num_layers*sizeof(double));
	valid_elevs_out = (bool *) malloc(num_layers*sizeof(bool));
	range_out = (double *) malloc(num_layers*sizeof(double));
	amp_co = (double *) malloc(num_layers*sizeof(double));
	amp_cross = (double *) malloc(num_layers*sizeof(double));
	wav_I = (double *) malloc(wav_lags*sizeof(double));
	wav_Q = (double *) malloc(wav_lags*sizeof(double));
	for(i=0; i<size_elevs; i++){
		memset(elevs_out, 0, sizeof(double)*num_layers);
		memset(valid_elevs_out, 0, sizeof(bool)*num_layers);
		memset(range_out, 0, sizeof(double)*num_layers);
		memset(amp_co, 0, sizeof(double)*num_layers);
		memset(amp_cross, 0, sizeof(double)*num_layers);
		memset(wav_I, 0, sizeof(double)*wav_lags);
		memset(wav_Q, 0, sizeof(double)*wav_lags);
		for(j=0; j<num_layers; j++){
			alpha_plane[j] = atan(tan(alpha_x[j]*PI_NUM/180.0)*cos(yaws[i]*PI_NUM/180.0) + tan(alpha_y[j]*PI_NUM/180.0)*sin(yaws[i]*PI_NUM/180.0))*180.0/PI_NUM;
		}
		num_valid_layers = compute_elevs_out(elevations[i], elevs_out, valid_elevs_out, depth_layer, alpha_plane, epsilon_r, epsilon_i, num_layers);
		if(num_valid_layers == 0){
			printf("ERROR! Valid layers = 0 at sample %d\n", i);
			free(alpha_plane);
			free(elevs_out);
			free(valid_elevs_out);
			free(range_out);
			free(amp_co);
			free(amp_cross);
			free(wav_I);
			free(wav_Q);
			free(lambda_func);
			return;
		}
		num_crossed_layers = get_layered_powrange(elevations[i], height_z, lambda, depth_layer, alpha_plane, epsilon_r, epsilon_i, elevs_out, valid_elevs_out, num_layers, range_out, amp_co, amp_cross, false);
		if((num_valid_layers - num_crossed_layers) == 0){
			printf("ERROR! Crossed layers = Valid layers at sample %d. (%d/%d total layers)\n", i, num_crossed_layers, num_layers);
			free(alpha_plane);
			free(elevs_out);
			free(valid_elevs_out);
			free(range_out);
			free(amp_co);
			free(amp_cross);
			free(wav_I);
			free(wav_Q);
			free(lambda_func);
			return;
		}
		accum_wav_I.set_waveform(wav_I, wav_lags);
		accum_wav_Q.set_waveform(wav_Q, wav_lags);
		accum_wavs = 0;
		for(j=0; j<num_layers; j++){
			if((range_out[j] > 0.0)&&(amp_cross[j] > 0.0)){
				rec_inc_vec[0] = -height_z*cos(yaws[i]*PI_NUM/180.0)/tan(elevs_out[i]*PI_NUM/180.0);
				rec_inc_vec[1] = -height_z*sin(yaws[i]*PI_NUM/180.0)/tan(elevs_out[i]*PI_NUM/180.0);
				rec_inc_vec[2] = height_z;
				gain_ant = pow(10.0, receiver.get_incvector_gain_dB(rec_inc_vec)/20.0);
				phase = -2.0*PI_NUM*range_out[j]/lambda;
				for(k=0; k<gnss_signal.lambda_size; k++){
					wav_I[k] = amp_cross[j]*sqrt(lambda_func[k])*gain_ant*cos(phase);
					wav_Q[k] = amp_cross[j]*sqrt(lambda_func[k])*gain_ant*sin(phase);
				}
				accum_wavs ++;
				accum_wav_I.add_waveform_retracking(wav_I, wav_lags, range_direct + range_out[j], 1.0/float(accum_wavs), true);
				accum_wav_Q.add_waveform_retracking(wav_Q, wav_lags, range_direct + range_out[j], 1.0/float(accum_wavs), true);
			}
		}
		accum_wav_I.get_waveform(wav_I, wav_lags);
		accum_wav_Q.get_waveform(wav_Q, wav_lags);
		waveforms.add_waveform_scale(wav_I, wav_lags, wav_Q, wav_lags, i, double(accum_wavs));
	}
	free(alpha_plane);
	free(elevs_out);
	free(valid_elevs_out);
	free(range_out);
	free(amp_co);
	free(amp_cross);
	free(wav_I);
	free(wav_Q);
	free(lambda_func);
	return;
}

void MRSR_Model::compute_LH_freqs_and_depths( double elev_range[2], double azim_range[2], double time_range[2], double* freq_LH, int samples_freq_LH, double* depth_LH, int samples_depth_LH ){
	int i, j, ind, num_valid_layers, num_crossed_layers;
	double lambda, mid_elev, elev_rate, surface_freq, cycles_degElev, chisq;
	double *alpha_plane, *elevs_out, *amp_co, *amp_cross, *range_out, *range_diff;
	bool *valid_elevs_out;
	gsl_matrix *X, *cov;
	gsl_vector *y, *w, *c;
	int poly_order = 7;
	if(num_layers == 0){
		printf("ERROR! There are no layers stored.\n");
		return;
	}
	if(samples_freq_LH != samples_depth_LH){
		printf("ERROR! Frequencies and depths must have the same length.\n");
		return;
	}
	if(samples_freq_LH <= 0){
		printf("ERROR! Frequencies and depths must have length > 0.\n");
		return;
	}
	if(elev_range[0] == elev_range[1]){
		printf("ERROR! There is not elevation variation.\n");
		return;
	}
	if(time_range[0] == time_range[1]){
		printf("ERROR! There is not time variation.\n");
		return;
	}
	mid_elev = (elev_range[0] + elev_range[1])/2.0;
	elev_rate = (elev_range[0] - elev_range[1])/(time_range[0] - time_range[1]);
	lambda = C_LIGHT/gnss_signal.frequency;
	surface_freq = -2.0*(height_z/lambda)*cos(mid_elev*PI_NUM/180.0)*PI_NUM/180.0;
	alpha_plane = (double *) malloc(num_layers*sizeof(double));
	elevs_out = (double *) malloc(num_layers*sizeof(double));
	valid_elevs_out = (bool *) malloc(num_layers*sizeof(bool));
	amp_co = (double *) malloc(num_layers*sizeof(double));
	amp_cross = (double *) malloc(num_layers*sizeof(double));
	range_out = (double*) calloc (num_layers, sizeof(double));
	range_diff = (double*) calloc (num_layers, sizeof(double));
	//Init range case
	memset(elevs_out, 0, sizeof(double)*num_layers);
	memset(amp_co, 0, sizeof(double)*num_layers);
	memset(amp_cross, 0, sizeof(double)*num_layers);
	for(i=0; i<num_layers; i++){
		alpha_plane[i] = atan(tan(alpha_x[i]*PI_NUM/180.0)*cos(azim_range[0]*PI_NUM/180.0) + tan(alpha_y[i]*PI_NUM/180.0)*sin(azim_range[0]*PI_NUM/180.0))*180.0/PI_NUM;
	}
	num_valid_layers = compute_elevs_out(elev_range[0], elevs_out, valid_elevs_out, depth_layer, alpha_plane, epsilon_r, epsilon_i, num_layers);
	if(num_valid_layers == 0){
		printf("ERROR-1! Valid layers = 0\n");
		free(alpha_plane);
		free(elevs_out);
		free(valid_elevs_out);
		free(range_out);
		free(range_diff);
		free(amp_co);
		free(amp_cross);
		return;
	}
	num_crossed_layers = get_layered_powrange(elev_range[0], height_z, lambda, depth_layer, alpha_plane, epsilon_r, epsilon_i, elevs_out, valid_elevs_out, num_layers, range_diff, amp_co, amp_cross, false);
	if((num_valid_layers - num_crossed_layers) == 0){
		printf("ERROR-1! Crossed layers = Valid layers. (%d/%d total layers)\n", num_crossed_layers, num_layers);
		free(alpha_plane);
		free(elevs_out);
		free(valid_elevs_out);
		free(range_out);
		free(range_diff);
		free(amp_co);
		free(amp_cross);
		return;
	}
	//End range case
	memset(elevs_out, 0, sizeof(double)*num_layers);
	memset(amp_co, 0, sizeof(double)*num_layers);
	memset(amp_cross, 0, sizeof(double)*num_layers);
	for(i=0; i<num_layers; i++){
		alpha_plane[i] = atan(tan(alpha_x[i]*PI_NUM/180.0)*cos(azim_range[1]*PI_NUM/180.0) + tan(alpha_y[i]*PI_NUM/180.0)*sin(azim_range[1]*PI_NUM/180.0))*180.0/PI_NUM;
	}
	num_valid_layers = compute_elevs_out(elev_range[1], elevs_out, valid_elevs_out, depth_layer, alpha_plane, epsilon_r, epsilon_i, num_layers);
	if(num_valid_layers == 0){
		printf("ERROR-2! Valid layers = 0\n");
		free(alpha_plane);
		free(elevs_out);
		free(valid_elevs_out);
		free(range_out);
		free(range_diff);
		free(amp_co);
		free(amp_cross);
		return;
	}
	num_crossed_layers = get_layered_powrange(elev_range[1], height_z, lambda, depth_layer, alpha_plane, epsilon_r, epsilon_i, elevs_out, valid_elevs_out, num_layers, range_out, amp_co, amp_cross, false);
	free(alpha_plane);
	free(elevs_out);
	free(valid_elevs_out);
	free(amp_co);
	free(amp_cross);
	if((num_valid_layers - num_crossed_layers) == 0){
		printf("ERROR-2! Crossed layers = Valid layers. (%d/%d total layers)\n", num_crossed_layers, num_layers);
		free(range_out);
		free(range_diff);
		return;
	}
	//Compute range difference
	num_valid_layers = 0;
	for(i=0; i<num_layers; i++){
		if((range_diff[i] > 0.0)&&(range_out[i] > 0.0)){
			range_diff[i] = range_diff[i] - range_out[i];
			num_valid_layers ++;
		}else{
			range_diff[i] = 0.0;
		}
	}
	free(range_out);
	if(num_valid_layers < (poly_order + 1)){
		printf("ERROR! Valid layers < polynomial fit order\n");
		free(range_diff);
		return;
	}
	//Compute least-squares polynomial fit
	X = gsl_matrix_alloc(num_valid_layers, (poly_order + 1));
	y = gsl_vector_alloc(num_valid_layers);
	w = gsl_vector_alloc(num_valid_layers);
	c = gsl_vector_alloc((poly_order + 1));
	cov = gsl_matrix_alloc((poly_order + 1), (poly_order + 1));
	ind = 0;
	for(i=0; i<num_layers; i++){
		if(fabs(range_diff[i]) > 0.0){
			cycles_degElev = -1.0*range_diff[i]/(lambda*(elev_range[0] - elev_range[1]));
			gsl_matrix_set(X, ind, 0, 1.0);
			for(j=1; j<=poly_order; j++){
				gsl_matrix_set(X, ind, j, pow(cycles_degElev, double(j)));
			}
			gsl_vector_set(y, ind, depth_layer[i]);
			gsl_vector_set(w, ind, 1.0);
			ind ++;
		}
	}
	free(range_diff);
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(num_valid_layers, (poly_order + 1));
	gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
	gsl_multifit_linear_free(work);
	gsl_matrix_free(X);
	gsl_vector_free(y);
	gsl_vector_free(w);
	gsl_matrix_free(cov);
	for(i=0; i<samples_freq_LH; i++){
		freq_LH[i] = double(i - samples_freq_LH/2 + 1)/(double(samples_freq_LH)*elev_rate);
		if(freq_LH[i] > 0.0){
			depth_LH[i] = -height_z - 1.0;
		}else{
			if(freq_LH[i] > surface_freq){
				depth_LH[i] = -height_z + (freq_LH[i]*lambda/(-2.0*cos(mid_elev*PI_NUM/180.0)*PI_NUM/180.0));
			}else{
				depth_LH[i] = gsl_vector_get(c, 0);
				for(j=1; j<=poly_order; j++){
					depth_LH[i] = depth_LH[i] + gsl_vector_get(c, j)*pow(freq_LH[i], double(j));
				}
			}
		}
	}
	gsl_vector_free(c);
	return;
}

void MRSR_Model::compute_pow_linearPol( double elevation, double yaw, double freq, double pow_HV[2] ){
	int i, num_valid_layers, num_crossed_layers;
	double lambda, phase, yaw_rad, gain_amp_lin;
	double rec_inc_vec[3], complex_H[2], complex_V[2];
	double *alpha_plane, *elevs_out, *range_out, *amp_v, *amp_h;
	bool *valid_elevs_out;
	lambda = C_LIGHT/freq;
	yaw_rad = yaw*PI_NUM/180.0;
	alpha_plane = (double*) calloc (num_layers, sizeof(double));
	elevs_out = (double*) calloc (num_layers, sizeof(double));
	valid_elevs_out = (bool*) calloc (num_layers, sizeof(bool));
	if(num_layers == 0){
		printf("ERROR! There are no layers stored.\n");
		return;
	}
	for(i=0; i<num_layers; i++){
		alpha_plane[i] = atan(tan(alpha_x[i]*PI_NUM/180.0)*cos(yaw_rad) + tan(alpha_y[i]*PI_NUM/180.0)*sin(yaw_rad))*180.0/PI_NUM;
	}
	num_valid_layers = compute_elevs_out(elevation, elevs_out, valid_elevs_out, depth_layer, alpha_plane, epsilon_r, epsilon_i, num_layers);
	if(num_valid_layers == 0){
		printf("ERROR! Valid layers = 0\n");
		free(alpha_plane);
		free(elevs_out);
		free(valid_elevs_out);
		return;
	}
	range_out = (double*) calloc (num_layers, sizeof(double));
	amp_v = (double*) calloc (num_layers, sizeof(double));
	amp_h = (double*) calloc (num_layers, sizeof(double));
	num_crossed_layers = get_layered_powrange(elevation, height_z, lambda, depth_layer, alpha_plane, epsilon_r, epsilon_i, elevs_out, valid_elevs_out, num_layers, range_out, amp_v, amp_h, true);
	if((num_valid_layers - num_crossed_layers) == 0){
		printf("ERROR! Crossed layers = Valid layers. (%d/%d total layers)\n", num_crossed_layers, num_layers);
		free(alpha_plane);
		free(elevs_out);
		free(valid_elevs_out);
		free(range_out);
		free(amp_v);
		free(amp_h);
		return;
	}
	pow_HV[0] = 0.0;
	pow_HV[1] = 0.0;
	complex_H[0] = 0.0;
	complex_H[1] = 0.0;
	complex_V[0] = 0.0;
	complex_V[1] = 0.0;
	for(i=0; i<num_layers; i++){
		if(range_out[i] > 0.0){
			rec_inc_vec[0] = -height_z*cos(yaw_rad)/tan(elevs_out[i]*PI_NUM/180.0);
			rec_inc_vec[1] = -height_z*sin(yaw_rad)/tan(elevs_out[i]*PI_NUM/180.0);
			rec_inc_vec[2] = height_z;
			gain_amp_lin = pow(10.0, receiver.get_incvector_gain_dB(rec_inc_vec)/20.0);
			phase = -2.0*PI_NUM*range_out[i]/lambda;
			complex_H[0] = complex_H[0] + amp_h[i]*cos(phase)*gain_amp_lin;
			complex_H[1] = complex_H[1] + amp_h[i]*sin(phase)*gain_amp_lin;
			complex_V[0] = complex_V[0] + amp_v[i]*cos(phase)*gain_amp_lin;
			complex_V[1] = complex_V[1] + amp_v[i]*sin(phase)*gain_amp_lin;
		}    
	}
	pow_HV[0] = complex_H[0]*complex_H[0] + complex_H[1]*complex_H[1];
	pow_HV[1] = complex_V[0]*complex_V[0] + complex_V[1]*complex_V[1];
	if(num_valid_layers < num_layers){
		printf("WARNING! Valid layers: %d/%d\n", num_valid_layers, num_layers);
	}
	if(num_crossed_layers > 0){
		printf("WARNING! Crossed layers: %d/%d\n", num_crossed_layers, num_valid_layers);
	}
	free(alpha_plane);
	free(elevs_out);
	free(valid_elevs_out);
	free(range_out);
	free(amp_v);
	free(amp_h);
	return;
}


int compute_elevs_out( double elev_in, double *elevs_out, bool *valid_elevs_out, double *depth_layer, double *alpha_plane, double *epsilon_r, double *epsilon_i, int num_layers )
{
	int i, j, outs_ko;
	bool out_ok;
	double theta_U_in, theta_D_in, theta_U_out, theta_D_out;
	outs_ko = 0;
	//Air to first layer -> IN direction
	theta_U_in = 90.0 - elev_in - alpha_plane[0];
	if((theta_U_in >= 90.0)||(theta_U_in < 0.0)){ //No reflection in the proper direction
		return 0;
	}
	//================================Reflection -> OUT direction
	theta_U_out = theta_U_in;
	elevs_out[0] = 90.0 + alpha_plane[0] - theta_U_out;
	valid_elevs_out[0] = true;
	//================================Transmission to first layer (air to 0) -> IN direction
	theta_D_in = theta_layer2_Snell(theta_U_in, 1.0, 0.0, epsilon_r[0], epsilon_i[0]);
	if(theta_D_in == 90.0){ //No further transmission
		return 1;
	}
	// Loop layers
	for(i=1; i<num_layers; i++){
		valid_elevs_out[i] = false;
		// Arrival to intersection between layers i-1 and i from previous intersection
		theta_U_in = theta_D_in + alpha_plane[i - 1] - alpha_plane[i];
		if((theta_U_in >= 90.0)||(theta_U_in < 0.0)){ //No reflection in the proper direction
			return i;
		}
		//=======================================Reflection -> OUT direction
		theta_U_out = theta_U_in;
		out_ok = true;
		for(j=(i-1); j>0; j--){
			if(out_ok){
				theta_D_out = theta_U_out - alpha_plane[j + 1] + alpha_plane[j];
				if((theta_D_out >= 90.0)||(theta_D_out < 0.0)){ //No arrival in the proper angle
					out_ok = false;
					outs_ko ++;
				}
			}
			if(out_ok){
				theta_U_out = theta_layer2_Snell(theta_D_out, epsilon_r[j], epsilon_i[j], epsilon_r[j - 1], epsilon_i[j - 1]);
				if(theta_U_out == 90.0){ //No further transmission
					out_ok = false;
					outs_ko ++;
				}
			}
		}
		//First layer to air -> OUT direction
		if(out_ok){
			theta_D_out = theta_U_out - alpha_plane[1] + alpha_plane[0];
			if((theta_D_out >= 90.0)||(theta_D_out < 0.0)){ //No arrival in the proper angle
				out_ok = false;
				outs_ko ++;
			}
		}
		if(out_ok){
			theta_U_out = theta_layer2_Snell(theta_D_out, epsilon_r[0], epsilon_i[0], 1.0, 0.0);
			if(theta_U_out == 90.0){ //No further transmission
				outs_ko ++;
			}else{
				elevs_out[i] = 90.0 + alpha_plane[0] - theta_U_out;
				valid_elevs_out[i] = true;
			}
		}
		//=========================================Transmission from layer i-1 to layer i -> IN direction
		theta_D_in = theta_layer2_Snell(theta_U_in, epsilon_r[i - 1], epsilon_i[i - 1], epsilon_r[i], epsilon_i[i]);
		if(theta_D_in == 90.0){ //No further transmission
			return (i + 1);
		}
	}
	return (num_layers - outs_ko);
}

double theta_layer2_Snell( double theta_layer1, double eps_r_layer1, double eps_i_layer1, double eps_r_layer2, double eps_i_layer2 )
{
	double n_real_layer1, n_real_layer2, critical_angle;
	n_real_layer1 = sqrt((sqrt(eps_r_layer1*eps_r_layer1 + eps_i_layer1*eps_i_layer1) + eps_r_layer1)/2.0);
	n_real_layer2 = sqrt((sqrt(eps_r_layer2*eps_r_layer2 + eps_i_layer2*eps_i_layer2) + eps_r_layer2)/2.0);
	if((n_real_layer2/n_real_layer1) < 1.0){
		critical_angle = asin(n_real_layer2/n_real_layer1)*180.0/PI_NUM;
		if(theta_layer1 >= critical_angle){
			return 90.0; //There is no transmission
		}
	}
	return (asin(n_real_layer1*sin(theta_layer1*PI_NUM/180.0)/n_real_layer2)*180.0/PI_NUM);
}

int get_layered_powrange( double elev_in, double height_rec, double lambda, double *depth_layer, double *alpha_plane, double *epsilon_r, double *epsilon_i, double *elevs_out, bool *valid_elevs_out, int num_layers, double *range_out, double *amp_pol1, double *amp_pol2, bool lin_pol )
{
	int i, j, num_crossed_layers;
	double epsilon_incLayer[2], tvv[2], thh[2], tco[2], tcross[2], rvv[2], rhh[2], rco[2], rcross[2];
	double atten_amp_v, atten_amp_h, atten_amp_co, atten_amp_cross, accum_delay, horizontal_pos;
	double delta_range, amp_atten;
	double theta_U_in, theta_D_in, theta_U_out, theta_D_out, theta_normal;
	double alpha, theta, x_surf, z_rec, delta_range_direct;
	bool valid_path;
	Reflecting_surface layer;
	num_crossed_layers = 0;
	for(i=0; i<num_layers; i++){
		atten_amp_v = 1.0;
		atten_amp_h = 1.0;
		atten_amp_co = 1.0;
		atten_amp_cross = 1.0;
		accum_delay = 0.0;
		horizontal_pos = 0.0;
		if(valid_elevs_out[i]){
			//Receiver to first layer (surface level)
			theta_normal = 90.0 - elevs_out[i];
			valid_path = traverse_single_layer(theta_normal, height_rec, 0.0, alpha_plane[0], 1.0, 0.0, lambda, true, horizontal_pos, delta_range, amp_atten);
			accum_delay = accum_delay + delta_range;
			theta_U_out = 90.0 + alpha_plane[0] - elevs_out[i];
			layer.epsilon_real = 1.0;
			layer.epsilon_imag = 0.0;
			//Transmission till layer i
			for(j=0; j<i; j++){
				if(valid_path){
					//Layer j+1 to layer j intersection
					theta_D_out = theta_layer2_Snell(theta_U_out, layer.epsilon_real, layer.epsilon_imag, epsilon_r[j], epsilon_i[j]);
					epsilon_incLayer[0] = epsilon_r[j];
					epsilon_incLayer[1] = epsilon_i[j];
					if(lin_pol){
						layer.compute_Tfresnel_linear(theta_D_out, epsilon_incLayer, tvv, thh);
					}else{
						layer.compute_Tfresnel_circular(theta_D_out, epsilon_incLayer, tco, tcross);
					}
					//Ray propagation through layer j
					theta_normal = theta_D_out - alpha_plane[j];
					valid_path = traverse_single_layer(theta_normal, (depth_layer[j + 1] - depth_layer[j]), alpha_plane[j], alpha_plane[j + 1], epsilon_r[j], epsilon_i[j], lambda, true, horizontal_pos, delta_range, amp_atten);
					accum_delay = accum_delay + delta_range;
					if(lin_pol){
						atten_amp_v = atten_amp_v*sqrt(tvv[0]*tvv[0] + tvv[1]*tvv[1])*amp_atten;
						atten_amp_h = atten_amp_h*sqrt(thh[0]*thh[0] + thh[1]*thh[1])*amp_atten;
					}else{
						atten_amp_co = atten_amp_co*sqrt(tco[0]*tco[0] + tco[1]*tco[1])*amp_atten;
						atten_amp_cross = atten_amp_cross*sqrt(tco[0]*tco[0] + tco[1]*tco[1])*amp_atten; //co-polar transmission of the cross-polar component
					}
				}
				if(valid_path){
					layer.epsilon_real = epsilon_r[j];
					layer.epsilon_imag = epsilon_i[j];
					theta_U_out = theta_D_out - alpha_plane[j] + alpha_plane[j + 1];
				}
			}
			//Reflection at layer i
			if(valid_path){
				epsilon_incLayer[0] = layer.epsilon_real;
				epsilon_incLayer[1] = layer.epsilon_imag;
				layer.epsilon_real = epsilon_r[i];
				layer.epsilon_imag = epsilon_i[i];
				if(lin_pol){
					layer.compute_Rfresnel_linear(theta_U_out, epsilon_incLayer, rvv, rhh);
					atten_amp_v = atten_amp_v*sqrt(rvv[0]*rvv[0] + rvv[1]*rvv[1]);
					atten_amp_h = atten_amp_h*sqrt(rhh[0]*rhh[0] + rhh[1]*rhh[1]);
				}else{
					layer.compute_Rfresnel_circular(theta_U_out, epsilon_incLayer, rco, rcross);
					atten_amp_co = atten_amp_co*sqrt(rco[0]*rco[0] + rco[1]*rco[1]);
					atten_amp_cross = atten_amp_cross*sqrt(rcross[0]*rcross[0] + rcross[1]*rcross[1]);
				}
				theta_U_in = theta_U_out;
			}
			//Back to surface
			for(j=i-1; j>=0; j--){
				if(valid_path){
					//Ray propagation through layer j
					theta_normal = theta_U_in + alpha_plane[j + 1];
					valid_path = traverse_single_layer(theta_normal, (depth_layer[j + 1] - depth_layer[j]), alpha_plane[j], alpha_plane[j + 1], epsilon_r[j], epsilon_i[j], lambda, false, horizontal_pos, delta_range, amp_atten);
					accum_delay = accum_delay + delta_range;
				}
				if(valid_path){
					//Layer j to layer j-1 intersection
					theta_D_in = theta_U_in + alpha_plane[j + 1] - alpha_plane[j];
					if(j==0){
						epsilon_incLayer[0] = 1.0;
						epsilon_incLayer[1] = 0.0;
					}else{
						epsilon_incLayer[0] = epsilon_r[j - 1];
						epsilon_incLayer[1] = epsilon_i[j - 1];
					}
					layer.epsilon_real = epsilon_r[j];
					layer.epsilon_imag = epsilon_i[j];
					theta_U_in = theta_layer2_Snell(theta_D_in, epsilon_r[j], epsilon_i[j], epsilon_incLayer[0], epsilon_incLayer[1]);
					if(lin_pol){
						layer.compute_Tfresnel_linear(theta_U_in, epsilon_incLayer, tvv, thh);
						atten_amp_v = atten_amp_v*sqrt(tvv[0]*tvv[0] + tvv[1]*tvv[1])*amp_atten;
						atten_amp_h = atten_amp_h*sqrt(thh[0]*thh[0] + thh[1]*thh[1])*amp_atten;
					}else{
						layer.compute_Tfresnel_circular(theta_U_in, epsilon_incLayer, tco, tcross);
						atten_amp_co = atten_amp_co*sqrt(tco[0]*tco[0] + tco[1]*tco[1])*amp_atten;
						atten_amp_cross = atten_amp_cross*sqrt(tco[0]*tco[0] + tco[1]*tco[1])*amp_atten; //co-polar transmission of the cross-polar component
					}
				}
			}
			//Surface level
			if(valid_path){
				if(alpha_plane[0] > 0.0){
					alpha = alpha_plane[0]*PI_NUM/180.0;
					theta = theta_U_in*PI_NUM/180.0;
					delta_range = horizontal_pos*tan(alpha)/cos(alpha + theta);
					x_surf = horizontal_pos*(1.0 + tan(alpha)*tan(alpha + theta));
					z_rec = height_rec;
				}else{
					alpha = -alpha_plane[0]*PI_NUM/180.0;
					delta_range = 0.0;
					x_surf = horizontal_pos;
					z_rec = height_rec - horizontal_pos*tan(alpha);
				}
				delta_range_direct = sqrt(x_surf*x_surf + z_rec*z_rec)*cos((elev_in*PI_NUM/180.0) + atan2(z_rec, x_surf));
				//================== Assign output values
				range_out[i] = accum_delay + delta_range - delta_range_direct;
				if(lin_pol){
					amp_pol1[i] = atten_amp_v;
					amp_pol2[i] = atten_amp_h;
				}else{
					amp_pol1[i] = atten_amp_co;
					amp_pol2[i] = atten_amp_cross;
				}
			}else{
				num_crossed_layers ++;
			}
		}
	}
	return num_crossed_layers;
}

bool traverse_single_layer( double theta_normal, double orig_thickness, double alpha_layer, double alpha_beneath_layer, double epsilon_r, double epsilon_i, double lambda, bool up2down, double &horizontal_pos, double &delta_range, double &amp_atten )
{
	double h, n_real, theta, alpha, hypo_inc, alpha_atten, delta_H_pos, delta_V_pos;
	h = orig_thickness + horizontal_pos*(tan(alpha_beneath_layer*PI_NUM/180.0) - tan(alpha_layer*PI_NUM/180.0));
	n_real = sqrt((sqrt(epsilon_r*epsilon_r + epsilon_i*epsilon_i) + epsilon_r)/2.0);
	if(h <= 0.0){ //crossed layers
		return false;
	}
	theta = theta_normal*PI_NUM/180.0;
	if(up2down){
		alpha = alpha_beneath_layer*PI_NUM/180.0;
	}else{
		alpha = alpha_layer*PI_NUM/180.0;
	}
	if((up2down)&&(alpha > 0.0)||(!up2down)&&(alpha < 0.0)){ //Increase delta range wrt planar layers
		alpha = fabs(alpha);
		hypo_inc = h*tan(theta)*tan(alpha)*(cos(theta) + sin(theta)*tan(theta + alpha));
		delta_H_pos = h*tan(theta) + sin(theta)*hypo_inc;
		delta_V_pos = h + cos(theta)*hypo_inc;
	}else{ //Decrease delta range wrt planar layers
		alpha = fabs(alpha);
		hypo_inc = h*((1.0/cos(theta)) - cos(theta) - (sin(theta)/tan(alpha - theta + PI_NUM/2.0)));
		delta_H_pos = h*tan(theta) - sin(theta)*hypo_inc;
		delta_V_pos = h - cos(theta)*hypo_inc;
	}
	horizontal_pos = horizontal_pos + delta_H_pos;
	delta_range = sqrt(delta_H_pos*delta_H_pos + delta_V_pos*delta_V_pos)*n_real;
	alpha_atten = 2.0*PI_NUM*sqrt((sqrt(epsilon_r*epsilon_r + epsilon_i*epsilon_i) - epsilon_r)/2.0)/lambda;
	amp_atten = exp(-alpha_atten*delta_range/n_real);
	return true;
}
