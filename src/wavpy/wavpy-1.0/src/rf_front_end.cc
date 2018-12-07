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
#include "rf_front_end.h"
#include "ancillary_functions.h"

RF_FrontEnd::RF_FrontEnd( void ){
	int i, j;
	//Set default values: omnidirectional GPS L1 antenna looking down
	frequency = FREQ_GPS_L1;
	for(i=0; i<181; i++){
		for(j=0; j<360; j++){
			antenna_pattern_dB[i][j] = 0.0;
		}
	}
	isotropic = true;
	antenna_vector_BF_E[0] = 1.0;
	antenna_vector_BF_E[1] = 0.0;
	antenna_vector_BF_E[2] = 0.0;
	antenna_vector_BF_H[0] = 0.0;
	antenna_vector_BF_H[1] = 1.0;
	antenna_vector_BF_H[2] = 0.0;
	vector3Prod(antenna_vector_BF_E, antenna_vector_BF_H, antenna_vector_BF_k);
	antenna_Gain_dB = 3.0;
	antenna_Aeff = (C_LIGHT/frequency)*(C_LIGHT/frequency)*pow(10.0, (antenna_Gain_dB/10.0))/(4.0*PI_NUM);
	antenna_T = 200.0;
	noise_F_dB = 3.0;
	noise_T = antenna_T + 290.0*(pow(10.0, (noise_F_dB/10.0)) - 1.0);
	filter_BB_BW = 12000000.0;
	noise_pow_dBW = 10.0*log10(K_BOLTZMANN*noise_T*filter_BB_BW);
	array_num_elements = 1;
	array_factor_ready = false;
	return;
}

void RF_FrontEnd::dump_parameters( void ){
	int i;
	printf("======================= ANTENNA =======================\n");
	printf("- E-plane (Body-Frame)   : %f %f %f\n", antenna_vector_BF_E[0], antenna_vector_BF_E[1], antenna_vector_BF_E[2]);
	printf("- H-plane (Body-Frame)   : %f %f %f\n", antenna_vector_BF_H[0], antenna_vector_BF_H[1], antenna_vector_BF_H[2]);
	printf("- 'Pointing' vector (BF) : %f %f %f\n", antenna_vector_BF_k[0], antenna_vector_BF_k[1], antenna_vector_BF_k[2]);
	printf("- Gain              [dB] : %f\n", antenna_Gain_dB);
	printf("- Effective area   [m^2] : %f\n", antenna_Aeff);
	printf("- Antenna Temp       [K] : %f\n", antenna_T);
	if(isotropic){
		printf("- Isotropic antenna      : YES\n");
	}else{
		printf("- Isotropic antenna      : NO\n");
	}
	printf("=======================================================\n");
	printf("==================== RF PARAMETERS ====================\n");
	printf("- Frequency         [Hz] : %f\n", frequency);
	printf("- B-B Bandwidth     [Hz] : %f\n", filter_BB_BW);
	printf("- Noise Temperature  [K] : %f\n", noise_T);
	printf("- Noise Figure      [dB] : %f\n", noise_F_dB);
	printf("- Noise Power      [dBW] : %f\n", noise_pow_dBW);
	printf("=======================================================\n");
	printf("======================= ARRAY =========================\n");
	printf("- # Antenna elements     : %d\n", array_num_elements);
	if(array_num_elements > 1){
		printf("-------------------------------------------------------\n");
		for(i=0; i<array_num_elements; i++){
			printf("- %d) Pos [m]: [%f, %f] | Phase delay [rad]: %f\n", i, element_pos_AF[i][0], element_pos_AF[i][1], phase_delay[i]);
		}
	}
	printf("=======================================================\n");
	return;
}

void RF_FrontEnd::set_antenna_orientation_BF_EH( double antenna_vector_BF_E_in[3], double antenna_vector_BF_H_in[3] ){
	int i;
	double norm_vec_E, norm_vec_H;
	norm_vec_E = norm3vec(antenna_vector_BF_E_in);
	norm_vec_H = norm3vec(antenna_vector_BF_H_in);
	for(i=0; i<3; i++){
		antenna_vector_BF_E[i] = antenna_vector_BF_E_in[i]/norm_vec_E;
		antenna_vector_BF_H[i] = antenna_vector_BF_H_in[i]/norm_vec_H;
	}
	vector3Prod(antenna_vector_BF_E, antenna_vector_BF_H, antenna_vector_BF_k);
	return;
}

void RF_FrontEnd::set_antenna_orientation_BF_k( double antenna_vector_BF_k_in[3] ){
	int i;
	double norm_vec_k, norm_vec_E;
	printf("WARNING! With this option, antenna vectors from planes E and H are arbitrary selected (but still perpendicular to vector k). This option is only valid if the antenna pattern is equal in both planes. If not, set_antenna_orientation_BF_EH must be called.\n");
	norm_vec_k = norm3vec(antenna_vector_BF_k_in);
	for(i=0; i<3; i++){
		antenna_vector_BF_k[i] = antenna_vector_BF_k_in[i]/norm_vec_k;
	}
	if(antenna_vector_BF_k[0] == 0.0){ //We construct an ARBITRARY antenna_vector_BF_E perpendicular to antenna_vector_BF_k
		antenna_vector_BF_E[0] = 1.0;
		antenna_vector_BF_E[1] = 0.0;
		antenna_vector_BF_E[2] = 0.0;
	}else{
		norm_vec_E = sqrt(((antenna_vector_BF_k[1]/antenna_vector_BF_k[0])*(antenna_vector_BF_k[1]/antenna_vector_BF_k[0])) + 1.0);
		antenna_vector_BF_E[0] = (-antenna_vector_BF_k[1]/antenna_vector_BF_k[0])/norm_vec_E;
		antenna_vector_BF_E[1] = 1.0/norm_vec_E;
		antenna_vector_BF_E[2] = 0.0;
	}
	vector3Prod(antenna_vector_BF_k, antenna_vector_BF_E, antenna_vector_BF_H);
	return;
}

void RF_FrontEnd::get_antenna_orientation_BF( double antenna_vector_BF_E_out[3], double antenna_vector_BF_H_out[3], double antenna_vector_BF_k_out[3] ){
	int i;
	for(i=0; i<3; i++){
		antenna_vector_BF_E_out[i] = antenna_vector_BF_E[i];
		antenna_vector_BF_H_out[i] = antenna_vector_BF_H[i];
		antenna_vector_BF_k_out[i] = antenna_vector_BF_k[i];
	}
	return;
}

void RF_FrontEnd::set_antenna_whole_pattern( double antenna_whole_pattern_dB_in[181][360] ){
	int i, j;
	for(i=0; i<181; i++){
		for(j=0; j<360; j++){
			antenna_pattern_dB[i][j] = antenna_whole_pattern_dB_in[i][j];
		}
	}
	isotropic = false;
	return;
}

void RF_FrontEnd::set_val_antenna_pattern( int phi_AntFrame, int theta_AntFrame, double pattern_dB_value ){
	if((theta_AntFrame > 180)||(theta_AntFrame < 0)){
		printf("ERROR! theta_AntFrame must lay within [0,180].\n");
		return;
	}
	if((phi_AntFrame > 360)||(phi_AntFrame < -180)){
		printf("ERROR! phi_AntFrame must lay within [-180,180] or [0,360] deg.\n");
		return;
	}
	if(phi_AntFrame < 0){
		phi_AntFrame = phi_AntFrame + 360;
	}
	if(phi_AntFrame == 360){
		phi_AntFrame = 0;
	}
	antenna_pattern_dB[theta_AntFrame][phi_AntFrame] = pattern_dB_value;
	return;
}

void RF_FrontEnd::set_antenna_pattern_FF( double antenna_full_pattern_dB_in[360] ){
	int i;
	double antenna_pattern_E_dB[360];
	double antenna_pattern_H_dB[360];
	double directivity_dB = -999.9;
	for(i=0; i<360; i++){
		if(antenna_full_pattern_dB_in[i] > directivity_dB){
			directivity_dB = antenna_full_pattern_dB_in[i];
		}
	}
	for(i=0; i<360; i++){
		antenna_pattern_E_dB[i] = antenna_full_pattern_dB_in[i] - directivity_dB;
		antenna_pattern_H_dB[i] = antenna_full_pattern_dB_in[i] - directivity_dB;
	}
	antenna_Gain_dB = antenna_Gain_dB + directivity_dB;
	Set_whole_pattern_from_EH_planes_RevLinInterp(antenna_pattern_E_dB, antenna_pattern_H_dB, antenna_pattern_dB);
	isotropic = false;
	return;
}

void RF_FrontEnd::set_antenna_pattern_FH( double antenna_half_pattern_dB_in[181] ){
	int i;
	double antenna_pattern_E_dB[360];
	double antenna_pattern_H_dB[360];
	double directivity_dB = -999.9;
	for(i=0; i<181; i++){
		if(antenna_half_pattern_dB_in[i] > directivity_dB){
			directivity_dB = antenna_half_pattern_dB_in[i];
		}
	}
	for(i=0; i<181; i++){
		antenna_pattern_E_dB[i] = antenna_half_pattern_dB_in[i] - directivity_dB;
		antenna_pattern_H_dB[i] = antenna_half_pattern_dB_in[i] - directivity_dB;
	}
	for(i=181; i<360; i++){
		antenna_pattern_E_dB[i] = antenna_pattern_E_dB[360 - i] - directivity_dB;
		antenna_pattern_H_dB[i] = antenna_pattern_H_dB[360 - i] - directivity_dB;
	}
	antenna_Gain_dB = antenna_Gain_dB + directivity_dB;
	Set_whole_pattern_from_EH_planes_RevLinInterp(antenna_pattern_E_dB, antenna_pattern_H_dB, antenna_pattern_dB);
	isotropic = false;
	return;
}

void RF_FrontEnd::set_antenna_pattern_interp( double* antenna_angles_deg, int angles_length, double* antenna_pattern_dB_in, int pattern_length, double min_level_dB ){
	int i;
	double antenna_pattern_E_dB[360];
	double antenna_pattern_H_dB[360];
	double directivity_dB = -999.9;
	if(angles_length != pattern_length){
		printf("ERROR! antenna_angles_deg and antenna_pattern_dB_in must have equal size.\n");
		return;
	}
	interpolate_antenna_pattern(antenna_angles_deg, antenna_pattern_dB_in, pattern_length, min_level_dB, antenna_pattern_E_dB);
	for(i=0; i<360; i++){
		if(antenna_pattern_E_dB[i] > directivity_dB){
			directivity_dB = antenna_pattern_E_dB[i];
		}
	}
	for(i=0; i<360; i++){
		antenna_pattern_E_dB[i] = antenna_pattern_E_dB[i] - directivity_dB;
		antenna_pattern_H_dB[i] = antenna_pattern_E_dB[i];
	}
	antenna_Gain_dB = antenna_Gain_dB + directivity_dB;
	Set_whole_pattern_from_EH_planes_RevLinInterp(antenna_pattern_E_dB, antenna_pattern_H_dB, antenna_pattern_dB);
	isotropic = false;
	return;
}

void RF_FrontEnd::set_antenna_patterns_FF( double antenna_full_pattern_E_dB_in[360], double antenna_full_pattern_H_dB_in[360] ){
	int i;
	double antenna_pattern_E_dB[360];
	double antenna_pattern_H_dB[360];
	double directivity_dB = -999.9;
	for(i=0; i<360; i++){
		if(antenna_full_pattern_E_dB_in[i] > directivity_dB){
			directivity_dB = antenna_full_pattern_E_dB_in[i];
		}
		if(antenna_full_pattern_H_dB_in[i] > directivity_dB){
			directivity_dB = antenna_full_pattern_H_dB_in[i];
		}
	}
	for(i=0; i<360; i++){
		antenna_pattern_E_dB[i] = antenna_full_pattern_E_dB_in[i] - directivity_dB;
		antenna_pattern_H_dB[i] = antenna_full_pattern_H_dB_in[i] - directivity_dB;
	}
	antenna_Gain_dB = antenna_Gain_dB + directivity_dB;
	Set_whole_pattern_from_EH_planes_RevLinInterp(antenna_pattern_E_dB, antenna_pattern_H_dB, antenna_pattern_dB);
	isotropic = false;
	return;
}

void RF_FrontEnd::set_antenna_patterns_FH( double antenna_half_pattern_E_dB_in[181], double antenna_half_pattern_H_dB_in[181] ){
	int i;
	double antenna_pattern_E_dB[360];
	double antenna_pattern_H_dB[360];
	double directivity_dB = -999.9;
	for(i=0; i<181; i++){
		if(antenna_half_pattern_E_dB_in[i] > directivity_dB){
			directivity_dB = antenna_half_pattern_E_dB_in[i];
		}
		if(antenna_half_pattern_H_dB_in[i] > directivity_dB){
			directivity_dB = antenna_half_pattern_H_dB_in[i];
		}
	}
	for(i=0; i<181; i++){
		antenna_pattern_E_dB[i] = antenna_half_pattern_E_dB_in[i] - directivity_dB;
		antenna_pattern_H_dB[i] = antenna_half_pattern_H_dB_in[i] - directivity_dB;
	}
	for(i=181; i<360; i++){
		antenna_pattern_E_dB[i] = antenna_pattern_E_dB[360 - i] - directivity_dB;
		antenna_pattern_H_dB[i] = antenna_pattern_H_dB[360 - i] - directivity_dB;
	}
	antenna_Gain_dB = antenna_Gain_dB + directivity_dB;
	Set_whole_pattern_from_EH_planes_RevLinInterp(antenna_pattern_E_dB, antenna_pattern_H_dB, antenna_pattern_dB);
	isotropic = false;
	return;
}

void RF_FrontEnd::set_antenna_patterns_interp( double* antenna_angles_E_deg, int angles_E_length, double* antenna_pattern_E_dB_in, int pattern_E_length, double* antenna_angles_H_deg, int angles_H_length, double* antenna_pattern_H_dB_in, int pattern_H_length, double min_level_dB ){
	int i;
	double antenna_pattern_E_dB[360];
	double antenna_pattern_H_dB[360];
	double directivity_dB = -999.9;
	if(angles_E_length != pattern_E_length){
		printf("ERROR! antenna_angles_E_deg and antenna_pattern_E_dB_in must have equal size.\n");
		return;
	}
	if(angles_H_length != pattern_H_length){
		printf("ERROR! antenna_angles_H_deg and antenna_pattern_H_dB_in must have equal size.\n");
		return;
	}
	interpolate_antenna_pattern(antenna_angles_E_deg, antenna_pattern_E_dB_in, pattern_E_length, min_level_dB, antenna_pattern_E_dB);
	interpolate_antenna_pattern(antenna_angles_H_deg, antenna_pattern_H_dB_in, pattern_H_length, min_level_dB, antenna_pattern_H_dB);
	for(i=0; i<360; i++){
		if(antenna_pattern_E_dB[i] > directivity_dB){
			directivity_dB = antenna_pattern_E_dB[i];
		}
		if(antenna_pattern_H_dB[i] > directivity_dB){
			directivity_dB = antenna_pattern_H_dB[i];
		}
	}
	for(i=0; i<360; i++){
		antenna_pattern_E_dB[i] = antenna_pattern_E_dB[i] - directivity_dB;
		antenna_pattern_H_dB[i] = antenna_pattern_H_dB[i] - directivity_dB;
	}
	antenna_Gain_dB = antenna_Gain_dB + directivity_dB;
	Set_whole_pattern_from_EH_planes_RevLinInterp(antenna_pattern_E_dB, antenna_pattern_H_dB, antenna_pattern_dB);
	isotropic = false;
	return;
}

void RF_FrontEnd::get_antenna_whole_pattern( double antenna_pattern_dB_out[181][360] ){
	int i, j;
	for(i=0; i<181; i++){
		for(j=0; j<360; j++){
			antenna_pattern_dB_out[i][j] = antenna_pattern_dB[i][j] + antenna_Gain_dB;
		}
	}
	return;
}

void RF_FrontEnd::get_antenna_patterns( double antenna_pattern_E_dB_out[360], double antenna_pattern_H_dB_out[360] ){
	int i;
	for(i=0; i<181; i++){
		antenna_pattern_E_dB_out[i] = antenna_pattern_dB[i][0] + antenna_Gain_dB;
		antenna_pattern_H_dB_out[i] = antenna_pattern_dB[i][90] + antenna_Gain_dB;
	}
	for(i=181; i<360; i++){
		antenna_pattern_E_dB_out[i] = antenna_pattern_dB[360 - i][180] + antenna_Gain_dB;
		antenna_pattern_H_dB_out[i] = antenna_pattern_dB[360 - i][270] + antenna_Gain_dB;
	}
	return;
}

void RF_FrontEnd::set_receiver_params( double antenna_Gain_dB_in, double antenna_T_in, double noise_F_dB_in, double filter_BB_BW_in, signed char isotropic_antenna ){
	if(filter_BB_BW_in <= 0.0){
		printf("ERROR! filter_BB_BW must be > 0.\n");
		return;
	}
	antenna_Gain_dB = antenna_Gain_dB_in;
	antenna_Aeff = (C_LIGHT/frequency)*(C_LIGHT/frequency)*pow(10.0, (antenna_Gain_dB/10.0))/(4.0*PI_NUM);
	antenna_T = antenna_T_in;
	noise_F_dB = noise_F_dB_in;
	noise_T = antenna_T + 290.0*(pow(10.0, (noise_F_dB/10.0)) - 1.0);
	filter_BB_BW = filter_BB_BW_in;
	noise_pow_dBW = 10.0*log10(K_BOLTZMANN*noise_T*filter_BB_BW);
	if(isotropic_antenna == 1){
		isotropic = true;
	}else{
		isotropic = false;
	}
	return;
}

void RF_FrontEnd::set_antenna_eff_area( double antenna_Aeff_in ){
	if(antenna_Aeff_in <= 0.0){
		printf("ERROR! antenna_Aeff must be > 0.\n");
		return;
	}
	antenna_Aeff = antenna_Aeff_in;
	antenna_Gain_dB = 10.0*log10(antenna_Aeff*4.0*PI_NUM/((C_LIGHT/frequency)*(C_LIGHT/frequency)));
	return;
}

void RF_FrontEnd::set_noise_T( double noise_T_in ){
	if(noise_T_in <= 0.0){
		printf("ERROR! noise_T must be > 0.\n");
		return;
	}
	noise_T = noise_T_in;
	antenna_T = noise_T - 290.0*(pow(10.0, (noise_F_dB/10.0)) - 1.0); //We keep previous noise_F_dB
	noise_pow_dBW = 10.0*log10(K_BOLTZMANN*noise_T*filter_BB_BW);
	return;
}

void RF_FrontEnd::set_noise_pow_dBW( double noise_pow_dBW_in ){
	noise_pow_dBW = noise_pow_dBW_in;
	noise_T = pow(10.0, (noise_pow_dBW/10.0))/(K_BOLTZMANN*filter_BB_BW);
	antenna_T = noise_T - 290.0*(pow(10.0, (noise_F_dB/10.0)) - 1.0); //We keep previous noise_F_dB
	return;
}

double RF_FrontEnd::get_anglesEH_gain_dB( double angle_E_plane_wrt_k, double angle_H_plane_wrt_k ){
	double total_gain, gain_pattern, phi_AntFrame, theta_AntFrame, array_factor_gain;
	if((angle_E_plane_wrt_k > 360.0)||(angle_E_plane_wrt_k < -180.0)){
		printf("ERROR! angle_E_plane_wrt_k must lay within [-180,180] or [0,360] deg.\n");
		return -999.999;
	}
	if((angle_H_plane_wrt_k > 360.0)||(angle_H_plane_wrt_k < -180.0)){
		printf("ERROR! angle_H_plane_wrt_k must lay within [-180,180] or [0,360] deg.\n");
		return -999.999;
	}
	//////////////////////WARNING: this function does not work properly
	Compute_PhiTheta_from_anglesEH(angle_E_plane_wrt_k, angle_H_plane_wrt_k, &phi_AntFrame, &theta_AntFrame);
	//////////////////////
	gain_pattern = Get_gain_pattern(theta_AntFrame, phi_AntFrame, antenna_pattern_dB);
	if(array_factor_ready){
		array_factor_gain = Get_gain_pattern(theta_AntFrame, phi_AntFrame, array_factor_dB);
	}else{
		array_factor_gain = 0.0;
	}
	total_gain = antenna_Gain_dB + gain_pattern + array_factor_gain;
	return total_gain;
}

double RF_FrontEnd::get_PhiTheta_gain_dB( double phi_AntFrame, double theta_AntFrame ){
	double total_gain, gain_pattern, array_factor_gain;
	if((theta_AntFrame > 180.0)||(theta_AntFrame < 0.0)){
		printf("ERROR! theta_AntFrame must lay within [0,180] deg.\n");
		return -999.999;
	}
	if((phi_AntFrame > 360.0)||(phi_AntFrame < -180.0)){
		printf("ERROR! phi_AntFrame must lay within [-180,180] or [0,360] deg.\n");
		return -999.999;
	}
	if(isotropic){
		return antenna_Gain_dB;
	}
	//2D angular interpolation
	if(phi_AntFrame < 0.0){
		phi_AntFrame = phi_AntFrame + 360.0;
	}
	if(phi_AntFrame == 360.0){
		phi_AntFrame = 0.0;
	}
	gain_pattern = Get_gain_pattern(theta_AntFrame, phi_AntFrame, antenna_pattern_dB);
	if(array_factor_ready){
		array_factor_gain = Get_gain_pattern(theta_AntFrame, phi_AntFrame, array_factor_dB);
	}else{
		array_factor_gain = 0.0;
	}
	//Result
	total_gain = antenna_Gain_dB + gain_pattern + array_factor_gain;
	return total_gain;
}

double RF_FrontEnd::get_incvector_gain_dB( double incvector[3] ){
	double norm_incvec, total_gain, gain_pattern, phi_AntFrame, theta_AntFrame, array_factor_gain;
	double point_vec[3];
	norm_incvec = norm3vec(incvector);
	if(norm_incvec == 0.0){
		printf("ERROR! Incoming vector with 0.0 norm.\n");
		return -999.999;
	}
	if(isotropic){
		return antenna_Gain_dB;
	}
	//Reverse sense of incoming vector
	point_vec[0] = -incvector[0]/norm_incvec;
	point_vec[1] = -incvector[1]/norm_incvec;
	point_vec[2] = -incvector[2]/norm_incvec;
	//Computation of phi and theta angles
	theta_AntFrame = acos(scalar3Prod(antenna_vector_BF_k, point_vec))*180.0/PI_NUM;
	phi_AntFrame = atan2(scalar3Prod(antenna_vector_BF_H, point_vec), scalar3Prod(antenna_vector_BF_E, point_vec))*180.0/PI_NUM;
	if(phi_AntFrame < 0.0){
		phi_AntFrame = phi_AntFrame + 360.0;
	}
	if(phi_AntFrame == 360.0){
		phi_AntFrame = 0.0;
	}
	//Get antenna gain
	gain_pattern = Get_gain_pattern(theta_AntFrame, phi_AntFrame, antenna_pattern_dB);
	//Get array factor gain
	if(array_factor_ready){
		array_factor_gain = Get_gain_pattern(theta_AntFrame, phi_AntFrame, array_factor_dB);
	}else{
		array_factor_gain = 0.0;
	}
	total_gain = antenna_Gain_dB + gain_pattern + array_factor_gain + array_factor_gain;
	return total_gain;
}

double RF_FrontEnd::get_frequency( void ){
	return frequency;
}

double RF_FrontEnd::get_antenna_Gain_dB( void ){
	return antenna_Gain_dB;
}

double RF_FrontEnd::get_antenna_Aeff( void ){
	return antenna_Aeff;
}

double RF_FrontEnd::get_antenna_T( void ){
	return antenna_T;
}

double RF_FrontEnd::get_noise_T( void ){
	return noise_T;
}

double RF_FrontEnd::get_noise_pow_dBW( void ){
	return noise_pow_dBW;
}

double RF_FrontEnd::get_noise_F_dB( void ){
	return noise_F_dB;
}

double RF_FrontEnd::get_filter_BB_BW( void ){
	return filter_BB_BW;
}

void RF_FrontEnd::set_antenna_elements_pos_AF( double* element_pos_in, int num_elem_in, int plane_dim, char lambda_units ){
	int i;
	double length_factor;
	if((plane_dim != 2)||(num_elem_in < 2)){
		printf("Dimensions not valid: minimum number of elements is 2 (planar array) and 2-D coordinates\n");
		return;
	}
	if(lambda_units == 1){
		length_factor = C_LIGHT/frequency; //lambda
	}else{
		length_factor = 1.0;
	}
	if(num_elem_in != array_num_elements){
		if(array_num_elements > 1){
			free(phase_delay);
			for(i=0;i<array_num_elements;i++)
				free(element_pos_AF[i]);
			free(element_pos_AF);
		}
		array_num_elements = num_elem_in;
		phase_delay = (double*) malloc(array_num_elements*sizeof(double));
		element_pos_AF = (double **) malloc(array_num_elements*sizeof(double *));
		for(i=0;i<array_num_elements;i++){
			element_pos_AF[i] = (double *) malloc(2*sizeof(double));
			phase_delay[i] = 0.0;
		}
	}
	for(i=0;i<array_num_elements;i++){ 
		element_pos_AF[i][0] = element_pos_in[2*i]*length_factor;
		element_pos_AF[i][1] = element_pos_in[2*i + 1]*length_factor;
	}
	array_factor_ready = false;
	return;
}

void RF_FrontEnd::set_phase_delays( double* phase_delay_in, int num_elems_in ){
	int i;
	if(array_num_elements < 2){
		printf("ERROR! Empty array.\n");
		return;
	}
	if(num_elems_in != array_num_elements){
		printf("ERROR! Number of elements not valid: %d (it has to be %d)\n", num_elems_in, array_num_elements);
		return;
	}
	for(i=0; i<array_num_elements; i++){
		phase_delay[i] = phase_delay_in[i];
	}
	array_factor_ready = false;
	return;
}

void RF_FrontEnd::get_phase_delays( double* phase_delay_out, int num_elems_out ){
	int i;
	if(array_num_elements < 2){
		printf("ERROR! Empty array.\n");
		return;
	}
	if(num_elems_out != array_num_elements){
		printf("ERROR! Number of elements not valid: %d (it has to be %d)\n", num_elems_out, array_num_elements);
		return;
	}
	for(i=0; i<array_num_elements; i++){
		phase_delay_out[i] = phase_delay[i];
	}
	return;
}

void RF_FrontEnd::compute_array_factor( void ){
	int i, j, k;
	double theta_rad, phi_rad, kx, ky, accumAF_real, accumAF_imag, arg;
	if(array_num_elements < 2){
		printf("ERROR! Empty array.\n");
		return;
	}
	for(i=0; i<181; i++){
		for(j=0; j<360; j++){
			theta_rad = double(i)*PI_NUM/180.0;
			phi_rad = double(j)*PI_NUM/180.0;
			kx = 2.0*PI_NUM*sin(theta_rad)*cos(phi_rad)*frequency/C_LIGHT;
			ky = 2.0*PI_NUM*sin(theta_rad)*sin(phi_rad)*frequency/C_LIGHT;
			accumAF_real = 0.0;
			accumAF_imag = 0.0;
			for(k=0; k<array_num_elements; k++){
				arg = phase_delay[k] - (kx*element_pos_AF[k][0] + ky*element_pos_AF[k][1]);
				accumAF_real = accumAF_real + cos(arg);
				accumAF_imag = accumAF_imag + sin(arg);
			}
			array_factor_dB[i][j] = 10.0*log10(sqrt(accumAF_real*accumAF_real + accumAF_imag*accumAF_imag));
		}
	}
	array_factor_ready = true;
	return;
}

void RF_FrontEnd::get_array_factor( double array_factor_dB_out[181][360] ){
	int i, j;
	if(!array_factor_ready){
		printf("ERROR! Array factor is not computed for current configuration.\n");
		return;
	}
	for(i=0; i<181; i++){
		for(j=0; j<360; j++){
			array_factor_dB_out[i][j] = array_factor_dB[i][j];
		}
	}
	return;
}

void RF_FrontEnd::compute_phase_delays_UPA( double theta_max, double phi_max ){
	int i, x_dim, y_dim;
	double theta_rad, phi_rad, k, dist_x, dist_y, alpha_x, alpha_y;
	if(!Check_if_UPA_distribution(x_dim, y_dim, dist_x, dist_y, array_num_elements, element_pos_AF)){
		printf("ERROR! Current distribution is not UPA (Uniformly-distributed Planar Array) with min (x,y) at first element.\n");
		return;
	}
	k = 2.0*PI_NUM*frequency/C_LIGHT;
	theta_rad = theta_max*PI_NUM/180.0;
	phi_rad = range180(phi_max)*PI_NUM/180.0;
	if(y_dim == 1){
		alpha_y = 0.0;
		alpha_x = -k*dist_x*sin(theta_rad)*cos(phi_rad);
	}
	if(x_dim == 1){
		alpha_x = 0.0;
		alpha_y = -k*dist_y*sin(theta_rad)*sin(phi_rad);
	}
	if((y_dim > 1)&&(x_dim > 1)){
		alpha_x = k*dist_x*sin(theta_rad)/sqrt(1.0 + tan(phi_rad)*tan(phi_rad));
		if(fabs(phi_rad) > (PI_NUM/2.0)){ //This requirement was found after comparison with vectorial results
			alpha_x = -alpha_x;
		}
		alpha_y = dist_y*alpha_x*tan(phi_rad)/dist_x;
	}
	phase_delay[0] = 0.0;
	for(i=1; i<array_num_elements; i++){
		phase_delay[i] = alpha_x*((element_pos_AF[i][0] - element_pos_AF[0][0])/dist_x) + alpha_y*((element_pos_AF[i][1] - element_pos_AF[0][1])/dist_y);
	}
	array_factor_ready = false;
	return;
}

void RF_FrontEnd::compute_phase_delays_pos_ECEF_RT( double inertials[3], double posR_km[3], double posT_km[3] ){
	int i;
	double ref_range, range_element;
	double pos_vec_BF[3], pos_vec_ECEF[3];
	for(i=0; i<array_num_elements; i++){
		pos_vec_BF[0] = antenna_vector_BF_E[0]*element_pos_AF[i][0] + antenna_vector_BF_H[0]*element_pos_AF[i][1];
		pos_vec_BF[1] = antenna_vector_BF_E[1]*element_pos_AF[i][0] + antenna_vector_BF_H[1]*element_pos_AF[i][1];
		pos_vec_BF[2] = antenna_vector_BF_E[2]*element_pos_AF[i][0] + antenna_vector_BF_H[2]*element_pos_AF[i][1];
		BF2ECEF(inertials[0], inertials[1], inertials[2], posR_km, pos_vec_BF, pos_vec_ECEF);
		range_element = sqrt(pow((posR_km[0]*1000.0 + pos_vec_ECEF[0] - posT_km[0]*1000.0),2.0) + pow((posR_km[1]*1000.0 + pos_vec_ECEF[1] - posT_km[1]*1000.0),2.0) + pow((posR_km[2]*1000.0 + pos_vec_ECEF[2] - posT_km[2]*1000.0),2.0));
		if(i==0){
			ref_range = range_element;
		}
		phase_delay[i] = (ref_range - range_element)*2.0*PI_NUM*frequency/C_LIGHT; //Phase that has to be applied in order to be aligned with first array element
	}
	array_factor_ready = false;
	return;
}



void Set_whole_pattern_from_EH_planes( double ant_pattern_E_dB[360], double ant_pattern_H_dB[360], double ant_pattern_dB[181][360] )
{
	int phi, theta;
	int index_E_dw, index_E_up, index_H_dw, index_H_up;
	double phi_rad, theta_rad, angle_E_plane_wrt_k, angle_H_plane_wrt_k, coeff_E_up, coeff_E_dw, coeff_H_up, coeff_H_dw;
	for(theta=0; theta<181; theta++){
		for(phi=0; phi<360; phi++){
			theta_rad = double(theta)*PI_NUM/180.0;
			phi_rad = double(phi)*PI_NUM/180.0;
			angle_E_plane_wrt_k = atan2(sin(theta_rad)*cos(phi_rad), cos(theta_rad))*180.0/PI_NUM;
			if(angle_E_plane_wrt_k < 0.0){
				angle_E_plane_wrt_k = angle_E_plane_wrt_k + 360.0;
			}
			angle_H_plane_wrt_k = atan2(sin(theta_rad)*sin(phi_rad), cos(theta_rad))*180.0/PI_NUM;
			if(angle_H_plane_wrt_k < 0.0){
				angle_H_plane_wrt_k = angle_H_plane_wrt_k + 360.0;
			}
			index_E_dw = int(floor(angle_E_plane_wrt_k));
			coeff_E_dw = 1.0 - (angle_E_plane_wrt_k - double(index_E_dw)); 
			index_E_up = (index_E_dw + 1)%360;
			coeff_E_up = 1.0 - coeff_E_dw;
			index_H_dw = int(floor(angle_H_plane_wrt_k));
			coeff_H_dw = 1.0 - (angle_H_plane_wrt_k - double(index_H_dw)); 
			index_H_up = (index_H_dw + 1)%360;
			coeff_H_up = 1.0 - coeff_H_dw;
			ant_pattern_dB[theta][phi] = coeff_E_dw*ant_pattern_E_dB[index_E_dw] + coeff_E_up*ant_pattern_E_dB[index_E_up] + coeff_H_dw*ant_pattern_H_dB[index_H_dw] + coeff_H_up*ant_pattern_H_dB[index_H_up];
		}
	}
	return;
}

void Set_whole_pattern_from_EH_planes_RevLinInterp( double ant_pattern_E_dB[360], double ant_pattern_H_dB[360], double ant_pattern_dB[181][360] )
{
	int phi, theta;
	int index_E, index_H;
	double coeff_E, coeff_H;
	for(phi=0; phi<90; phi++){ //First quadrant in EH plane
		coeff_E = double(90 - phi)/90.0;
		coeff_H = double(phi)/90.0;
		for(theta=0; theta<181; theta++){
			index_E = theta;
			index_H = theta;
			ant_pattern_dB[theta][phi] = coeff_E*ant_pattern_E_dB[index_E] + coeff_H*ant_pattern_H_dB[index_H];
		}
	}
	for(phi=90; phi<180; phi++){ //Second quadrant in EH plane
		coeff_E = double(phi - 90)/90.0;
		coeff_H = double(180 - phi)/90.0;
		for(theta=0; theta<181; theta++){
			index_E = (360 - theta)%360;
			index_H = theta;
			ant_pattern_dB[theta][phi] = coeff_E*ant_pattern_E_dB[index_E] + coeff_H*ant_pattern_H_dB[index_H];
		}
	}
	for(phi=180; phi<270; phi++){ //Third quadrant in EH plane
		coeff_E = double(270 - phi)/90.0;
		coeff_H = double(phi - 180)/90.0;
		for(theta=0; theta<181; theta++){
			index_E = (360 - theta)%360;
			index_H = (360 - theta)%360;
			ant_pattern_dB[theta][phi] = coeff_E*ant_pattern_E_dB[index_E] + coeff_H*ant_pattern_H_dB[index_H];
		}
	}
	for(phi=270; phi<360; phi++){ //Fourth quadrant in EH plane
		coeff_E = double(phi - 270)/90.0;
		coeff_H = double(360 - phi)/90.0;
		for(theta=0; theta<181; theta++){
			index_E = theta;
			index_H = (360 - theta)%360;
			ant_pattern_dB[theta][phi] = coeff_E*ant_pattern_E_dB[index_E] + coeff_H*ant_pattern_H_dB[index_H];
		}
	}
	return;
}

//////////////////////WARNING: this function does not work properly
void Compute_PhiTheta_from_anglesEH( double angleE, double angleH, double* phi, double* theta )
{
	if(angleE == 90.0){
		*theta = 90.0;
		*phi = 0.0;
		return;
	}
	if((angleE == -90.0)||(angleE == 270.0)){
		*theta = 90.0;
		*phi = 180.0;
		return;
	}
	if(angleH == 90.0){
		*theta = 90.0;
		*phi = 90.0;
		return;
	}
	if((angleH == -90.0)||(angleH == 270.0)){
		*theta = 90.0;
		*phi = 270.0;
		return;
	}
	*phi = atan2(tan(angleH*PI_NUM/180.0), tan(angleE*PI_NUM/180.0))*180.0/PI_NUM;
	*theta = atan2(tan(angleE*PI_NUM/180.0), cos((*phi)*PI_NUM/180.0))*180.0/PI_NUM;
	return;
}
//////////////////////

double Get_gain_pattern( double theta, double phi, double ant_pattern_dB[181][360] )
{
	double gain_pattern, coeff_theta_up, coeff_theta_dw, coeff_phi_up, coeff_phi_dw;
	int index_theta_dw, index_theta_up, index_phi_dw, index_phi_up;
	index_theta_dw = int(floor(theta));
	coeff_theta_dw = 1.0 - (theta - double(index_theta_dw));
	index_theta_up = index_theta_dw + 1;
	coeff_theta_up = 1.0 - coeff_theta_dw;
	if(index_theta_up == 181){
		index_theta_up = 180;
		coeff_theta_up = 0.0;
	}
	index_phi_dw = int(floor(phi));
	coeff_phi_dw = 1.0 - (phi - double(index_phi_dw)); 
	index_phi_up = (index_phi_dw + 1)%360;
	coeff_phi_up = 1.0 - coeff_phi_dw;
	gain_pattern = coeff_theta_dw*coeff_phi_dw*ant_pattern_dB[index_theta_dw][index_phi_dw] + coeff_theta_up*coeff_phi_dw*ant_pattern_dB[index_theta_up][index_phi_dw] + coeff_theta_dw*coeff_phi_up*ant_pattern_dB[index_theta_dw][index_phi_up] + coeff_theta_up*coeff_phi_up*ant_pattern_dB[index_theta_up][index_phi_up];
	return gain_pattern;
}

bool Check_if_UPA_distribution( int &dimX, int &dimY, double &distX, double &distY, int num, double** positions )
{
	int i, j;
	bool diff_x, diff_y;
	double min_distX, min_distY, dist_error;
	//Check if: 1) first element is at position with min (x,y)
	//          2) num_elements = x_dim*y_dim
	//          3) uniform dist_x and dist_y (full grid)
	dimX = 1;
	dimY = 1;
	distX = 0.0;
	distY = 0.0;
	dist_error = 0.0000000001; //10^(-10) or 1 Angström
	for(i=1; i<num; i++){
		if(((positions[i][0] - positions[0][0]) < 0.0) || ((positions[i][1] - positions[0][1]) < 0.0)){
			printf("UPA-check: first element is not at position with min (x,y).\n");
			return false;
		}
		if((positions[i][0] - positions[0][0]) > 0.0){
			if(distX == 0.0) distX = (positions[i][0] - positions[0][0]);
			distX = std::min(distX, (positions[i][0] - positions[0][0]));
		}
		if((positions[i][1] - positions[0][1]) > 0.0){
			if(distY == 0.0) distY = (positions[i][1] - positions[0][1]);
			distY = std::min(distY, (positions[i][1] - positions[0][1]));
		}
		diff_x = true;
		diff_y = true;
		for(j=(i-1); j>=0; j--){
			if(positions[i][0] == positions[j][0]){
				diff_x = false;
			}
			if(positions[i][1] == positions[j][1]){
				diff_y = false;
			}
		}
		if(diff_x) dimX++;
		if(diff_y) dimY++;
	}
	if(dimX*dimY != num){
		printf("UPA-check: dimX*dimY must be equal to num_elements (rectangular grid).\n");
		return false;
	}
	for(i=1; i<num; i++){
		min_distX = distX*dimX;
		min_distY = distY*dimY;
		for(j=0; j<num; j++){
			if((fabs(positions[i][0] - positions[j][0]) > dist_error)&&(fabs(positions[i][0] - positions[j][0]) < min_distX)){
				min_distX = fabs(positions[i][0] - positions[j][0]);
			}
			if((fabs(positions[i][1] - positions[j][1]) > dist_error)&&(fabs(positions[i][1] - positions[j][1]) < min_distY)){
				min_distY = fabs(positions[i][1] - positions[j][1]);
			}
		}
		if((fabs(min_distX - distX) > dist_error)||(fabs(min_distY - distY) > dist_error)){
			printf("UPA-check: not uniform distances.\n");
			return false;
		}
	}
	return true;
}