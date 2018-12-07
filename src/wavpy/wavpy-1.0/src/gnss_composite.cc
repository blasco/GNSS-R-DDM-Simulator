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
#include "gnss_composite.h"
#include "ancillary_functions.h"

GNSS_composite::GNSS_composite( void ){
	frequency = FREQ_GPS_L1;
	sampling_rate = 80000000.0;
	filter_BB_BW = 12000000.0;
	weight_CA = 1.0;
	weight_PY = 1.0;
	weight_M = 1.0;
	weight_IM = 1.0;
	weight_E1A = 0.0;
	weight_E1B = 0.0;
	weight_E1C = 0.0;
	weight_B1I = 0.0;
	weight_L1C = 0.0;
	lambda_size = 2*int(round(sampling_rate/GPS_CA_CHIP_RATE)) + 1;
	lambda_func = (double*) malloc (lambda_size*sizeof(double));
	Compute_GPS_L1_composite(sampling_rate, filter_BB_BW, weight_CA, weight_PY, weight_M, weight_IM, weight_L1C, lambda_func);
	lambda_updated = true;
	return;
}

void GNSS_composite::dump_parameters( void ){
	printf("===================== GPS L1 CODES ====================\n");
	printf("- C/A   [Tx. dBW][weight] : %f %f\n", POW_CA_Trans_dBW, weight_CA);
	printf("- P(Y)  [Tx. dBW][weight] : %f %f\n", POW_PY_Trans_dBW, weight_PY);
	printf("- M     [Tx. dBW][weight] : %f %f\n", POW_M_Trans_dBW, weight_M);
	printf("- IM    [Tx. dBW][weight] : %f %f\n", POW_IM_Trans_dBW, weight_IM);
	printf("- L1C   [Tx. dBW][weight] : %f %f\n", POW_L1C_Trans_dBW, weight_L1C);
	printf("=======================================================\n");
	printf("=================== GALILEO E1 CODES ===================\n");
	printf("- E1-A  [Rx. dBW][weight] : %f %f\n", POW_E1A_RecEarth_dBW, weight_E1A);
	printf("- E1-B  [Rx. dBW][weight] : %f %f\n", POW_E1B_RecEarth_dBW, weight_E1B);
	printf("- E1-C  [Rx. dBW][weight] : %f %f\n", POW_E1C_RecEarth_dBW, weight_E1C);
	printf("=======================================================\n");
	printf("=================== BEIDOU B1 CODES ===================\n");
	printf("- B1I   [Rx. dBW][weight] : %f %f\n", POW_B1I_RecEarth_dBW, weight_B1I);
	printf("=======================================================\n");
	printf("====================== QZSS CODES =====================\n");
	printf("- C/A   [Rx. dBW][weight] : %f %f\n", POW_CA_RecEarth_dBW, weight_CA);
	printf("- L1C   [Rx. dBW][weight] : %f %f\n", POW_L1C_RecEarth_dBW, weight_L1C);
	printf("=======================================================\n");
	printf("================= INSTRUMENTAL PARAMS =================\n");
	printf("- Sampling rate      [Hz] : %f\n", sampling_rate);
	printf("- Filter B-B BW      [Hz] : %f\n", filter_BB_BW);
	printf("=======================================================\n");
	return;
}

void GNSS_composite::set_instrumental_params( double input_sampling_rate, double input_filter_BB_BW, char computeLambda ){
	if(sampling_rate != input_sampling_rate){
		sampling_rate = input_sampling_rate;
		lambda_size = 2*int(round(sampling_rate/GPS_CA_CHIP_RATE)) + 1;
		free (lambda_func);
		lambda_func = (double*) malloc (lambda_size*sizeof(double));
		lambda_updated = false;
	}
	if(filter_BB_BW != input_filter_BB_BW){
		filter_BB_BW = input_filter_BB_BW;
		lambda_updated = false;
	}
	if(computeLambda == 1){
		compute_lambda_func();
	}
	return;
}

void GNSS_composite::compute_lambda_func( void ){
	if((weight_CA > 0.0)||(weight_PY > 0.0)||(weight_M > 0.0)||(weight_IM > 0.0)||(weight_L1C > 0.0)){ //GPS or QZSS
		Compute_GPS_L1_composite(sampling_rate, filter_BB_BW, weight_CA, weight_PY, weight_M, weight_IM, weight_L1C, lambda_func);
		frequency = FREQ_GPS_L1;
		weight_E1A = 0.0;
		weight_E1B = 0.0;
		weight_E1C = 0.0;
		weight_B1I = 0.0;
	}else{
		if((weight_E1A > 0.0)||(weight_E1B > 0.0)||(weight_E1C > 0.0)){ //Galileo
			Compute_Galileo_E1_composite(sampling_rate, filter_BB_BW, weight_E1A, weight_E1B, weight_E1C, lambda_func);
			frequency = FREQ_GPS_L1;
			weight_CA = 0.0;
			weight_PY = 0.0;
			weight_M = 0.0;
			weight_IM = 0.0;
			weight_B1I = 0.0;
			weight_L1C = 0.0;
		}else{
			if(weight_B1I > 0.0){ //BeiDou
				Compute_BeiDou_B1(sampling_rate, filter_BB_BW, weight_B1I, lambda_func);
				frequency = FREQ_BEIDOU_B1;
				weight_CA = 0.0;
				weight_PY = 0.0;
				weight_M = 0.0;
				weight_IM = 0.0;
				weight_E1A = 0.0;
				weight_E1B = 0.0;
				weight_E1C = 0.0;
				weight_L1C = 0.0;
			}else{
				printf("ERROR! All code weights are <= 0.0! Lambda is not computed.\n");
			}
		}
	}
	lambda_updated = true;
	return;
}

void GNSS_composite::get_lambda_func( double* range_vec, int range_length, double* lambda_out, int lambda_length ){
	int i;
	if((lambda_length != lambda_size)||(range_length != lambda_size)){
		printf("ERROR! lambda_length and range_length must be equal to lambda_size");
		return;
	}
	if(!lambda_updated){
		compute_lambda_func();
	}
	for(i=0; i<lambda_size; i++){
		lambda_out[i] = lambda_func[i];
		range_vec[i] = double(i - (lambda_size/2))*C_LIGHT/sampling_rate;
	}
	return;
}
