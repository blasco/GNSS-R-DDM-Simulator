#include "wavpy_global_variables.h"
#include "ancillary_functions.h"
#include "waveform_power.h"
#include "compute_delays_wav.h"

Waveform_power::Waveform_power( void ){
	int i;
	wav_length = 0;
	scale_factor = 1.0;
	positionMax = 0.0;
	posSampleMax = 0;
	sigma_posMax = 999999.999999;
	powerMax = 0.0;
	positionDer = 0.0;
	posSampleDer = 0;
	sigma_posDer = 999999.999999;
	power_posDer = 0.0;
	powerDer_posDer = 0.0;
	positionRel = 0.0;
	posSampleRel = 0;
	sigma_posRel = 999999.999999;
	floorNoise = 0.0;
	slope_normTail = 0.0;
	sigma_slope_normTail = 999999.999999;
	sampling_rate = 20000000.0;
	rel_factor = 0.5;
	init_range = 0.0;
	min_resolution_fft_interp = 0.15;
	fit_length = 10.0;
	normTail_length = 50.0; //50 meters is OK for GPS. For Galileo, 100 meters.
	return;
}

void Waveform_power::set_waveform( double* waveform_in, int wav_in_length ){
	int i;
	if(wav_in_length <= 0){
		printf("ERROR! Waveform size not valid (< 0)\n");
		return;
	}
	if(wav_length != wav_in_length){
		if(wav_length > 0){
			free (waveform);
		}
		waveform = (double*) malloc (wav_in_length*sizeof(double));
	}
	wav_length = wav_in_length;
	scale_factor = 0.0;
	for(i=0; i<wav_length; i++){
		waveform[i] = waveform_in[i];
		if(fabs(waveform[i]) > fabs(scale_factor)){
			scale_factor = waveform[i];
		}
	}
	if(scale_factor == 0.0){
		scale_factor = 1.0;
	}else{
		for(i=0; i<wav_length; i++){
			waveform[i] = waveform[i]/scale_factor;
		}
	}
	return;
}

void Waveform_power::set_float_waveform( float* float_waveform_in, int wav_in_length ){
	int i;
	if(wav_in_length <= 0){
		printf("ERROR! Waveform size not valid (< 0)\n");
		return;
	}
	if(wav_length != wav_in_length){
		if(wav_length > 0){
			free (waveform);
		}
		waveform = (double*) malloc (wav_in_length*sizeof(double));
	}
	wav_length = wav_in_length;
	scale_factor = 0.0;
	for(i=0; i<wav_length; i++){
		waveform[i] = double(float_waveform_in[i]);
		if(fabs(waveform[i]) > fabs(scale_factor)){
			scale_factor = waveform[i];
		}
	}
	if(scale_factor == 0.0){
		scale_factor = 1.0;
	}else{
		for(i=0; i<wav_length; i++){
			waveform[i] = waveform[i]/scale_factor;
		}
	}
	return;
}

void Waveform_power::set_norm_waveform( float* norm_waveform_in, int norm_wav_in_length, double scale_factor_in ){
	int i;
	if(norm_wav_in_length <= 0){
		printf("ERROR! Waveform size not valid (< 0)\n");
		return;
	}
	if(wav_length != norm_wav_in_length){
		if(wav_length > 0){
			free (waveform);
		}
		waveform = (double*) malloc (norm_wav_in_length*sizeof(double));
	}
	wav_length = norm_wav_in_length;
	for(i=0; i<wav_length; i++){
		waveform[i] = double(norm_waveform_in[i]);
	}
	scale_factor = scale_factor_in;
	return;
}

void Waveform_power::get_waveform( double* waveform_out, int wav_out_length ){
	int i;
	if(wav_length != wav_out_length){
		printf("ERROR! Waveform size = %d\n", wav_length);
		return;
	}
	for(i=0; i<wav_length; i++){
		waveform_out[i] = waveform[i]*scale_factor;
	}
	return;
}

void Waveform_power::add_waveform_retracking( double* waveform_in, int wav_in_length, double retrack_range, float wav_weight, bool apply_safety_margin ){
	int i;
	double retrack_samples, max_new_wav;
	double *wavdata;
	if(wav_length != wav_in_length){
		printf("ERROR! Waveform size = %d\n", wav_length);
		return;
	}
	if((wav_weight > 1.0)||(wav_weight <= 0.0)){
		printf("ERROR! Incorrect weight (it must be between 0 and 1)\n");
		return; 
	}
	wavdata = (double*) malloc (wav_length*sizeof(double));
	for(i=0;i<wav_length;i++){
		wavdata[i] = waveform_in[i];
	}
	retrack_samples = retrack_range*sampling_rate/C_LIGHT;
	retrack_real_waveform(wavdata, wav_length, retrack_samples);
	max_new_wav = 0.0;
	for(i=0; i<wav_length; i++){
		if((apply_safety_margin)&&(((retrack_samples > 0.0)&&(i <= int(retrack_samples)))||((retrack_samples < 0.0)&&(i >= (wav_length - 1 + int(trunc(retrack_samples))))))){
			waveform[i] = waveform[i]*scale_factor*(1.0 - wav_weight);
		}else{
			waveform[i] = waveform[i]*scale_factor*(1.0 - wav_weight) + wavdata[i]*wav_weight;
		}
		if(fabs(waveform[i]) > fabs(max_new_wav)){
			max_new_wav = waveform[i];
		}
	}
	scale_factor = max_new_wav;
	free (wavdata);
	for(i=0; i<wav_length; i++){
		waveform[i] = waveform[i]/scale_factor;
	}
	return;
}

void Waveform_power::set_sampling_rate( double samplingRate ){
	if(samplingRate > 0.0){
		sampling_rate = samplingRate;
	}else{
		printf("ERROR! Sampling-rate must be >0.0. Change not applied.\n");
	}
	return;
}

void Waveform_power::set_rel_factor( double relFactor ){
	if((relFactor < 1.0)&&(relFactor > 0.0)){
		rel_factor = relFactor;
	}else{
		printf("ERROR! Rel-Factor must be <1.0 and >0.0. Change not applied.\n");
	}
	return;
}

double Waveform_power::get_rel_factor( void ){
	return rel_factor;
}

void Waveform_power::set_min_resolution_fft_interp( double resolution_in ){
	if(resolution_in > 0.0){
		min_resolution_fft_interp = resolution_in;
	}else{
		printf("ERROR! Min_resolution_fft_interp must be >0.0. Change not applied.\n");
	}
	return;
}

void Waveform_power::set_fit_length( double fit_length_in ){
	if(fit_length_in > 0.0){
		fit_length = fit_length_in;
	}else{
		printf("ERROR! Fit_length must be greater than zero. Change not applied.\n");
	}
	return;
}

void Waveform_power::set_normtail_length( double normtail_length_in ){
	if(normtail_length_in > 0.0){
		normTail_length = normtail_length_in;
	}else{
		printf("ERROR! normTail_length must be >0.0. Change not applied.\n");
	}
	return;
}

void Waveform_power::compute_delays( void ){
	if(wav_length == 0){
		printf("ERROR! Waveform not stored.\n");
		return;
	}
	compute_delays_wav(waveform, wav_length, sampling_rate, rel_factor, positionMax, posSampleMax, sigma_posMax, powerMax, positionDer, posSampleDer, sigma_posDer, power_posDer, powerDer_posDer, positionRel, posSampleRel, sigma_posRel, floorNoise, slope_normTail, sigma_slope_normTail, normTail_length, min_resolution_fft_interp, fit_length, false, 0, false, 0.0, 0.0, 0.0);
	if(scale_factor != 1.0){
		powerMax = powerMax*scale_factor;
		power_posDer = power_posDer*scale_factor;
		powerDer_posDer = powerDer_posDer*scale_factor;
		floorNoise = floorNoise*scale_factor;
	}
	positionMax = positionMax + init_range;
	posSampleMax = posSampleMax + init_range;
	positionDer = positionDer + init_range;
	posSampleDer = posSampleDer + init_range;
	positionRel = positionRel + init_range;
	posSampleRel = posSampleRel + init_range;
	return;
}

void Waveform_power::compute_delays_wspeckle( int num_incoh ){
	if(wav_length == 0){
		printf("ERROR! Waveform not stored.\n");
		return;
	}
	compute_delays_wav(waveform, wav_length, sampling_rate, rel_factor, positionMax, posSampleMax, sigma_posMax, powerMax, positionDer, posSampleDer, sigma_posDer, power_posDer, powerDer_posDer, positionRel, posSampleRel, sigma_posRel, floorNoise, slope_normTail, sigma_slope_normTail, normTail_length, min_resolution_fft_interp, fit_length, true, num_incoh, false, 0.0, 0.0, 0.0);
	if(scale_factor != 1.0){
		powerMax = powerMax*scale_factor;
		power_posDer = power_posDer*scale_factor;
		powerDer_posDer = powerDer_posDer*scale_factor;
		floorNoise = floorNoise*scale_factor;
	}
	positionMax = positionMax + init_range;
	posSampleMax = posSampleMax + init_range;
	positionDer = positionDer + init_range;
	posSampleDer = posSampleDer + init_range;
	positionRel = positionRel + init_range;
	posSampleRel = posSampleRel + init_range;
	return;
}

void Waveform_power::compute_delays_wlimits( double limits_center, double limits_width, double apriori_scattdel ){
	double end_range, limits_init, limits_end;
	if(wav_length == 0){
		printf("ERROR! Waveform not stored.\n");
		return;
	}
	end_range = init_range + double(wav_length)*C_LIGHT/sampling_rate;
	limits_init = limits_center - limits_width/2;
	limits_end = limits_center + limits_width/2;
	if((limits_init < init_range)||(limits_end > end_range)){
		printf("ERROR! Limits [%f, %f] outside waveform's range [%f, %f].\n", limits_init, limits_end, init_range, end_range);
		return;
	}
	compute_delays_wav(waveform, wav_length, sampling_rate, rel_factor, positionMax, posSampleMax, sigma_posMax, powerMax, positionDer, posSampleDer, sigma_posDer, power_posDer, powerDer_posDer, positionRel, posSampleRel, sigma_posRel, floorNoise, slope_normTail, sigma_slope_normTail, normTail_length, min_resolution_fft_interp, fit_length, false, 0, true, (limits_center - init_range), limits_width, apriori_scattdel);
	if(scale_factor != 1.0){
		powerMax = powerMax*scale_factor;
		power_posDer = power_posDer*scale_factor;
		powerDer_posDer = powerDer_posDer*scale_factor;
		floorNoise = floorNoise*scale_factor;
	}
	positionMax = positionMax + init_range;
	posSampleMax = posSampleMax + init_range;
	positionDer = positionDer + init_range;
	posSampleDer = posSampleDer + init_range;
	positionRel = positionRel + init_range;
	posSampleRel = posSampleRel + init_range;
	return;
}

void Waveform_power::set_init_range( double init_range_in ){
	positionMax = positionMax - init_range + init_range_in;
	posSampleMax = posSampleMax - init_range + init_range_in;
	positionDer = positionDer - init_range + init_range_in;
	posSampleDer = posSampleDer - init_range + init_range_in;
	positionRel = positionRel - init_range + init_range_in;
	posSampleRel = posSampleRel - init_range + init_range_in;
	init_range = init_range_in;
	return;
}

void Waveform_power::get_range_waveform( double* range, int size_range ){
	int i;
	if(wav_length != size_range){
		printf("ERROR! Waveform size = %d\n", wav_length);
		return;
	}
	for(i=0; i<wav_length; i++){
		range[i] = init_range + double(i)*C_LIGHT/sampling_rate;
	}
	return;
}

void Waveform_power::dump_norm_waveform( void ){
	int i;
	double dt = C_LIGHT/sampling_rate;
	if(wav_length == 0){
		printf("ERROR! Waveform not stored.\n");
		return;
	}
	for(i=0; i<wav_length; i++){
		printf("WAV %f %f\n", double(i)*dt, waveform[i]);
	}
	return;
}

void Waveform_power::dump_delays( void ){
	if(wav_length == 0){
		printf("ERROR! Waveform not stored.\n");
		return;
	}
	printf("MAX %f %f %f %f\n", positionMax, sigma_posMax, powerMax, (10.0*log10(powerMax)));
	printf("DER %f %f %f %f\n", positionDer, sigma_posDer, power_posDer, powerDer_posDer);
	printf("REL %f %f %f %f\n", positionRel, sigma_posRel, powerMax*rel_factor, (10.0*log10(powerMax*rel_factor)));
	return;
}

int Waveform_power::get_wav_length( void ){
	return wav_length;
}
