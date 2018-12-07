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
#include "gsl/gsl_errno.h"
#include "gsl/gsl_fft_real.h"
#include "gsl/gsl_fft_halfcomplex.h"
#include <gsl/gsl_fit.h>

void compute_delays_wav( double* waveform_in, int wav_length, double sampling_rate, double relFactor, double &positionMax, double &posSampleMax, double &sigma_posMax, double &power_posMax, double &positionDer, double &posSampleDer, double &sigma_posDer, double &power_posDer, double &powerDer_posDer, double &positionRel, double &posSampleRel, double &sigma_posRel, double &noiseFloor, double &slope_NormTail, double &sigma_slope_NormTail, double normTail_length, double min_resolution_fft_interp, double fit_length, bool apply_speckle_weights, int num_incoh, bool apply_limits, double limits_center, double limits_width, double apriori_scattdel )
{
//============================== VARIABLES
	double *tau_ext, *waveform, *waveform_ext, *spectrum_ext, *waveform_first_deriv, *waveform_second_deriv, *x_lfit, *y_lfit, *w_lfit;
	double c0, c1, cov00, cov01, cov11, chisq, dt_ext, dt, coef_left, coef_right, twopifreq, sigma_y;
	int i, imin, iposmax, iposder, iposrel, fit_samples, left_lags, right_lags, interpfactor, wav_ext_length, wav_tail_length;
//===========================================================
// Notes on computation of sigma_y when applying speckle weigths:
// Waveform's sigma_y[i] is given by y[i]/sqrt(num incoherent samples).
// At first derivative, sigma_y is obtained from the approximation of y[i]' = (y[i+1] - y[i-1])/(x[i+1] - x[i-1]), then 
// sigma_y[i]' = sqrt(sigma_y[i+1]^2 + sigma_y[i-1]^2 - 2·r_y[i+1]y[i-1]·sigma_y[i+1]·sigma_y[i-1])/(x[i+1] - x[i-1]). In addition,
// from Cardellach et al 2013 (PIRA campaign), we model the speckle noise as the correlation of GNSS-code Autocorr^2 with
// white noise, then r_y[i+1]y[i-1] = Rnn[2i] = 1, then sigma_y[i]' = |(sigma_y[i+1] - sigma_y[i-1])|/(x[i+1] - x[i-1]).
// The same type of approach is followed when computing sigma_y[i]''
//============================= Add "tail" and connect both waveform's ends
	wav_tail_length = 2;
	while(wav_tail_length <= wav_length){
		wav_tail_length = wav_tail_length*2;
	}
	waveform = (double*) malloc (wav_tail_length*sizeof(double));
	for(i=0; i<wav_length; i++){
		waveform[i] = waveform_in[i];
	}
	left_lags = 6;
	right_lags = 6;
	noiseFloor = 0.0;
	for(i=0; i<left_lags; i++){
		noiseFloor = noiseFloor + waveform[i];
	}
	noiseFloor = noiseFloor/double(left_lags);
	for(i=0; i<left_lags; i++){
		coef_left = double(left_lags - i)/double(left_lags);
		coef_right = double(i)/double(left_lags);
		waveform[i] = noiseFloor*coef_left + waveform[left_lags]*coef_right;
	}
	x_lfit = (double*) malloc (right_lags*sizeof(double));
	y_lfit = (double*) malloc (right_lags*sizeof(double));
	for(i=0; i<right_lags; i++){
		x_lfit[i] = i + wav_length - right_lags;
		y_lfit[i] = waveform[int(x_lfit[i])];
	}
	gsl_fit_linear(x_lfit, 1, y_lfit, 1, right_lags, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
	if((c1<0.0)&&((c0 + double(wav_tail_length - 1)*c1)<=noiseFloor)){
		for(i=wav_length; i<wav_tail_length; i++){
			waveform[i] = c0 + double(i)*c1;
			if(waveform[i] < noiseFloor){
				waveform[i] = noiseFloor;
			}
		}
	}else{
		for(i=wav_length; i<wav_tail_length; i++){
			coef_left = double(wav_tail_length - 1 - i)/double(wav_tail_length - wav_length);
			coef_right = double(i - wav_length + 1)/double(wav_tail_length - wav_length);
			waveform[i] = waveform[wav_length - 1]*coef_left + noiseFloor*coef_right;
		}
	}
	free(x_lfit);
	free(y_lfit);
//===========================================================

//============================== PREPARE ARRAYS
	dt = C_LIGHT/sampling_rate;
	interpfactor = 1;
	dt_ext = dt/double(interpfactor);
	while(dt_ext > min_resolution_fft_interp){
		interpfactor = interpfactor*2;
		dt_ext = dt/double(interpfactor);
	}
	wav_ext_length = interpfactor*wav_tail_length;
	tau_ext = (double*) malloc (wav_ext_length*sizeof(double));
	waveform_ext = (double*) malloc (wav_ext_length*sizeof(double));
	spectrum_ext = (double*) malloc (wav_ext_length*sizeof(double));
	waveform_first_deriv = (double*) malloc (wav_ext_length*sizeof(double));
	waveform_second_deriv = (double*) malloc (wav_ext_length*sizeof(double));
	for(i=0; i<wav_ext_length; i++){
		tau_ext[i] = double(i)*dt_ext;
	}
	fit_samples = int(fit_length/dt_ext);
	if(fit_samples > wav_ext_length/2){
		fit_samples = wav_ext_length/2;
		printf("WARNING! fit_samples set to wav_ext_length/2\n");
	}
	x_lfit = (double*) malloc (fit_samples*sizeof(double));
	y_lfit = (double*) malloc (fit_samples*sizeof(double));
	if(apply_speckle_weights){
		w_lfit = (double*) malloc (fit_samples*sizeof(double));
	}
//===========================================================

//============================== Interpolation: zero padding in spectrum
	gsl_fft_real_radix2_transform(waveform, 1, wav_tail_length);
	memset(spectrum_ext, 0, sizeof(double)*wav_ext_length);
	memset(waveform_ext, 0, sizeof(double)*wav_ext_length);
	for(i=0; i<=wav_tail_length/2; i++){
		spectrum_ext[i] = waveform[i]*double(interpfactor);
		waveform_ext[i] = spectrum_ext[i];
	}
	for(i=1; i<wav_tail_length/2; i++){
		spectrum_ext[wav_ext_length - i] = waveform[wav_tail_length - i]*double(interpfactor);
		waveform_ext[wav_ext_length - i] = spectrum_ext[wav_ext_length - i];
	}
	gsl_fft_halfcomplex_radix2_inverse(waveform_ext, 1, wav_ext_length);
//===========================================================

//============================== First derivative: spectrum*(j2pi*freq)
	waveform_first_deriv[0] = 0.0;
	waveform_first_deriv[wav_ext_length/2] = 0.0;
	for(i=1; i<wav_ext_length/2; i++){
		twopifreq = (2.0*PI_NUM*double(i)/double(wav_ext_length/2))/dt_ext;
		waveform_first_deriv[i] = spectrum_ext[wav_ext_length - i]*(-1.0)*twopifreq;
		waveform_first_deriv[wav_ext_length - i] = spectrum_ext[i]*twopifreq;
	}
	gsl_fft_halfcomplex_radix2_inverse(waveform_first_deriv, 1, wav_ext_length);
//===========================================================

//============================== Second Derivative: spectrum*(j2pi*freq)^2
	waveform_second_deriv[0] = 0.0;
	waveform_second_deriv[wav_ext_length/2] = spectrum_ext[wav_ext_length/2]*4.0*PI_NUM*PI_NUM;
	for(i=1; i<wav_ext_length/2; i++){
		twopifreq = (2.0*PI_NUM*double(i)/double(wav_ext_length/2))/dt_ext;
		waveform_second_deriv[i] = spectrum_ext[i]*twopifreq*twopifreq;
		waveform_second_deriv[wav_ext_length - i] = spectrum_ext[wav_ext_length - i]*twopifreq*twopifreq;
	}
	gsl_fft_halfcomplex_radix2_inverse(waveform_second_deriv, 1, wav_ext_length);
//============================================================

//============================== Compute MaxWav delay from first derivative
	power_posMax = 0.0;
	iposmax = 0;
	if(apply_limits){
		powerDer_posDer = 0.0;
		iposder = 0;
		for(i=int(round((limits_center - limits_width/2)/dt_ext)); i<=int(round((limits_center + limits_width/2)/dt_ext)); i++){
			if((i >= 0)&&(i < wav_ext_length)){
				if(waveform_first_deriv[i] > powerDer_posDer){
					powerDer_posDer = waveform_first_deriv[i];
					iposder = i;
				}
			}
		}
		for(i=int(round((dt_ext*double(iposder) + apriori_scattdel - limits_width/4)/dt_ext)); i<=int(round((dt_ext*double(iposder) + apriori_scattdel + limits_width/4)/dt_ext)); i++){
			if((i >= 0)&&(i < wav_ext_length)){
				if(waveform_ext[i] > power_posMax){
					power_posMax = waveform_ext[i];
					iposmax = i;
				}
			}
		}
	}else{
		for(i=0; i<wav_ext_length; i++){
			if(waveform_ext[i] > power_posMax){
				power_posMax = waveform_ext[i];
				iposmax = i;
			}
		}
	}
	imin = iposmax - fit_samples/2;
	if(imin < 1){
		imin = 1;
	}
	if(imin > (wav_ext_length - fit_samples - 3)){
		imin = wav_ext_length - fit_samples - 3;
	}
	for(i=0; i<fit_samples; i++){
		x_lfit[i] = tau_ext[imin + i];
		y_lfit[i] = waveform_first_deriv[imin + i];
		if(apply_speckle_weights){
			sigma_y = fabs(waveform_ext[imin + i + 1] - waveform_ext[imin + i - 1])/(sqrt(double(num_incoh))*2.0*dt_ext);
			w_lfit[i] = 1.0/(sigma_y*sigma_y);
		}
	}
	if(apply_speckle_weights){
		gsl_fit_wlinear(x_lfit, 1, w_lfit, 1, y_lfit, 1, fit_samples, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
	}else{
		gsl_fit_linear(x_lfit, 1, y_lfit, 1, fit_samples, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
	}
	posSampleMax = tau_ext[iposmax];
	positionMax = -c0/c1;
	sigma_posMax = sqrt(cov00/(c1*c1) + cov11*c0*c0/(c1*c1*c1*c1) - 2.0*cov01*c0/(c1*c1*c1));
//============================================================

//============================== Compute MaxDer delay from second derivative
	if(!apply_limits){
		imin = 0;
		while(waveform_ext[imin] < (noiseFloor + 0.3*(power_posMax - noiseFloor))){
			imin = imin + 1;
		}
		powerDer_posDer = 0.0;
		iposder = 0;
		for(i=imin; i<iposmax; i++){
			if(waveform_first_deriv[i] > powerDer_posDer){
				powerDer_posDer = waveform_first_deriv[i];
				iposder = i;
			}
		}
	}
	power_posDer = waveform_ext[iposder];
	imin = iposder - fit_samples/2;
	if(imin < 2){
		imin = 2;
	}
	if(imin > (wav_ext_length - fit_samples - 4)){
		imin = wav_ext_length - fit_samples - 4;
	}
	for(i=0; i<fit_samples; i++){
		x_lfit[i] = tau_ext[imin + i];
		y_lfit[i] = waveform_second_deriv[imin + i];
		if(apply_speckle_weights){
			sigma_y = fabs(fabs(waveform_ext[imin + i + 2] - waveform_ext[imin + i]) - fabs(waveform_ext[imin + i] - waveform_ext[imin + i -2]))/(sqrt(double(num_incoh))*4.0*dt_ext*dt_ext);
			w_lfit[i] = 1.0/(sigma_y*sigma_y);
		}
	}
	if(apply_speckle_weights){
		gsl_fit_wlinear(x_lfit, 1, w_lfit, 1, y_lfit, 1, fit_samples, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
	}else{
		gsl_fit_linear(x_lfit, 1, y_lfit, 1, fit_samples, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
	}
	posSampleDer = tau_ext[iposder];
	positionDer = -c0/c1;
	sigma_posDer = sqrt(cov00/(c1*c1) + cov11*c0*c0/(c1*c1*c1*c1) - 2.0*cov01*c0/(c1*c1*c1));
//============================================================

//============================== Compute RelPos delay from waveform
	sigma_posRel = -1.0;
	positionRel = -1.0;
	posSampleRel = -1.0;
	if((relFactor*power_posMax > noiseFloor)&&(relFactor < 1.0)&&(relFactor > 0.0)){
		//iposrel = 0;
		//while(waveform_ext[iposrel] < relFactor*power_posMax){
		//	iposrel = iposrel + 1;
		//}
		iposrel = iposmax;
		while((waveform_ext[iposrel] > relFactor*power_posMax)&&(iposrel > 0)){
			iposrel = iposrel - 1;
		}
		imin = iposrel - fit_samples/2;
		if(imin < 0){
			imin = 0;
		}
		if(imin > (wav_ext_length - fit_samples - 1)){
			imin = wav_ext_length - fit_samples - 1;
		}
		for(i=0; i<fit_samples; i++){
			x_lfit[i] = tau_ext[imin + i];
			y_lfit[i] = waveform_ext[imin + i] - relFactor*power_posMax;
			if(apply_speckle_weights){
				sigma_y = waveform_ext[imin + i]/sqrt(double(num_incoh));
				w_lfit[i] = 1.0/(sigma_y*sigma_y);
			}
		}
		if(apply_speckle_weights){
			gsl_fit_wlinear(x_lfit, 1, w_lfit, 1, y_lfit, 1, fit_samples, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
		}else{
			gsl_fit_linear(x_lfit, 1, y_lfit, 1, fit_samples, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
		}
		positionRel = -c0/c1;
		posSampleRel = tau_ext[iposrel];
		sigma_posRel = sqrt(cov00/(c1*c1) + cov11*c0*c0/(c1*c1*c1*c1) - 2.0*cov01*c0/(c1*c1*c1));
	}else{
		positionRel = 0.0;
		posSampleRel = 0;
		sigma_posRel = 999999.999999;
	}
//===========================================================

//==============================
	free(x_lfit);
	free(y_lfit);
	if(apply_speckle_weights){
		free(w_lfit);
	}
	fit_samples = int(normTail_length/dt_ext);
	x_lfit = (double*) malloc (fit_samples*sizeof(double));
	y_lfit = (double*) malloc (fit_samples*sizeof(double));
	if(apply_speckle_weights){
		w_lfit = (double*) malloc (fit_samples*sizeof(double));
	}
	if(iposmax + fit_samples < wav_ext_length){
		for(i=0; i<fit_samples; i++){
			x_lfit[i] = tau_ext[iposmax + i];
			y_lfit[i] = waveform_ext[iposmax + i]/power_posMax;
			if(apply_speckle_weights){
				sigma_y = waveform_ext[iposmax + i]/(sqrt(double(num_incoh))*power_posMax);
				w_lfit[i] = 1.0/(sigma_y*sigma_y);
			}
		}
		if(apply_speckle_weights){
			gsl_fit_wlinear(x_lfit, 1, w_lfit, 1, y_lfit, 1, fit_samples, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
		}else{
			gsl_fit_linear(x_lfit, 1, y_lfit, 1, fit_samples, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
		}
		slope_NormTail = c1;
		sigma_slope_NormTail = cov11;
	}else{
		slope_NormTail = 0.0;
		sigma_slope_NormTail = 999999.999999;
	}
//============================================================

	free(waveform);
	free(waveform_ext);
	free(spectrum_ext);
	free(tau_ext);
	free(waveform_first_deriv);
	free(waveform_second_deriv);
	free(x_lfit);
	free(y_lfit);
	if(apply_speckle_weights){
		free(w_lfit);
	}

	return;
}