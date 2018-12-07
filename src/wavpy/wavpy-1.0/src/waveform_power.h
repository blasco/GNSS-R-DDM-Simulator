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

class Waveform_power
{
	double *waveform;
	int wav_length;
	double sampling_rate, scale_factor, rel_factor, init_range, min_resolution_fft_interp, fit_length, normTail_length;
  public:
	double positionMax, posSampleMax, sigma_posMax, powerMax, positionDer, posSampleDer, sigma_posDer, power_posDer, powerDer_posDer, floorNoise, positionRel, posSampleRel, sigma_posRel, slope_normTail, sigma_slope_normTail;
	Waveform_power( void );
	void set_waveform( double* waveform_in, int wav_in_length );
	void set_float_waveform( float* float_waveform_in, int wav_in_length );
	void set_norm_waveform( float* norm_waveform_in, int norm_wav_in_length, double scale_factor_in );
	void get_waveform( double* waveform_out, int wav_out_length );
	void add_waveform_retracking( double* waveform_in, int wav_in_length, double retrack_delay, float wav_weight, bool apply_safety_margin );
	void set_sampling_rate( double samplingRate );
	void set_rel_factor( double relFactor );
	double get_rel_factor( void );
	void set_min_resolution_fft_interp( double resolution_in );
	void set_fit_length( double fit_length_in );
	void set_normtail_length( double normtail_length_in );
	void compute_delays( void );
	void compute_delays_wspeckle( int num_incoh );
	void compute_delays_wlimits( double limits_center, double limits_width, double apriori_scattdel );
	void set_init_range( double init_range_in );
	void get_range_waveform( double* range, int size_range );
	void dump_norm_waveform( void );
	void dump_delays( void );
	int get_wav_length( void );
};
