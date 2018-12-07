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

class Waveform_complex_cluster
{
	double **Icomponents;
	double **Qcomponents;
	bool *valid_wavs;
	double *phasorI;
	double *phasorQ;
	bool *valid_phasor;
	int wav_length;
	int cluster_length;
  public:
	int num_valid_wavs;
	int num_phasor_iter;
	Waveform_complex_cluster( void );
	void initialize( int in_cluster_length, int wav_in_length );
	void add_waveform( double* Icomp_in, int len_I, double* Qcomp_in, int len_Q, int cluster_pos );
	void add_waveform_scale( double* Icomp_in, int len_I, double* Qcomp_in, int len_Q, int cluster_pos, double scale_factor );
	void add_waveform_GOLD( signed char* Icomponents_in, int wav_in_length_I, signed char* Qcomponents_in, int wav_in_length_Q, int cluster_pos );
	void add_waveform_PIR( short* XiYi, int wav_in_length_II, short* XqYq, int wav_in_length_QQ, short* XiYq, int wav_in_length_IQ, short* XqYi, int wav_in_length_QI, int cluster_pos );
	double load_ITF_waveforms_SPIR( const char* namefile, double peak_delay_estimate, double BF_phases_UP[8], double BF_phases_DW[8], int filter_num );
	void integrate_waveforms( int coherent_int, double* wav_out, int wav_out_length );
	void integrate_waveforms_remdir( int coherent_int, int coherent_int_dir, double* wav_out, int wav_out_length );
	void integrate_waveforms_retracking( int coherent_int, double sampling_rate, double* retracking_meters, int retracking_meters_length, double* wav_out, int wav_out_length );
	void dump_phase( int lag_pos );
	void dump_phase_peak( void );
	void store_phasor_wavs( int lag_pos );
	void get_phasor( double* phasorI_out, int phasor_length_I, double* phasorQ_out, int phasor_length_Q, signed char* valid_phasor_out, int phasor_length_bool );
	double get_sigma_phase_phasor( int min_valid_samples );
	double get_sigma_phase_phasor_interv( int init_sample, int interv_samples, int min_valid_samples );
	void counterrot_phasor( double* phases_rad, int phases_length, signed char* valid_phases, int phases_length_bool );
	void counterrot_waveforms( double* phases_rad, int phases_length, signed char* valid_phases, int phases_length_bool );
	void correct_navigation_bit( int lag_pos, int store_navbit_phasorI );
	double compute_coherence_time( int lag_peak, int store_acf_phasorQ );
	void compute_singlefreq_DDM( int coherent_int, double doppler_freq, double* freq_ddm, int delay_samples );
	void compute_singlefreq_DDM_remdir( int coherent_int, int coherent_int_dir, double doppler_freq, double* freq_ddm, int delay_samples );
	void compute_singlelag_DDM( int coherent_int, int lag_pos, double delta_freq, double* lag_ddm, int freq_samples );
	void compute_singlelag_DDM_remdir( int coherent_int, int coherent_int_dir, int lag_pos, double delta_freq, double* lag_ddm, int freq_samples );
	double compute_DopplerMap_BW( int coherent_int, int lag_pos, int freq_samples, double delta_freq, double pos_pow_Max[2] );
	double compute_DopplerMap_BW_remdir( int coherent_int, int coherent_int_dir, int lag_pos, int freq_samples, double delta_freq, double pos_pow_Max[2] );
	void compute_whole_DDM( int coherent_int, double ddm_delta_freq, double** ddm, int ddm_delay_samples, int ddm_freq_samples );
	void compute_whole_DDM_remdir( int coherent_int, int coherent_int_dir, double ddm_delta_freq, double** ddm, int ddm_delay_samples, int ddm_freq_samples );
	void compute_whole_DDM_retracking( int coherent_int, double ddm_delta_freq, double** ddm, int ddm_delay_samples, int ddm_freq_samples, double sampling_rate, double* retracking_meters, int retracking_meters_length );
	void compute_whole_LagHologram( double** powerLagHolo, int lagHolo_lags, int lagHolo_freqs );
	void compute_LagHologram( int lag, double* powerLagHolo_singleLag, int fft_samples );
	int get_wav_length( void );
	int get_cluster_length( void );
};

void integrate_wavs( int coh_int, int coh_int_dir, int cl_length, int w_length, double freq, bool* valid, double** Icomp, double** Qcomp, double* wav, bool retracking, double* retrack_samples );
double integrate_wav_lag( int coh_int, int coh_int_dir, int cl_length, int lag, double freq, bool* valid, double** Icomp, double** Qcomp );
double compute_sigma_phase_phasor( int init, int interval, int min_samples, double* phasI, double* phasQ, bool* valid );
bool check_SPIR_blocks( int* spirdata );
void compute_RF_filter( float* filter );
void compute_RF_filter_GalileoE1( float* filter );
void compute_RF_filter_GPS_L5( float* filter );
void compute_ITFwav_SPIR( int* spirdata, double interm_freq, float* filter, double* wav_real, double* wav_imag, int init_sample, int wav_size, double* phases_UP, double* phases_DW, int msec );
double getSPIRdatabit( int word, int pos );
