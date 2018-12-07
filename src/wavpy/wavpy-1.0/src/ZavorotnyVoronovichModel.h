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

class ZaVoModel_GNSSR
{
	bool dump_isolines_data;
	std::string isolines_data_namefile;
	int size_ddm_stored[2];
	double **ddm;
	int len_cov_stored;
	double **cov;
	double **chol;
  public:
	char polarization;
	int exponent_wav_model_length;
	int num_angles;
	int wav_length;
	int ddm_half_dopplers;
	double sampling_rate;
	double delta_doppler;
	double delta_freq;
	double coherent_integration;
	Specular_geometry geometry;
	Reflecting_surface surface;
	RF_FrontEnd receiver_Up;
	RF_FrontEnd receiver_Down;
	GNSS_composite gnss_signal;
	Waveform_power waveform_POW;
	ZaVoModel_GNSSR( void );
	void enable_isolines_data_dump( const char* namefile );
	void disable_isolines_data_dump( void );
	void compute_waveform( int interferometric, int apply_curvature, int add_coherent_pow, int compute_wav_cov, int compute_ddm_cov );
	void get_DDM_doppler_slice( int doppler_index, double* dm_out, int dm_out_length );
	void get_cov_slice( int cov_index, double* cov_out, int cov_out_length );
	void get_noisy_waveform( double* wav_out, int wav_out_length, unsigned long int seed_in );
	void get_noisy_DDM( double* ddm_row_out, int ddm_row_out_length, unsigned long int seed_in );
};

double Compute_power_trans( GNSS_composite signal, double elev, double atm_loss, double posTnorm );
double Doppler_func( double k_scat[3], double velR[3], double k_inc[3], double velT[3], double frequency );
double Sigma0_func( Reflecting_surface surf, double n_scat[3], double n_inc[3], char pol, double azimuthT );
double PDFunction( Reflecting_surface surf, double q[3], double azimuthT );
void AntGain_EH_angles( double ant_k[3], double ant_E[3], double ant_H[3], double n_target_to_ant[3], double out_angles[2] );
void AntGain_PhiTheta_angles( double ant_k[3], double ant_E[3], double ant_H[3], double n_target_to_ant[3], double out_angles[2] );
void Compute_Covariance_DDM( double power_direct_surf, double power_direct_rec, double sampling_rate, double coh_time, double BW, double pow_nd, double pow_nr, double delta_doppler, double** power_surf, int doppler_len, int delay_len, double* lambda_func, int lambda_len, double* wav_out, double** ddm_out, double** cov_out, double** chol_out, int wav_len, int interferometric, int ddm_freq_factor, bool compute_ddm );
void Get_noise_string( double* out_string, double* mean_string, int len_string, double** chol_mat, int seed_rn );