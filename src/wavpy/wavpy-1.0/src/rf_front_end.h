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

class RF_FrontEnd
{
	double antenna_pattern_dB[181][360]; //Antenna pattern as a function of theta and phi in the antenna frame
	double antenna_vector_BF_E[3]; //X from antenna frame in body frame
	double antenna_vector_BF_H[3]; //Y from antenna frame in body frame
	double antenna_vector_BF_k[3]; //Z from antenna frame in body frame
	bool isotropic;
	double frequency;
	double antenna_Gain_dB;
	double antenna_Aeff;
	double antenna_T;
	double noise_T;
	double noise_pow_dBW;
	double noise_F_dB;
	double filter_BB_BW;
	//2D PLANAR ARRAY
	int array_num_elements;
	double **element_pos_AF; //In meters
	double *phase_delay;  //In radians
	double array_factor_dB[181][360];
	bool array_factor_ready;
  public:
	RF_FrontEnd( void );
	void dump_parameters( void );
	void set_antenna_orientation_BF_EH( double antenna_vector_BF_E_in[3], double antenna_vector_BF_H_in[3] );
	void set_antenna_orientation_BF_k( double antenna_vector_BF_k_in[3] );
	void get_antenna_orientation_BF( double antenna_vector_BF_E_out[3], double antenna_vector_BF_H_out[3], double antenna_vector_BF_k_out[3] );
	void set_antenna_whole_pattern( double antenna_whole_pattern_dB_in[181][360] );
	void set_val_antenna_pattern( int phi_AntFrame, int theta_AntFrame, double pattern_dB_value );
	void set_antenna_pattern_FF( double antenna_full_pattern_dB_in[360] );
	void set_antenna_pattern_FH( double antenna_half_pattern_dB_in[181] );
	void set_antenna_pattern_interp( double* antenna_angles_deg, int angles_length, double* antenna_pattern_dB_in, int pattern_length, double min_level_dB );
	void set_antenna_patterns_FF( double antenna_full_pattern_E_dB_in[360], double antenna_full_pattern_H_dB_in[360] );
	void set_antenna_patterns_FH( double antenna_half_pattern_E_dB_in[181], double antenna_half_pattern_H_dB_in[181] );
	void set_antenna_patterns_interp( double* antenna_angles_E_deg, int angles_E_length, double* antenna_pattern_E_dB_in, int pattern_E_length, double* antenna_angles_H_deg, int angles_H_length, double* antenna_pattern_H_dB_in, int pattern_H_length, double min_level_dB );
	void get_antenna_whole_pattern( double antenna_pattern_dB_out[181][360] );
	void get_antenna_patterns( double antenna_pattern_E_dB_out[360], double antenna_pattern_H_dB_out[360] );
	void set_receiver_params( double antenna_Gain_dB_in, double antenna_T_in, double noise_F_dB_in, double filter_BB_BW_in, signed char isotropic_antenna );
	void set_antenna_eff_area( double antenna_Aeff_in );
	void set_noise_T( double noise_T_in );
	void set_noise_pow_dBW( double noise_pow_dBW_in );
	double get_anglesEH_gain_dB( double angle_E_plane, double angle_H_plane );
	double get_PhiTheta_gain_dB( double phi_AntFrame, double theta_AntFrame );
	double get_incvector_gain_dB( double incvector[3] );
	double get_frequency( void );
	double get_antenna_Gain_dB( void );
	double get_antenna_Aeff( void );
	double get_antenna_T( void );
	double get_noise_T( void );
	double get_noise_pow_dBW( void );
	double get_noise_F_dB( void );
	double get_filter_BB_BW( void );
	//2D PLANAR ARRAY
	void set_antenna_elements_pos_AF( double* element_pos_in, int num_elem_in, int plane_dim, char lambda_units );
	void set_phase_delays( double* phase_delay_in, int num_elems_in );
	void get_phase_delays( double* phase_delay_out, int num_elems_out );
	void compute_array_factor( void );
	void get_array_factor( double array_factor_dB_out[181][360] );
	void compute_phase_delays_UPA( double theta_max, double phi_max );
	void compute_phase_delays_pos_ECEF_RT( double inertials[3], double posR_km[3], double posT_km[3] );
};

void Set_whole_pattern_from_EH_planes( double ant_pattern_E_dB[360], double ant_pattern_H_dB[360], double ant_pattern_dB[181][360] );
void Set_whole_pattern_from_EH_planes_RevLinInterp( double ant_pattern_E_dB[360], double ant_pattern_H_dB[360], double ant_pattern_dB[181][360] );
void Compute_PhiTheta_from_anglesEH( double angleE, double angleH, double* phi, double* theta );
double Get_gain_pattern( double theta, double phi, double ant_pattern_dB[181][360] );
bool Check_if_UPA_distribution( int &dimX, int &dimY, double &distX, double &distY, int num, double** positions );
