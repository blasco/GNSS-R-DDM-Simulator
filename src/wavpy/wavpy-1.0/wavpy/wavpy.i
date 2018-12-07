%module wavpy
%{
#define SWIG_FILE_WITH_INIT
#include "../src/waveform_complex.h"
#include "../src/waveform_power.h"
#include "../src/rf_front_end.h"
#include "../src/gnss_composite.h"
#include "../src/specular_geometry.h"
#include "../src/reflecting_surface.h"
#include "../src/ZavorotnyVoronovichModel.h"
#include "../src/MultiRaySingleRefl_model.h"
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (double* IN_ARRAY1, int DIM1) {(double* Icomp_in, int len_I)}
%apply (double* IN_ARRAY1, int DIM1) {(double* Qcomp_in, int len_Q)}
%apply (signed char* IN_ARRAY1, int DIM1) {(signed char* Icomponents_in, int wav_in_length_I)}
%apply (signed char* IN_ARRAY1, int DIM1) {(signed char* Qcomponents_in, int wav_in_length_Q)}
%apply (short* IN_ARRAY1, int DIM1) {(short* XiYi, int wav_in_length_II)}
%apply (short* IN_ARRAY1, int DIM1) {(short* XqYq, int wav_in_length_QQ)}
%apply (short* IN_ARRAY1, int DIM1) {(short* XiYq, int wav_in_length_IQ)}
%apply (short* IN_ARRAY1, int DIM1) {(short* XqYi, int wav_in_length_QI)}
%apply (double IN_ARRAY1[ANY]) {(double BF_phases_UP[8])}
%apply (double IN_ARRAY1[ANY]) {(double BF_phases_DW[8])}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* wav_out, int wav_out_length)}
%apply (double* IN_ARRAY1, int DIM1) {(double* retracking_meters, int retracking_meters_length)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* phasorI_out, int phasor_length_I)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* phasorQ_out, int phasor_length_Q)}
%apply (signed char* ARGOUT_ARRAY1, int DIM1) {(signed char* valid_phasor_out, int phasor_length_bool)}
%apply (double* IN_ARRAY1, int DIM1) {(double* phases_rad, int phases_length)}
%apply (signed char* IN_ARRAY1, int DIM1) {(signed char* valid_phases, int phases_length_bool)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* freq_ddm, int delay_samples)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* lag_ddm, int freq_samples)}
%apply (double ARGOUT_ARRAY1[ANY]) {(double pos_pow_Max[2])}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* powerLagHolo_singleLag, int fft_samples)}
%include "../src/waveform_complex.h"
%clear (double* Icomp_in, int len_I);
%clear (double* Qcomp_in, int len_Q);
%clear (signed char* Icomponents_in, int wav_in_length_I);
%clear (signed char* Qcomponents_in, int wav_in_length_Q);
%clear (short* XiYi, int wav_in_length_II);
%clear (short* XqYq, int wav_in_length_QQ);
%clear (short* XiYq, int wav_in_length_IQ);
%clear (short* XqYi, int wav_in_length_QI);
%clear (double BF_phases_UP[8]);
%clear (double BF_phases_DW[8]);
%clear (double* wav_out, int wav_out_length);
%clear (double* retracking_meters, int retracking_meters_length);
%clear (double* phasorI_out, int phasor_length_I);
%clear (double* phasorQ_out, int phasor_length_Q);
%clear (signed char* valid_phasor_out, int phasor_length_bool);
%clear (double* phases_rad, int phases_length);
%clear (signed char* valid_phases, int phases_length_bool);
%clear (double* freq_ddm, int delay_samples);
%clear (double* lag_ddm, int freq_samples);
%clear (double pos_pow_Max[2]);
%clear (double* powerLagHolo_singleLag, int fft_samples);

%apply (double* IN_ARRAY1, int DIM1) {(double* waveform_in, int wav_in_length)}
%apply (float* IN_ARRAY1, int DIM1) {(float* float_waveform_in, int wav_in_length)}
%apply (float* IN_ARRAY1, int DIM1) {(float* norm_waveform_in, int norm_wav_in_length)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* waveform_out, int wav_out_length)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* range, int size_range)}
%include "../src/waveform_power.h"
%clear (double* waveform_in, int wav_in_length);
%clear (float* float_waveform_in, int wav_in_length);
%clear (float* norm_waveform_in, int norm_wav_in_length);
%clear (double* waveform_out, int wav_out_length);
%clear (double* range, int size_range);

%apply (double IN_ARRAY1[ANY]) {(double antenna_vector_BF_E_in[3])}
%apply (double IN_ARRAY1[ANY]) {(double antenna_vector_BF_H_in[3])}
%apply (double IN_ARRAY1[ANY]) {(double antenna_vector_BF_k_in[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double antenna_vector_BF_E_out[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double antenna_vector_BF_H_out[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double antenna_vector_BF_k_out[3])}
%apply (double IN_ARRAY2[ANY][ANY]) {(double antenna_whole_pattern_dB_in[181][360])}
%apply (double IN_ARRAY1[ANY]) {(double antenna_full_pattern_dB_in[360])}
%apply (double IN_ARRAY1[ANY]) {(double antenna_half_pattern_dB_in[181])}
%apply (double* IN_ARRAY1, int DIM1) {(double* antenna_angles_deg, int angles_length)}
%apply (double* IN_ARRAY1, int DIM1) {(double* antenna_pattern_dB_in, int pattern_length)}
%apply (double IN_ARRAY1[ANY]) {(double antenna_full_pattern_E_dB_in[360])}
%apply (double IN_ARRAY1[ANY]) {(double antenna_full_pattern_H_dB_in[360])}
%apply (double IN_ARRAY1[ANY]) {(double antenna_half_pattern_E_dB_in[181])}
%apply (double IN_ARRAY1[ANY]) {(double antenna_half_pattern_H_dB_in[181])}
%apply (double* IN_ARRAY1, int DIM1) {(double* antenna_angles_E_deg, int angles_E_length)}
%apply (double* IN_ARRAY1, int DIM1) {(double* antenna_pattern_E_dB_in, int pattern_E_length)}
%apply (double* IN_ARRAY1, int DIM1) {(double* antenna_angles_H_deg, int angles_H_length)}
%apply (double* IN_ARRAY1, int DIM1) {(double* antenna_pattern_H_dB_in, int pattern_H_length)}
%apply (double ARGOUT_ARRAY2[ANY][ANY]) {(double antenna_pattern_dB_out[181][360])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double antenna_pattern_E_dB_out[360])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double antenna_pattern_H_dB_out[360])}
%apply (double IN_ARRAY1[ANY]) {(double incvector[3])}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* element_pos_in, int num_elem_in, int plane_dim)}
%apply (double* IN_ARRAY1, int DIM1) {(double* phase_delay_in, int num_elems_in)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* phase_delay_out, int num_elems_out)}
%apply (double ARGOUT_ARRAY2[ANY][ANY]) {(double array_factor_dB_out[181][360])}
%apply (double IN_ARRAY1[ANY]) {(double inertials[3])}
%apply (double IN_ARRAY1[ANY]) {(double posR_km[3])}
%apply (double IN_ARRAY1[ANY]) {(double posT_km[3])}
%include "../src/rf_front_end.h"
%clear (double antenna_vector_BF_E_in[3]);
%clear (double antenna_vector_BF_H_in[3]);
%clear (double antenna_vector_BF_k_in[3]);
%clear (double antenna_vector_BF_E_out[3]);
%clear (double antenna_vector_BF_H_out[3]);
%clear (double antenna_vector_BF_k_out[3]);
%clear (double antenna_whole_pattern_dB_in[181][360]);
%clear (double antenna_full_pattern_dB_in[360]);
%clear (double antenna_half_pattern_dB_in[181]);
%clear (double* antenna_angles_deg, int angles_length);
%clear (double* antenna_pattern_dB_in, int pattern_length);
%clear (double antenna_full_pattern_E_dB_in[360]);
%clear (double antenna_full_pattern_H_dB_in[360]);
%clear (double antenna_half_pattern_E_dB_in[181]);
%clear (double antenna_half_pattern_H_dB_in[181]);
%clear (double* antenna_angles_E_deg, int angles_E_length);
%clear (double* antenna_pattern_E_dB_in, int pattern_E_length);
%clear (double* antenna_angles_H_deg, int angles_H_length);
%clear (double* antenna_pattern_H_dB_in, int pattern_H_length);
%clear (double antenna_pattern_dB_out[181][360]);
%clear (double antenna_pattern_E_dB_out[360]);
%clear (double antenna_pattern_H_dB_out[360]);
%clear (double incvector[3]);
%clear (double* element_pos_in, int num_elem_in, int plane_dim);
%clear (double* phase_delay_in, int num_elems_in);
%clear (double* phase_delay_out, int num_elems_out);
%clear (double array_factor_dB_out[181][360]);
%clear (double inertials[3]);
%clear (double posR_km[3]);
%clear (double posT_km[3]);

%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* range_vec, int range_length)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* lambda_out, int lambda_length)}
%include "../src/gnss_composite.h"
%clear (double* range_vec, int range_length);
%clear (double* lambda_out, int lambda_length);

%apply (double IN_ARRAY1[ANY]) {(double posR_in[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double posR_out[3])}
%apply (double IN_ARRAY1[ANY]) {(double velR_in[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double velR_out[3])}
%apply (double IN_ARRAY1[ANY]) {(double posT_in[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double posT_out[3])}
%apply (double IN_ARRAY1[ANY]) {(double velT_in[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double velT_out[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double posS_out[3])}
%apply (double IN_ARRAY1[ANY]) {(double LonLatHeight_R_in[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double LonLatHeight_R_out[3])}
%apply (double IN_ARRAY1[ANY]) {(double LonLatHeight_T_in[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double LonLatHeight_T_out[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double elevAzimT_out[2])}
%apply (double IN_ARRAY1[ANY]) {(double vector_BF_in[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double vector_local_out[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double vector_ECEF_out[3])}
%apply (double IN_ARRAY1[ANY]) {(double vector_r_a_BF[3])}
%apply (double IN_ARRAY1[ANY]) {(double vector_r_t_BF[3])}
%apply (double IN_ARRAY1[ANY]) {(double rvv[2])}
%apply (double IN_ARRAY1[ANY]) {(double rhh[2])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double windup_phase_R_L[2])}
%include "../src/specular_geometry.h"
%clear (double posR_in[3]);
%clear (double posR_out[3]);
%clear (double velR_in[3]);
%clear (double velR_out[3]);
%clear (double posT_in[3]);
%clear (double posT_out[3]);
%clear (double velT_in[3]);
%clear (double velT_out[3]);
%clear (double posS_out[3]);
%clear (double LonLatHeight_R_in[3]);
%clear (double LonLatHeight_R_out[3]);
%clear (double LonLatHeight_T_in[3]);
%clear (double LonLatHeight_T_out[3]);
%clear (double elevAzimT_out[2]);
%clear (double vector_BF_in[3]);
%clear (double vector_local_out[3]);
%clear (double vector_ECEF_out[3]);
%clear (double vector_r_a_BF[3]);
%clear (double vector_r_t_BF[3]);
%clear (double rvv[2]);
%clear (double rhh[2]);
%clear (double windup_phase_R_L[2]);

%apply (double IN_ARRAY1[ANY]) {(double epsilon_upLayer[2])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double rvv[2])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double rhh[2])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double rco[2])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double rcross[2])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double tvv[2])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double thh[2])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double tco[2])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double tcross[2])}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* spectrum_in, int dim_x, int dim_y)}
%apply (double* IN_ARRAY1, int DIM1) {(double* kx_spec_in, int dim_kx)}
%apply (double* IN_ARRAY1, int DIM1) {(double* ky_spec_in, int dim_ky)}
%apply (double* IN_ARRAY1, int DIM1) {(double* spectrum_in_omnidir, int dim_x_omnidir)}
%apply (double* IN_ARRAY1, int DIM1) {(double* k_spec_in_omnidir, int dim_k_omnidir)}
%apply (double ARGOUT_ARRAY1[ANY]) {(double kx_ky_spectrum[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double k_spectrum[2])}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* wspeed_grid_in, int dim_wsx, int dim_wsy)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* wazim_grid_in, int dim_wax, int dim_way)}
%apply (double* IN_ARRAY1, int DIM1) {(double* lon_in, int dim_lon)}
%apply (double* IN_ARRAY1, int DIM1) {(double* lat_in, int dim_lat)}
%include "../src/reflecting_surface.h"
%clear (double epsilon_upLayer[2]);
%clear (double rvv[2]);
%clear (double rhh[2]);
%clear (double rco[2]);
%clear (double rcross[2]);
%clear (double tvv[2]);
%clear (double thh[2]);
%clear (double tco[2]);
%clear (double tcross[2]);
%clear (double* spectrum_in, int dim_x, int dim_y);
%clear (double* kx_spec_in, int dim_kx);
%clear (double* ky_spec_in, int dim_ky);
%clear (double* spectrum_in_omnidir, int dim_x_omnidir);
%clear (double* k_spec_in_omnidir, int dim_k_omnidir);
%clear (double kx_ky_spectrum[3]);
%clear (double k_spectrum[2]);
%clear (double* wspeed_grid_in, int dim_wsx, int dim_wsy);
%clear (double* wazim_grid_in, int dim_wax, int dim_way);
%clear (double* lon_in, int dim_lon);
%clear (double* lat_in, int dim_lat);

%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* dm_out, int dm_out_length)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* cov_out, int cov_out_length)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* wav_out, int wav_out_length)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* ddm_row_out, int ddm_row_out_length)}
%include "../src/ZavorotnyVoronovichModel.h"
%clear (double* dm_out, int dm_out_length);
%clear (double* cov_out, int cov_out_length);
%clear (double* wav_out, int wav_out_length);
%clear (double* ddm_row_out, int ddm_row_out_length);

%apply (double* IN_ARRAY1, int DIM1) {(double* depths_in, int num_depths)}
%apply (double* IN_ARRAY1, int DIM1) {(double* alpha_x_in, int num_alpha_x)}
%apply (double* IN_ARRAY1, int DIM1) {(double* alpha_y_in, int num_alpha_y)}
%apply (double* IN_ARRAY1, int DIM1) {(double* epsilon_r_in, int num_epsilon_r)}
%apply (double* IN_ARRAY1, int DIM1) {(double* epsilon_i_in, int num_epsilon_i)}
%apply (double* IN_ARRAY1, int DIM1) {(double* snow_dens_in, int num_snow_dens)}
%apply (double* IN_ARRAY1, int DIM1) {(double* elevations, int size_elevs)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yaws, int size_yaws)}
%apply (double IN_ARRAY1[ANY]) {(double elev_range[2])}
%apply (double IN_ARRAY1[ANY]) {(double azim_range[2])}
%apply (double IN_ARRAY1[ANY]) {(double time_range[2])}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* freq_LH, int samples_freq_LH)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* depth_LH, int samples_depth_LH)}
%apply (double ARGOUT_ARRAY1[ANY]) {(double pow_HV[2])}
%include "../src/MultiRaySingleRefl_model.h"
%clear (double* depths_in, int num_depths);
%clear (double* alpha_x_in, int num_alpha_x);
%clear (double* alpha_y_in, int num_alpha_y);
%clear (double* epsilon_r_in, int num_epsilon_r);
%clear (double* epsilon_i_in, int num_epsilon_i);
%clear (double* snow_dens_in, int num_snow_dens);
%clear (double* elevations, int size_elevs);
%clear (double* yaws, int size_yaws);
%clear (double elev_range[2]);
%clear (double azim_range[2]);
%clear (double time_range[2]);
%clear (double* freq_LH, int samples_freq_LH);
%clear (double* depth_LH, int samples_depth_LH);
%clear (double pow_HV[2]);
