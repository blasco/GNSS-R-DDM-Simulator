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

double range180( double val );
int Check_week_sow_3values_file( const char* namefile );
void Read_week_sow_3values_file( const char* namefile, int n_lines, double** values );
bool Interpol_ECEFpos_file( int week, double sow, int nlines, double* results, double** txyz );
void Read_Interpol_ECEFpos_file( const char* namefile, int week, double sow, double* pos_ECEF_out, double* vel_ECEF_out );
void Read_Interpol_Inertials_file( const char* namefile, int week, double sow, double* roll, double* pitch, double* heading );
int GPSday(int day, int month, int year);
void CheckSP3File( const char* namefile, int &nsats, int &nlines, int &sow_diff );
bool ReadSP3File( const char* namefile, double* tt, double** xx, double** yy, double** zz, double** vxx, double** vyy, double** vzz, int nlines, int sow_diff, int &weekGPS_ref, char gnss_identifier );
bool Interpol_sat_sp3( double sow, int prn, int nlines, double* results, double* tt, double** xx, double** yy, double** zz, bool satvel, double** vxx, double** vyy, double** vzz );
void Read_Interpol_ECEFpos_SP3file( const char* namefile, int week, double sow, int prn, char gnss_identifier, double* pos_ECEF_out, double* vel_ECEF_out );
void ReadUndulationFile( double* undMat, double* undLats, double* undLongs );
double Interpol_und( double latitude, double longitude, double* undMat, double* undLats, double* undLongs );
void XYZ2GEOID( double x[3], double PhiLambdaH[3], double ECEF2ENU[3][3] );
void XYZ2AZELH( double x[3], double r[3], double ECEF2ENU[3][3], double &Azimuth, double &Elevation, double &height_local );
void Compute_ECEF_from_LatLonH( double lon, double lat, double height, double* pos_ECEF );
void Compute_LatLonH_from_ECEF( double* pos_ECEF, double* longitude, double* latitude, double* height );
double get_undulation( double longitude, double latitude );
void Compute_SpecPoint( double* posR_ECEF, double* posT_ECEF, double* posS_ECEF, double* undu, double* elev, double* azimR, double* azimT, double* local_H_R, double* local_H_T, double* longitudeS, double* latitudeS, bool computeUndu );
void Rot3DaxisX( double vector_rot[3], double angle_deg );
void Rot3DaxisY( double vector_rot[3], double angle_deg );
void Rot3DaxisZ( double vector_rot[3], double angle_deg );
void ECEF2LocalRot( double satAzim, double vecIn[3], double vecOut[3], double ECEF2ENU[3][3] );
void Local2ECEF_Rot( double satAzim, double vecIn[3], double vecOut[3], double ECEF2ENU[3][3] );
void BF2ECEF( double roll_deg, double pitch_deg, double heading_deg, double pos_ref_ECEF[3], double vector_BF_in[3], double vector_ECEF_out[3] );
double InertialDelayComputation( double roll_deg, double pitch_deg, double heading_deg, double elevation_deg, double azimuth_deg, double posR_ECEF[3], double posS_ECEF[3], double vector_BF_in[3], double vector_local_out[3] );
void get_local_geometry_vectors( double azimuthT, double posS_ECEF_km[3], double posR_ECEF_km[3], double velR_ECEF_kms[3], double posT_ECEF_km[3], double velT_ECEF_kms[3], double posR_local[3], double velR_local[3], double posT_local[3], double velT_local[3] );
void vector3Prod( double vec1[3], double vec2[3], double resultvec[3] );
double scalar3Prod( double vec1[3], double vec2[3] );
double norm3vec( double vec[3] );
void Compute_Beyerle_windup( double vector_r_a_BF[3], double vector_r_t_BF[3], double posR_ECEF[3], double posT_ECEF[3], double posS_ECEF[3], double roll_deg, double pitch_deg, double heading_deg, double elevationT_deg, double azimuthT_deg, double rvv[2], double rhh[2], int week, double sow, bool computeDirectLink, double windup_phase_R_L[2] );
void Compute_GPS_L1_composite( double sampling_rate, double filter_BW, float weight_CA, float weight_PY, float weight_M, float weight_IM, float weight_L1C, double lambda[] );
void Compute_Galileo_E1_composite( double sampling_rate, double filter_BW, float weight_E1A, float weight_E1B, float weight_E1C, double lambda[] );
void Compute_BeiDou_B1( double sampling_rate, double filter_BW, float weight_B1I, double lambda[] );
bool check_input_antenna_pattern( double antenna_angles[], double antenna_pattern[], int pattern_length, double zero_level );
//void interpolate_Spline_antenna_pattern( double input_angles[], double input_pattern[], int pattern_length, double antenna_pattern[] );
void interpolate_antenna_pattern( double input_angles[], double input_pattern[], int pattern_length, double zero_level, double antenna_pattern[] );
void convolution_gsl( double vec1[], int size_vec1, double vec2[], int size_vec2, double vec_out[] );
void correlation_gsl( double vec1[], double vec2[], double vec_out[], int size_vec );
void Compute_power_spectrum( int samples, double* real_part, double* imag_part, double* power_spectrum );
double Compute_Doppler_bandwidth( int samples, double* doppler_freqs, double* doppler_map, double pos_pow_Max[2] );
void retrack_real_waveform( double* waveform, int wav_length, double delay );