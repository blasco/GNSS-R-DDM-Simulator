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

class Specular_geometry
{
	//IMPORTANT: Distances in km, angles in deg
	double posR_ECEF[3];
	double posT_ECEF[3];
	double velR_ECEF[3];
	double velT_ECEF[3];
	double posS_ECEF[3];
	double longitudeR, latitudeR;
	double longitudeT, latitudeT;
	double heightR, heightT;
	double local_heightR, local_heightT;
	double undulation;
	double roll, pitch, heading;
  public:
	double longitudeS, latitudeS;
	double elevation;
	double azimuthR, azimuthT;
	double geometric_delay;
	Specular_geometry( void );
	void dump_parameters( void );
	void set_ECEFpos_Receiver( double posR_in[3] );
	void get_ECEFpos_Receiver( double posR_out[3] );
	void set_ECEFvel_Receiver( double velR_in[3] );
	void get_ECEFvel_Receiver( double velR_out[3] );
	void set_ECEFpos_Transmitter( double posT_in[3] );
	void get_ECEFpos_Transmitter( double posT_out[3] );
	void set_ECEFvel_Transmitter( double velT_in[3] );
	void get_ECEFvel_Transmitter( double velT_out[3] );
	void get_ECEFpos_Specular( double posS_out[3] );
	void set_LongLatHeight_Receiver( double LonLatHeight_R_in[3] );
	void get_LongLatHeight_Receiver( double LonLatHeight_R_out[3] );
	void set_LongLatHeight_Transmitter( double LonLatHeight_T_in[3] );
	void get_LongLatHeight_Transmitter( double LonLatHeight_T_out[3] );
	void set_geometry_from_ElevHeightsSpec( double elev_in, double heightR_in, double heightT_in, double lonS_in, double latS_in, double azimT_in, double heightS_in, char computeUndu );
	void set_tangEarthVel_Receiver( double velocity, double specAzim_deg );
	void set_tangEarthVel_Transmitter( double velocity, double specAzim_deg );
	void set_Undulation( double undu );
	double get_Undulation( void );
	void read_ECEFpos_Receiver( const char* namefile, int week, double sow );
	void read_ECEFpos_Transmitter( const char* namefile, int week, double sow );
	void read_ECEFpos_GNSS_Transmitter( const char* namefile, int week, double sow, int prn, char gnss_identifier );
	void compute_specular_point( char computeUndu );
	void compute_ElevAzimT_from_receiver( double elevAzimT_out[2] );
	void set_inertials( double roll_in, double pitch_in, double heading_in );
	void rotate_vector_BF_to_local( double vector_BF_in[3], double vector_local_out[3] );
	void rotate_vector_BF_to_ECEF( double vector_BF_in[3], double vector_ECEF_out[3] );
	double compute_inertial_delay( double vector_BF_in[3] );
	void read_Inertials_Receiver( const char* namefile, int week, double sow );
	void compute_Beyerle_windup_direct( double vector_r_a_BF[3], double vector_r_t_BF[3], int week, double sow, double windup_phase_R_L[2] );
	void compute_Beyerle_windup_reflected( double vector_r_a_BF[3], double vector_r_t_BF[3], double rvv[2], double rhh[2], int week, double sow, double windup_phase_R_L[2] );
};
