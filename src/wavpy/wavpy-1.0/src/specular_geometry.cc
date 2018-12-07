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
#include "specular_geometry.h"
#include "ancillary_functions.h"

Specular_geometry::Specular_geometry( void ){
	double norm_vec_TR, norm_vec_TS, norm_vec_SR;
	int i;
	//Set default values: Receiver at 3 km, Transmitter at 20000 km (GPS), both at lat/lon=(0,0) and static
	longitudeR = 0.0;
	latitudeR  = 0.0;
	heightR    = 3.0;
	Compute_ECEF_from_LatLonH(longitudeR, latitudeR, heightR, posR_ECEF);
	longitudeT = 0.0;
	latitudeT  = 0.0;
	heightT    = 20000.0;
	Compute_ECEF_from_LatLonH(longitudeT, latitudeT, heightT, posT_ECEF);
	undulation = 0.0;
	Compute_SpecPoint(posR_ECEF, posT_ECEF, posS_ECEF, &undulation, &elevation, &azimuthR, &azimuthT, &local_heightR, &local_heightT, &longitudeS, &latitudeS, false);
	norm_vec_TR = sqrt(pow((posR_ECEF[0]-posT_ECEF[0]),2.0) + pow((posR_ECEF[1]-posT_ECEF[1]),2.0) + pow((posR_ECEF[2]-posT_ECEF[2]),2.0));
	norm_vec_TS = sqrt(pow((posS_ECEF[0]-posT_ECEF[0]),2.0) + pow((posS_ECEF[1]-posT_ECEF[1]),2.0) + pow((posS_ECEF[2]-posT_ECEF[2]),2.0));
	norm_vec_SR = sqrt(pow((posR_ECEF[0]-posS_ECEF[0]),2.0) + pow((posR_ECEF[1]-posS_ECEF[1]),2.0) + pow((posR_ECEF[2]-posS_ECEF[2]),2.0));
	geometric_delay = norm_vec_TS + norm_vec_SR - norm_vec_TR;
	for(i=0; i<3; i++){
		velR_ECEF[i] = 0.0;
		velT_ECEF[i] = 0.0;
	}
	roll    = 0.0;
	pitch   = 0.0;
	heading = 0.0;
	return;
}

void Specular_geometry::dump_parameters( void ){
	printf("====================== RECEIVER =======================\n");
	printf("- ECEF position    [km] : %f %f %f\n", posR_ECEF[0], posR_ECEF[1], posR_ECEF[2]);
	printf("- ECEF velocity  [km/s] : %f %f %f\n", velR_ECEF[0], velR_ECEF[1], velR_ECEF[2]);
	printf("- Lon, Lat        [deg] : %f %f\n", longitudeR, latitudeR);
	printf("- Height, local-H  [km] : %f %f\n", heightR, local_heightR);
	printf("- Azim, Elevation [deg] : %f %f\n", azimuthR, elevation);
	printf("=======================================================\n");
	printf("===================== TRANSMITTER =====================\n");
	printf("- ECEF position    [km] : %f %f %f\n", posT_ECEF[0], posT_ECEF[1], posT_ECEF[2]);
	printf("- ECEF velocity  [km/s] : %f %f %f\n", velT_ECEF[0], velT_ECEF[1], velT_ECEF[2]);
	printf("- Lon, Lat        [deg] : %f %f\n", longitudeT, latitudeT);
	printf("- Height, local-H  [km] : %f %f\n", heightT, local_heightT);
	printf("- Azim, Elevation [deg] : %f %f\n", azimuthT, elevation);
	printf("=======================================================\n");
	printf("====================== SPECULAR =======================\n");
	printf("- ECEF position    [km] : %f %f %f\n", posS_ECEF[0], posS_ECEF[1], posS_ECEF[2]);
	printf("- Lon, Lat        [deg] : %f %f\n", longitudeS, latitudeS);
	printf("- Geometric delay  [km] : %f\n", geometric_delay);
	printf("- Undulation        [m] : %f\n", undulation*1000.0);
	printf("=======================================================\n");
	printf("====================== INERTIALS ======================\n");
	printf("- Roll, Pitch     [deg] : %f %f\n", roll, pitch);
	printf("- Heading         [deg] : %f\n", heading);
	printf("=======================================================\n");
	return;
}

void Specular_geometry::set_ECEFpos_Receiver( double posR_in[3] ){
	posR_ECEF[0] = posR_in[0];
	posR_ECEF[1] = posR_in[1];
	posR_ECEF[2] = posR_in[2];
	Compute_LatLonH_from_ECEF(posR_ECEF, &longitudeR, &latitudeR, &heightR);
	return;
}

void Specular_geometry::get_ECEFpos_Receiver( double posR_out[3] ){
	posR_out[0] = posR_ECEF[0];
	posR_out[1] = posR_ECEF[1];
	posR_out[2] = posR_ECEF[2];
	return;
}

void Specular_geometry::set_ECEFvel_Receiver( double velR_in[3] ){
	velR_ECEF[0] = velR_in[0];
	velR_ECEF[1] = velR_in[1];
	velR_ECEF[2] = velR_in[2];
	return;
}

void Specular_geometry::get_ECEFvel_Receiver( double velR_out[3] ){
	velR_out[0] = velR_ECEF[0];
	velR_out[1] = velR_ECEF[1];
	velR_out[2] = velR_ECEF[2];
	return;
}

void Specular_geometry::set_ECEFpos_Transmitter( double posT_in[3] ){
	posT_ECEF[0] = posT_in[0];
	posT_ECEF[1] = posT_in[1];
	posT_ECEF[2] = posT_in[2];
	Compute_LatLonH_from_ECEF(posT_ECEF, &longitudeT, &latitudeT, &heightT);
	return;
}

void Specular_geometry::get_ECEFpos_Transmitter( double posT_out[3] ){
	posT_out[0] = posT_ECEF[0];
	posT_out[1] = posT_ECEF[1];
	posT_out[2] = posT_ECEF[2];
	return;
}

void Specular_geometry::set_ECEFvel_Transmitter( double velT_in[3] ){
	velT_ECEF[0] = velT_in[0];
	velT_ECEF[1] = velT_in[1];
	velT_ECEF[2] = velT_in[2];
	return;
}

void Specular_geometry::get_ECEFvel_Transmitter( double velT_out[3] ){
	velT_out[0] = velT_ECEF[0];
	velT_out[1] = velT_ECEF[1];
	velT_out[2] = velT_ECEF[2];
	return;
}

void Specular_geometry::get_ECEFpos_Specular( double posS_out[3] ){
	posS_out[0] = posS_ECEF[0];
	posS_out[1] = posS_ECEF[1];
	posS_out[2] = posS_ECEF[2];
	return;
}

void Specular_geometry::set_LongLatHeight_Receiver( double LonLatHeight_R_in[3] ){
	longitudeR = LonLatHeight_R_in[0];
	latitudeR  = LonLatHeight_R_in[1];
	heightR    = LonLatHeight_R_in[2];
	Compute_ECEF_from_LatLonH(longitudeR, latitudeR, heightR, posR_ECEF);
	return;
}

void Specular_geometry::get_LongLatHeight_Receiver( double LonLatHeight_R_out[3] ){
	LonLatHeight_R_out[0] = longitudeR;
	LonLatHeight_R_out[1] = latitudeR;
	LonLatHeight_R_out[2] = heightR;
	return;
}

void Specular_geometry::set_LongLatHeight_Transmitter( double LonLatHeight_T_in[3] ){
	longitudeT = LonLatHeight_T_in[0];
	latitudeT  = LonLatHeight_T_in[1];
	heightT    = LonLatHeight_T_in[2];
	Compute_ECEF_from_LatLonH(longitudeT, latitudeT, heightT, posT_ECEF);
	return;
}

void Specular_geometry::get_LongLatHeight_Transmitter( double LonLatHeight_T_out[3] ){
	LonLatHeight_T_out[0] = longitudeT;
	LonLatHeight_T_out[1] = latitudeT;
	LonLatHeight_T_out[2] = heightT;
	return;
}

void Specular_geometry::set_geometry_from_ElevHeightsSpec( double elev_in, double heightR_in, double heightT_in, double lonS_in, double latS_in, double azimT_in, double heightS_in, char computeUndu ){
	double ECEF2ENU[3][3];
	double posT_local[3], posR_local[3], PhiLambdaH[3];
	double alpha, tge, z_specular, radial_dist, dist_ratio;
	double norm_vec_TR, norm_vec_TS, norm_vec_SR;
	if((elev_in <= 0.0)||(elev_in > 90.0)){
		printf("ERROR! Incorrect elevation.\n");
		return;
	}
	if((heightR_in <= 0.0)||(heightT_in <= 0.0)){
		printf("ERROR! Incorrect receiver and/or transmitter height.\n");
		return;
	}
	if(computeUndu == 1){
		undulation = get_undulation(lonS_in, latS_in);
	}else{
		undulation = 0.0;
	}
	undulation = undulation + heightS_in;
	Compute_ECEF_from_LatLonH(lonS_in, latS_in, undulation, posS_ECEF);
	//Earth centered system where Z-axis points towards the specular point and Y-axis towards the receiver
	z_specular = sqrt(posS_ECEF[0]*posS_ECEF[0] + posS_ECEF[1]*posS_ECEF[1] + posS_ECEF[2]*posS_ECEF[2]);
	tge = tan(elev_in*PI_NUM/180.0);
	//Transmitter case
	radial_dist = heightT_in + z_specular - undulation;
	dist_ratio = z_specular/radial_dist;
	alpha = asin((-dist_ratio*tge + sqrt(tge*tge + 1.0 - dist_ratio*dist_ratio))/(tge*tge + 1.0));
	posT_local[0] = 0.0;
	posT_local[1] = radial_dist*sin(alpha);
	posT_local[2] = radial_dist*cos(alpha);
	//Receiver case
	radial_dist = heightR_in + z_specular - undulation;
	dist_ratio = z_specular/radial_dist;
	alpha = asin((-dist_ratio*tge + sqrt(tge*tge + 1.0 - dist_ratio*dist_ratio))/(tge*tge + 1.0));
	posR_local[0] = 0.0;
	posR_local[1] = -radial_dist*sin(alpha);
	posR_local[2] = radial_dist*cos(alpha);
	//Z-axis rotation at the local system, where Y-axis points towards the PRN sat, to convert to ENU-specular
	Rot3DaxisZ(posT_local, -azimT_in);
	Rot3DaxisZ(posR_local, -azimT_in);
	//ENU-specular to ECEF (using ECEF2ENU^(-1) <transpose equals inverse in this case>)
	PhiLambdaH[0] = lonS_in;
	PhiLambdaH[1] = latS_in;
	PhiLambdaH[2] = undulation;
	XYZ2GEOID(posS_ECEF, PhiLambdaH, ECEF2ENU);
	posT_ECEF[0] = ECEF2ENU[0][0]*posT_local[0] + ECEF2ENU[1][0]*posT_local[1] + ECEF2ENU[2][0]*posT_local[2];
	posT_ECEF[1] = ECEF2ENU[0][1]*posT_local[0] + ECEF2ENU[1][1]*posT_local[1] + ECEF2ENU[2][1]*posT_local[2];
	posT_ECEF[2] = ECEF2ENU[0][2]*posT_local[0] + ECEF2ENU[1][2]*posT_local[1] + ECEF2ENU[2][2]*posT_local[2];
	Compute_LatLonH_from_ECEF(posT_ECEF, &longitudeT, &latitudeT, &heightT);
	posR_ECEF[0] = ECEF2ENU[0][0]*posR_local[0] + ECEF2ENU[1][0]*posR_local[1] + ECEF2ENU[2][0]*posR_local[2];
	posR_ECEF[1] = ECEF2ENU[0][1]*posR_local[0] + ECEF2ENU[1][1]*posR_local[1] + ECEF2ENU[2][1]*posR_local[2];
	posR_ECEF[2] = ECEF2ENU[0][2]*posR_local[0] + ECEF2ENU[1][2]*posR_local[1] + ECEF2ENU[2][2]*posR_local[2];
	Compute_LatLonH_from_ECEF(posR_ECEF, &longitudeR, &latitudeR, &heightR);
	//Compute specular point with the results obtained
	Compute_SpecPoint(posR_ECEF, posT_ECEF, posS_ECEF, &undulation, &elevation, &azimuthR, &azimuthT, &local_heightR, &local_heightT, &longitudeS, &latitudeS, false);
	norm_vec_TR = sqrt(pow((posR_ECEF[0]-posT_ECEF[0]),2.0) + pow((posR_ECEF[1]-posT_ECEF[1]),2.0) + pow((posR_ECEF[2]-posT_ECEF[2]),2.0));
	norm_vec_TS = sqrt(pow((posS_ECEF[0]-posT_ECEF[0]),2.0) + pow((posS_ECEF[1]-posT_ECEF[1]),2.0) + pow((posS_ECEF[2]-posT_ECEF[2]),2.0));
	norm_vec_SR = sqrt(pow((posR_ECEF[0]-posS_ECEF[0]),2.0) + pow((posR_ECEF[1]-posS_ECEF[1]),2.0) + pow((posR_ECEF[2]-posS_ECEF[2]),2.0));
	geometric_delay = norm_vec_TS + norm_vec_SR - norm_vec_TR;
	return;
}

void Specular_geometry::set_tangEarthVel_Receiver( double velocity, double specAzim_deg ){
	int i;
	double x_unitary[3], y_unitary[3], aux[3];
	double norm_aux, specAzim_rad;
	//Compute unitary vectors in a local system centered at the receiver, plane XY parallel to tangecial Earth, X pointing towards the specular point and Z pointing towards the Earth's center
	vector3Prod(posS_ECEF, posR_ECEF, aux);
	norm_aux = norm3vec(aux);
	for(i=0; i<3; i++){
		y_unitary[i] = aux[i]/norm_aux;
	}
	vector3Prod(posR_ECEF, y_unitary, aux);
	norm_aux = norm3vec(aux);
	for(i=0; i<3; i++){
		x_unitary[i] = aux[i]/norm_aux;
	}
	//Compute velocity vector with pointing direction given by specular azimuth
	specAzim_rad = specAzim_deg*PI_NUM/180.0;
	for(i=0; i<3; i++){
		velR_ECEF[i] = (x_unitary[i]*cos(specAzim_rad) + y_unitary[i]*sin(specAzim_rad))*velocity;
	}
	return;
}

void Specular_geometry::set_tangEarthVel_Transmitter( double velocity, double specAzim_deg ){
	int i;
	double x_unitary[3], y_unitary[3], aux[3];
	double norm_aux, specAzim_rad;
	//Compute unitary vectors in a local system centered at the transmitter, plane XY parallel to tangecial Earth, X pointing towards the specular point and Z pointing towards the Earth's center
	vector3Prod(posS_ECEF, posT_ECEF, aux);
	norm_aux = norm3vec(aux);
	for(i=0; i<3; i++){
		y_unitary[i] = aux[i]/norm_aux;
	}
	vector3Prod(posT_ECEF, y_unitary, aux);
	norm_aux = norm3vec(aux);
	for(i=0; i<3; i++){
		x_unitary[i] = aux[i]/norm_aux;
	}
	//Compute velocity vector with pointing direction given by specular azimuth
	specAzim_rad = specAzim_deg*PI_NUM/180.0;
	for(i=0; i<3; i++){
		velT_ECEF[i] = (x_unitary[i]*cos(specAzim_rad) + y_unitary[i]*sin(specAzim_rad))*velocity;
	}
	return;
}

void Specular_geometry::set_Undulation( double undu ){
	undulation = undu;
	return;
}

double Specular_geometry::get_Undulation( void ){
	return undulation;
}

void Specular_geometry::read_ECEFpos_Receiver( const char* namefile, int week, double sow ){
	Read_Interpol_ECEFpos_file(namefile, week, sow, posR_ECEF, velR_ECEF);
	if((posR_ECEF[0]==0.0)&&(posR_ECEF[1]==0.0)&&(posR_ECEF[2]==0.0)){ //posR_ECEF at [0,0,0] means that there is no valid data
		longitudeR = -1;
		latitudeR = -1;
		heightR = -1;
	}else{
		Compute_LatLonH_from_ECEF(posR_ECEF, &longitudeR, &latitudeR, &heightR);
	}
	return;
}

void Specular_geometry::read_ECEFpos_Transmitter( const char* namefile, int week, double sow ){
	Read_Interpol_ECEFpos_file(namefile, week, sow, posT_ECEF, velT_ECEF);
	if((posT_ECEF[0]==0.0)&&(posT_ECEF[1]==0.0)&&(posT_ECEF[2]==0.0)){ //posT_ECEF at [0,0,0] means that there is no valid data
		longitudeT = -1;
		latitudeT = -1;
		heightT = -1;
	}else{
		Compute_LatLonH_from_ECEF(posT_ECEF, &longitudeT, &latitudeT, &heightT);
	}
	return;
}

void Specular_geometry::read_ECEFpos_GNSS_Transmitter( const char* namefile, int week, double sow, int prn, char gnss_identifier ){
	Read_Interpol_ECEFpos_SP3file(namefile, week, sow, prn, gnss_identifier, posT_ECEF, velT_ECEF);
	if((posT_ECEF[0]==0.0)&&(posT_ECEF[1]==0.0)&&(posT_ECEF[2]==0.0)){ //posT_ECEF at [0,0,0] means that there is no valid data
		longitudeT = -1;
		latitudeT = -1;
		heightT = -1;
	}else{
		Compute_LatLonH_from_ECEF(posT_ECEF, &longitudeT, &latitudeT, &heightT);
	}
	return;
}

void Specular_geometry::compute_specular_point( char computeUndu ){
	double norm_vec_TR, norm_vec_TS, norm_vec_SR;
	bool computeUndu_bool = false;
	if((posR_ECEF[0]==0.0)&&(posR_ECEF[1]==0.0)&&(posR_ECEF[2]==0.0)){
		printf("ERROR! Missing receiver's position.\n");
		return;
	}
	if((posT_ECEF[0]==0.0)&&(posT_ECEF[1]==0.0)&&(posT_ECEF[2]==0.0)){
		printf("ERROR! Missing transmitter's position.\n");
		return;
	}
	if(computeUndu==1){
		computeUndu_bool = true;
	}
	Compute_SpecPoint(posR_ECEF, posT_ECEF, posS_ECEF, &undulation, &elevation, &azimuthR, &azimuthT, &local_heightR, &local_heightT, &longitudeS, &latitudeS, computeUndu_bool);
	norm_vec_TR = sqrt(pow((posR_ECEF[0]-posT_ECEF[0]),2.0) + pow((posR_ECEF[1]-posT_ECEF[1]),2.0) + pow((posR_ECEF[2]-posT_ECEF[2]),2.0));
	norm_vec_TS = sqrt(pow((posS_ECEF[0]-posT_ECEF[0]),2.0) + pow((posS_ECEF[1]-posT_ECEF[1]),2.0) + pow((posS_ECEF[2]-posT_ECEF[2]),2.0));
	norm_vec_SR = sqrt(pow((posR_ECEF[0]-posS_ECEF[0]),2.0) + pow((posR_ECEF[1]-posS_ECEF[1]),2.0) + pow((posR_ECEF[2]-posS_ECEF[2]),2.0));
	geometric_delay = norm_vec_TS + norm_vec_SR - norm_vec_TR;
	return;
}

void Specular_geometry::compute_ElevAzimT_from_receiver( double elevAzimT_out[2] ){
	double PhiLambdaH_R[3];
	double ECEF2ENU_R[3][3];
	double azimT_R, elevT_R, heightT_R;
	//Coordinate conversion
	XYZ2GEOID(posR_ECEF, PhiLambdaH_R, ECEF2ENU_R);
	//Azimuth, elevation and local height of satellite with respect to the receiver
	XYZ2AZELH(posR_ECEF, posT_ECEF, ECEF2ENU_R, azimT_R, elevT_R, heightT_R);
	elevAzimT_out[0] = elevT_R;
	elevAzimT_out[1] = azimT_R;
	return;
}

void Specular_geometry::set_inertials( double roll_in, double pitch_in, double heading_in ){
	roll = roll_in;
	pitch = pitch_in;
	heading = heading_in;
	return;
}

void Specular_geometry::rotate_vector_BF_to_local( double vector_BF_in[3], double vector_local_out[3] ){
	double inertialDelay;
	inertialDelay = InertialDelayComputation(roll, pitch, heading, elevation, azimuthT, posR_ECEF, posS_ECEF, vector_BF_in, vector_local_out);
	return;
}

void Specular_geometry::rotate_vector_BF_to_ECEF( double vector_BF_in[3], double vector_ECEF_out[3] ){
	BF2ECEF(roll, pitch, heading, posR_ECEF, vector_BF_in, vector_ECEF_out);
	return;
}

double Specular_geometry::compute_inertial_delay( double vector_BF_in[3] ){
	double inertialDelay;
	double vector_local_out[3];
	inertialDelay = InertialDelayComputation(roll, pitch, heading, elevation, azimuthT, posR_ECEF, posS_ECEF, vector_BF_in, vector_local_out);
	return inertialDelay;
}

void Specular_geometry::read_Inertials_Receiver( const char* namefile, int week, double sow ){
	Read_Interpol_Inertials_file(namefile, week, sow, &roll, &pitch, &heading);
	return;
}

void Specular_geometry::compute_Beyerle_windup_direct( double vector_r_a_BF[3], double vector_r_t_BF[3], int week, double sow, double windup_phase_R_L[2] ){
	double rvv[2], rhh[2];
	rvv[0] = 0.0;
	rvv[1] = 0.0;
	rhh[0] = 0.0;
	rhh[1] = 0.0;
	Compute_Beyerle_windup(vector_r_a_BF, vector_r_t_BF, posR_ECEF, posT_ECEF, posS_ECEF, 0.0, 0.0, 0.0, 0.0, azimuthT, rvv, rhh, week, sow, true, windup_phase_R_L);
	return;
}

void Specular_geometry::compute_Beyerle_windup_reflected( double vector_r_a_BF[3], double vector_r_t_BF[3], double rvv[2], double rhh[2], int week, double sow, double windup_phase_R_L[2] ){
	Compute_Beyerle_windup(vector_r_a_BF, vector_r_t_BF, posR_ECEF, posT_ECEF, posS_ECEF, roll, pitch, heading, elevation, azimuthT, rvv, rhh, week, sow, false, windup_phase_R_L);
	return;
}

