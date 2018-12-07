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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include "wavpy_global_variables.h"
#include "ancillary_functions.h"

# define EPS  1.e-25      // epsilon for a determinant
# define PRECISIO  1.e-11 // Precision for 0 (norma 1)
# define MAXITER   50     // Max number of iteracions for s


using namespace std;


double range180( double val )
{
	double val_r360;
	if(fabs(val) <= 180.0){
		return val;
	}
	val_r360 = fmod(val, 360.0);
	if(val_r360 > 180.0){
		return (val_r360 - 360.0);
	}
	if(val_r360 < -180.0){
		return (val_r360 + 360.0);
	}
	return val_r360;
}


int Check_week_sow_3values_file( const char* namefile )
{
	stringstream lineRead;
	string line;
	int n_lines = 0;
	int i;
	ifstream fpi_pos (namefile);
	if (!fpi_pos){
		cout << "ERROR! Unable to read week_sow_3values file: " << namefile << endl ;
		return 0;
	}
	if(fpi_pos.is_open()){
		while(!fpi_pos.eof()){
			getline(fpi_pos,line);
			n_lines ++;
		}
		fpi_pos.close();
	}
	n_lines --;
	return n_lines;
}


void Read_week_sow_3values_file( const char* namefile, int n_lines, double** values )
{
	stringstream lineRead;
	string line;
	int i;
	ifstream fpi_pos (namefile);
	if (!fpi_pos){
		cout << "ERROR! Unable to read week_sow_3values file: " << namefile << endl ;
		return;
	}
	if(fpi_pos.is_open()){
		for(i=0;i<n_lines;i++){
			getline(fpi_pos,line);
			lineRead << line;
			lineRead.clear();
			lineRead.str(line);
			lineRead >> values[0][i] >> values[1][i] >> values[2][i] >> values[3][i] >> values[4][i];
		}
		fpi_pos.close();
	}
	return;
}


bool Interpol_ECEFpos_file( int week, double sow, int nlines, double* results, double** txyz )
{
	int i;
	int ind = 0;
	double second_GPS, ref_sec_GPS, prev_ref_sec_GPS;
	second_GPS = sow + week*604800.0;
	if(nlines <= 0){
		return false;
	}else{
		if((second_GPS <= (txyz[1][0] + txyz[0][0]*604800.0))||(second_GPS > (txyz[1][nlines-1] + txyz[0][nlines-1]*604800.0))){
			cout << "ERROR! Time out of week_sow_ECEFxyz file." << endl;
			return false;
		}else{
			while(second_GPS > (txyz[1][ind] + txyz[0][ind]*604800.0))
				ind++;
			ref_sec_GPS = txyz[1][ind] + txyz[0][ind]*604800.0;
			prev_ref_sec_GPS = txyz[1][ind-1] + txyz[0][ind-1]*604800.0;
			for(i=0;i<3;i++)
				results[i+3] = (txyz[i+2][ind] - txyz[i+2][ind-1])/(ref_sec_GPS - prev_ref_sec_GPS); //Velocity results in km/s
			for(i=0;i<3;i++)
				results[i] = txyz[i+2][ind-1] + results[i+3]*(second_GPS - prev_ref_sec_GPS);   //Position results
			return true;
		}
	}
}


void Read_Interpol_ECEFpos_file( const char* namefile, int week, double sow, double* pos_ECEF_out, double* vel_ECEF_out )
{
	double **week_sow_XYZ;
	double resultsPos[6];
	int nlines, i;
	bool correctData;
	nlines = Check_week_sow_3values_file(namefile);
	if(nlines <= 0){
		cout << "ERROR! Empty or wrong week_sow_ECEFxyz file: " << namefile << endl;
		pos_ECEF_out[0] = 0.0;
		pos_ECEF_out[1] = 0.0;
		pos_ECEF_out[2] = 0.0;
		vel_ECEF_out[0] = 0.0;
		vel_ECEF_out[1] = 0.0;
		vel_ECEF_out[2] = 0.0;
		return;
	}
	week_sow_XYZ = (double **) malloc(5*sizeof(double *));
	for(i=0; i<5; i++){
		week_sow_XYZ[i] = (double *) malloc(nlines*sizeof(double));
	}
	Read_week_sow_3values_file(namefile, nlines, week_sow_XYZ);
	correctData = Interpol_ECEFpos_file(week, sow, nlines, resultsPos, week_sow_XYZ);
	if(correctData){
		pos_ECEF_out[0] = resultsPos[0];
		pos_ECEF_out[1] = resultsPos[1];
		pos_ECEF_out[2] = resultsPos[2];
		vel_ECEF_out[0] = resultsPos[3];
		vel_ECEF_out[1] = resultsPos[4];
		vel_ECEF_out[2] = resultsPos[5];
	}else{
		pos_ECEF_out[0] = 0.0;
		pos_ECEF_out[1] = 0.0;
		pos_ECEF_out[2] = 0.0;
		vel_ECEF_out[0] = 0.0;
		vel_ECEF_out[1] = 0.0;
		vel_ECEF_out[2] = 0.0;
	}
	for(i=0; i<5; i++){
		free(week_sow_XYZ[i]);
	}
	free(week_sow_XYZ);
	return;
}


void Read_Interpol_Inertials_file( const char* namefile, int week, double sow, double* roll, double* pitch, double* heading )
{
	double **week_sow_RollPitchYaw;
	int nlines, i;
	int ind = 1;
	double second_GPS, ref_sec_GPS, prev_ref_sec_GPS, coef_low, coef_grt;
	second_GPS = sow + week*604800.0;
	nlines = Check_week_sow_3values_file(namefile);
	if(nlines <= 0){
		cout << "ERROR! Empty or wrong week_sow_RollPitchYaw file: " << namefile << endl;
		return;
	}
	week_sow_RollPitchYaw = (double **) malloc(5*sizeof(double *));
	for(i=0; i<5; i++){
		week_sow_RollPitchYaw[i] = (double *) malloc(nlines*sizeof(double));
	}
	Read_week_sow_3values_file(namefile, nlines, week_sow_RollPitchYaw);
	if((second_GPS < (week_sow_RollPitchYaw[1][0] + week_sow_RollPitchYaw[0][0]*604800.0))||(second_GPS > (week_sow_RollPitchYaw[1][nlines-1] + week_sow_RollPitchYaw[0][nlines-1]*604800.0))){
		cout << "ERROR! Time out of week_sow_RollPitchYaw file." << endl;
		return;
	}
	while(second_GPS > (week_sow_RollPitchYaw[1][ind] + week_sow_RollPitchYaw[0][ind]*604800.0)) ind++;
	ref_sec_GPS = week_sow_RollPitchYaw[1][ind] + week_sow_RollPitchYaw[0][ind]*604800.0;
	prev_ref_sec_GPS = week_sow_RollPitchYaw[1][ind-1] + week_sow_RollPitchYaw[0][ind-1]*604800.0;
	coef_low = (second_GPS - prev_ref_sec_GPS)/(ref_sec_GPS - prev_ref_sec_GPS);
	coef_grt = (ref_sec_GPS - second_GPS)/(ref_sec_GPS - prev_ref_sec_GPS);
	*roll = week_sow_RollPitchYaw[2][ind]*coef_low + week_sow_RollPitchYaw[2][ind-1]*coef_grt;
	*pitch = week_sow_RollPitchYaw[3][ind]*coef_low + week_sow_RollPitchYaw[3][ind-1]*coef_grt;
	*heading = week_sow_RollPitchYaw[4][ind]*coef_low + week_sow_RollPitchYaw[4][ind-1]*coef_grt;
	for(i=0; i<5; i++){
		free(week_sow_RollPitchYaw[i]);
	}
	free(week_sow_RollPitchYaw);
	return;
}


int GPSday(int day, int month, int year)
{
	int julday, julday_GPS;
	if((year < 1980)||(year > 2099)){
		printf("ERROR! Wrong date (%d/%d/%d) when computing Julian Day.\n", day, month, year);
		return 0;
	}
	if((day < 1)||(day > 31)||(month < 1)||(month > 12)){
		printf("ERROR! Wrong date (%d/%d/%d) when computing Julian Day.\n", day, month, year);
		return 0;
	}
	if(((month == 4)||(month == 6)||(month == 9)||(month == 11))&&(day>30)){
		printf("ERROR! Wrong date (%d/%d/%d) when computing Julian Day.\n", day, month, year);
		return 0;
	}
	if((month == 2)&&(((year%4 != 0)&&(day > 28))||(day > 29))){
		printf("ERROR! Wrong date (%d/%d/%d) when computing Julian Day.\n", day, month, year);
		return 0;
	}
	julday_GPS = 367*1980 - int((7*(1980 + int((1 + 9)/12)))/4) + int((275*1)/9) + 6 + 1721013;
	julday = 367*year - int((7*(year + int((month + 9)/12)))/4) + int((275*month)/9) + day + 1721013;
	return (julday - julday_GPS);
}


void CheckSP3File( const char* namefile, int &nsats, int &nlines, int &sow_diff )
{
	nsats = 0;
	nlines = 0;
	ifstream sp3file(namefile);
	if(!sp3file){
		cout << "ERROR! Unable to read SP3 file: " << namefile << endl;
		return;
	}
	int year, month, day, hour, minute;
	int gpsDay, gpsWeek, gpsDoW;
	int gpsWeek_0 = 0;
	double second, gpsSoW;
	double gpsSoW_0 = 0.;
	char mySym;
	bool n_lines_done = false;
	stringstream lineRead;
	string line, blank;
	int k = 0;
	if(sp3file.is_open()){
		while(!sp3file.eof()){
			getline(sp3file,line);
			k++;
			if(k == 3){
				blank=line.substr(4,2);
				nsats=atoi(blank.c_str());
			}
			if(!n_lines_done){
				if(line[0] == '*'){ //time data
					nlines++;
					lineRead.clear();
					lineRead.str(line);
					lineRead >>  mySym >> year >> month >> day >> hour >> minute >> second;
					gpsDay = GPSday(day, month, year);
					gpsWeek = int(gpsDay/7);
					gpsDoW  = gpsDay%7;
					gpsSoW  = gpsDoW*86400. + hour*3600. + minute*60. + second;
					if((gpsWeek > gpsWeek_0)||((gpsWeek == gpsWeek_0)&&(gpsSoW > gpsSoW_0))){
						sow_diff = int(gpsSoW - gpsSoW_0) + (gpsWeek - gpsWeek_0)*604800;
						gpsSoW_0 = gpsSoW;
						gpsWeek_0 = gpsWeek;
					}else{
						n_lines_done = true;
						cout << "WARNING! Non-consecutive time in sp3 file. Some data won't be taken into account." << endl;
					}
				}
			}
		}
		sp3file.close();
	}
	return;
}


bool ReadSP3File( const char* namefile, double* tt, double** xx, double** yy, double** zz, double** vxx, double** vyy, double** vzz, int nlines, int sow_diff, int &weekGPS_ref, char gnss_identifier )
{
	ifstream sp3file(namefile);
	if(!sp3file){
		cout << "ERROR! Unable to read SP3 file: " << namefile << endl;
		return false;
	}
	double x, y, z, vx, vy, vz;
	int i;
	int n = 0; 
	int year, month, day, hour, minute, prn;
	int gpsDay, gpsWeek, gpsDoW;
	int gpsWeek_0 = 0;
	double second, gpsSoW;
	double gpsSoW_0 = 0.;
	char mySym;
	bool read_vel = false;
	bool first_time_set = true;
	stringstream lineRead;
	string line, blank;
	int k = 0;
	for(i=0; i<32; i++)
		xx[0][i] = 999999.;  //xx[0][prn-1] is initialized because it will be used to know if there is data for a determined prn
	if(sp3file.is_open()){
		while(!sp3file.eof()){
			getline(sp3file,line);
			lineRead << line;
			if(line[0] == '*'){ //time data
				lineRead.clear();
				lineRead.str(line);
				lineRead >>  mySym >> year >> month >> day >> hour >> minute >> second;
				gpsDay = GPSday(day, month, year);
				gpsWeek = int(gpsDay/7);
				gpsDoW  = gpsDay%7;
				gpsSoW  = gpsDoW*86400. + hour*3600. + minute*60. + second;
				if(first_time_set){
					gpsSoW_0 = gpsSoW;
					gpsWeek_0 = gpsWeek;
					weekGPS_ref = gpsWeek_0;
					first_time_set = false;
				}
			}
			if(line[0] == 'P'&&line[1] == gnss_identifier){  //Position data
				line[0] = ' ';
				line[1] = ' ';
				lineRead.clear();
				lineRead.str(line);
				lineRead  >> prn >> x >> y >> z;
				n = int(0.5+(gpsSoW - gpsSoW_0 + double((gpsWeek-gpsWeek_0)*604800.))/sow_diff);
				if(n>=0 && n<nlines && (prn-1)<32){
					tt[n]        = gpsSoW + double((gpsWeek-gpsWeek_0)*604800.);
					xx[n][prn-1] = x;
					yy[n][prn-1] = y;
					zz[n][prn-1] = z;
				}
			}
			if(line[0] == 'V'&&line[1] == gnss_identifier){  //velocity data
				line[0] = ' ';
				line[1] = ' ';
				lineRead.clear();
				lineRead.str(line);
				lineRead  >> prn >> vx >> vy >> vz;
				n = int(0.5+(gpsSoW - gpsSoW_0 + double((gpsWeek-gpsWeek_0)*604800.))/sow_diff);
				if(n>=0 && n<nlines && (prn-1)<32){
					vxx[n][prn-1] = vx/10000; //velocity in dm/s to km/s
					vyy[n][prn-1] = vy/10000;
					vzz[n][prn-1] = vz/10000;
				}
				if(!read_vel){
					read_vel = true;
				}
			}
		}
		sp3file.close();
	}
	return read_vel;
}


bool Interpol_sat_sp3( double sow, int prn, int nlines, double* results, double* tt, double** xx, double** yy, double** zz, bool satvel, double** vxx, double** vyy, double** vzz )
{
	double tt_prn[10], xx_prn[10], yy_prn[10], zz_prn[10], vxx_prn[10], vyy_prn[10], vzz_prn[10];
	double pos;
	int ind = 0;
	int i, i_min;
	gsl_interp *workspace;
	gsl_interp_accel *accel;
	if(xx[0][prn-1]==999999.){
		cout << "ERROR! No data with PRN" << prn << " in .sp3 file" << endl;
		return false;
	}
	if((sow < tt[0]) || (sow > tt[nlines-1])){
		cout << "ERROR! Time out of .sp3 file." << endl;
		return false;
	}
	if(nlines < 10){
		cout << "ERROR! Not enough points in .sp3 file for polynomial-10 fitting." << endl;
		return false;
	}
	workspace = gsl_interp_alloc(gsl_interp_polynomial, 10);
	accel = gsl_interp_accel_alloc();
	while(sow > tt[ind]) ind++;
	i_min = max(0, (ind - 5));
	i_min = min(i_min, (nlines - 10));
	for(i=0; i<10; i++){
		tt_prn[i] = tt[i_min + i];
		xx_prn[i] = xx[i_min + i][prn-1];
		yy_prn[i] = yy[i_min + i][prn-1];
		zz_prn[i] = zz[i_min + i][prn-1];
	}
	gsl_interp_init(workspace, tt_prn, xx_prn, 10);
	results[0] = gsl_interp_eval(workspace, tt_prn, xx_prn, sow, accel);
	gsl_interp_init(workspace, tt_prn, yy_prn, 10);
	results[1] = gsl_interp_eval(workspace, tt_prn, yy_prn, sow, accel);
	gsl_interp_init(workspace, tt_prn, zz_prn, 10);
	results[2] = gsl_interp_eval(workspace, tt_prn, zz_prn, sow, accel);
	if(satvel){
		for(i=0; i<10; i++){
			vxx_prn[i] = vxx[i_min + i][prn-1];
			vyy_prn[i] = vyy[i_min + i][prn-1];
			vzz_prn[i] = vzz[i_min + i][prn-1];
		}
		gsl_interp_init(workspace, tt_prn, vxx_prn, 10);
		results[3] = gsl_interp_eval(workspace, tt_prn, vxx_prn, sow, accel);
		gsl_interp_init(workspace, tt_prn, vyy_prn, 10);
		results[4] = gsl_interp_eval(workspace, tt_prn, vyy_prn, sow, accel);
		gsl_interp_init(workspace, tt_prn, vzz_prn, 10);
		results[5] = gsl_interp_eval(workspace, tt_prn, vzz_prn, sow, accel);
	}else{
		if(sow<(tt[0] + 1.0)){   //Velocity results in km/s
			gsl_interp_init(workspace, tt_prn, xx_prn, 10);
			pos = gsl_interp_eval(workspace, tt_prn, xx_prn, (sow + 1.0), accel);
			results[3] = pos - results[0];
			gsl_interp_init(workspace, tt_prn, yy_prn, 10);
			pos = gsl_interp_eval(workspace, tt_prn, yy_prn, (sow + 1.0), accel);
			results[4] = pos - results[1];
			gsl_interp_init(workspace, tt_prn, zz_prn, 10);
			pos = gsl_interp_eval(workspace, tt_prn, zz_prn, (sow + 1.0), accel);
			results[5] = pos - results[2];
		}else{
			gsl_interp_init(workspace, tt_prn, xx_prn, 10);
			pos = gsl_interp_eval(workspace, tt_prn, xx_prn, (sow - 1.0), accel);
			results[3] = results[0] - pos;
			gsl_interp_init(workspace, tt_prn, yy_prn, 10);
			pos = gsl_interp_eval(workspace, tt_prn, yy_prn, (sow - 1.0), accel);
			results[4] = results[1] - pos;
			gsl_interp_init(workspace, tt_prn, zz_prn, 10);
			pos = gsl_interp_eval(workspace, tt_prn, zz_prn, (sow - 1.0), accel);
			results[5] = results[2] - pos;
		}
	}
	gsl_interp_free(workspace);
	gsl_interp_accel_free(accel);
	return true;
}


void Read_Interpol_ECEFpos_SP3file( const char* namefile, int week, double sow, int prn, char gnss_identifier, double* pos_ECEF_out, double* vel_ECEF_out )
{
	double **xx, **yy, **zz, **vxx, **vyy, **vzz;
	double *tt;
	double resultsT[6];  //[x, y, z, vx, vy, vz]
	int i, nsats, nlines, sow_diff, ref_GPSweek;
	bool correctData, sat_vel;
	double sow_diffWeek;
	CheckSP3File(namefile, nsats, nlines, sow_diff);
	if((nsats == 0)||(nlines == 0)){
		cout << "ERROR! Empty or wrong GPS.sp3 file: " << namefile << endl;
		pos_ECEF_out[0] = 0.0;
		pos_ECEF_out[1] = 0.0;
		pos_ECEF_out[2] = 0.0;
		vel_ECEF_out[0] = 0.0;
		vel_ECEF_out[1] = 0.0;
		vel_ECEF_out[2] = 0.0;
		return;
	}
	tt = (double *) malloc(nlines*sizeof(double));
	xx = (double **) malloc(nlines*sizeof(double *));
	yy = (double **) malloc(nlines*sizeof(double *));
	zz = (double **) malloc(nlines*sizeof(double *));
	vxx = (double **) malloc(nlines*sizeof(double *));
	vyy = (double **) malloc(nlines*sizeof(double *));
	vzz = (double **) malloc(nlines*sizeof(double *));
	for(i=0; i<nlines; i++){
		xx[i] = (double *) malloc(32*sizeof(double));
		yy[i] = (double *) malloc(32*sizeof(double));
		zz[i] = (double *) malloc(32*sizeof(double));
		vxx[i] = (double *) malloc(32*sizeof(double));
		vyy[i] = (double *) malloc(32*sizeof(double));
		vzz[i] = (double *) malloc(32*sizeof(double));
	}
	sat_vel = ReadSP3File(namefile, tt, xx, yy, zz, vxx, vyy, vzz, nlines, sow_diff, ref_GPSweek, gnss_identifier);
	sow_diffWeek = sow + (week - ref_GPSweek)*604800.0;
	correctData = Interpol_sat_sp3(sow_diffWeek, prn, nlines, resultsT, tt, xx, yy, zz, sat_vel, vxx, vyy, vzz);
	free(tt);
	for(i=0; i<nlines; i++){
		free(xx[i]);
		free(yy[i]);
		free(zz[i]);
		free(vxx[i]);
		free(vyy[i]);
		free(vzz[i]);
	}
	free(xx);
	free(yy);
	free(zz);
	free(vxx);
	free(vyy);
	free(vzz);
	if(correctData){
		pos_ECEF_out[0] = resultsT[0];
		pos_ECEF_out[1] = resultsT[1];
		pos_ECEF_out[2] = resultsT[2];
		vel_ECEF_out[0] = resultsT[3];
		vel_ECEF_out[1] = resultsT[4];
		vel_ECEF_out[2] = resultsT[5];
	}else{
		pos_ECEF_out[0] = 0.0;
		pos_ECEF_out[1] = 0.0;
		pos_ECEF_out[2] = 0.0;
		vel_ECEF_out[0] = 0.0;
		vel_ECEF_out[1] = 0.0;
		vel_ECEF_out[2] = 0.0;
	}
	return;
}


void ReadUndulationFile( double* undMat, double* undLats, double* undLongs )
{
	stringstream lineRead;
	string line;
	int i, j, k;
	double val0;
	string data_path = LOCAL_DATA_PATH;
	string namefile = data_path + "p85und_qrtdeg_egm96_to360.mean_tide";
	ifstream fpi (namefile.c_str());
	if(!fpi){
		cout << "ERROR! Unable to read undulation file: " << namefile << endl ;
		return;
	}
	if(fpi.is_open()){
		while(!fpi.eof()){
			for(i=0; i<681; i++){
				for(j=0; j<144; j++){
					getline(fpi,line);
					lineRead << line;
					lineRead.clear();
					lineRead.str(line);
					lineRead >> undMat[i*1440 + j*10] >> undMat[i*1440 + j*10 + 1] >> undMat[i*1440 + j*10 + 2] >> undMat[i*1440 + j*10 + 3] >> undMat[i*1440 + j*10 + 4] >> undMat[i*1440 + j*10 + 5] >> undMat[i*1440 + j*10 + 6] >> undMat[i*1440 + j*10 + 7] >> undMat[i*1440 + j*10 + 8] >> undMat[i*1440 + j*10 + 9];
					for(k=0; k<10; k++){
						undLongs[j*10 + k] = j*2.5 + k*0.25;
					}
				}
				undLats[i] = -(85. - i*0.25);
				getline(fpi,line);
				lineRead << line;
				lineRead.clear();
				lineRead.str(line);
				lineRead >> val0;
			}
		}
		fpi.close();
	}
	return;
}


double Interpol_und( double latitude, double longitude, double* undMat, double* undLats, double* undLongs )
{
	double resultm;
	double result;
	const gsl_interp2d_type *T = gsl_interp2d_bicubic;
	gsl_spline2d *spline = gsl_spline2d_alloc(T, 1440, 681);
	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
	gsl_interp_accel *yacc = gsl_interp_accel_alloc();
	gsl_spline2d_init(spline, undLongs, undLats, undMat, 1440, 681);
	resultm = gsl_spline2d_eval(spline, longitude, -latitude, xacc, yacc);
	gsl_spline2d_free(spline);
	gsl_interp_accel_free(xacc);
	gsl_interp_accel_free(yacc);
	result = resultm/1000.;
	return result;
}


void XYZ2GEOID( double x[3], double PhiLambdaH[3], double ECEF2ENU[3][3] )
{
	double  eccsq, equ_rad, lat_i, lat_p, longitude, h_i, h_p, rad_lat, rad_curve, tolerance, earth_flat, earth_rad;
	double A_[3];
	int niter, i;
	double ScaleElipsoide = 1.;
	A_[0] = A_EARTH_SEMIAXIS_KM/ScaleElipsoide;
	A_[1] = B_EARTH_SEMIAXIS_KM/ScaleElipsoide;
	A_[2] = C_EARTH_SEMIAXIS_KM/ScaleElipsoide;
	earth_rad = A_[0];
	earth_flat = (A_[0] - A_[2])/A_[1];
	tolerance = 0.000000001;
	equ_rad = sqrt(x[0]*x[0] + x[1]*x[1]) ;
	eccsq = 2.*earth_flat - earth_flat*earth_flat;
	lat_p = atan2(x[2], equ_rad);
	h_p = 0.; 
	int converged;
	converged = 0;
	niter = 0;
	double rad2deg = 180./PI_NUM;
	while(converged == 0){
		rad_curve = earth_rad/sqrt(1. - eccsq*sin(lat_p)*sin(lat_p));
		rad_lat   = equ_rad*(1. - eccsq*rad_curve/(rad_curve + h_p));
		lat_i = atan2(x[2], rad_lat);
		if(fabs(lat_i) < PI_NUM/4){
			h_i = equ_rad/cos(lat_i) - rad_curve;
		}else{
			h_i = x[2]/sin(lat_i) - (1. - eccsq)*rad_curve;
		}
		if((fabs(h_i-h_p) < tolerance)&& (fabs(lat_i-lat_p) < tolerance)){
			converged = 1;
		}
		niter = niter+1;
		if (niter > 50){
			printf("Failure to converge \n");
			printf("                  x \n");
			for (i=0;i <3;i++){
				printf("Input vectors   %12.6f  \n",x[i]);
			}
			converged=1;
		}
		h_p = h_i;
		lat_p = lat_i;
	}
	PhiLambdaH[0] = lat_i;
	PhiLambdaH[1] = atan2(x[1],x[0]);
	longitude = PhiLambdaH[1]; 
	if(PhiLambdaH[1] < 0.){
		PhiLambdaH[1] = PhiLambdaH[1] + 2.*PI_NUM;
	}
	PhiLambdaH[0] = PhiLambdaH[0]*rad2deg;
	PhiLambdaH[1] = PhiLambdaH[1]*rad2deg;
	PhiLambdaH[2] = h_i;
	//EAST component
	ECEF2ENU[0][0] = -sin(longitude);
	ECEF2ENU[0][1] = cos(longitude);
	ECEF2ENU[0][2] = 0.;
	//NORTH component
	ECEF2ENU[1][0] = -sin(lat_i)*cos(longitude);
	ECEF2ENU[1][1] = -sin(lat_i)*sin(longitude);
	ECEF2ENU[1][2] = cos(lat_i);
	//UP component
	ECEF2ENU[2][0] = cos(lat_i)*cos(longitude);
	ECEF2ENU[2][1] = cos(lat_i)*sin(longitude);
	ECEF2ENU[2][2] = sin(lat_i);
	return;
}


void XYZ2AZELH( double x[3], double r[3], double ECEF2ENU[3][3], double &Azimuth, double &Elevation, double &height_local )
{
	double b[3];
	double Horizontal, Norm_b;
	int i, j;
	double rad2deg = 180./PI_NUM;
	for(i=0; i<3; i++){
		b[i] = 0.;
		for(j=0; j<3; j++){
			b[i] = b[i] + ECEF2ENU[i][j]*(r[j] - x[j]);
		}
	}
	height_local = b[2];
	Norm_b = 0.;
	for (i=0;i <3;i++){
		Norm_b = Norm_b + b[i]*b[i];
	}
	Norm_b = sqrt(Norm_b);
	if(Norm_b > 0.){
		for (i=0;i<3;i++){b[i]=b[i]/Norm_b;}
		Horizontal = sqrt(b[0]*b[0]+b[1]*b[1]);
		Azimuth = 0.;
		Elevation = 90.;
		if(Horizontal > 0.){
			Azimuth = (PI_NUM/2. - atan2(b[1],b[0]))*rad2deg;
			if(Azimuth < 0.){
				Azimuth = Azimuth + 360.;
			}
			Elevation = asin(b[2])*rad2deg;
		}
	}
	else{
		Azimuth = 0.0;
		Elevation = 90.;
	}
	return;
}


void axnf3( double v[], double x[], double *xn, double dxn[] )
{
	*xn = sqrt((v[0] - x[0])*(v[0] - x[0]) + (v[1] - x[1])*(v[1] - x[1]) + (v[2] - x[2])*(v[2] - x[2]));
	dxn[0] = (x[0] - v[0])/(*xn);
	dxn[1] = (x[1] - v[1])/(*xn); 
	dxn[2] = (x[2] - v[2])/(*xn);
	return;
}


void auxf3( double v[], double x[], double *esc, double desc[], double A, double B, double C )
{
	*esc = (v[0] - x[0])*x[0]/(A*A) + (v[1] - x[1])*x[1]/(B*B) + (v[2] - x[2])*x[2]/(C*C);
	desc[0] = (v[0] - 2.e0*x[0])/(A*A);
	desc[1] = (v[1] - 2.e0*x[1])/(B*B);
	desc[2] = (v[2] - 2.e0*x[2])/(C*C);
	return;
}


double deter( double a[3][3] )
{
	return(a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[0][2]*a[1][1]*a[2][0] - a[0][1]*a[1][0]*a[2][2] - a[0][0]*a[1][2]*a[2][1]);
}
  

void F1( double x[], double *f1, double df1[], double A, double B, double C )
{
	/* Entra punt x. Calcula f1 i diferencial en punt */
	*f1 = (x[0]*x[0])/(A*A) + (x[1]*x[1])/(B*B) + (x[2]*x[2])/(C*C) - 1.e0;
	df1[0] = 2.e0*x[0]/(A*A);
	df1[1] = 2.e0*x[1]/(B*B);
	df1[2] = 2.e0*x[2]/(C*C);
	return;
}


void F2( double r[], double t[], double x[], double *f2, double df2[], int ind, double A, double B, double C )
{
	/* Entra punt x i ind. Calcula f2 i diferencial en punt.
	  Si ind=0 refa constants internes a partir de r i t, en cas 
	  contrari no */
	double rt0, rt1, rt2;
	static double cx, cy, cz, cxy, cxz, cyz;
	if(ind==0){
		rt0 = t[0] - r[0];
		rt1 = t[1] - r[1];
		rt2 = t[2] - r[2];
		cx = (rt1*r[2] - rt2*r[1])/(A*A);
		cy = (rt2*r[0] - rt0*r[2])/(B*B);
		cz = (rt0*r[1] - rt1*r[0])/(C*C);
		cxy = rt2*(1.e0/(A*A) - 1.e0/(B*B));
		cxz = rt1*(1.e0/(C*C) - 1.e0/(A*A));
		cyz = rt0*(1.e0/(B*B) - 1.e0/(C*C));
	}
	*f2 = (cx + cxy*x[1])*x[0] + (cy + cyz*x[2])*x[1] + (cz + cxz*x[0])*x[2];
	df2[0] = cx + cxy*x[1] + cxz*x[2];
	df2[1] = cy + cxy*x[0] + cyz*x[2];
	df2[2] = cz + cxz*x[0] + cyz*x[1];
	return;
}


void F3( double r[], double t[], double x[], double *f3, double df3[], double A, double B, double C )
{
	/* Entra punt x i vect r, t. Calcula f3 i diferencial en punt. */
	double xnr, xnt, dxnr[3], dxnt[3], er, et, der[3], det[3];
	axnf3(r, x, &xnr, dxnr);
	axnf3(t, x, &xnt, dxnt);
	auxf3(r, x, &er, der, A, B, C);
	auxf3(t, x, &et, det, A, B, C);
	*f3 = xnt*er - xnr*et;
	df3[0] = dxnt[0]*er + xnt*der[0] - dxnr[0]*et - xnr*det[0];
	df3[1] = dxnt[1]*er + xnt*der[1] - dxnr[1]*et - xnr*det[1];
	df3[2] = dxnt[2]*er + xnt*der[2] - dxnr[2]*et - xnr*det[2];
	return;
}


void CiRflx( double r[], double t[], double x[], double A, double B, double C )
{
	/* Troba una c.i. senzilla pel problema feta a partir de suposits plans */
	double xnt, xnr, xlp, ra[3], ta[3];
	int i;
	xnr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
	xnt = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
	for(i=0; i<3; i++) ra[i] = r[i]/xnr;
	for(i=0; i<3; i++) ta[i] = t[i]/xnt;
	xlp = xnt/(xnt+xnr);
	for(i=0; i<3; i++) x[i] = (A + B)*0.5e0*(xlp*ra[i] + (1.e0 - xlp)*ta[i]);
	return;
}


void IterRflx( double r[],double t[],double x[],int iter, int *isol, double A, double B, double C )
{
	/* Fa una iteracio de Newton per trobar la solucio.
	  Entra punt x i surt punt x modificat. 
	  La primera crida s'ha de fer amb iter=0 per inicialitzar ctants.
	  A la sortida dona isol=0 (no es te prec PRECISIO) o isol=1 en cas de tenir la solucio amb la prec PRECISIO. */
	double f[3], df[3][3], xda[3][3], det, xh[3];
	int i, j;
	F1(x, f, df[0], A, B, C);
	F2(r, t, x, (f + 1), df[1], iter, A, B, C);
	F3(r, t, x, (f + 2), df[2], A, B, C);
	*isol=1;
	for(i=0; i<3; i++) for(j=0; j<3; j++) xda[i][j]=df[i][j];
	det = deter(xda);
	if(fabs(f[0]) + fabs(f[1]) + fabs(f[2]) >= PRECISIO){
		for(i=0; i<3; i++) for(j=0; j<3; j++) xda[i][j] = df[i][j];
		det = deter(xda); 
		for(j=0; j<3; j++){
			if(j>0) for(i=0; i<3; i++) xda[i][j-1] = df[i][j-1];
			for(i=0; i<3; i++) xda[i][j] = -f[i];
			xh[j] = deter(xda)/det;
		}
		x[0] += xh[0];
		x[1] += xh[1];
		x[2] += xh[2];
		*isol=0;
	}
	if(fabs(det)<=EPS) *isol=2;
	return;
}
  

void TestS( double r[], double t[], double s[], double A, double B, double C )
{
	/*  Donat  s  calcula el valor de f1, f2 i f3.
	*  f1  val zero si el punt  s  esta sobre l'el.lipsoide
	*  f2  val zero si els angles que forma la normal a l'el.lipsoide al punt s amb els vector  sr  i  st  son iguals
	*  f3  val zero si els vectors r, t i n son coplanaris */
	double f1, f2, f3, n[3], rs[3], ts[3], nrs, nts;
	f1 = (s[0]*s[0])/(A*A) + (s[1]*s[1])/(B*B) + (s[2]*s[2])/(C*C) - 1.e0;
	n[0] = s[0]/(A*A);
	n[1] = s[1]/(B*B);
	n[2] = s[2]/(C*C);
	rs[0] = r[0] - s[0];
	rs[1] = r[1] - s[1];
	rs[2] = r[2] - s[2];
	ts[0] = t[0] - s[0];
	ts[1] = t[1] - s[1];
	ts[2] = t[2] - s[2];
	nrs = sqrt(rs[0]*rs[0] + rs[1]*rs[1] + rs[2]*rs[2]);
	nts = sqrt(ts[0]*ts[0] + ts[1]*ts[1] + ts[2]*ts[2]);
	rs[0] = rs[0]/nrs;
	rs[1] = rs[1]/nrs;
	rs[2] = rs[2]/nrs;
	ts[0] = ts[0]/nts;
	ts[1] = ts[1]/nts;
	ts[2] = ts[2]/nts;
	f2 = n[0]*rs[0] + n[1]*rs[1] + n[2]*rs[2] - n[0]*ts[0] - n[1]*ts[1] - n[2]*ts[2];
	f3 = rs[0]*ts[1]*n[2] + rs[1]*ts[2]*n[0] + rs[2]*ts[0]*n[1] - rs[2]*ts[1]*n[0] - rs[1]*ts[0]*n[2] - rs[0]*ts[2]*n[1];
	return;
}


void CalcS( double r[], double t[], double s[], int iout, int *iter, double A, double B, double C )
{
	/* Donats r i t mira de calcular s damunt de l'elipsoid de manera que el raig llencat per t es reflexa a s
	fins a r. si iout=0 tot O.K. si iout=1 no es troba solucio prou aprox o amb el nombre d'iteracions */
	int isol, itera;
	itera = 0;
	iout = 1;
	CiRflx(r, t, s, A, B, C);
	do{
		IterRflx(r, t, s, itera, &isol, A, B, C);
		itera++;
	}while(isol==0 && itera<MAXITER);
	if(isol==1) iout = 0;
	if(isol==2) iout = 1;
	*iter = itera;
	return; 
}


void rt2sp( double t[3], double r[3], double s[3], double e[4], double geo[4], double &undulation, int &iout )
{
	double sr[3], rt[3], Dst, Dsr, Ds, Dd, cainc, cainc2, cainc3;
	double rlocnev[3], rloccart[3], tlocnev[3], tloccart[3], rr[3][3];
	double LatGcentr, LatGgrph, Lon, epsi, rad, elix, eliy, eliz, ond;
	double A = A_EARTH_SEMIAXIS_KM;
	double B = B_EARTH_SEMIAXIS_KM;
	double C = C_EARTH_SEMIAXIS_KM;
	int i;
	int iter;
	iout = 0;
	ond = undulation;
	s[0] = 0.0e0;
	s[1] = 0.0e0;
	s[2] = 0.0e0;
	geo[0] = 0.0e0;
	geo[1] = 0.0e0;
	geo[2] = 0.0e0;
	geo[3] = 0.0e0;
	e[0] = 0.0e0;
	e[1] = 0.0e0;
	e[2] = 0.0e0;
	e[3] = 0.0e0;
	rad = 0.0e0;
	/* STEP 1: SCALE FACTOR TO APPLY TO THE REFERENCE ELLIPSOID TO ACCOUNT FOR GEOID ONDULATION
	 * The factor should be ond/Geograph_Radial where ond is the geoid ondulation and Geograph_Radial is the radial distance
	 * along the geographical latitude angle (NOT toward center of the Earth)
	 * LatGrph is related to LatGcentr as:
	 * tan(LatGcentr) = (1-f)^2*tan(lat_Grph)
	 * where f= (A-C)/A
	*/
	LatGcentr = atan2(r[2], sqrt(r[0]*r[0] + r[1]*r[1]));
	LatGgrph = atan(tan(LatGcentr)/((1.0e0 - ((A - C)/A))*(1.0e0 - ((A - C)/A))));
	/* then the Geograph_Radial is
	 *          r=A(1-e^2)/(1+e cos(LatGrph))
	 */
	epsi = sqrt(1.0e0 - (C*C/(A*A)));
	rad = A*(1.0e0 - epsi*epsi)/(1.0e0 + epsi*cos(LatGgrph));
	A = A*(1.0e0 + ond/rad);
	B = A;
	C = C*(1.0e0 + ond/rad);
	/* STEP 2: complementary information:
	 *  Longitude, Latitude--ellipsoidal, Height */
	Lon = atan2(r[1], r[0]);
	/* location of sub-rcv point on geoid (ellipsoid+ondul) in XYZ ecef 
	 * we need rad wrt center of Earth (it's wrt ell+ondul) because A,B,C values
	 * now are modified to account for geoid) */
	rad  = A*(1.0e0 - (A - C)/A*sin(LatGcentr)*sin(LatGcentr));
	elix = rad*cos(LatGcentr)*cos(Lon);
	eliy = rad*cos(LatGcentr)*sin(Lon);
	eliz = rad*sin(LatGcentr);
	/* cartesian */
	rloccart[0] = r[0] - elix;
	rloccart[1] = r[1] - eliy;
	rloccart[2] = r[2] - eliz;
	tloccart[0] = t[0] - r[0];
	tloccart[1] = t[1] - r[1];
	tloccart[2] = t[2] - r[2];
	/* rotation matrix to express rloccart in NEV */
	rr[0][0] = -sin(LatGcentr)*cos(Lon);
	rr[0][1] = -sin(LatGcentr)*sin(Lon);
	rr[0][2] = cos(LatGcentr);
	rr[1][0] = -sin(Lon);
	rr[1][1] = cos(Lon);
	rr[1][2] = 0.0e0;
	rr[2][0] = cos(LatGcentr)*cos(Lon);
	rr[2][1] = cos(LatGcentr)*sin(Lon);
	rr[2][2] = sin(LatGcentr);
	/* rlocnev: */
	rlocnev[0] = rr[0][0]*rloccart[0] + rr[0][1]*rloccart[1] + rr[0][2]*rloccart[2];
	rlocnev[1] = rr[1][0]*rloccart[0] + rr[1][1]*rloccart[1] + rr[1][2]*rloccart[2];
	rlocnev[2] = rr[2][0]*rloccart[0] + rr[2][1]*rloccart[1] + rr[2][2]*rloccart[2];
	tlocnev[0] = rr[0][0]*tloccart[0] + rr[0][1]*tloccart[1] + rr[0][2]*tloccart[2];
	tlocnev[1] = rr[1][0]*tloccart[0] + rr[1][1]*tloccart[1] + rr[1][2]*tloccart[2];
	tlocnev[2] = rr[2][0]*tloccart[0] + rr[2][1]*tloccart[1] + rr[2][2]*tloccart[2];
	/* prepare output vector including AZ from receiver local NEV */
	geo[0] = Lon*180.0e0/PI_NUM;
	geo[1] = LatGgrph*180.0e0/PI_NUM;
	geo[2] = rlocnev[2];
	geo[3] = 90 - atan2(tlocnev[0], tlocnev[1])*180.0e0/PI_NUM;
	/* STEP 3: estimation of specular point */
	iout=0;
	iter=0;
	CalcS(r, t, s, iout, &iter, A, B, C);
	if(iout==0){ 
		TestS(r, t, s, A, B, C);
		Dst = sqrt((s[0] - t[0])*(s[0] - t[0]) + (s[1] - t[1])*(s[1] - t[1]) + (s[2] - t[2])*(s[2] - t[2]));
		sr[0] = r[0] - s[0];
		sr[1] = r[1] - s[1];
		sr[2] = r[2] - s[2];
		Dsr = sqrt((s[0] - r[0])*(s[0] - r[0]) + (s[1] - r[1])*(s[1] - r[1]) + (s[2] - r[2])*(s[2] - r[2]));
		Ds = Dst + Dsr;
		rt[0] = t[0] - r[0];
		rt[1] = t[1] - r[1];
		rt[2] = t[2] - r[2];  
		Dd = sqrt((r[0]-t[0])*(r[0]-t[0])+(r[1]-t[1])*(r[1]-t[1])+(r[2]-t[2])*(r[2]-t[2]));
		cainc  = (s[0]*sr[0]+s[1]*sr[1]+s[2]*sr[2])/(sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2])*sqrt(sr[0]*sr[0]+sr[1]*sr[1]+sr[2]*sr[2]));
		cainc2 = (r[0]*sr[0]+r[1]*sr[1]+r[2]*sr[2])/(sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])*sqrt(sr[0]*sr[0]+sr[1]*sr[1]+sr[2]*sr[2]));
		cainc3 = (r[0]*rt[0]+r[1]*rt[1]+r[2]*rt[2])/(sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])*sqrt(rt[0]*rt[0]+rt[1]*rt[1]+rt[2]*rt[2]));
		e[3] = Ds - Dd;
		/* Elevations: wrt-SP, wrt-Rcv-SP, wrt-R-T (deg) */
		e[0] = 90.0e0 - atan2(sqrt(1.0e0 - cainc*cainc), cainc)*180.0e0/PI_NUM;
		e[1] = 90.0e0 - atan2(sqrt(1.0e0 - cainc2*cainc2), cainc2)*180.0e0/PI_NUM;
		e[2] = 90.0e0 - atan2(sqrt(1.0e0 - cainc3*cainc3), cainc3)*180.0e0/PI_NUM;
	}
	else{
		printf("ERROR! Problems in rt2sp.\n");
		iout = 1;
	}
	return;
}


void Compute_ECEF_from_LatLonH( double lon, double lat, double height, double* pos_ECEF )
{
	double A, f, C, S;
	A = A_EARTH_SEMIAXIS_KM;
	f = 1.0/298.257223563;
	C = 1.0/(sqrt((cos(lat*PI_NUM/180.)*cos(lat*PI_NUM/180.)) + ((1.0 - f)*(1.0 - f)*sin(lat*PI_NUM/180.)*sin(lat*PI_NUM/180.))));
	S = (1.0 - f)*(1.0 - f)*C;
	pos_ECEF[0] = (A*C + height)*cos(lat*PI_NUM/180.)*cos(lon*PI_NUM/180.);
	pos_ECEF[1] = (A*C + height)*cos(lat*PI_NUM/180.)*sin(lon*PI_NUM/180.);
	pos_ECEF[2] = (A*S + height)*sin(lat*PI_NUM/180.);
	//From http://mathforum.org/library/drmath/view/51832.html
	return;
}
  
  
void Compute_LatLonH_from_ECEF( double* pos_ECEF, double* longitude, double* latitude, double* height )
{
	double PhiLambdaH[3];
	double pos_vec[3];
	double ECEF2ENU[3][3];
	int i;
	if((pos_ECEF[0]==0.0)&&(pos_ECEF[1]==0.0)&&(pos_ECEF[2]==0.0)) return;
	for(i=0; i<3; i++){
		pos_vec[i] = pos_ECEF[i];
	}
	XYZ2GEOID(pos_vec, PhiLambdaH, ECEF2ENU);
	*longitude = PhiLambdaH[1];
	*latitude  = PhiLambdaH[0];
	*height    = PhiLambdaH[2];
	return;
}


double get_undulation( double longitude, double latitude )
{
	double *undMat, *undLats, *undLongs;
	double undulation;
	undMat = (double*) malloc (681*1440*sizeof(double));
	undLats = (double*) malloc (681*sizeof(double));
	undLongs = (double*) malloc (1440*sizeof(double));
	ReadUndulationFile(undMat, undLats, undLongs);
	undulation = Interpol_und(longitude, latitude, undMat, undLats, undLongs);
	free(undMat);
	free(undLats);
	free(undLongs);
	return undulation;
}


void Compute_SpecPoint( double* posR_ECEF, double* posT_ECEF, double* posS_ECEF, double* undu, double* elev, double* azimR, double* azimT, double* local_H_R, double* local_H_T, double* longitudeS, double* latitudeS, bool computeUndu )
{
	double PhiLambdaH[3];
	double pos_vec_S[3], pos_vec_R[3], pos_vec_T[3];
	double E[4], GEO[4];
	double ECEF2ENU[3][3];
	double local_height, azimuth, elevation, undulation;
	int iout = 0;
	int i;
	for(i=0; i<3; i++){
		pos_vec_R[i] = posR_ECEF[i];
		pos_vec_T[i] = posT_ECEF[i];
	}
	undulation = *undu;
	if(computeUndu){
		undulation = 0.0;
		rt2sp(pos_vec_T, pos_vec_R, pos_vec_S, E, GEO, undulation, iout);
		XYZ2GEOID(pos_vec_S, PhiLambdaH, ECEF2ENU);
		undulation = get_undulation(PhiLambdaH[0], PhiLambdaH[1]);
		*undu = undulation;
	}
	rt2sp(pos_vec_T, pos_vec_R, pos_vec_S, E, GEO, undulation, iout);
	for(i=0; i<3; i++){
		posS_ECEF[i] = pos_vec_S[i];
	}
	//Coordinate conversion
	XYZ2GEOID(pos_vec_S, PhiLambdaH, ECEF2ENU);
	*longitudeS = PhiLambdaH[1];
	*latitudeS  = PhiLambdaH[0];
	//Azimuth, elevation and local height of aircraft and satellite with respect to the specular point
	XYZ2AZELH(pos_vec_S, pos_vec_R, ECEF2ENU, azimuth, elevation, local_height);
	*azimR = azimuth;
	*local_H_R = local_height;
	XYZ2AZELH(pos_vec_S, pos_vec_T, ECEF2ENU, azimuth, elevation, local_height);
	*azimT = azimuth;
	*local_H_T = local_height;
	*elev = elevation;
	return;
}


void Rot3DaxisX( double vector_rot[3], double angle_deg )
{
	double angle_rad;
	double rotMat[2][2];
	double aux[2];
	angle_rad  = angle_deg*PI_NUM/180.;
	rotMat[0][0] = cos(angle_rad);
	rotMat[0][1] = -sin(angle_rad);
	rotMat[1][0] = sin(angle_rad);
	rotMat[1][1] = cos(angle_rad);
	aux[0] = rotMat[0][0]*vector_rot[1] + rotMat[0][1]*vector_rot[2];
	aux[1] = rotMat[1][0]*vector_rot[1] + rotMat[1][1]*vector_rot[2];
	vector_rot[1] = aux[0];
	vector_rot[2] = aux[1];
	return;
}


void Rot3DaxisY( double vector_rot[3], double angle_deg )
{
	double angle_rad;
	double rotMat[2][2];
	double aux[2];
	angle_rad  = angle_deg*PI_NUM/180.;
	rotMat[0][0] = cos(angle_rad);
	rotMat[0][1] = sin(angle_rad);
	rotMat[1][0] = -sin(angle_rad);
	rotMat[1][1] = cos(angle_rad);
	aux[0] = rotMat[0][0]*vector_rot[0] + rotMat[0][1]*vector_rot[2];
	aux[1] = rotMat[1][0]*vector_rot[0] + rotMat[1][1]*vector_rot[2];
	vector_rot[0] = aux[0];
	vector_rot[2] = aux[1];
	return;
}


void Rot3DaxisZ( double vector_rot[3], double angle_deg )
{
	double angle_rad;
	double rotMat[2][2];
	double aux[2];
	angle_rad  = angle_deg*PI_NUM/180.;
	rotMat[0][0] = cos(angle_rad);
	rotMat[0][1] = -sin(angle_rad);
	rotMat[1][0] = sin(angle_rad);
	rotMat[1][1] = cos(angle_rad);
	aux[0] = rotMat[0][0]*vector_rot[0] + rotMat[0][1]*vector_rot[1];
	aux[1] = rotMat[1][0]*vector_rot[0] + rotMat[1][1]*vector_rot[1];
	vector_rot[0] = aux[0];
	vector_rot[1] = aux[1];
	return;
}


void ECEF2LocalRot( double satAzim, double vecIn[3], double vecOut[3], double ECEF2ENU[3][3] )
{
	double aux[3];
	//first, convert vecIn from ECEF to ENU
	aux[0] = ECEF2ENU[0][0]*vecIn[0] + ECEF2ENU[0][1]*vecIn[1] + ECEF2ENU[0][2]*vecIn[2];
	aux[1] = ECEF2ENU[1][0]*vecIn[0] + ECEF2ENU[1][1]*vecIn[1] + ECEF2ENU[1][2]*vecIn[2];
	aux[2] = ECEF2ENU[2][0]*vecIn[0] + ECEF2ENU[2][1]*vecIn[1] + ECEF2ENU[2][2]*vecIn[2];
	//finally, the ENU-vector aux is rotated to be expressed in the local Specular Point system, where Y-axis points towards the PRN sat
	Rot3DaxisZ(aux, satAzim);
	vecOut[0] = aux[0];
	vecOut[1] = aux[1];
	vecOut[2] = aux[2];
	return;
}


void Local2ECEF_Rot( double satAzim, double vecIn[3], double vecOut[3], double ECEF2ENU[3][3] )
{
	double aux[3];
	aux[0] = vecIn[0];
	aux[1] = vecIn[1];
	aux[2] = vecIn[2];
	//First, counter-rotated Local vector with satellite's Azimuth to convert to ENU
	Rot3DaxisZ(aux, -satAzim);
	//Then, ENU to ECEF (using ECEF2ENU^(-1) <transpose equals inverse in this case>)
	vecOut[0] = ECEF2ENU[0][0]*aux[0] + ECEF2ENU[1][0]*aux[1] + ECEF2ENU[2][0]*aux[2];
	vecOut[1] = ECEF2ENU[0][1]*aux[0] + ECEF2ENU[1][1]*aux[1] + ECEF2ENU[2][1]*aux[2];
	vecOut[2] = ECEF2ENU[0][2]*aux[0] + ECEF2ENU[1][2]*aux[1] + ECEF2ENU[2][2]*aux[2];
	return;
}


void BF2ECEF( double roll_deg, double pitch_deg, double heading_deg, double pos_ref_ECEF[3], double vector_BF_in[3], double vector_ECEF_out[3] )
{
	double ECEF2ENU[3][3];
	double aux[3], aux2[3], PhiLambdaH[3];
	aux[0] = vector_BF_in[0];
	aux[1] = vector_BF_in[1];
	aux[2] = vector_BF_in[2];
	//Rotation of roll around craft's X-axis (along craft)
	Rot3DaxisX(aux, roll_deg);
	//Rotation of pitch around craft's Y-axis (across craft)
	Rot3DaxisY(aux, pitch_deg);
	//Rotation of yaw around craft's Z-axis (Z-axis down positive)
	Rot3DaxisZ(aux, heading_deg);
	//NED-receiver to ENU-receiver
	aux2[0] = aux[1];
	aux2[1] = aux[0];
	aux2[2] = -aux[2];
	//ENU-receiver to ECEF (using ECEF2ENU^(-1) <transpose equals inverse in this case>)
	XYZ2GEOID(pos_ref_ECEF, PhiLambdaH, ECEF2ENU);
	vector_ECEF_out[0] = ECEF2ENU[0][0]*aux2[0] + ECEF2ENU[1][0]*aux2[1] + ECEF2ENU[2][0]*aux2[2];
	vector_ECEF_out[1] = ECEF2ENU[0][1]*aux2[0] + ECEF2ENU[1][1]*aux2[1] + ECEF2ENU[2][1]*aux2[2];
	vector_ECEF_out[2] = ECEF2ENU[0][2]*aux2[0] + ECEF2ENU[1][2]*aux2[1] + ECEF2ENU[2][2]*aux2[2];
	return;
}


double InertialDelayComputation( double roll_deg, double pitch_deg, double heading_deg, double elevation_deg, double azimuth_deg, double posR_ECEF[3], double posS_ECEF[3], double vector_BF_in[3], double vector_local_out[3] )
{
	double ECEF2ENU[3][3];
	double vec_ECEF[3], PhiLambdaH[3];
	//Body-Frame to ECEF: Rotation of roll, pitch and yaw at local BF; NED to ENU (receiver); and ENU to ECEF
	BF2ECEF(roll_deg, pitch_deg, heading_deg, posR_ECEF, vector_BF_in, vec_ECEF);
	//ECEF to ENU-specular and then to local Specular Point system
	XYZ2GEOID(posS_ECEF, PhiLambdaH, ECEF2ENU);
	ECEF2LocalRot(azimuth_deg, vec_ECEF, vector_local_out, ECEF2ENU);
	//The result is the projection of the antenna vector over the reflected ray that would receive the upper antenna
	return (vector_local_out[1]*(-cos(elevation_deg*PI_NUM/180.)) + vector_local_out[2]*sin(elevation_deg*PI_NUM/180.));
}


void get_local_geometry_vectors( double azimuthT, double posS_ECEF_km[3], double posR_ECEF_km[3], double velR_ECEF_kms[3], double posT_ECEF_km[3], double velT_ECEF_kms[3], double posR_local[3], double velR_local[3], double posT_local[3], double velT_local[3] )
{
	int i;
	double ECEF2ENU_S[3][3];
	double posR[3], posT[3], PhiLambdaH_S[3];
	for(i=0; i<3; i++){
		posR[i] = posR_ECEF_km[i] - posS_ECEF_km[i];
		posT[i] = posT_ECEF_km[i] - posS_ECEF_km[i];
	}
        XYZ2GEOID(posS_ECEF_km, PhiLambdaH_S, ECEF2ENU_S);
	ECEF2LocalRot(azimuthT, posR, posR_local, ECEF2ENU_S);
	ECEF2LocalRot(azimuthT, velR_ECEF_kms, velR_local, ECEF2ENU_S);
	ECEF2LocalRot(azimuthT, posT, posT_local, ECEF2ENU_S);
	ECEF2LocalRot(azimuthT, velT_ECEF_kms, velT_local, ECEF2ENU_S);
	for(i=0; i<3; i++){
		posR_local[i] = posR_local[i]*1000.0;
		velR_local[i] = velR_local[i]*1000.0;
		posT_local[i] = posT_local[i]*1000.0;
		velT_local[i] = velT_local[i]*1000.0;
	}
	return;
}


void Compute_Sun_Position( int week, int sow, double posSun[3] )
{
	double g, q, EclipticLongitude, EclipticLatitude, Obliquity, TimeAfter2000Century, GST_h, GST0_h, UT_h, JulianDateAfter2000, JulianDate, sinDec, cosRAcosDec, sinRAcosDec, Dec, RA, EastLongitude;
	double AU_km = 149597871.0;
	double P_ECEF[3];
	//Universal time after 00:00
	UT_h = double(sow%86400)/3600.;
	//Julian Date  page 110 Misra
	JulianDate = 2444244.5 + week*7.0 + sow/86400.0;
	JulianDateAfter2000 = JulianDate - 2451545.0;
	TimeAfter2000Century = JulianDateAfter2000/36525.0;
	//Greenwich sidereal time
	//http://www.astro.uio.no/~bgranslo/aares/calculate.html
	GST0_h = 6.6974 + 2400.0515*TimeAfter2000Century;
	GST_h  = fmod((GST0_h + (366.2422/365.2422)*UT_h),24.0);
	//Sun ecliptic coordinates
	//http://aa.usno.navy.mil/faq/docs/SunApprox.php
	g = (357.529 + 0.98560028*JulianDateAfter2000);
	q = (280.459 + 0.98564736*JulianDateAfter2000);
	g = fmod(g,360.0);
	q = fmod(q,360.0);
	EclipticLongitude = (q + 1.915 * sin(g*PI_NUM/180.0) + 0.020*sin(2.0*g*PI_NUM/180.0))*PI_NUM/180.0;
	EclipticLatitude  = 0.;
	//Ecliptic obliquity
	//http://www.astro.uio.no/~bgranslo/aares/calculate.html
	Obliquity = (23.439 - 0.013*TimeAfter2000Century)*PI_NUM/180.0;
	//we assume that EclipticLatitude = 0
	//Formulas from   http://en.wikipedia.org/wiki/Ecliptic_longitude
	cosRAcosDec = cos(EclipticLongitude);
	sinRAcosDec = cos(Obliquity)*sin(EclipticLongitude);
	RA          = atan2(sinRAcosDec, cosRAcosDec);
	sinDec      = sin(Obliquity)*sin(EclipticLongitude);
	Dec         = asin(sinDec);
	EastLongitude = RA - GST_h*15.0*PI_NUM/180.0;
	P_ECEF[0] = cos(Dec)*cos(EastLongitude);
	P_ECEF[1] = cos(Dec)*sin(EastLongitude);
	P_ECEF[2] = sin(Dec);
	posSun[0] = AU_km*P_ECEF[0];
	posSun[1] = AU_km*P_ECEF[1];
	posSun[2] = AU_km*P_ECEF[2];
	return;
}


void vector3Prod( double vec1[3], double vec2[3], double resultvec[3] )
{
	resultvec[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
	resultvec[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
	resultvec[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
	return;
}


double scalar3Prod( double vec1[3], double vec2[3] )
{
	double res;
	res = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
	return res;
}


double norm3vec( double vec[3] )
{
	double res;
	res = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	return res;
}


void Compute_Beyerle_windup( double gb_r_a_BF[3], double gb_r_t_BF[3], double posR_ECEF[3], double posT_ECEF[3], double posS_ECEF[3], double roll_deg, double pitch_deg, double heading_deg, double elevationT_deg, double azimuthT_deg, double rvv[2], double rhh[2], int week, double sow, bool computeDirectLink, double windup_phase_R_L[2] )
{
	//Variables used for emulating the computations shown by Georg Beyerle, 2008 "Carrier phase wind-up in GPS reflectometry" (GB08)
	double temp_vec1[3], temp_vec2[3];
	double gb_r_a[3], gb_r_t[3]; //Receiver antenna vectors
	double gb_t_b[3], gb_t_b_Local[3], gb_t_t[3], gb_t_a[3], gb_k[3], gb_T_a_k_ECEF[3], gb_T_t_k_ECEF[3], gb_T_a_k[3], gb_T_t_k[3]; //Transmitter antenna vectors
	double gb_s_perp[3], gb_s_i_parl[3], gb_s_o_parl[3]; //Incidence plane vectors
	double gb_r_parl, gb_r_perp; //Reflection coefficients
	double gb_S_a[3], gb_S_t[3];
	double x_R, y_R, x_L, y_L;
	int i;
	double inertialD;
	double posSun_ECEF[3];
	double ECEF2ENU_S[3][3];
	double PhiLambdaH_S[3];
	// Coordinate conversion
        XYZ2GEOID(posS_ECEF, PhiLambdaH_S, ECEF2ENU_S);
	// Sun ECEF position
	Compute_Sun_Position(week, sow, posSun_ECEF);
	// Transmitter antenna vectors from GB08 with ECEF orientation
	//t_b points towards the Earth's center
	for(i=0 ;i<3; i++){
		gb_t_b[i] = -posT_ECEF[i]/norm3vec(posT_ECEF);
	}
	//t_t is perpendicular to the plane defined by t_b and the Earth-Sun line
	for(i=0; i<3; i++){
		temp_vec1[i] = posSun_ECEF[i] - posT_ECEF[i];
	}
	for(i=0; i<3; i++){
		temp_vec2[i] = temp_vec1[i]/norm3vec(temp_vec1);
	}
	vector3Prod(gb_t_b, temp_vec2, temp_vec1);
	for(i=0; i<3; i++){
		gb_t_t[i] = temp_vec1[i]/norm3vec(temp_vec1);
	}
	//t_a = t_t x t_b
	vector3Prod(gb_t_t, gb_t_b, temp_vec1);
	for(i=0; i<3; i++){
		gb_t_a[i] = temp_vec1[i]/norm3vec(temp_vec1);
	}
	//k points towards the specular point (direction of propagation)
	for(i=0; i<3; i++){
		temp_vec1[i] = posS_ECEF[i] - posT_ECEF[i];
	}
	for(i=0; i<3; i++){
		gb_k[i] = temp_vec1[i]/norm3vec(temp_vec1);
	}
	// Get T_a(k) and T_t(k) vectors from GB08 with ECEF orientation
	vector3Prod(gb_k, gb_t_a, temp_vec1);
	vector3Prod(temp_vec1, gb_k, gb_T_a_k_ECEF);
	vector3Prod(gb_k, gb_t_t, temp_vec1);
	vector3Prod(temp_vec1, gb_k, gb_T_t_k_ECEF);
	// Convert T_a(k), T_t(k) and t_b from ECEF to local frame centered at specular point, Y-axis pointing towards the satellite, Z up
	ECEF2LocalRot(azimuthT_deg, gb_T_a_k_ECEF, gb_T_a_k, ECEF2ENU_S);
	ECEF2LocalRot(azimuthT_deg, gb_T_t_k_ECEF, gb_T_t_k, ECEF2ENU_S);
	ECEF2LocalRot(azimuthT_deg, gb_t_b, gb_t_b_Local, ECEF2ENU_S);
	// Convert receiver antenna vectors from Aircraft Body Frame to local frame
	inertialD = InertialDelayComputation(roll_deg, pitch_deg, heading_deg, 0.0, azimuthT_deg, posR_ECEF, posS_ECEF, gb_r_a_BF, gb_r_a);
	inertialD = InertialDelayComputation(roll_deg, pitch_deg, heading_deg, 0.0, azimuthT_deg, posR_ECEF, posS_ECEF, gb_r_t_BF, gb_r_t);
	// Next steps are 
	if(computeDirectLink){
		// Compute windup phases components
		y_R = scalar3Prod(gb_T_t_k, gb_r_a) + scalar3Prod(gb_T_a_k, gb_r_t);
		x_R = scalar3Prod(gb_T_a_k, gb_r_a) - scalar3Prod(gb_T_t_k, gb_r_t);
		y_L = scalar3Prod(gb_T_t_k, gb_r_a) - scalar3Prod(gb_T_a_k, gb_r_t);
		x_L = scalar3Prod(gb_T_a_k, gb_r_a) + scalar3Prod(gb_T_t_k, gb_r_t);
	}else{
		// Incidence plane vectors (parallel to it)
		gb_s_i_parl[0] = 0.0;
		gb_s_i_parl[1] = -sin(elevationT_deg*PI_NUM/180.);
		gb_s_i_parl[2] = cos(elevationT_deg*PI_NUM/180.);
		gb_s_o_parl[0] = 0.0;
		gb_s_o_parl[1] = sin(elevationT_deg*PI_NUM/180.);
		gb_s_o_parl[2] = cos(elevationT_deg*PI_NUM/180.);
		// Compute reflection coefficients
		gb_r_parl = sqrt(rvv[0]*rvv[0] + rvv[1]*rvv[1]);
		gb_r_perp = sqrt(rhh[0]*rhh[0] + rhh[1]*rhh[1]);
		// Get S_a and S_t vectors (gb_s* components with fixed 0 values are avoided)
		gb_S_a[0] = gb_r_perp*scalar3Prod(gb_T_a_k, gb_s_perp)*gb_s_perp[0];
		gb_S_a[1] = gb_r_parl*scalar3Prod(gb_T_a_k, gb_s_i_parl)*gb_s_o_parl[1];
		gb_S_a[2] = gb_r_parl*scalar3Prod(gb_T_a_k, gb_s_i_parl)*gb_s_o_parl[2];
		gb_S_t[0] = gb_r_perp*scalar3Prod(gb_T_t_k, gb_s_perp)*gb_s_perp[0];
		gb_S_t[1] = gb_r_parl*scalar3Prod(gb_T_t_k, gb_s_i_parl)*gb_s_o_parl[1];
		gb_S_t[2] = gb_r_parl*scalar3Prod(gb_T_t_k, gb_s_i_parl)*gb_s_o_parl[2];
		// Compute windup phases components
		y_R = scalar3Prod(gb_S_t, gb_r_a) + scalar3Prod(gb_S_a, gb_r_t);
		x_R = scalar3Prod(gb_S_a, gb_r_a) - scalar3Prod(gb_S_t, gb_r_t);
		y_L = scalar3Prod(gb_S_t, gb_r_a) - scalar3Prod(gb_S_a, gb_r_t);
		x_L = scalar3Prod(gb_S_a, gb_r_a) + scalar3Prod(gb_S_t, gb_r_t);
	}
	windup_phase_R_L[0] = atan2(y_R, x_R);
	windup_phase_R_L[1] = atan2(y_L, x_L);
	return;
}


void positive_normalization( double vec_in[], int size_vec )
{
	double max_val = 0.0;
	int i;
	for(i=0; i<size_vec; i++){
		if(vec_in[i] > max_val){
			max_val = vec_in[i];
		}
	}
	if(max_val > 0.0){
		for(i=0; i<size_vec; i++){
			vec_in[i] = vec_in[i]/max_val;
		}
	}
	return;
}


void Compute_GPS_L1_composite( double sampling_rate, double filter_BW, float weight_CA, float weight_PY, float weight_M, float weight_IM, float weight_L1C, double lambda[] )
{
	int i, j;
	int lambda_size = 2*int(round(sampling_rate/GPS_CA_CHIP_RATE)) + 1;
	int samples_CA_chip = 16384; //2^14
	int samples_PY_chip = int(round(double(samples_CA_chip)/10.0));
	int samples_M_chip = int(round(double(samples_CA_chip)/5.0));
	int autocorr_size = 4*samples_CA_chip; //It has to be a power of 2
	int samples_filter = 2*(int(round(double(samples_CA_chip)*GPS_CA_CHIP_RATE/(2.0*filter_BW)))/2) + 1; //It has to be an odd integer
	double *aux_vec, *aux_lambda, *aux_vec_CA, *aux_vec_PY, *aux_vec_M, *aux_vec_IM, *aux_vec_L1C_BOC1, *aux_vec_L1C_BOC6, *filter;
	double amplitude_CA, amplitude_PY, amplitude_M, amplitude_IM, amplitude_L1C;
	if(lambda_size > autocorr_size){
		printf("ERROR! lambda_size > autocorr_size! Decrease sampling rate!\n");
		return;
	}
	if((weight_CA == 0.0)&&(weight_PY == 0.0)&&(weight_M == 0.0)&&(weight_IM == 0.0)&&(weight_L1C == 0.0)){
		for(i=0; i<lambda_size; i++){
			lambda[i] = 0.0;
		}
		return;
	}
	aux_vec = (double*) malloc (autocorr_size*sizeof(double));
	aux_lambda = (double*) malloc (autocorr_size*sizeof(double)); 
	aux_vec_CA = (double*) malloc (autocorr_size*sizeof(double)); 
	aux_vec_PY = (double*) malloc (autocorr_size*sizeof(double)); 
	aux_vec_M  = (double*) malloc (autocorr_size*sizeof(double));
	aux_vec_IM = (double*) malloc (autocorr_size*sizeof(double));
	aux_vec_L1C_BOC1  = (double*) malloc (autocorr_size*sizeof(double));
	aux_vec_L1C_BOC6  = (double*) malloc (autocorr_size*sizeof(double));
	filter = (double*) malloc (samples_filter*sizeof(double));
	//C/A code contribution
	if(weight_CA == 0.0){
		for(i=0; i<autocorr_size; i++){
			aux_vec_CA[i] = 0.0;
		}
	}else{
		for(i=0; i<autocorr_size; i++){
			if(i < samples_CA_chip){
				aux_vec[i] = 1.0/sqrt(double(samples_CA_chip));
			}else{
				aux_vec[i] = 0.0;
			}
		}
		//correl(aux_vec, aux_vec, aux_vec_CA);
		correlation_gsl(aux_vec, aux_vec, aux_vec_CA, autocorr_size);
		positive_normalization(aux_vec_CA, autocorr_size);
	}
	//P(Y) code contribution
	if(weight_PY == 0.0){
		for(i=0; i<autocorr_size; i++){
			aux_vec_PY[i] = 0.0;
		}
	}else{
		for(i=0; i<autocorr_size; i++){
			if(i < samples_PY_chip){
				aux_vec[i] = 1.0/sqrt(double(samples_PY_chip));
			}else{
				aux_vec[i] = 0.0;
			}
		}
		//correl(aux_vec, aux_vec, aux_vec_PY);
		correlation_gsl(aux_vec, aux_vec, aux_vec_PY, autocorr_size);
		positive_normalization(aux_vec_PY, autocorr_size);
	}
	//M code contribution BOCs(10, 5)
	if(weight_M == 0.0){
		for(i=0; i<autocorr_size; i++){
			aux_vec_M[i] = 0.0;
		}
	}else{
		for(i=0; i<autocorr_size; i++){
			if(i < samples_M_chip){
				if(sin(2.0*PI_NUM*double(i)/double(samples_PY_chip)) >= 0.0){
					aux_vec[i] = 1.0/sqrt(double(samples_M_chip));
				}else{
					aux_vec[i] = -1.0/sqrt(double(samples_M_chip));
				}
			}else{
				aux_vec[i] = 0.0;
			}
		}
		//correl(aux_vec, aux_vec, aux_vec_M);
		correlation_gsl(aux_vec, aux_vec, aux_vec_M, autocorr_size);
		positive_normalization(aux_vec_M, autocorr_size);
	}
	//IM code contribution BOCs(10, 10)
	if(weight_IM == 0.0){
		for(i=0; i<autocorr_size; i++){
			aux_vec_IM[i] = 0.0;
		}
	}else{
		for(i=0; i<autocorr_size; i++){
			if(i < samples_PY_chip){
				if(sin(2.0*PI_NUM*double(i)/double(samples_PY_chip)) >= 0.0){
					aux_vec[i] = 1.0/sqrt(double(samples_PY_chip));
				}else{
					aux_vec[i] = -1.0/sqrt(double(samples_PY_chip));
				}
			}else{
				aux_vec[i] = 0.0;
			}
		}
		//correl(aux_vec, aux_vec, aux_vec_IM);
		correlation_gsl(aux_vec, aux_vec, aux_vec_IM, autocorr_size);
		positive_normalization(aux_vec_IM, autocorr_size);
	}
	//L1C codes: 25%-data BOCs(1,1), 75%-pilot TMBOC(6, 1, 1/11) -> 29/33 Autocorr[BOCs(1,1)] + 4/33 Autocorr[BOCs(6,1)]
	if(weight_L1C == 0.0){
		for(i=0; i<autocorr_size; i++){
			aux_vec_L1C_BOC1[i] = 0.0;
			aux_vec_L1C_BOC6[i] = 0.0;
		}
	}else{
		for(i=0; i<autocorr_size; i++){
			if(i < samples_CA_chip){
				if(sin(2.0*PI_NUM*double(i)/double(samples_CA_chip)) >= 0.0){
					aux_vec[i] = 1.0/sqrt(double(samples_CA_chip));
				}else{
					aux_vec[i] = -1.0/sqrt(double(samples_CA_chip));
				}
			}else{
				aux_vec[i] = 0.0;
			}
		}
		//correl(aux_vec, aux_vec, aux_vec_L1C_BOC1);
		correlation_gsl(aux_vec, aux_vec, aux_vec_L1C_BOC1, autocorr_size);
		positive_normalization(aux_vec_L1C_BOC1, autocorr_size);
		for(i=0; i<autocorr_size; i++){
			if(i < samples_CA_chip){
				if(sin(12.0*PI_NUM*double(i)/double(samples_CA_chip)) >= 0.0){
					aux_vec[i] = 1.0/sqrt(double(samples_CA_chip));
				}else{
					aux_vec[i] = -1.0/sqrt(double(samples_CA_chip));
				}
			}else{
				aux_vec[i] = 0.0;
			}
		}
		//correl(aux_vec, aux_vec, aux_vec_L1C_BOC6);
		correlation_gsl(aux_vec, aux_vec, aux_vec_L1C_BOC6, autocorr_size);
		positive_normalization(aux_vec_L1C_BOC6, autocorr_size);
	}
	//Sum all code autocorrelations
	amplitude_CA = sqrt(pow(10.0, (POW_CA_Trans_dBW/10.0)))*double(weight_CA);
	amplitude_PY = sqrt(pow(10.0, (POW_PY_Trans_dBW/10.0)))*double(weight_PY);
	amplitude_M = sqrt(pow(10.0, (POW_M_Trans_dBW/10.0)))*double(weight_M);
	amplitude_IM = sqrt(pow(10.0, (POW_IM_Trans_dBW/10.0)))*double(weight_IM);
	amplitude_L1C = sqrt(pow(10.0, (POW_L1C_Trans_dBW/10.0)))*double(weight_L1C);
	for(i=0; i<(autocorr_size/2); i++){
		aux_vec[i + (autocorr_size/2)] = amplitude_CA*aux_vec_CA[i] + amplitude_PY*aux_vec_PY[i] + amplitude_M*aux_vec_M[i] + amplitude_IM*aux_vec_IM[i] + amplitude_L1C*(0.25*aux_vec_L1C_BOC1[i] + 0.75*((29.0/33.0)*aux_vec_L1C_BOC1[i] + (4.0/33.0)*aux_vec_L1C_BOC6[i]));
	}
	for(i=(autocorr_size/2); i<autocorr_size; i++){
		aux_vec[i - (autocorr_size/2)] = amplitude_CA*aux_vec_CA[i] + amplitude_PY*aux_vec_PY[i] + amplitude_M*aux_vec_M[i] + amplitude_IM*aux_vec_IM[i] + amplitude_L1C*(0.25*aux_vec_L1C_BOC1[i] + 0.75*((29.0/33.0)*aux_vec_L1C_BOC1[i] + (4.0/33.0)*aux_vec_L1C_BOC6[i]));
	}
	//Convolution with box filter
	for(i=0; i<samples_filter; i++){
		filter[i] = 1.0/double(samples_filter);
	}
	convolution_gsl(aux_vec, autocorr_size, filter, samples_filter, aux_lambda);
	for(i=0; i<lambda_size; i++){
		//lambda's center is at (lambda_size/2)
		//aux_lambda's center is at (autocorr_size/2)
		//the sample resolution is corrected multiplying by (samples_CA_chip*GPS_CA_CHIP_RATE/sampling_rate)
		j = int(round(double(i - (lambda_size/2))*(double(samples_CA_chip)*GPS_CA_CHIP_RATE/sampling_rate))) + (autocorr_size/2);
		lambda[i] = aux_lambda[j]*aux_lambda[j]/(aux_lambda[autocorr_size/2]*aux_lambda[autocorr_size/2]); //Normalized power
	}
	free(aux_vec);
	free(aux_lambda); 
	free(aux_vec_CA); 
	free(aux_vec_PY); 
	free(aux_vec_M);
	free(aux_vec_IM);
	free(aux_vec_L1C_BOC1);
	free(aux_vec_L1C_BOC6);
	return;
}


void Compute_Galileo_E1_composite(double sampling_rate, double filter_BW, float weight_E1A, float weight_E1B, float weight_E1C, double lambda[])
{
	int i, j;
	int lambda_size = 2*int(round(sampling_rate/GPS_CA_CHIP_RATE)) + 1;
	int samples_CA_chip = 16384; //2^14
	int samples_E1A_chip = int(round(double(samples_CA_chip)/2.5));
	int samples_E1A_cos = int(round(double(samples_CA_chip)/15.0));
	int autocorr_size = 4*samples_CA_chip; //It has to be a power of 2
	int samples_filter = 2*(int(round(double(samples_CA_chip)*GPS_CA_CHIP_RATE/(2.0*filter_BW)))/2) + 1; //It has to be an odd integer
	double aux_vec[autocorr_size], aux_lambda[autocorr_size], aux_vec_E1A[autocorr_size], aux_vec_E1B[autocorr_size], aux_vec_E1C[autocorr_size], filter[samples_filter];
	double amplitude_E1A, amplitude_E1B, amplitude_E1C;
	double alpha = sqrt(10.0/11.0);
	double beta = sqrt(1.0/11.0);
	if(lambda_size > autocorr_size){
		printf("ERROR! lambda_size > autocorr_size! Decrease sampling rate!\n");
		return;
	}
	if((weight_E1A == 0.0)&&(weight_E1B == 0.0)&&(weight_E1C == 0.0)){
		for(i=0; i<lambda_size; i++){
			lambda[i] = 0.0;
		}
		return;
	}
	//E1A code contribution BOCc(15, 2.5)
	if(weight_E1A == 0.0){
		for(i=0; i<autocorr_size; i++){
			aux_vec_E1A[i] = 0.0;
		}
	}else{
		for(i=0; i<autocorr_size; i++){
			if(i < samples_E1A_chip){
				if(cos(2.0*PI_NUM*double(i)/double(samples_E1A_cos)) >= 0.0){
					aux_vec[i] = 1.0;
				}else{
					aux_vec[i] = -1.0;
				}
			}else{
				aux_vec[i] = 0.0;
			}
		}
		correlation_gsl(aux_vec, aux_vec, aux_vec_E1A, autocorr_size);
		positive_normalization(aux_vec_E1A, autocorr_size);
	}
	//E1B code contribution CBOC(6, 1, 1/11)
	if(weight_E1B == 0.0){
		for(i=0; i<autocorr_size; i++){
			aux_vec_E1B[i] = 0.0;
		}
	}else{
		for(i=0; i<autocorr_size; i++){
			if(i < samples_CA_chip){
				if(sin(2.0*PI_NUM*double(i)/double(samples_CA_chip)) >= 0.0){
					aux_vec[i] = alpha;
				}else{
					aux_vec[i] = -alpha;
				}
				if(sin(12.0*PI_NUM*double(i)/double(samples_CA_chip)) >= 0.0){
					aux_vec[i] = aux_vec[i] + beta;
				}else{
					aux_vec[i] = aux_vec[i] - beta;
				}
			}else{
				aux_vec[i] = 0.0;
			}
		}
		correlation_gsl(aux_vec, aux_vec, aux_vec_E1B, autocorr_size);
		positive_normalization(aux_vec_E1B, autocorr_size);
	}
	//E1C code contribution CBOC(6, 1, 1/11)
	if(weight_E1C == 0.0){
		for(i=0; i<autocorr_size; i++){
			aux_vec_E1C[i] = 0.0;
		}
	}else{
		for(i=0; i<autocorr_size; i++){
			if(i < samples_CA_chip){
				if(sin(2.0*PI_NUM*double(i)/double(samples_CA_chip)) >= 0.0){
					aux_vec[i] = alpha;
				}else{
					aux_vec[i] = -alpha;
				}
				if(sin(12.0*PI_NUM*double(i)/double(samples_CA_chip)) >= 0.0){
					aux_vec[i] = aux_vec[i] - beta;
				}else{
					aux_vec[i] = aux_vec[i] + beta;
				}
			}else{
				aux_vec[i] = 0.0;
			}
		}
		correlation_gsl(aux_vec, aux_vec, aux_vec_E1C, autocorr_size);
		positive_normalization(aux_vec_E1C, autocorr_size);
	}
	//Sum all code autocorrelations
	amplitude_E1A = sqrt(pow(10.0, (POW_E1A_RecEarth_dBW/10.0)))*double(weight_E1A);
	amplitude_E1B = sqrt(pow(10.0, (POW_E1B_RecEarth_dBW/10.0)))*double(weight_E1B);
	amplitude_E1C = sqrt(pow(10.0, (POW_E1C_RecEarth_dBW/10.0)))*double(weight_E1C);
	for(i=0; i<(autocorr_size/2); i++){
		aux_vec[i + (autocorr_size/2)] = amplitude_E1A*aux_vec_E1A[i] + amplitude_E1B*aux_vec_E1B[i] + amplitude_E1C*aux_vec_E1C[i];
	}
	for(i=(autocorr_size/2); i<autocorr_size; i++){
		aux_vec[i - (autocorr_size/2)] = amplitude_E1A*aux_vec_E1A[i] + amplitude_E1B*aux_vec_E1B[i] + amplitude_E1C*aux_vec_E1C[i];
	}
	//Convolution with box filter
	for(i=0; i<samples_filter; i++){
		filter[i] = 1.0/double(samples_filter);
	}
	convolution_gsl(aux_vec, autocorr_size, filter, samples_filter, aux_lambda);
	for(i=0; i<lambda_size; i++){
		//lambda's center is at (lambda_size/2)
		//aux_lambda's center is at (autocorr_size/2)
		//the sample resolution is corrected multiplying by (samples_CA_chip*GPS_CA_CHIP_RATE/sampling_rate)
		j = int(round(double(i - (lambda_size/2))*(double(samples_CA_chip)*GPS_CA_CHIP_RATE/sampling_rate))) + (autocorr_size/2);
		lambda[i] = aux_lambda[j]*aux_lambda[j]/(aux_lambda[autocorr_size/2]*aux_lambda[autocorr_size/2]); //Normalized power
	}
	return;
}


void Compute_BeiDou_B1( double sampling_rate, double filter_BW, float weight_B1I, double lambda[] )
{
	int i, j;
	int lambda_size = 2*int(round(sampling_rate/GPS_CA_CHIP_RATE)) + 1;
	int samples_CA_chip = 16384; //2^14
	int samples_B1I_chip = int(round(double(samples_CA_chip)/2.0));
	int autocorr_size = 4*samples_CA_chip; //It has to be a power of 2
	int samples_filter = 2*(int(round(double(samples_CA_chip)*GPS_CA_CHIP_RATE/(2.0*filter_BW)))/2) + 1; //It has to be an odd integer
	double aux_vec[autocorr_size], aux_lambda[autocorr_size], aux_vec_B1I[autocorr_size], filter[samples_filter];
	if(lambda_size > autocorr_size){
		printf("ERROR! lambda_size > autocorr_size! Decrease sampling rate!\n");
		return;
	}
	//B1I code contribution
	if(weight_B1I == 0.0){
		for(i=0; i<lambda_size; i++){
			lambda[i] = 0.0;
		}
		return;
	}else{
		for(i=0; i<autocorr_size; i++){
			if(i < samples_B1I_chip){
				aux_vec[i] = 1.0/sqrt(double(samples_B1I_chip));
			}else{
				aux_vec[i] = 0.0;
			}
		}
		correlation_gsl(aux_vec, aux_vec, aux_vec_B1I, autocorr_size);
		positive_normalization(aux_vec_B1I, autocorr_size);
	}
	for(i=0; i<(autocorr_size/2); i++){
		aux_vec[i + (autocorr_size/2)] = aux_vec_B1I[i];
	}
	for(i=(autocorr_size/2); i<autocorr_size; i++){
		aux_vec[i - (autocorr_size/2)] = aux_vec_B1I[i];
	}
	//Convolution with box filter
	for(i=0; i<samples_filter; i++){
		filter[i] = 1.0/double(samples_filter);
	}
	convolution_gsl(aux_vec, autocorr_size, filter, samples_filter, aux_lambda);
	for(i=0; i<lambda_size; i++){
		//lambda's center is at (lambda_size/2)
		//aux_lambda's center is at (autocorr_size/2)
		//the sample resolution is corrected multiplying by (samples_CA_chip*GPS_CA_CHIP_RATE/sampling_rate)
		j = int(round(double(i - (lambda_size/2))*(double(samples_CA_chip)*GPS_CA_CHIP_RATE/sampling_rate))) + (autocorr_size/2);
		lambda[i] = aux_lambda[j]*aux_lambda[j]/(aux_lambda[autocorr_size/2]*aux_lambda[autocorr_size/2]); //Normalized power
	}
	return; 
}


bool check_input_antenna_pattern( double antenna_angles[], double antenna_pattern[], int pattern_length, double zero_level )
{
	int i;
	int min_length = 2;
	bool check_zero_level = false;
	double prev_angle;
	if(pattern_length < min_length){
		printf("ERROR! Input pattern too short (minimum = %d).\n", min_length);
		return false;
	}
	if(antenna_pattern[0] <= zero_level){
		check_zero_level = true;
	}
	if((antenna_angles[0] <= -180.0) || (antenna_angles[pattern_length - 1] >= 360.0)){
		printf("ERROR! Antenna-angles must lay (in degrees) within (-180,180], [0,360) or [0,180] (symmetric case).\n");
		return false;
	}
	if((antenna_angles[0] < 0.0) && (antenna_angles[pattern_length - 1] > 180.0)){
		printf("ERROR! Antenna-angles must lay (in degrees) within (-180,180], [0,360) or [0,180] (symmetric case).\n");
		return false;
	}
	prev_angle = antenna_angles[0];
	for(i=1; i<pattern_length; i++){
		if(antenna_pattern[i] <= zero_level){
			check_zero_level = true;
		}
		if(antenna_angles[i] <= prev_angle){
			printf("ERROR! Antenna-angles has to be a monotonic sequence.\n");
			return false;
		}
		prev_angle = antenna_angles[i];
	}
	if(check_zero_level){
		return true;
	}else{
		printf("ERROR! Pattern min-level has to be reached at least one time.\n");
		return false;
	}
}


void interpolate_antenna_pattern( double input_angles[], double input_pattern[], int pattern_length, double zero_level, double antenna_pattern[] )
{
	int i, j, first_zero_pos, last_zero_pos, offset_left, offset_right;
	int size_reordered_pattern, size_interval;
	double *interval_angles, *interval_pattern;
	double *reordered_angles, *reordered_pattern, *interm_angles, *interm_pattern;
	gsl_interp_accel *accel;
	gsl_spline *spline;
	if(!check_input_antenna_pattern(input_angles, input_pattern, pattern_length, zero_level)) return;
	i = 0;
	first_zero_pos = -1;
	while(first_zero_pos == -1){
		if(input_pattern[i] <= zero_level){
			first_zero_pos = i;
		}
		i ++;
	}
	if((input_angles[0] >= 0.0)&&(input_angles[pattern_length - 1] <= 180.0)){
		//==========SYMMETRIC CASE
		offset_left = 0;
		offset_right = 0;
		if(input_angles[0] == 0.0) offset_left = 1;
		if(input_angles[pattern_length - 1] == 180.0) offset_right = 1;
		size_reordered_pattern = 2*pattern_length + 1 - offset_left - offset_right;
		reordered_angles = (double*) malloc (size_reordered_pattern*sizeof(double));
		reordered_pattern = (double*) malloc (size_reordered_pattern*sizeof(double));
		//Construct intermediate pattern by applying symmetry around 180. From [0,180] to [0,360)
		interm_angles = (double*) malloc ((size_reordered_pattern - 1)*sizeof(double));
		interm_pattern = (double*) malloc ((size_reordered_pattern - 1)*sizeof(double));
		for(i=0; i<pattern_length; i++){
			interm_angles[i] = input_angles[i];
			interm_pattern[i] = input_pattern[i];
			j = size_reordered_pattern - 2 - i + offset_left;
			if(j<(size_reordered_pattern - 1)){
				interm_angles[j] = 360.0 - input_angles[i];
				interm_pattern[j] = input_pattern[i];
			}
		}
		//Reorder input vectors to start from first "zero" and finish at first "zero" + 360
		for(i=first_zero_pos; i<(size_reordered_pattern - 1); i++){
			reordered_angles[i - first_zero_pos] = interm_angles[i];
			reordered_pattern[i - first_zero_pos] = interm_pattern[i];
		}
		for(i=0; i<=first_zero_pos; i++){
			reordered_angles[i - first_zero_pos + size_reordered_pattern - 1] = interm_angles[i] + 360.0;
			reordered_pattern[i - first_zero_pos + size_reordered_pattern - 1] = interm_pattern[i];
		}
		free(interm_angles);
		free(interm_pattern);
	}else{
		//==========ASYMMETRIC CASE
		size_reordered_pattern = pattern_length + 1;
		reordered_angles = (double*) malloc (size_reordered_pattern*sizeof(double));
		reordered_pattern = (double*) malloc (size_reordered_pattern*sizeof(double));
		//Reorder input vectors to start from first "zero" and finish at first "zero" + 360
		for(i=first_zero_pos; i<pattern_length; i++){
			reordered_angles[i - first_zero_pos] = input_angles[i];
			reordered_pattern[i - first_zero_pos] = input_pattern[i];
		}
		for(i=0; i<=first_zero_pos; i++){
			reordered_angles[i - first_zero_pos + size_reordered_pattern - 1] = input_angles[i] + 360.0;
			reordered_pattern[i - first_zero_pos + size_reordered_pattern - 1] = input_pattern[i];
		}
	}
	last_zero_pos = 0;
	for(i=1; i<size_reordered_pattern; i++){
		if(reordered_pattern[i] <= zero_level){
			size_interval = i - last_zero_pos + 1;
			interval_angles = (double*) malloc (size_interval*sizeof(double));
			interval_pattern = (double*) malloc (size_interval*sizeof(double));
			for(j=0; j<size_interval; j++){
				interval_angles[j] = reordered_angles[j + last_zero_pos];
				interval_pattern[j] = reordered_pattern[j + last_zero_pos];
			}
			if(size_interval < 3){ //3 points needed for cubic spline
				for(j=ceil(reordered_angles[last_zero_pos]); j<=floor(reordered_angles[i]); j++){
					antenna_pattern[(j + 360)%360] = zero_level;
				}
			}else{
				accel = gsl_interp_accel_alloc();
				spline = gsl_spline_alloc(gsl_interp_cspline, size_interval);
				gsl_spline_init(spline, interval_angles, interval_pattern, size_interval);
				for(j=ceil(reordered_angles[last_zero_pos]); j<=floor(reordered_angles[i]); j++){
					antenna_pattern[(j + 360)%360] = max(zero_level, gsl_spline_eval(spline, double(j), accel));
				}
				gsl_spline_free(spline);
				gsl_interp_accel_free(accel);
			}
			free(interval_angles);
			free(interval_pattern);
			last_zero_pos = i;
		}
	}
	free(reordered_angles);
	free(reordered_pattern);
	return;
}


void convolution_gsl( double vec1[], int size_vec1, double vec2[], int size_vec2, double vec_out[] )
{
	int i;
	double *aux1, *aux2;
	aux1 = (double*) malloc (size_vec1*sizeof(double));
	aux2 = (double*) malloc (size_vec1*sizeof(double));
	for(i=0; i<size_vec1; i++){
		aux1[i] = vec1[i];
	}
	memset(aux2, 0, sizeof(double)*size_vec1);
	aux2[0] = vec2[0];
	for(i=1; i<((size_vec2 + 1)/2); i++){
		aux2[i] = vec2[i];
		aux2[size_vec1 - i] = vec2[size_vec2 - i];
	}
	gsl_fft_real_radix2_transform(aux1, 1, size_vec1);
	gsl_fft_real_radix2_transform(aux2, 1, size_vec1);
	vec_out[0] = aux1[0]*aux2[0];
	vec_out[size_vec1/2] = aux1[size_vec1/2]*aux2[size_vec1/2];
	for(i=1; i<(size_vec1/2); i++){
		vec_out[i] = aux1[i]*aux2[i] - aux1[size_vec1 - i]*aux2[size_vec1 - i];
		vec_out[size_vec1 - i] = aux1[size_vec1 - i]*aux2[i] + aux1[i]*aux2[size_vec1 - i];
	}
	gsl_fft_halfcomplex_radix2_inverse(vec_out, 1, size_vec1);
	free(aux1);
	free(aux2);
	return;
}


void correlation_gsl( double vec1[], double vec2[], double vec_out[], int size_vec )
{
	int i;
	double *aux1, *aux2;
	aux1 = (double*) malloc (size_vec*sizeof(double));
	aux2 = (double*) malloc (size_vec*sizeof(double));
	for(i=0; i<size_vec; i++){
		aux1[i] = vec1[i];
		aux2[i] = vec2[i];
	}
	gsl_fft_real_radix2_transform(aux1, 1, size_vec);
	gsl_fft_real_radix2_transform(aux2, 1, size_vec);
	vec_out[0] = aux1[0]*aux2[0];
	vec_out[size_vec/2] = aux1[size_vec/2]*aux2[size_vec/2];
	for(i=1; i<(size_vec/2); i++){
		vec_out[i] = aux1[i]*aux2[i] + aux1[size_vec - i]*aux2[size_vec - i];
		vec_out[size_vec - i] = aux1[size_vec - i]*aux2[i] - aux1[i]*aux2[size_vec - i];
	}
	gsl_fft_halfcomplex_radix2_inverse(vec_out, 1, size_vec);
	free(aux1);
	free(aux2);
	return;
}


void Compute_power_spectrum( int samples, double* real_part, double* imag_part, double* power_spectrum )
{
	int i;
	bool max_val_found;
	double max_val;
	double *data_fft, *power_fft;
	data_fft = (double*) malloc (2*samples*sizeof(double));
	power_fft = (double*) malloc (samples*sizeof(double));
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc(samples);
	workspace = gsl_fft_complex_workspace_alloc(samples);
	for(i=0; i<samples; i++){
		data_fft[2*i] = double(real_part[i]);
		data_fft[2*i + 1] = double(imag_part[i]);
	}
	gsl_fft_complex_forward(data_fft, 1, samples, wavetable, workspace);
	gsl_fft_complex_wavetable_free(wavetable);
	gsl_fft_complex_workspace_free(workspace);
	max_val = 0.0;
	max_val_found = false;
	for(i=0; i<samples; i++){
		power_fft[i] = sqrt(data_fft[2*i]*data_fft[2*i] + data_fft[2*i + 1]*data_fft[2*i + 1]);
		if(power_fft[i] > max_val){
			max_val = power_fft[i];
			max_val_found = true;
		}
	}
	if(!max_val_found){
		max_val = 1.0;
	}
	for(i=0; i<(samples/2 - 1); i++){
		power_spectrum[i] = power_fft[i + samples/2 + 1]/max_val;
	}
	for(i=(samples/2 - 1); i<samples; i++){
		power_spectrum[i] = power_fft[i - samples/2 + 1]/max_val;
	}
	free(data_fft);
	free(power_fft);
	return;
}


double Compute_Doppler_bandwidth( int samples, double* doppler_freqs, double* doppler_map, double pos_pow_Max[2] )
{
	int i, pos_max;
	double max_val, x_interp, delta_freq, ref_val, left_freq, right_freq, coef_left, coef_right;
	double *x_segment, *y_segment;
	if(samples < 22){ //11 samples needed for 10-order polynomial fit
		cout << "ERROR! Not enough samples at Compute_Doppler_bandwidth" << endl;
		return 0.0;
	}
	x_segment = (double*) malloc ((samples/2)*sizeof(double));
	y_segment = (double*) malloc ((samples/2)*sizeof(double));
	gsl_interp *workspace;
	gsl_interp_accel *accel;
	for(i=0; i<(samples/2); i++){
		x_segment[i] = doppler_freqs[i + (samples/4)];
		y_segment[i] = doppler_map[i + (samples/4)];
	}
	gsl_interp_init(workspace, x_segment, y_segment, 10);
	max_val = 0.0;
	for(i=0; i<(samples/2); i++){
		if(gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i], accel) > max_val){
			max_val = gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i], accel);
			pos_max = i;
		}
	}
	if((pos_max == 0)||(pos_max == ((samples/2)-1))){
		free(x_segment);
		free(y_segment);
		gsl_interp_free(workspace);
		gsl_interp_accel_free(accel);
		return 0.0;
	}
	pos_pow_Max[1] = 0.0;
	delta_freq = x_segment[pos_max] - x_segment[pos_max - 1];
	for(i=0; i<=samples; i++){
		x_interp = x_segment[pos_max] - (delta_freq/2.0) + delta_freq*double(i)/double(samples);
		if(gsl_interp_eval(workspace, x_segment, y_segment, x_interp, accel) > pos_pow_Max[1]){
			pos_pow_Max[1] = gsl_interp_eval(workspace, x_segment, y_segment, x_interp, accel);
			pos_pow_Max[0] = x_interp;
		}
	}
	ref_val = pos_pow_Max[1]/2.0;
	if((gsl_interp_eval(workspace, x_segment, y_segment, x_segment[0], accel) > ref_val)||(gsl_interp_eval(workspace, x_segment, y_segment, x_segment[(samples/2)-1], accel) > ref_val)){
		free(x_segment);
		free(y_segment);
		gsl_interp_free(workspace);
		gsl_interp_accel_free(accel);
		return 0.0;
	}
	i = pos_max - 1;
	while((gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i], accel) > ref_val)&&(i>=0)) i--;
	coef_left =  (gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i + 1], accel) - ref_val)/(gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i + 1], accel) - gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i], accel));
	coef_right = (ref_val - gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i], accel))/(gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i + 1], accel) - gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i], accel));
	left_freq = x_segment[i]*coef_left + x_segment[i + 1]*coef_right;
	i = pos_max + 1;
	while((gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i], accel) > ref_val)&&(i<(samples/2))) i++;
	coef_left =  (gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i], accel) - ref_val)/(gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i], accel) - gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i - 1], accel));
	coef_right = (ref_val - gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i - 1], accel))/(gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i], accel) - gsl_interp_eval(workspace, x_segment, y_segment, x_segment[i - 1], accel));
	right_freq = x_segment[i - 1]*coef_left + x_segment[i]*coef_right;
	free(x_segment);
	free(y_segment);
	gsl_interp_free(workspace);
	gsl_interp_accel_free(accel);
	return (right_freq - left_freq);
}


void retrack_real_waveform( double* waveform, int wav_length, double delay )
{
	int i;
	double real_part, imag_part, cos_del, sin_del;
	//GSL FFT variables
	gsl_fft_real_wavetable * real;
	gsl_fft_halfcomplex_wavetable * hc;
	gsl_fft_real_workspace * work;
	real = gsl_fft_real_wavetable_alloc(wav_length);
	hc = gsl_fft_halfcomplex_wavetable_alloc(wav_length);
	work = gsl_fft_real_workspace_alloc(wav_length);
	//FFT transform of sequence of real numbers (power waveform)
	gsl_fft_real_transform(waveform, 1, wav_length, real, work);
	gsl_fft_real_wavetable_free(real);
	//Delay shift by means of *e^(-i2PIf*delay) in the frequency domain
	for(i=1; i<=((wav_length - 1)/2); i++){
		cos_del = cos(2.0*PI_NUM*delay*double(i)/double(wav_length));
		sin_del = sin(2.0*PI_NUM*delay*double(i)/double(wav_length));
		real_part = waveform[2*i - 1]*cos_del + waveform[2*i]*sin_del;
		imag_part = waveform[2*i]*cos_del - waveform[2*i - 1]*sin_del;
		waveform[2*i - 1] = real_part;
		waveform[2*i] = imag_part;
	}
	//Inverse FFT
	gsl_fft_halfcomplex_inverse(waveform, 1, wav_length, hc, work);
	gsl_fft_halfcomplex_wavetable_free(hc);
	gsl_fft_real_workspace_free(work);
	return;
}

