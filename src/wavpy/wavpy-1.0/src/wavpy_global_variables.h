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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <fstream>

// Global variables (units in MKS if not otherwise specified)
const double C_LIGHT = 299792458;
const double PI_NUM = 3.141592653589793238462643383279502884197;
const double FREQ_GPS_L1 = 1575420000;
const double FREQ_GPS_L5 = 1176450000;
const double FREQ_BEIDOU_B1 = 1561098000;
const double A_EARTH_SEMIAXIS_KM = 6378.137;
const double B_EARTH_SEMIAXIS_KM = 6378.137;
const double C_EARTH_SEMIAXIS_KM = 6356.7523142;
const double GPS_CA_CHIP_RATE = 1023000.0;
const double POW_CA_Trans_dBW = 14.3;
const double POW_L1C_Trans_dBW = 15.8;
const double POW_PY_Trans_dBW = 13.05; //11.6;
const double POW_M_Trans_dBW = 15.27; //15.6;
const double POW_IM_Trans_dBW = 14.3;
const double POW_B1I_RecEarth_dBW = -163.0;
const double POW_CA_RecEarth_dBW = -158.5;
const double POW_L1C_RecEarth_dBW = -157.0;
const double POW_E1A_RecEarth_dBW = -157.0;
const double POW_E1B_RecEarth_dBW = -160.0;
const double POW_E1C_RecEarth_dBW = -160.0;
const double K_BOLTZMANN = 1.3806488e-23;

