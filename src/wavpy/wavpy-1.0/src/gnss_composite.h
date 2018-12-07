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

class GNSS_composite
{
	double *lambda_func;
	double sampling_rate;
	double filter_BB_BW;
	bool lambda_updated;
  public:
	int lambda_size;
	float weight_CA;
	float weight_PY;
	float weight_M;
	float weight_IM;
	float weight_E1A;
	float weight_E1B;
	float weight_E1C;
	float weight_B1I;
	float weight_L1C;
	double frequency;
	GNSS_composite( void );
	void dump_parameters( void );
	void set_instrumental_params( double input_sampling_rate, double input_filter_BW, char computeLambda );
	void compute_lambda_func( void );
	void get_lambda_func( double* range_vec, int range_length, double* lambda_out, int lambda_length );
};