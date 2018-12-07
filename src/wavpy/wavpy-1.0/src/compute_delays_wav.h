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

void compute_delays_wav( double* waveform_in, int wav_length, double sampling_rate, double relFactor, double &positionMax, double &posSampleMax, double &sigma_posMax, double &power_posMax, double &positionDer, double &posSampleDer, double &sigma_posDer, double &power_posDer, double &powerDer_posDer, double &positionRel, double &posSampleRel, double &sigma_posRel, double &noiseFloor, double &slope_NormTail, double &sigma_slope_NormTail, double normTail_length, double min_resolution_fft_interp, double fit_length, bool apply_speckle_weights, int num_incoh, bool apply_limits, double limits_center, double limits_width, double apriori_scattdel );