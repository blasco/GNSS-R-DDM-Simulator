
AM_CPPFLAGS= -I/usr/local/gsl-2.1/include -DLOCAL_DATA_PATH=\"${PWD}/data_share/\" -fPIC

# create static and dynamic libraries
lib_LTLIBRARIES = libwaveform.la
libwaveform_la_SOURCES = reflecting_surface.cc \
                 Surface_routines.f90 \
                 ancillary_functions.cc \
                 compute_delays_wav.cc \
                 specular_geometry.cc \
                 gnss_composite.cc \
                 rf_front_end.cc \
                 waveform_complex.cc \
                 waveform_power.cc \
                 ZavorotnyVoronovichModel.cc \
                 MultiRaySingleRefl_model.cc

libwaveform_la_LIBADD =  -L/usr/local/gsl-2.1/lib -lgsl -lgslcblas -lm -lfftw3 -lgfortran -lstdc++
