BUILT_SOURCES = wavpy_wrap.cpp
SWIG_SOURCES = wavpy.i

pkgpython_PYTHON = wavpy.py
pkgpyexec_LTLIBRARIES = _wavpy.la
_wavpy_la_SOURCES = wavpy_wrap.cpp $(SWIG_SOURCES)
_wavpy_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS)
_wavpy_la_LDFLAGS = -module
_wavpy_la_LIBADD = ../src/libwaveform.la

wavpy_wrap.cpp : $(SWIG_SOURCES)
	$(SWIG) $(SWIG_PYTHON_OPT) -o $@ $<