# Prelude
AC_INIT([waveform_pylib],[1.0])
AM_INIT_AUTOMAKE([subdir-objects])

# Checks programs
AC_PROG_CC
AC_PROG_FC([gfortran])
AC_PROG_CXX
AC_PROG_LIBTOOL

AC_CONFIG_MACRO_DIR([m4])

AC_DISABLE_STATIC
AM_PATH_PYTHON(2.7)
AC_PROG_SWIG
SWIG_ENABLE_CXX
SWIG_PYTHON

# Checks libraries
# Checks header files
AC_CONFIG_HEADERS([config.h])
# Checks typedefs, structures, and compiler characteristics
# Checks library functions

# Output files
AC_PROG_INSTALL
AC_CONFIG_FILES([Makefile src/Makefile wavpy/Makefile])
AC_OUTPUT

