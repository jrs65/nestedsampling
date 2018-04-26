dnl Nested sampler.

AC_INIT(nested, 0.1)
AC_CONFIG_SRCDIR(src/test1d.F90)
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE

dnl Set to Fortran 95
AC_PROG_FC([], 95)

dnl Look for lapack
AC_LANG(Fortran 77)
AC_CHECK_LIB(lapack, dsyev)

AC_PROG_LIBTOOL

dnl Output Makefile
AC_OUTPUT(Makefile lib/Makefile src/Makefile)