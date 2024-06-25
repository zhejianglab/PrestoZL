AC_INIT([tempo], [13.103])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])


AC_PROG_CC

# Copied the list of fortran compilers from autoconf, and moved gfortran
# to the front so that it is preferred if available.  Would be nice if 
# there were a less kludgey way to do this..
AC_PROG_F77([gfortran g77 xlf f77 frt pgf77 cf77 fort77 fl32 af77 \
	     xlf90 f90 pgf90 pghpf epcf90 g95 xlf95 f95 \
	     fort ifort ifc efc pgfortran pgf95 lf95 ftn nagfor])
case "$F77" in 
    *gfortran*) FFLAGS+=" -ffixed-form -ffixed-line-length-none -std=legacy"
     ;;
esac

AC_C_BIGENDIAN([BIGENDIAN=".true."], [BIGENDIAN=".false."])
AC_SUBST(BIGENDIAN)

AC_CHECK_HEADERS([malloc.h])

AX_LAPACK
AM_CONDITIONAL([HAVE_LAPACK],[test x"$ax_lapack_ok" = xyes])

AC_CONFIG_FILES([
	Makefile
	src/Makefile
	tempo.cfg
	src/config.h
	clock/Makefile
	ephem/Makefile
	test/Makefile
	tzpar/Makefile
	util/print_resid/Makefile
])
AC_OUTPUT

# Misc messages
echo 
echo "Tempo is ready to be compiled.  Please run 'make'."
echo 
echo "'make install' will copy 'tempo' to $prefix/bin."
echo
