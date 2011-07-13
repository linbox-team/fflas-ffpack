dnl  Check for BLAS
dnl  Copyright 2011 Brice Boyer <bboyer@imag.fr>
dnl  This file is part of Fflas-Fpack (and comes from LinBox)
dnl  See COPYING for licence information.


dnl FF_CHECK_BLAS ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for BLAS and define BLAS_LIBS

AC_DEFUN([FF_CHECK_GSL],
		[
		AC_ARG_WITH(gsl,
			[AC_HELP_STRING([--with-gsl=<path|yes>],
				[Use GSL blas library. BLAS are mandatory for FFLAS-FFPACK
				compilation. If argument is  <yes> that means
				the library is reachable with the standard search path
				(/usr or /usr/local). Otherwise you give the <path> to
				the directory which contains the library. If empty, GSL is not
				searched for.  ])
			])
		dnl  echo $with_gsl
		dnl  echo $withval

		AS_IF([ test -n "$with_gsl" ],
			[ BLAS_HOME_PATH="${DEFAULT_CHECKING_PATH}"
			AS_IF([ test "$with_gsl" != "yes" ],
				[ BLAS_HOME_PATH="$with_gsl ${DEFAULT_CHECKING_PATH}" ])

			BACKUP_CXXFLAGS=${CXXFLAGS}
			BACKUP_LIBS=${LIBS}

			AC_MSG_CHECKING(for C interface to BLAS with -lgsl -lgslcblas)


			for BLAS_HOME in ${BLAS_HOME_PATH} ; do
				CBLAS="yes"
				CBLAS_FLAG="-D__FFLASFFPACK_HAVE_CBLAS"

				AS_IF([ test -r "$BLAS_HOME/lib/libgsl.a" -o -r "$BLAS_HOME/lib/libgsl.so"  ],
				[BLAS_LIBS="-lgsl -lgslcblas -lm"
				BLAS_PATH="${BLAS_HOME}/lib"
				AS_IF([ test "x$BLAS_HOME" != "x/usr" -a "x$BLAS_HOME" != "x/usr/local"],
					[BLAS_LIBS="-L${BLAS_HOME}/lib ${BLAS_LIBS}"])
				],
				[test -r "$BLAS_HOME/libgsl.a" -o -r "$BLAS_HOME/libgsl.so" ],
				[ BLAS_LIBS="-lgsl -lgslcblas -lm"
				BLAS_PATH="${BLAS_HOME}"
				AS_IF([ test "x$BLAS_HOME" != "x/usr" -a "x$BLAS_HOME" != "x/usr/local"],
					[BLAS_LIBS="-L${BLAS_HOME} ${BLAS_LIBS}"])
				])

				CXXFLAGS="${BACKUP_CXXFLAGS} ${CBLAS_FLAG}"
				LIBS="${BACKUP_LIBS} ${BLAS_LIBS}"

				AC_TRY_LINK(
						[#define __FFLASFFPACK_CONFIGURATION
						#include "fflas-ffpack/config-blas.h"],
						[double a;],
						[
						AC_TRY_RUN(
							[#define __FFLASFFPACK_CONFIGURATION
							#include "fflas-ffpack/config-blas.h"
							int main ()
							{  double a[4] = {1.,2.,3.,4.}; double b[4]= {4.,3.,2.,1.}; double c[4];
							cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans,2,2,2,1., a,2,b,2,0.,c,2);
							if ( (c[0]!=8.) && (c[1]!=5.) && (c[2]!=20.) && (c[3]!=13))
								return -1;
							else
								return 0;
							} ],
							[ blas_found="yes"
							break ],
							[ blas_problem="$problem $BLAS_HOME"
							unset BLAS_LIBS ],
							[ blas_found="yes"
							blas_cross="yes"
							break ])
						],
						[
						blas_found="no"
						blas_checked="$checked $BLAS_HOME"
						unset BLAS_LIBS ]
				)
			done

	AS_IF([ test "x$blas_found" = "xyes" ],[
			BLAS_VENDOR="GSL"
			AC_SUBST(BLAS_VENDOR)
			AC_SUBST(BLAS_LIBS)
			AC_SUBST(CBLAS_FLAG)
			AC_SUBST(BLAS_PATH)
			AC_DEFINE(HAVE_BLAS,1,[Define if BLAS is installed])
			AC_DEFINE(HAVE_CBLAS,1,[Define if C interface to BLAS is available])
			dnl  AC_DEFINE(BLAS_AVAILABLE,,[Define if BLAS routines are available])
			BLAS_FOUND=true
			AC_SUBST(BLAS_FOUND)
			HAVE_BLAS=yes
			AS_IF([ test "x$blas_cross" != "xyes" ], [
				AC_MSG_RESULT(found)],
				[AC_MSG_RESULT(unknown)
				echo "WARNING: You appear to be cross compiling, so there is no way to determine"
				echo "whether your BLAS are good. I am assuming it is."
				])
			ifelse([$2], , :, [$2])
			],
			[ test -n "$blas_problem" ],
			[ AC_MSG_RESULT(not working) ],
			dnl  echo "Sorry, your BLAS are not working. Disabling."
			[test "x$blas_found" = "xno" ],
			[AC_MSG_RESULT(not found)]
	)


	CXXFLAGS=${BACKUP_CXXFLAGS}
	LIBS=${BACKUP_LIBS}
	])

])


