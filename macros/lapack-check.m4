dnl  Check for LAPACK
dnl  Copyright 2011 Brice Boyer <bboyer@imag.fr>
dnl  This file is part of Fflas-Fpack
dnl  See COPYING for licence information.

dnl **********************************
dnl *              TODO              *
dnl **********************************
dnl no support yet to MKL
dnl AS_IF([test -r "$BLAS_VAL/include/mkl_cblas.h"],
dnl [ BLAS_LIBS="-L${BLAS_VAL}/lib/${MKL_ARCH}/ -lmkl_lapack64 -lmkl -lvml -lguide" ])
dnl **********************************

AC_DEFUN([FF_CHECK_LAPACK], [

		BACKUP_CXXFLAGS=${CXXFLAGS}
		BACKUP_LIBS=${LIBS}


		AC_ARG_WITH(lapack,
			[AC_HELP_STRING([--with-lapack=<blas|path>],
				[Use LAPACK functions. This library is mandatory for LinBox
				compilation. If argument is <empty> that means
				the library is reachable with the standard search path
				(/usr or /usr/local). Or, you can give the <path> to
				the directory which contains the library. If the argument
				is 'blas', then we look in the BLAS vendor library.
				We look for a C interface (clapack_), and if not present,
				look for standard functions (as dgetrf_). First one available
				in order in '$path /usr /usr/local', first chosen, even if it is
                not clapack_ (example: clapack_ in /usr but dgetrf_ in $path : dgetrf_ chosen,
					$path not even looked into).
				])
			])


		AC_MSG_CHECKING(for LAPACK)


		AS_IF([ test "$with_lapack" = "blas"], [
				dnl  echo "vendor ${BLAS_VENDOR} in ${BLAS_PATH}"
			dnl check for lapack function in vendor lib
			AS_CASE([${BLAS_VENDOR}],
				["ATLAS"],[
				dnl atlas provides a liblapack next to its libcblas
				LAPACK_LIBS="-llapack"
				dnl why would we need lapack_atlas when llapack is enough ?
				dnl could llapack not provide the symbols ?
				AS_IF([test -r "${BLAS_PATH}/liblapack_atlas.a" -o -r "${BLAS_PATH}/liblapack_atlas.so"],
					[LAPACK_LIBS="${LAPACK_LIBS} -llapack_atlas"])
				dnl  AS_IF([ test "x$BLAS_PATH" != "x/usr/lib" -a "x$BLAS_PATH" != "x/usr/local/lib"],
					dnl  [LAPACK_LIBS="-L${BLAS_PATH} ${BLAS_LIBS}"])

				],
				dnl GSL provides no lapack ! why would you use GSL ?
				["GSL"],
				[LAPACK_LIBS=""],
				dnl lapack is in libgoto2
				["GOTO2"],
				[LAPACK_LIBS=""],
				dnl maybe lapack is in libblas ?
				["OTHER"],
				[LAPACK_LIBS=""],
				dnl defaulting somewhere...
				[LAPACK_LIBS=""])

			CXXFLAGS="${BACKUP_CXXFLAGS} ${CBLAS_FLAG} "
			LIBS="${BACKUP_LIBS} ${BLAS_LIBS} ${LAPACK_LIBS}"

			dnl  echo ${LAPACK_LIBS}

			AC_TRY_RUN(
					[#define __FFLAFLAS_CONFIGURATION
					#define __FFLAFLAS_HAVE_LAPACK 1
					#define __FFLAFLAS_HAVE_CLAPACK 1
					#include "fflas-ffpack/config-blas.h"
					int main () {  double a[4] = {1.,2.,3.,4.};
					int ipiv[2];
					clapack_dgetrf(CblasRowMajor, 2, 2, a, 2, ipiv);
					if ( (a[0]!=2.) && (a[1]!=0.5) && (a[2]!=4.) && (a[3]!=1.))
					return -1;
					else
					return 0;
					} ],
					[ dgetrf_found="yes" ],
					[ dgetrf_problem="problem" ],
					[ dgetrf_found="" ])

				AS_IF( [test "${dgetrf_found}" = "yes"],
						[ AC_SUBST(LAPACK_LIBS)
						AC_MSG_RESULT( yes (clapack))
						AC_DEFINE(HAVE_LAPACK,1,[Define if LAPACK is installed])
						AC_DEFINE(HAVE_CLAPACK,1,[Define if C interface to LAPACK is available])
						],
                        [dnl not found : trying only lapack
						CXXFLAGS="${BACKUP_CXXFLAGS} ${CBLAS_FLAG} "
						LIBS="${BACKUP_LIBS} ${BLAS_LIBS} ${LAPACK_LIBS}"


						AC_TRY_RUN(
							[#define __FFLAFLAS_CONFIGURATION
							#define __FFLAFLAS_HAVE_LAPACK 1
							//#define __FFLAFLAS_HAVE_CLAPACK 1
							#include "fflas-ffpack/config-blas.h"
							int main () {  double a[4] = {1.,2.,3.,4.};
							int ipiv[2];
							clapack_dgetrf(CblasRowMajor, 2, 2, a, 2, ipiv);
							if ( (a[0]!=2.) && (a[1]!=0.5) && (a[2]!=4.) && (a[3]!=1.))
							return -1;
							else
							return 0;
							} ],
							[ dgetrf_found="yes" ],
							[ dgetrf_problem="problem" ],
							[ dgetrf_found="" ])
						AS_IF([test "${dgetrf_found}" = "yes"],
								[ AC_SUBST(LAPACK_LIBS)
								AC_MSG_RESULT( yes (lapack))
								AC_DEFINE(HAVE_LAPACK,1,[Define if LAPACK is installed])
								],
								[ AC_MSG_RESULT(no) ])
						])
			],[ dnl not BLAS vendor asked, so looking in DEFAULT_CHECKING_PATH
			dnl  echo "path"

			AS_IF([test "x$BLAS_VENDOR" = "xUSER"], [ dnl this is temporary  -- because the user supplies everything in --with-blas.
					CXXFLAGS="${BACKUP_CXXFLAGS} ${CBLAS_FLAG} "
					LIBS="${BACKUP_LIBS} ${BLAS_LIBS}"


					AC_TRY_RUN(
						[#define __FFLAFLAS_CONFIGURATION
						#define __FFLAFLAS_HAVE_LAPACK 1
						#define __FFLAFLAS_HAVE_CLAPACK 1
						#include "fflas-ffpack/config-blas.h"
						int main () {  double a[4] = {1.,2.,3.,4.};
						int ipiv[2];
						clapack_dgetrf(CblasRowMajor, 2, 2, a, 2, ipiv);
						if ( (a[0]!=2.) && (a[1]!=0.5) && (a[2]!=4.) && (a[3]!=1.))
						return -1;
						else
						return 0;
						} ],
						[ dgetrf_found="yes"
						dnl  echo "yes"
						],
						[ dgetrf_problem="problem"
						dnl  echo "no"
						],
						[ ])

					AS_IF([ test "${dgetrf_found}" = "yes"],
							[	 AC_SUBST(LAPACK_LIBS)
							AC_MSG_RESULT( yes (clapack))
							AC_DEFINE(HAVE_LAPACK,1,[Define if LAPACK is installed])
							AC_DEFINE(HAVE_CLAPACK,1,[Define if C interface to LAPACK is available])
							], dnl clapack not found. looking for lapack
							[

							AC_TRY_RUN(
								[#define __FFLAFLAS_CONFIGURATION
								#define __FFLAFLAS_HAVE_LAPACK 1
								//#define __FFLAFLAS_HAVE_CLAPACK 1
								#include "fflas-ffpack/config-blas.h"
								int main () {  double a[4] = {1.,2.,3.,4.};
								int ipiv[2];
								clapack_dgetrf(CblasRowMajor, 2, 2, a, 2, ipiv);
								if ( (a[0]!=2.) && (a[1]!=0.5) && (a[2]!=4.) && (a[3]!=1.))
								return -1;
								else
								return 0;
								} ],
								[ dgetrf_found="yes"
								 ],
								[ dgetrf_problem="$problem"
								],
								[  ])

							AS_IF([ test "x${dgetrf_found}" = "xyes"],
									[	 AC_SUBST(LAPACK_LIBS)
									AC_MSG_RESULT( yes (lapack))
									AC_DEFINE(HAVE_LAPACK,1,[Define if LAPACK is installed])
									], dnl clapack not found. looking for lapack
									[
									AC_MSG_RESULT( no )
									])
							])

							],[

			LAPACK_HOME_PATH="$with_lapack ${DEFAULT_CHECKING_PATH}"
			for LAPACK_HOME in ${LAPACK_HOME_PATH} ; do
				dnl  echo "in ${LAPACK_HOME} for clapack"
				AS_IF(
						[test -r "$LAPACK_HOME/lib/liblapack.a" -o -r "$LAPACK_HOME/lib/liblapack.so"  ],
						[LAPACK_LIBS="-llapack"
						LAPACK_PATH="${LAPACK_HOME}/lib"
						AS_IF([ test "x$LAPACK_HOME" != "x/usr" -a "x$LAPACK_HOME" != "x/usr/local"],
							[LAPACK_LIBS="-L${LAPACK_HOME}/lib  -llapack"])
						],
						[test -r "$LAPACK_HOME/liblapack.a" -o -r "$LAPACK_HOME/liblapack.so" ],
						[ LAPACK_LIBS="-llapack"
						LAPACK_PATH="${LAPACK_HOME}"
						AS_IF([ test "x$LAPACK_HOME" != "x/usr" -a "x$LAPACK_HOME" != "x/usr/local"],
							[LAPACK_LIBS="-L${LAPACK_HOME}  -llapack"])
						]
					 )
				dnl  echo "lapack libs : $LAPACK_LIBS"
				CXXFLAGS="${BACKUP_CXXFLAGS} ${CBLAS_FLAG} "
				LIBS="${BACKUP_LIBS} ${BLAS_LIBS} ${LAPACK_LIBS}"


				AC_TRY_RUN(
						[#define __FFLAFLAS_CONFIGURATION
						#define __FFLAFLAS_HAVE_LAPACK 1
						#define __FFLAFLAS_HAVE_CLAPACK 1
						#include "fflas-ffpack/config-blas.h"
						int main () {  double a[4] = {1.,2.,3.,4.};
						int ipiv[2];
						clapack_dgetrf(CblasRowMajor, 2, 2, a, 2, ipiv);
						if ( (a[0]!=2.) && (a[1]!=0.5) && (a[2]!=4.) && (a[3]!=1.))
						return -1;
						else
						return 0;
						} ],
						[ dgetrf_found="yes"
						dnl  echo "yes"
						break ],
						[ dgetrf_problem="problem"
						unset LAPACK_LIBS
						dnl  echo "no" ],
						[ break ])
			done ;
			AS_IF([ test "${dgetrf_found}" = "yes"],
					[	 AC_SUBST(LAPACK_LIBS)
					AC_MSG_RESULT( yes (clapack))
					AC_DEFINE(HAVE_LAPACK,1,[Define if LAPACK is installed])
					AC_DEFINE(HAVE_CLAPACK,1,[Define if C interface to LAPACK is available])
					], dnl clapack not found. looking for lapack
					[
					for LAPACK_HOME in ${LAPACK_HOME_PATH} ; do
						dnl  echo "in ${LAPACK_HOME}"
						AS_IF(
							[test -r "$LAPACK_HOME/lib/liblapack.a" -o -r "$LAPACK_HOME/lib/liblapack.so"  ],
							[LAPACK_LIBS="-llapack"
							LAPACK_PATH="${LAPACK_HOME}/lib"
							AS_IF([ test "x$LAPACK_HOME" != "x/usr" -a "x$LAPACK_HOME" != "x/usr/local"],
								[LAPACK_LIBS="-L${LAPACK_HOME}/lib  -llapack"])
							],
							[test -r "$LAPACK_HOME/liblapack.a" -o -r "$LAPACK_HOME/liblapack.so" ],
							[ LAPACK_LIBS="-llapack"
							LAPACK_PATH="${LAPACK_HOME}"
							AS_IF([ test "x$LAPACK_HOME" != "x/usr" -a "x$LAPACK_HOME" != "x/usr/local"],
								[LAPACK_LIBS="-L${LAPACK_HOME}  -llapack"])
							]
						 )
						CXXFLAGS="${BACKUP_CXXFLAGS} ${CBLAS_FLAG} "
						LIBS="${BACKUP_LIBS} ${BLAS_LIBS} ${LAPACK_LIBS}"


						AC_TRY_RUN(
								[#define __FFLAFLAS_CONFIGURATION
								#define __FFLAFLAS_HAVE_LAPACK 1
								//#define __FFLAFLAS_HAVE_CLAPACK 1
								#include "fflas-ffpack/config-blas.h"
								int main () {  double a[4] = {1.,2.,3.,4.};
								int ipiv[2];
								clapack_dgetrf(CblasRowMajor, 2, 2, a, 2, ipiv);
								if ( (a[0]!=2.) && (a[1]!=0.5) && (a[2]!=4.) && (a[3]!=1.))
								return -1;
								else
								return 0;
								} ],
								[ dgetrf_found="yes"
								break ],
								[ dgetrf_problem="$problem"
								unset LAPACK_LIBS
								],
								[ break ])
					done ;
					AS_IF([ test "${dgetrf_found}" = "yes"],
							[	 AC_SUBST(LAPACK_LIBS)
							AC_MSG_RESULT( yes (lapack))
							AC_DEFINE(HAVE_LAPACK,1,[Define if LAPACK is installed])
							], dnl clapack not found. looking for lapack
							[
							AC_MSG_RESULT( no )
							])
					])
					])
					])

					dnl  AM_CONDITIONAL(FFLAFFLAS_HAVE_LAPACK, test "x$HAVE_LAPACK" = "xyes")

					CXXFLAGS=${BACKUP_CXXFLAGS}
					LIBS=${BACKUP_LIBS}
					dnl  unset LD_LIBRARY_PATH


					])

