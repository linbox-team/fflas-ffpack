dnl  Check for BLAS
dnl  Copyright Pascal Giorgi 2005
dnl  Modified Brice Boyer 2011

dnl
dnl ========LICENCE========
dnl This file is part of the library FFLAS-FFPACK.
dnl
dnl FFLAS-FFPACK is free software: you can redistribute it and/or modify
dnl it under the terms of the  GNU Lesser General Public
dnl License as published by the Free Software Foundation; either
dnl version 2.1 of the License, or (at your option) any later version.
dnl
dnl This library is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl Lesser General Public License for more details.
dnl
dnl You should have received a copy of the GNU Lesser General Public
dnl License along with this library; if not, write to the Free Software
dnl Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
dnl ========LICENCE========
dnl/



dnl FF_CHECK_CBLAS ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for C interface to BLAS and define BLAS_LIBS

AC_DEFUN([FF_CHECK_CBLAS],
		[ AC_ARG_WITH(cblas,
			[AC_HELP_STRING([--with-cblas=<lib>], [Use BLAS library. This library is mandatory for FFLAS-FFPACK
				compilation. If argument is <empty> that means
				the library is reachable with the standard search path
				(/usr or /usr/local). Otherwise you give the <path> to
				the directory which contains the library.  ])
			])

		BLAS_HOME_PATH="$with_cblas ${DEFAULT_CHECKING_PATH}"

		dnl Check for existence

		BACKUP_CXXFLAGS=${CXXFLAGS}
		BACKUP_LIBS=${LIBS}

		AC_MSG_CHECKING(for C interface to BLAS with -lcblas)

		dnl **************************************
		dnl Check first for C interface to BLAS
		dnl **************************************


		for BLAS_HOME in ${BLAS_HOME_PATH} ; do
			dnl  echo looking in ${BLAS_HOME}
			CBLAS="yes"
			CBLAS_FLAG="-D__FFLASFFPACK_HAVE_CBLAS"
			ATLAS_LIBS="-lcblas"
			AS_IF(
					dnl obscure
					[ test -r "/System/Library/Frameworks/Accelerate.framework" ],
					[BLAS_LIBS="-Wl,-framework -Wl,Accelerate"],
					dnl lib/libcblas.* ?
					[ test -r "$BLAS_HOME/lib/libcblas.a" -o -r "$BLAS_HOME/lib/libcblas.so" ],
					[ ATLAS_NEEDED=`nm  -u $BLAS_HOME/lib/libcblas.a  | grep ATL`
					ATLAS_NEEDED2=`nm -Du $BLAS_HOME/lib/libcblas.so | grep ATL`
					AS_IF( [test -n "$ATLAS_NEEDED" -o -n "$ATLAS_NEEDED2"],
						[ATLAS_LIBS=" ${ATLAS_LIBS} -latlas"])

					BLAS_LIBS=" ${ATLAS_LIBS}"
					BLAS_PATH="${BLAS_HOME}/lib"

					AS_IF([  test "x$BLAS_HOME" != "x/usr" -a "x$BLAS_HOME" != "x/usr/local"],
						[BLAS_LIBS="-L${BLAS_HOME}/lib ${ATLAS_LIBS}"])
					],
					dnl libcblas.* ?
					[ test -r "$BLAS_HOME/libcblas.a" -o -r "$BLAS_HOME/libcblas.so" ],
					[ ATLAS_NEEDED=`nm  -u $BLAS_HOME/libcblas.a  | grep ATL`
					ATLAS_NEEDED2=`nm -Du $BLAS_HOME/libcblas.so | grep ATL`
					AS_IF( [test -n "$ATLAS_NEEDED" -o -n "$ATLAS_NEEDED2"],
						[ATLAS_LIBS=" ${ATLAS_LIBS} -latlas"])

					BLAS_LIBS=" ${ATLAS_LIBS}"
					BLAS_PATH="${BLAS_HOME}"

					AS_IF([  test "x$BLAS_HOME" != "x/usr" -a "x$BLAS_HOME" != "x/usr/local"],
							[BLAS_LIBS="-L${BLAS_HOME} ${ATLAS_LIBS}"])
					]
					)

				CXXFLAGS="${BACKUP_CXXFLAGS} ${CBLAS_FLAG}"
				LIBS="${BACKUP_LIBS} ${BLAS_LIBS}"

					dnl  echo $LIBS

				AC_TRY_LINK(
					[#define __FFLASFFPACK_CONFIGURATION
					#include "fflas-ffpack/config-blas.h"],
					[double a;],
					[
					AC_TRY_RUN(
						[#define __FFLASFFPACK_CONFIGURATION
						#include "fflas-ffpack/config-blas.h"
						int main () {  double a[4] = {1.,2.,3.,4.}; double b[4]= {4.,3.,2.,1.}; double c[4];
						cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans,2,2,2,1., a,2,b,2,0.,c,2);
						if ( (c[0]!=8.) && (c[1]!=5.) && (c[2]!=20.) && (c[3]!=13))
						return -1;
						else
						return 0;
						}
						],[
						blas_found="yes"
						break
						],[
						blas_problem="$problem $BLAS_HOME"
						unset BLAS_LIBS
						],[
						blas_found="yes"
						blas_cross="yes"
						break
						]
					)
					],
					[
					blas_found="no"
					blas_checked="$checked $BLAS_HOME"
					unset BLAS_LIBS
					]
			)
		done



		AS_IF([ test "x$blas_found" = "xyes" ],
				BLAS_VENDOR="ATLAS"
				AC_SUBST(BLAS_VENDOR)

				[ AC_SUBST(BLAS_LIBS)
				AC_SUBST(BLAS_PATH)
				AC_SUBST(CBLAS_FLAG)
				AC_DEFINE(HAVE_BLAS,1,[Define if BLAS is installed])
				AC_DEFINE(HAVE_CBLAS,1,[Define if C interface to BLAS is available])
				dnl  AC_DEFINE(BLAS_AVAILABLE,,[Define if BLAS routines are available])
				HAVE_BLAS=yes
				BLAS_FOUND=true
				AC_SUBST(BLAS_FOUND)

				AS_IF([ test "x$blas_cross" != "xyes" ],
					[AC_MSG_RESULT(found)],
					[AC_MSG_RESULT(unknown)
					echo "WARNING: You appear to be cross compiling, so there is no way to determine"
					echo "whether your BLAS are good. I am assuming it is."
					])

				ifelse([$2], , :, [$2])
				],
				[ test -n "$blas_problem" ],
				[ AC_MSG_RESULT(not working) ],
				[ test "x$blas_found" = "xno" ],
				[ AC_MSG_RESULT(not found)]
				)
				AS_IF([  test "x$blas_found" != "xyes" ],
						[
						AC_MSG_CHECKING(for C interface to BLAS with -lblas)

		dnl **************************************
		dnl Check first for C interface to BLAS
		dnl **************************************


		for BLAS_HOME in ${BLAS_HOME_PATH} ; do
			dnl  echo looking in ${BLAS_HOME}
			CBLAS="yes"
			CBLAS_FLAG="-D__FFLASFFPACK_HAVE_CBLAS"
			ATLAS_LIBS="-lblas"
			AS_IF(
					dnl lib/libblas.* ?
					[ test -r "$BLAS_HOME/lib/libblas.a" -o -r "$BLAS_HOME/lib/libblas.so" ],
					[
					dnl  CBLAS_SYM=`nm -Du $BLAS_HOME/lib/libblas.so | grep cblas_'
					BLAS_LIBS=" ${ATLAS_LIBS}"
					BLAS_PATH="${BLAS_HOME}/lib"
					AS_IF([  test "x$BLAS_HOME" != "x/usr" -a "x$BLAS_HOME" != "x/usr/local"],
						[BLAS_LIBS="-L${BLAS_HOME}/lib ${ATLAS_LIBS}"])
					],
					dnl libblas.* ?
					[ test -r "$BLAS_HOME/libblas.a" -o -r "$BLAS_HOME/libblas.so" ],
					[ BLAS_LIBS=" ${ATLAS_LIBS}"
					BLAS_PATH="${BLAS_HOME}"
					AS_IF([ test "x$BLAS_HOME" != "x/usr" -a "x$BLAS_HOME" != "x/usr/local"],
						[BLAS_LIBS="-L${BLAS_HOME} ${ATLAS_LIBS}"])
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
						int main () {  double a[4] = {1.,2.,3.,4.}; double b[4]= {4.,3.,2.,1.}; double c[4];
						cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans,2,2,2,1., a,2,b,2,0.,c,2);
						if ( (c[0]!=8.) && (c[1]!=5.) && (c[2]!=20.) && (c[3]!=13))
						return -1;
						else
						return 0;
						}
						],[
						blas_found="yes"
						break
						],[
						blas_problem="$problem $BLAS_HOME"
						unset BLAS_LIBS
						],[
						blas_found="yes"
						blas_cross="yes"
						break
						]
					)
					],
					[
					blas_found="no"
					blas_checked="$checked $BLAS_HOME"
					unset BLAS_LIBS
					]
			)
		done ;

		AS_IF([ test "x$blas_found" = "xyes" ],
				[ 	BLAS_VENDOR="OTHER"
				AC_SUBST(BLAS_VENDOR)
				AC_SUBST(BLAS_LIBS)
				AC_SUBST(BLAS_PATH)
				AC_SUBST(CBLAS_FLAG)
				AC_DEFINE(HAVE_BLAS,1,[Define if BLAS is installed])
				AC_DEFINE(HAVE_CBLAS,1,[Define if C interface to BLAS is available])
				dnl  AC_DEFINE(BLAS_AVAILABLE,,[Define if BLAS routines are available])
				HAVE_BLAS=yes
				BLAS_FOUND=true
				AC_SUBST(BLAS_FOUND)
				AS_IF([ test "x$blas_cross" != "xyes" ],
					[AC_MSG_RESULT(found)],
					[AC_MSG_RESULT(unknown)
					echo "WARNING: You appear to be cross compiling, so there is no way to determine"
					echo "whether your BLAS are good. I am assuming it is."
					])
				ifelse([$2], , :, [$2])
				],
				[ test -n "$blas_problem" ],
				[ AC_MSG_RESULT(not working) ],
				dnl  echo "Sorry, your BLAS are not working. Disabling."
				[ test "x$blas_found" = "xno" ],
				[ AC_MSG_RESULT(not found)]
			)
	])


	AM_CONDITIONAL(FFLASFFPACK_HAVE_BLAS, test "x$HAVE_BLAS" = "xyes")

	CXXFLAGS=${BACKUP_CXXFLAGS}
	LIBS=${BACKUP_LIBS}
	dnl  unset LD_LIBRARY_PATH


])
