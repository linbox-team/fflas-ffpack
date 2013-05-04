dnl  Check for BLAS
dnl  Copyright Pascal Giorgi 2005
dnl  Modified Brice Boyer 2011
dnl  This file is part of FFLAS-FFPACK (and comes from LinBox)

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



dnl **********************************
dnl *              TODO              *
dnl **********************************
dnl no support yet to MKL
dnl AS_IF([test -r "$BLAS_VAL/include/mkl_cblas.h"],
dnl [ BLAS_LIBS="-L${BLAS_VAL}/lib/${MKL_ARCH}/ -lmkl_lapack64 -lmkl -lvml -lguide" ])
dnl **********************************



dnl FF_CHECK_BLAS ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for BLAS and define BLAS_LIBS

AC_DEFUN([FF_CHECK_BLAS],
		[ AC_ARG_WITH(blas,
			[AC_HELP_STRING([--with-blas=<lflags>],
				[Use BLAS library. This library is mandatory for FFLAS-FFPACK
				compilation. The user has the responsability to
				provide library flags such that the compiler
				will find and use BLAS (and LAPACK). An example
				could be --with-blas=/path/to/blas or
				--with-blas="-L/path/to/blas/lib -lsomeblas".\n
				* Warning : we don't really handle .a archives alone...])
			])

		BACKUP_CXXFLAGS=${CXXFLAGS}
		BACKUP_LIBS=${LIBS}

		AS_IF([ test -n "$with_blas"],[

		AC_MSG_CHECKING("for BLAS ($with_blas)")

		BLAS_PATH="$with_blas"
		BLAS_LIBS="$with_blas"
		CBLAS_FLAG="-D__FFLASFFPACK_HAVE_CBLAS"
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
					],[
					blas_problem="$problem $BLAS_HOME"
					unset BLAS_LIBS
					],[
					blas_found="yes"
					blas_cross="yes"
					])
				],
				[
				blas_found="no"
				])


		AS_IF([ test "x$blas_found" = "xyes" ],
				[ 	BLAS_VENDOR="USER"
				AC_SUBST(BLAS_VENDOR)
				AC_SUBST(BLAS_LIBS)
				AC_SUBST(CBLAS_FLAG)
				AC_SUBST(BLAS_PATH)
				AC_DEFINE(HAVE_BLAS,1,[Define if BLAS is installed])
				AC_DEFINE(HAVE_CBLAS,1,[Define if C interface to BLAS is installed])
				BLAS_FOUND=true
				AC_SUBST(BLAS_FOUND)
				dnl  AC_DEFINE(BLAS_AVAILABLE,,[Define if BLAS routines are available])
				HAVE_BLAS=yes
				AS_IF([test "x$blas_cross" != "xyes"],
					[ AC_MSG_RESULT(found (cblas)) ] ,
					[AC_MSG_RESULT(unknown)
					echo "WARNING: You appear to be cross compiling, so there is no way to determine"
					echo "whether your BLAS are good. I am assuming it is."])
				], dnl CBLAS not found. Looking for BLAS
				[
				CBLAS_FLAG=""
				CXXFLAGS="${BACKUP_CXXFLAGS}  ${CBLAS_FLAG}"

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
						],[
						blas_problem="$problem $BLAS_HOME"
						unset BLAS_LIBS
						],[
						blas_found="yes"
						blas_cross="yes"
						])
					],
					[
					blas_found="no"
					])
					AS_IF([test "x$blas_found" = "xyes"],
							[	BLAS_VENDOR="USER"
							AC_SUBST(BLAS_VENDOR)
							AC_SUBST(BLAS_LIBS)
							AC_SUBST(CBLAS_FLAG)
							AC_SUBST(BLAS_PATH)
							AC_DEFINE(HAVE_BLAS,1,[Define if BLAS is installed])
							dnl  AC_DEFINE(HAVE_CBLAS,1,[Define if C interface to BLAS is installed])
							BLAS_FOUND=true
							AC_SUBST(BLAS_FOUND)
							dnl  AC_DEFINE(BLAS_AVAILABLE,,[Define if BLAS routines are available])
							HAVE_BLAS=yes
							AC_MSG_RESULT(found (cblas)) ] ,
							[ AC_MSG_RESULT(no) ])

					])


		])


		AM_CONDITIONAL(FFLASFFPACK_HAVE_BLAS, test "x$HAVE_BLAS" = "xyes")

		CXXFLAGS=${BACKUP_CXXFLAGS}
		LIBS=${BACKUP_LIBS}
		dnl  unset LD_LIBRARY_PATH


])


