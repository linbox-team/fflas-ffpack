dnl  Check for BLAS
dnl  Copyright 2014 Brice Boyer (briceboyer) <boyer.brice@gmail.com>
dnl This file is part of FFLAS-FFPACK
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

dnl Tests BLAS for and define CBLAS_FLAG and CBLAS_LIBS
dnl Defines HAVE_LAPACK,  HAVE_CLAPACK,  HAVE_BLAS,  HAVE_CBLAS if available

AC_DEFUN([FF_CHECK_BLAS_CFLAGS],
		[ AC_ARG_WITH(blas-cflags,
			[AC_HELP_STRING([--with-blas-cflags=<cflags>],
				[ CFLAGS for BLAS/LAPACK (i.e. -I/path/to/toto-blas) ])
			])
		CBLAS_FLAG="$with_blas_cflags -D__FFLASFFPACK_HAVE_CBLAS"
		AC_SUBST(CBLAS_FLAG)
		dnl  echo $CBLAS_FLAG;
		]
	)

dnl
AC_DEFUN([FF_CHECK_BLAS_LIBS],
		[ AC_ARG_WITH(blas-libs,
			[AC_HELP_STRING([--with-blas-libs=<libs>],
				[ LIBS for BLAS/LAPACK (i.e. -L/path/to/toto-blas -ltoto-blas) ])
			])
		CBLAS_LIBS="$with_blas_libs"
		AC_SUBST(CBLAS_LIBS)
		dnl  echo $CBLAS_LIBS;
		]
	)

dnl
AC_DEFUN([FF_CHECK_USER_BLAS],
		[
		BACKUP_CXXFLAGS=${CXXFLAGS}
		BACKUP_LIBS=${LIBS}
		saved_LD_RUN_PATH="$LD_RUN_PATH"
		blas_lib_path=`echo $CBLAS_LIBS | $EGREP '\-L' | $SED -e 's/-L//;s/ .*//'`
		LD_RUN_PATH="${LD_RUN_PATH:+$LD_RUN_PATH$PATH_SEPARATOR}$blas_lib_path"
		export LD_RUN_PATH
		CODE_CBLAS=`cat macros/CodeChunk/cblas.C`

		AC_MSG_CHECKING(for USER BLAS)

		CXXFLAGS="${BACKUP_CXXFLAGS} ${CBLAS_FLAG} -I. -I.. -I`pwd` -I`pwd`/fflas-ffpack ${GIVARO_CFLAGS}"
		LIBS="${BACKUP_LIBS} ${CBLAS_LIBS}"

		AC_TRY_LINK( [
#define __FFLASFFPACK_CONFIGURATION
#include "fflas-ffpack/config-blas.h"],
			[double a;],
			[
			AC_TRY_RUN(
				[ ${CODE_CBLAS} ],[
				blas_found="yes"
				],[
				blas_problem="$problem"
				],[
				blas_found="yes"
				blas_cross="yes"
				])
			],
			[
			blas_found="no"
			])

		AS_IF([ test "x$blas_found" = "xyes" ],
				[
				BLAS_VENDOR="USER"
				AC_SUBST(BLAS_VENDOR)
				dnl  AC_SUBST(CBLAS_FLAG)
				dnl  AC_SUBST(BLAS_PATH)
				AC_DEFINE(HAVE_BLAS,1,[Define if BLAS is installed])
				AC_DEFINE(HAVE_CBLAS,1,[Define if C interface to BLAS is installed])
				BLAS_FOUND=true
				AC_SUBST(BLAS_FOUND)
				dnl  AC_DEFINE(BLAS_AVAILABLE,,[Define if BLAS routines are available])
				#echo ${CBLAS_FLAG}
				#echo ${CBLAS_LIBS}
				HAVE_BLAS=yes
				AS_IF([test "x$blas_cross" != "xyes"],
					[ AC_MSG_RESULT(found (cblas)) ] ,
					[AC_MSG_RESULT(unknown)
					echo "WARNING: You appear to be cross compiling, so there is no way to determine"
					echo "whether your BLAS are good. I am assuming it is."])
				],
				[
				AC_MSG_RESULT(problem)
				]
		)


		AM_CONDITIONAL(FFLASFFPACK_HAVE_BLAS, test "x$HAVE_BLAS" = "xyes")

		CXXFLAGS=${BACKUP_CXXFLAGS}
		LIBS=${BACKUP_LIBS}
		LD_RUN_PATH="$saved_LD_RUN_PATH"
		export LD_RUN_PATH
		unset saved_LD_RUN_PATH
		dnl  unset LD_LIBRARY_PATH


		]
	)

dnl
AC_DEFUN([FF_CHECK_USER_LAPACK],
	[
		BACKUP_CXXFLAGS=${CXXFLAGS}
		BACKUP_LIBS=${LIBS}

		CODE_CLAPACK=`cat macros/CodeChunk/clapack.C`
		CODE_LAPACK=`cat macros/CodeChunk/lapack.C`

		AC_MSG_CHECKING(for USER LAPACK)

		CXXFLAGS="${BACKUP_CXXFLAGS} ${CBLAS_FLAG}  -I. -I.. -I`pwd` -I`pwd`/fflas-ffpack ${GIVARO_CFLAGS}"
		LIBS="${BACKUP_LIBS} ${CBLAS_LIBS}"

		AC_TRY_RUN(
			[ ${CODE_CLAPACK} ],
			[ dgetrf_found="yes" ],
			[ dgetrf_problem="problem" ],
			[ dgetrf_found="" ]
		)

		AS_IF([ test "${dgetrf_found}" = "yes"],
			[
				AC_MSG_RESULT( yes (clapack))
				AC_DEFINE(HAVE_LAPACK,1,[Define if LAPACK is installed])
				AC_DEFINE(HAVE_CLAPACK,1,[Define if C interface to LAPACK is available])
				HAVE_LAPACK=yes
			],
			[
				AC_TRY_RUN(
					[  ${CODE_LAPACK} ],
					[ dgetrf_found="yes"],
					[ dgetrf_problem="$problem"],
					[ dgetrf_found="" ]
				)
				AS_IF([ test "x${dgetrf_found}" = "xyes"],
					[
						AC_SUBST(LAPACK_LIBS)
						AC_MSG_RESULT( yes (lapack))
						AC_DEFINE(HAVE_LAPACK,1,[Define if LAPACK is installed])
						HAVE_LAPACK=yes
					], dnl clapack not found. looking for lapack
					[
					AC_MSG_RESULT( no )
					]
				)
			]
		)

		dnl
		AM_CONDITIONAL(FFLASFFPACK_HAVE_LAPACK, test "x$HAVE_LAPACK" = "xyes")
		CXXFLAGS=${BACKUP_CXXFLAGS}
		LIBS=${BACKUP_LIBS}
		dnl  unset LD_LIBRARY_PATH

	]
)

