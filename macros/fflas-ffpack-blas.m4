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

dnl Tests BLAS for and define BLAS_CFLAGS and BLAS_LIBS
dnl Defines HAVE_LAPACK,  HAVE_CLAPACK,  HAVE_BLAS,  HAVE_CBLAS if available

AC_DEFUN([FF_CHECK_BLAS_CFLAGS],
		[ AC_ARG_WITH(blas-cflags,
			[AC_HELP_STRING([--with-blas-cflags=<cflags>],
				[ CFLAGS for BLAS/LAPACK (i.e. -I/path/to/toto-blas) ])
			])
                BLAS_CFLAGS="$with_blas_cflags"
		AC_SUBST(BLAS_CFLAGS)
		]
	)

dnl
AC_DEFUN([FF_CHECK_BLAS_LIBS],
		[ AC_ARG_WITH(blas-libs,
			[AC_HELP_STRING([--with-blas-libs=<libs>],
				[ LIBS for BLAS/LAPACK (i.e. -L/path/to/toto-blas -ltoto-blas) ])
			])
		BLAS_LIBS="$with_blas_libs"
		AC_SUBST(BLAS_LIBS)
		]
	)

dnl
AC_DEFUN([FF_CHECK_USER_BLAS],
		[
		BACKUP_CXXFLAGS=${CXXFLAGS}
		BACKUP_LIBS=${LIBS}
		saved_LD_RUN_PATH="$LD_RUN_PATH"
		blas_lib_path=`echo $BLAS_LIBS | $EGREP '\-L' | $SED -e 's/-L//;s/ .*//'`
		LD_RUN_PATH="${LD_RUN_PATH:+$LD_RUN_PATH$PATH_SEPARATOR}$blas_lib_path"
		export LD_RUN_PATH
		CODE_CBLAS=`cat ${srcdir}/macros/CodeChunk/cblas.C`
		CODE_FBLAS=`cat ${srcdir}/macros/CodeChunk/fblas.C`

		CXXFLAGS="${BACKUP_CXXFLAGS} ${BLAS_CFLAGS} -I. -I.. -I${srcdir} -I${srcdir}/fflas-ffpack ${GIVARO_CFLAGS}"
		LIBS="${BACKUP_LIBS} ${BLAS_LIBS}"

                dnl is there user CBLAS accessible ?
		AC_MSG_CHECKING(for BLAS)

		AC_TRY_RUN([ ${CODE_CBLAS} ],[
				blas_found="yes"
                                is_cblas="yes"
                                AC_MSG_RESULT(found CBLAS)
				],[
                                dnl No, then checking for Fortran BLAS
                                AC_TRY_RUN(
                                  [ ${CODE_FBLAS} ],
                                   [ blas_found="yes"
                                     is_cblas="no"
                                     AC_MSG_RESULT(found Fortran BLAS)
                                   ],[
                                     dnl No, then checking for  OpenBLAS
                                     BLAS_LIBS="${BLAS_LIBS} -lopenblas -lpthread"
                                     AS_CASE([$CCNAM], [gcc*], [BLAS_LIBS="${BLAS_LIBS} -lgfortran"])
                                     LIBS="${BACKUP_LIBS} ${BLAS_LIBS}"
				     AC_TRY_RUN(
					[ ${CODE_CBLAS} ],[
                                        AC_MSG_RESULT(found OpenBLAS)
					blas_found="yes"
                                        is_cblas="yes"
					AC_SUBST(BLAS_LIBS)
					],[
					blas_problem="$problem"
			                AC_MSG_RESULT(problem)
					],[
					blas_found="yes"
                                        is_cblas="yes"
			                AC_MSG_RESULT(cross compiling)
					blas_cross="yes"
					AC_SUBST(BLAS_LIBS)
					])
                                        ],[
		                       blas_found="yes"
                                       is_cblas="no"
			               AC_MSG_RESULT(cross compiling)
				       blas_cross="yes"
                                   ])
				],[
				blas_found="yes"
		                AC_MSG_RESULT(cross compiling)
                                is_cblas="yes"
				blas_cross="yes"
				])

		AS_IF([ test "x$blas_found" = "xyes" ],
				[
				BLAS_VENDOR="USER"
				AC_SUBST(BLAS_VENDOR)
                                AC_DEFINE(HAVE_BLAS,1,[Define if BLAS is available])
                                AS_IF([test "x$is_cblas" = "xyes" ],[
                                            AC_DEFINE(HAVE_CBLAS,1,[Define if BLAS is a CBLAS])
                                      ],[])

				AS_IF([test "x$blas_cross" = "xyes"], [
					echo "WARNING: You appear to be cross compiling, so there is no way to determine"
					echo "whether your BLAS are good. I am assuming it is."],[])
				],
				[
                                echo ''
	echo '*******************************************************************************'
	echo ' ERROR: BLAS not found!'
	echo
	echo ' BLAS routines are required for this library to compile. Please'
	echo ' make sure BLAS are installed and specify its location with the option'
	echo ' --with-blas-libs=<libs> and if necessary --with-blas-cflags=<cflags>'
	echo ' when running configure.'
	echo '*******************************************************************************'
	exit 1
        ])


        dnl	AM_CONDITIONAL(FFLASFFPACK_HAVE_BLAS, test "x$blas_found" = "xyes")
        dnl     AM_CONDITIONAL(FFLASFFPACK_HAVE_CBLAS, test "x$is_cblas" = "xyes")
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

		CODE_CLAPACK=`cat ${srcdir}/macros/CodeChunk/clapack.C`
		CODE_LAPACK=`cat ${srcdir}/macros/CodeChunk/lapack.C`

		AC_MSG_CHECKING(for USER LAPACK)

		CXXFLAGS="${BACKUP_CXXFLAGS} ${BLAS_CFLAGS}  -I. -I.. -I${srcdir} -I${srcdir}/flas-ffpack ${GIVARO_CFLAGS}"
		LIBS="${BACKUP_LIBS} ${BLAS_LIBS}"

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

AC_DEFUN([FF_OPENBLAS_NUM_THREADS],
		[ AC_ARG_WITH(openblas-num-threads,
			[AC_HELP_STRING([--with-openblas-num-threads=<num-threads>],
				[ Set the number of threads given to OpenBLAS])
				])
		dnl testing if we are using openblas
		BACKUP_CXXFLAGS=${CXXFLAGS}
		BACKUP_LIBS=${LIBS}

		CODE_OPENBLAS='extern "C"{void openblas_set_num_threads(int num_threads);} int main(){openblas_set_num_threads(1);return 0;}'

		AC_MSG_CHECKING(if this is OpenBLAS)

		CXXFLAGS="${BACKUP_CXXFLAGS} ${BLAS_CFLAGS}"
		LIBS="${BACKUP_LIBS} ${BLAS_LIBS}"

		AC_TRY_RUN(
			[ ${CODE_OPENBLAS} ],
			[ openblas_found="yes" ],
			[ openblas_problem="problem" ],
			[ openblas_found="" ]
		)

		AS_IF([test "x$openblas_found" = "xyes"],
		      [
		       AC_MSG_RESULT(yes)
                       AC_MSG_CHECKING(for OPENBLAS numthreads)
		       AS_IF([test "x$with_openblas_num_threads" = "x"],
		       [
			AC_MSG_RESULT(none specified (using default value 1))
			numthreads="1"
			],
		       [AC_MSG_RESULT($with_openblas_num_threads)
		        numthreads=$with_openblas_num_threads
			])
		       AC_DEFINE_UNQUOTED(OPENBLAS_NUM_THREADS,$numthreads,[Sets the number of threads given to OpenBLAS (default is 1)])
		       ],
		       [AC_MSG_RESULT(no)]
		       )
	])

