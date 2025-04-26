dnl turn on OPENMP
dnl  Copyright (c) 2011 FFLAS-FFPACK
dnl Created by BB, 2014-07-01
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
dnl License along with this library; if not, see
dnl <https://www.gnu.org/licenses/>.
dnl ========LICENCE========
dnl

dnl FF_CHECK_OMP
dnl
dnl turn on OpenMP if available

AC_DEFUN([FF_CHECK_OMP],
	[ AC_ARG_ENABLE(openmp,
		[AC_HELP_STRING([--enable-openmp],
				[ Use OpenMP ])
		],
		[ avec_omp=$enable_openmp],
		[ avec_omp=yes ]
		)
	  AC_MSG_CHECKING(for OpenMP)
	  AS_IF([ test "x$avec_omp" != "xno" ],
		[
		BACKUP_CXXFLAGS=${CXXFLAGS}
		CXXFLAGS="${BACKUP_CXXFLAGS} ${OPENMP_CXXFLAGS}"
		AC_TRY_RUN([
#include <omp.h>
			int main() {
			int p = omp_get_num_threads();
			return 0;
			}
		],
		[ omp_found="yes" ],
		[ omp_found="no" ],
		[
			echo "cross compiling...disabling"
			omp_found="no"
		])
		AS_IF(	[ test "x$omp_found" = "xyes" ],
			[
				AC_DEFINE(USE_OPENMP,1,[Define if OMP is available])
				AC_SUBST(OPENMP_CXXFLAGS)
				AC_MSG_RESULT(yes)
				HAVE_OMP=yes
			],
			[
				OPENMP_CXXFLAGS=
				AC_SUBST(OPENMP_CXXFLAGS)
				AC_MSG_RESULT(no)
			]
		)
		CXXFLAGS=${BACKUP_CXXFLAGS}
		],
		[ AC_MSG_RESULT(no) ]
	)
	AM_CONDITIONAL(FFLASFFPACK_HAVE_OMP, test "x$HAVE_OMP" = "xyes")
]	
)
