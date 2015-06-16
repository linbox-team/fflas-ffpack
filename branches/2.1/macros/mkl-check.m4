dnl  Check for MKL
dnl  Brice Boyer 2014
dnl  This file is part of FFLAS-FFPACK

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



AC_DEFUN([FF_CHECK_MKL],
		[
		AC_MSG_CHECKING(for use of MKL)
		dnl  echo $CBLAS_LIBS
                USE_MKL="false"
		MKL_USED=`echo $CBLAS_LIBS | grep -i MKL`
		AS_IF( [test -n "$MKL_USED"] , [
			AC_DEFINE(HAVE_MKL,1,[Define if we use MKL for blas/lapack])
			USE_MKL="true"
			AC_SUBST(USE_MKL)
			AC_MSG_RESULT( yes )
			]
			,
			[
			AC_MSG_RESULT( no )
			]
		     )
		]
	)
