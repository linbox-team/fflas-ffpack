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

dnl **********************************
dnl *              TODO              *
dnl **********************************
dnl
dnl no support yet to MKL
dnl AS_IF([test -r "$BLAS_VAL/include/mkl_cblas.h"],
dnl [ BLAS_LIBS="-L${BLAS_VAL}/lib/${MKL_ARCH}/ -lmkl_lapack64 -lmkl -lvml -lguide" ])
dnl
dnl test in a chunk by using some intel type, possibly the standard one with define/ifndef TRY_MKL
dnl won't do/works for me.
dnl
dnl **********************************



AC_DEFUN([FF_CHECK_MKL],
		MKL_USED=`echo $BLAS_LIBS | grep -i LMKL`
		AS_IF( [test -n "$MKL_USED"] , [AC_DEFINE(HAVE_MLK,1,[Define if we use MKL for blas/lapack])])
	)


