dnl Copyright (c) 2012 FFLAS-FFPACK
dnl Written by Clément Pernet, Brice Boyer.
dnl This file was taken from LinBox linbox-opt.m4
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






AC_DEFUN([FF_PRECOMPILE],
[

AC_MSG_CHECKING([whether to compile the standard specializations])

AC_ARG_ENABLE(precompilation,
[AC_HELP_STRING([--enable-precompilation], [ Enable precompilation of the standard specializations])])
AM_CONDITIONAL(FFLASFFPACK_PRECOMPILED, test "x$enable_precompilation" = "xyes")
AS_IF([test "x$enable_precompilation" = "xyes"],
	    [
		AC_MSG_RESULT(yes)
		PRECOMPILE_FLAGS="-DFFLAS_COMPILED -DFFPACK_COMPILED"
		PRECOMPILE_LIBS="-L${libdir} -lfflas -lffpack"
		AC_SUBST(PRECOMPILE_FLAGS)
		AC_SUBST(PRECOMPILE_LIBS)
	    ],
	    [AC_MSG_RESULT(no)]
     )
])
