dnl Copyright(c)'2011 FFLAS-FFPACK
dnl Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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


AC_DEFUN([FF_DOC],
[

AC_MSG_CHECKING(whether to build documentation)


AC_ARG_WITH(docdir,
[AC_HELP_STRING([--with-docdir=<path>], [Where the FFLAS-FFPACK documentation should be installed])],
            [
		FFLASFFPACK_DOC_PATH="$withval"
	    ],
	    [
		eval FFLASFFPACK_DOC_PATH="${prefix}/docs"
	    ])

AC_SUBST(FFLASFFPACK_DOC_PATH)

AC_ARG_WITH(doxygen,
[AC_HELP_STRING([--with-doxygen=<path>], [Give the path to Doxygen. Note: --enable-doc needed])],
            [
		DOXYGEN_PATH="$PATH $withval"
	    ],
	    [
		DOXYGEN_PATH="$PATH"
	    ])

AC_ARG_ENABLE(doc,[AC_HELP_STRING([--enable-doc], [Enable building documentation])],
[
AC_MSG_RESULT(yes)
AC_MSG_CHECKING(whether doxygen works)
export PATH=$DOXYGEN_PATH
(doxygen --version) < /dev/null > /dev/null 2>&1 || {
	AC_MSG_RESULT(no)
	echo
	echo "You must have doxygen installed to create documentation for"
	echo "FFLAS-FFPACK. This error only happens if you use --enable-doc."
	echo "Download the appropriate package for your distribution, or get"
	echo "the source tarball from http://www.stack.nl/~dimitri/doxygen/"
	exit -1
}
AC_MSG_RESULT(yes)
AM_CONDITIONAL(FFLASFFPACK_BUILD_DOC, true)
],
[
AC_MSG_RESULT(no)
AM_CONDITIONAL(FFLASFFPACK_BUILD_DOC, false)
])
])
