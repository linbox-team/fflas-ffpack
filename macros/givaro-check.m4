dnl Check for GIVARO
dnl Copyright (c) the Givaro group
dnl This file is part of FFLAS-FFPACK

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



dnl adapted from LinBox by BB.

dnl FF_CHECK_GIVARO ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Tests for Givaro and define GIVARO_CFLAGS and GIVARO_LIBS
dnl Defines HAVE_GIVARO

AC_DEFUN([FF_CHECK_GIVARO],
[

AC_ARG_WITH(givaro,
[AC_HELP_STRING([--with-givaro=<path>|yes], [Use Givaro library. This library is mandatory for
                           LinBox compilation. If argument is yes or <empty>
			   that means the library is reachable with the standard
			   search path (/usr or /usr/local). Otherwise you give
			   the <path> to the directory which contains the
			   library.
])],
	     [if test "$withval" = yes ; then
			GIVARO_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	      elif test "$withval" != no ; then
			GIVARO_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH}"
	     fi],
	     [GIVARO_HOME_PATH="${DEFAULT_CHECKING_PATH}"])

dnl -------------- dnl
dnl GIVARO VERSION dnl
dnl -------------- dnl

dnl As we need Integer and Modular, should be updated on each interface changes
version_min=40001
version_max=40003

dnl Check for existence

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}
saved_LD_RUN_PATH="$LD_RUN_PATH"

AC_MSG_CHECKING(for GIVARO >= $version_min and < $version_max)

for GIVARO_HOME in ${GIVARO_HOME_PATH}
 do
if test -r "$GIVARO_HOME/include/givaro/givconfig.h"; then

	# Givaro Libs + CFlags contain GMP info - AB 2014-12-12
	GIVARO_LIBS=`$GIVARO_HOME/bin/givaro-config --libs`
	GIVARO_CFLAGS=`$GIVARO_HOME/bin/givaro-config --cflags`
	givaro_lib_path=`$GIVARO_HOME/bin/givaro-config --prefix`/lib
	CXXFLAGS="${BACKUP_CXXFLAGS} ${GIVARO_CFLAGS}"
	LIBS="${BACKUP_LIBS} ${GIVARO_LIBS}"
	LD_RUN_PATH="${LD_RUN_PATH:+$LD_RUN_PATH$PATH_SEPARATOR}$givaro_lib_path"
	export LD_RUN_PATH
	AC_TRY_LINK(
	[#include <givaro/givinteger.h>],
	[Givaro::Integer a;],
	[
	AC_TRY_RUN(
	[#include <givaro/givconfig.h>
	 int main () { if (GIVARO_VERSION >= $version_min && GIVARO_VERSION < $version_max) return 0; else return -1; /* old version of Givaro are defined as hexa 0x03yyzz*/ }
	],[
	givaro_found="yes"
	break
	],[
	givaro_problem="$problem $GIVARO_HOME"
	unset GIVARO_CFLAGS
	unset GIVARO_LIBS
	],[
	givaro_found="yes"
	givaro_cross="yes"

	break
	])
	],
	[
	givaro_found="yes"
	givaro_checked="$checked $GIVARO_HOME"
#unset GIVARO_CFLAGS
#unset GIVARO_LIBS
	break

	])
else
	givaro_found="no"
fi
done

if test "x$givaro_found" = "xyes" ; then
	AC_SUBST(GIVARO_CFLAGS)
	AC_SUBST(GIVARO_LIBS)
	dnl  echo $GIVARO_CFLAGS $GIVARO_LIBS
	AC_DEFINE(HAVE_GIVARO,1,[Define if GIVARO is installed])
	HAVE_GIVARO=yes
	
	if test "x$givaro_cross" != "xyes"; then
		AC_MSG_RESULT(found)
	else
		AC_MSG_RESULT(unknown)
		echo "WARNING: You appear to be cross compiling, so there is no way to determine"
		echo "whether your GIVARO version is new enough. I am assuming it is."
	fi
	
	ifelse([$2], , :, [$2])
elif test -n "$givaro_problem"; then
	AC_MSG_RESULT(problem)
	echo "Sorry, your GIVARO version is too old. Disabling."
	ifelse([$3], , :, [$3])
elif test "x$givaro_found" = "xno" ; then
	AC_MSG_RESULT(not found)
	ifelse([$3], , :, [$3])
fi

AM_CONDITIONAL(FFLASFFPACK_HAVE_GIVARO, test "x$HAVE_GIVARO" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
LD_RUN_PATH="$saved_LD_RUN_PATH"
export LD_RUN_PATH
unset saved_LD_RUN_PATH
#unset LD_LIBRARY_PATH

])

AC_DEFUN([FF_CHECK_GIVARO_USABILITY],
[
    dnl backup
    BACKUP_CXXFLAGS=${CXXFLAGS}
    BACKUP_LIBS=${LIBS}
    dnl add GIVARO flags
	CXXFLAGS="${BACKUP_CXXFLAGS} ${GIVARO_CFLAGS}"
	LIBS="${BACKUP_LIBS} ${GIVARO_LIBS}"
    dnl try to compile a small example
    AC_MSG_CHECKING([for GIVARO usability])
    AC_TRY_LINK([#include <givaro/givinteger.h>], [Givaro::Integer a;],
                [AC_MSG_RESULT(yes)],
                [AC_MSG_RESULT(no)
                 AC_MSG_ERROR(The Givaro library could not be used with the compiler and the flags set up by the configure script)
                ])
    dnl restore backu
    CXXFLAGS=${BACKUP_CXXFLAGS}
    LIBS=${BACKUP_LIBS}
])
