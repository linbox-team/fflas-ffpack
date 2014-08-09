dnl Check for CUDA
dnl Copyright(c)'1994-2009,2003,2013 by The Givaro group
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


dnl Modified by Pascal Giorgi, 2003-12-03
dnl Modified by BB, 2013-5-22 and other times

dnl Test for CUDA
dnl Sets CUDA_CFLAGS and CUDA_LIBS
dnl Defines HAVE_CUDA

AC_DEFUN([FF_CHECK_CUDA], [

		AC_ARG_WITH(cuda,
			[AC_HELP_STRING([--with-cuda=<path>|yes|no],[
				Use CUDA library.
				If argument is no, you do not have the library installed on your machine.
				If argument is yes or <empty> that means the library is reachable with the standard
				search path "/usr" or "/usr/local"  (set as default).
				Otherwise you give the <path> to the directory which contain the library.
				])],
			[if test "$withval" = yes ; then
			CUDA_HOME_PATH="${DEFAULT_CHECKING_PATH}"
			elif test "$withval" != no ; then
			CUDA_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH}"
			fi],
			[CUDA_HOME_PATH="${DEFAULT_CHECKING_PATH}"])

		min_cuda_version=ifelse([$1], ,5.5.0,$1)

		dnl Check for existence

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for CUDA >= $min_cuda_version )

dnl todo lib (32) and lib64.
CUDA_PATH=
for CUDA_HOME in ${CUDA_HOME_PATH}
do
	if test "x$CUDA_HOME" != "x/usr" -a "x$CUDA_HOME" != "x/usr/local"; then
		if test -r "$CUDA_HOME/include/cuda.h" ; then
			CUDA_CFLAGS="-I${CUDA_HOME}/include"
			CUDA_PATH="-L${CUDA_HOME}/lib64"
			CUDA_LIBS="-L${CUDA_HOME}/lib64 -lcusparse"
		else
			echo "($CUDA_HOME) seems an invalid CUDA prefix"
			echo "Searching CUDA in PATH"
			CUDA_CFLAGS=""
			CUDA_LIBS="-lcusparse"
		fi
	else
		CUDA_CFLAGS=""
		CUDA_LIBS="-lcusparse"
	fi

	CXXFLAGS="${CXXFLAGS} ${CUDA_CFLAGS}"
	LIBS="${LIBS} ${CUDA_LIBS}"
	CODE_CUDA=`cat macros/CodeChunk/cuda.C`

	AC_TRY_LINK(
		[
		#include <cuda.h>
		],
		[ CUresult a;],
		[
		dnl  # See if we are running CUDA 4.0 with --enable-cxx
		AC_TRY_RUN(
			[ ${CODE_CUDA} ],
			[
			AC_MSG_RESULT(found)
			AC_DEFINE(HAVE_CUDA,1,[Define if CUDA is installed])

			dnl  CUDA_VERSION="" dnl I could find it but why is it here ?
			CUDA_LIBS="${CUDA_PATH} -lcusparse"
			dnl  AC_SUBST(CUDA_VERSION)
			AC_SUBST(CUDA_LIBS)
			AC_SUBST(CUDA_CFLAGS)
			break;
			],[
			AC_MSG_RESULT(no : cuda is too old or not found)
			dnl  AC_SUBST(CUDA_VERSION)
			],[ dnl This should never happen
			AC_MSG_RESULT(no)
			])
		],[
			AC_MSG_RESULT(unknown)
			echo "WARNING: You appear to be cross compiling, so there is no way to determine"
			echo "whether your CUDA version is new enough. I am assuming it is."
			AC_SUBST(CUDA_CFLAGS)
			AC_SUBST(CUDA_LIBS)
			AC_DEFINE(HAVE_CUDA,1,[Define if CUDA is installed])
		])
	unset CUDA_CFLAGS
	unset CUDA_LIBS
done

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
