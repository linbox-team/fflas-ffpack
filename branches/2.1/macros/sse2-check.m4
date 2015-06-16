dnl Check for SSE
dnl  Copyright (c) 2011 FFLAS-FFPACK
dnl Created by BB, 2014-03-25
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
dnl

dnl FF_CHECK_SSE
dnl
dnl turn on  SSE4.1 extensions if available

AC_DEFUN([FF_CHECK_SSE],
		[
		AC_ARG_ENABLE(sse,
			[AC_HELP_STRING([--enable-sse],
				[ Use Intel(r) SSE 4.1])
			],
			[ avec_sse=$enable_sse ],
			[ avec_sse=yes ]
			)
		AC_MSG_CHECKING(for SSE 4.1)
		AS_IF([ test  "x$avec_sse" != "xno" ],
			[
			BACKUP_CXXFLAGS=${CXXFLAGS}
			dnl  SSEFLAGS="-msse2"
			SSEFLAGS="-msse4.1"
			CXXFLAGS="${BACKUP_CXXFLAGS} ${SSEFLAGS}"
			CODE_SSE=`cat macros/CodeChunk/sse.C`
			AC_TRY_RUN([
				${CODE_SSE}
				],
				[ sse_found="yes" ],
				[ sse_found="no" ],
				[
				echo "cross compiling...disabling"
				sse_found="no"
				])
			AS_IF([ test "x$sse_found" = "xyes" ],[
				AC_DEFINE(USE_SSE,1,[Define if SSE is available])
				AC_SUBST(SSEFLAGS)
				AC_MSG_RESULT(yes (SSE))
				],
				[
				SSEFLAGS=""
				AC_MSG_RESULT(no)
				]
				)
			CXXFLAGS=${BACKUP_CXXFLAGS}
			],
			[ AC_MSG_RESULT(no) ]
	)
	])
