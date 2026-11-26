dnl Copyright (c) 2026 FFLAS-FFPACK group
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

dnl FF_CHECK_METAL
dnl
dnl --enable-metal: on an Apple-silicon Mac, share the single precision matrix products
dnl with the GPU (Metal Performance Shaders); see fflas-ffpack/fflas/fflas_metal/fflas_metal.h.
dnl Defines HAVE_METAL, sets METAL_FRAMEWORKS and METAL_LIBS and the conditional
dnl FFLASFFPACK_HAVE_METAL, under which libfflas_metal.a is built.

AC_DEFUN([FF_CHECK_METAL],
	[ AC_ARG_ENABLE(metal,
		[AC_HELP_STRING([--enable-metal],
				[ Share single precision matrix products with the GPU of an Apple-silicon Mac (Metal Performance Shaders) ])
		],
		[ avec_metal=$enable_metal ],
		[ avec_metal=no ]
		)
	  AC_MSG_CHECKING(whether to share single precision products with the GPU (Metal))
	  metal_found=no
	  METAL_FRAMEWORKS=
	  METAL_LIBS=
	  AS_IF([ test "x$avec_metal" = "xyes" ],
		[
		BACKUP_LIBS=${LIBS}
		METAL_FRAMEWORKS="-framework Foundation -framework Metal -framework MetalPerformanceShaders"
		LIBS="${BACKUP_LIBS} ${METAL_FRAMEWORKS}"
		AC_LINK_IFELSE([AC_LANG_PROGRAM([[]], [[return 0;]])],
			[ metal_found=yes ],
			[ metal_found=no ])
		LIBS=${BACKUP_LIBS}
		AS_IF([ test "x$metal_found" = "xyes" ],
			[
				AC_DEFINE(HAVE_METAL,1,[Define if single precision products are shared with the GPU through Metal])
				METAL_LIBS="-L${libdir} -lfflas_metal ${METAL_FRAMEWORKS}"
				AC_MSG_RESULT(yes)
			],
			[
				METAL_FRAMEWORKS=
				AC_MSG_RESULT(no)
				AC_MSG_ERROR([--enable-metal needs the Foundation, Metal and MetalPerformanceShaders frameworks (macOS)])
			])
		],
		[ AC_MSG_RESULT(no) ]
	  )
	  AC_SUBST(METAL_FRAMEWORKS)
	  AC_SUBST(METAL_LIBS)
	  AM_CONDITIONAL(FFLASFFPACK_HAVE_METAL, test "x$metal_found" = "xyes")
	]
)
