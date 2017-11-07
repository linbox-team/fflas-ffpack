dnl Check for SIMD
dnl  Copyright (c) 2011 FFLAS-FFPACK
dnl Created by BB, 2014-03-25
dnl modified by CP, 2016-07-11
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

dnl FF_CHECK_SIMD
dnl
dnl turn on SSE4.1 AVX, AVX2 extensions if available

AC_DEFUN([FF_CHECK_SIMD],
[
	AC_ARG_ENABLE(simd,[AC_HELP_STRING([--disable-simd], [ Disable vectorized instructions: SSE4.1, AVX, AVX2])])
	AS_IF([ test  "x$enable_simd" != "xno" ],
	[
		AS_ECHO("SIMD enabled")
		arch=`echo $target | cut -d"-" -f1`
		# if we are on a x86 (32 or 64 bits) with gcc>=4.8 then run the AX_CHECK_X86_FEATURES macro
		AS_IF([test "x$arch" = "xx86_64" -o "x$arch" = "xi686"],
			    [archx86="yes"],
			    [archx86="no"]
		     )
		AS_IF([ test  "x${CCNAM:0:3}" != "xgcc" -o "x$archx86" = "xno" ],
		[
		   CUSTOM_SIMD="yes"
		   echo "Compiling with $CCNAM for a $arch target: running custom checks for SSE4.1 and AVX1,2"
		   AC_MSG_CHECKING(for SSE 4.1)
		   BACKUP_CXXFLAGS=${CXXFLAGS}
		   SSEFLAGS="-msse4.1"
		   CXXFLAGS="${BACKUP_CXXFLAGS} ${SSEFLAGS}"
		   CODE_SSE=`cat macros/CodeChunk/sse.C`
		   AC_TRY_RUN([ ${CODE_SSE} ],
			      [ sse_found="yes" ],
			       [ sse_found="no" ],
			       [ 
			       echo "cross compiling...disabling"
				 sse_found="no"
			       ])
	           AS_IF([ test "x$sse_found" = "xyes" ],
		   [
			AC_SUBST(SSEFLAGS)
			AC_MSG_RESULT(yes)
                   ],
		   [
			SSEFLAGS=""
			AC_MSG_RESULT(no)
		   ])
		   CXXFLAGS=${BACKUP_CXXFLAGS}
		   
		   dnl Check for AVX
		   AC_MSG_CHECKING(for AVX)
		   CODE_AVX=`cat macros/CodeChunk/avx.C`
		   dnl Intel compilers usually do not require option to enable avx
		   dnl Thus, we test with no option on
		   for switch_avxflags in "" "-mavx"; do
		       CXXFLAGS="${BACKUP_CXXFLAGS} -O0 ${switch_avxflags}"
		       AC_TRY_RUN([ ${CODE_AVX} ],
		       [
				avx_found="yes"
		        	AVXFLAGS=${switch_avxflags}
				break
		       ],
		       [ avx_found="no" ],
		       [
		        echo "cross compiling...disabling"
		        avx_found="no"
		        break
		       ])
		   done
			
		   dnl Is AVX found?
		   AS_IF([ test "x$avx_found" = "xyes" ],
		   [
			AC_MSG_RESULT(yes)
	                dnl Check for AVX2
			AC_MSG_CHECKING(for AVX2)
			for switch_avx2flags in "" "-mfma -mavx2"; do
			    CXXFLAGS="${BACKUP_CXXFLAGS} -O0 ${switch_avx2flags}"
			    AC_TRY_RUN(
			    [
			        #define __try_avx2
				${CODE_AVX}
			    ],
			    [
			        avx2_found="yes"
			        AVX2FLAGS="${switch_avx2flags}"
			        break
		            ],
			    [ avx2_found="no" ],
			    [
			        echo "cross compiling...disabling"
			        avx2_found = "no"
			        break
			    ])
			done
				
	                dnl Is AVX2 found?
			AS_IF([ test "x$avx2_found" = "xyes" ],
			[
				AC_MSG_RESULT(yes)
				AVXFLAGS=${AVX2FLAGS}
			],
			[ AC_MSG_RESULT(no) ]
			)
		    ],
		    [
			dnl No AVX
		    	AC_MSG_RESULT(no)
		    ])
		
		    CXXFLAGS=${BACKUP_CXXFLAGS}
		],
		[ ])
	],[ AS_ECHO("SIMD disabled")
	    CUSTOM_SIMD="yes" ])
])
