dnl Check for AVX
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

dnl FF_CHECK_AVX
dnl
dnl turn on AVX or AVX2 extensions if available

AC_DEFUN([FF_CHECK_AVX],
		[
		AC_ARG_ENABLE(avx,
			[AC_HELP_STRING([--enable-avx],
				[ Use Intel(r) AVX ])
			],
			[ avec_avx=$enable_avx ],
			[ avec_avx=yes ]
			)
		AC_MSG_CHECKING(for AVX)
		AS_IF([ test  "x$avec_avx" != "xno" ],
			[
			BACKUP_CXXFLAGS=${CXXFLAGS}
			AVXFLAGS="-mavx"
			CXXFLAGS="${BACKUP_CXXFLAGS} ${AVXFLAGS}"
			AC_TRY_RUN([
				#include <immintrin.h>
				int main() {
				register __m256d P ;
				double p = 0;
				P   = _mm256_set1_pd(p);
				return 0;
				}
				],
				[ avx_found="yes" ],
				[ avx_found="no" ],
				[
				echo "cross compiling...disabling"
				avx_found="no"
				])
			AS_IF([ test "x$avx_found" = "xyes" ],[
				AC_DEFINE(USE_AVX,1,[Define if AVX is available])
				AC_SUBST(AVXFLAGS)
				AC_MSG_RESULT(yes (AVX))
				AC_MSG_CHECKING(for AVX2)
				AVX2FLAGS="-mfma -mavx2"
				CXXFLAGS="${BACKUP_CXXFLAGS} ${AVX2FLAGS}"
				AC_TRY_RUN(
					[
					#include <immintrin.h>
					int main() {
					register __m256d P ;
					double p = 0;
					P = _mm256_set1_pd(p);
					P = _mm256_fnmadd_pd(P,P,P);
					return 0;
					}
					],
					[ avx2_found="yes" ],
					[ avx2_found="no" ],
					[
					echo "cross compiling...disabling"
					avx2_found = "no"
					])
				AS_IF([ test "x$avx2_found" = "xyes" ],[
					AC_DEFINE(USE_AVX2,1,[Define if AVX2 is available])
					AC_MSG_RESULT(yes (AVX2))
					AVXFLAGS=${AVX2FLAGS}
					AC_SUBST(AVXFLAGS)
					],
					[ AC_MSG_RESULT(no) ]
					)
					],
				[ AC_MSG_RESULT(no) ]
				)
				CXXFLAGS=${BACKUP_CXXFLAGS}
				],
				[ AC_MSG_RESULT(no) ]
				)
	])
