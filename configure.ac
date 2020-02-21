# Copyright (c) 2011 FFLAS-FFPACK
# written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
# adapted from LinBox configuration
#
# ========LICENCE========
# This file is part of the library FFLAS-FFPACK.
#
# FFLAS-FFPACK is free software: you can redistribute it and/or modify
# it under the terms of the  GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========
#/

AC_PREREQ([2.67])


AC_INIT([FFLAS-FFPACK], [2.4.3],[ffpack-devel@googlegroups.com],[fflas-ffpack],
		[https://github.com/linbox-team/fflas-ffpack])


AC_CONFIG_MACRO_DIR([macros])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADERS([config.h])

AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([1.10 -Wall -Wno-portability foreign])
AX_PREFIX_CONFIG_H(fflas-ffpack/config.h, __FFLASFFPACK)

AC_LANG([C++])

echo "-----------------------------------------------"

AC_DEBUG
AC_PROFILE
AC_WARNINGS

echo "-----------------------------------------------"

# CFLAGS=${CFLAGS:-$DEFAULT_CFLAGS}
# CXXFLAGS=${CXXFLAGS:-$DEFAULT_CXXFLAGS}

######################################################
# Try and pass different flags according to compiler #
######################################################


# disable default -g -O2 CXXFLAGS
: ${CXXFLAGS=""}

#set CXX
AC_PROG_CXX

AC_COMPILER_NAME

# We need a C++11 compiler now - AB 2014-12-12
# clang-3.8 does not support __float128 without explicitly passing it -std=c++11
if test "x${CCNAM}" = "xclang38" ; then
   AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])
else
   AX_CXX_COMPILE_STDCXX_11([ext],[mandatory])
fi
# appending CXX11FLAGS for the remaining checks, to avoid clang __float128 undefined bug
# duplicates will be removed later on
CXXFLAGS="$CXXFLAGS $CXX11FLAGS"

REQUIRED_FLAGS="$CXX11FLAGS"


DEBUG_CFLAGS="-g"
DEFAULT_CFLAGS=""
WARN_CFLAGS="-Wall"

#TODO use -fast for icc, -ipa for eko...
if test "x$USE_DEBUG" = "xyes" ; then
	DEFAULT_CFLAGS="-O0 ${DEFAULT_CFLAGS} ${DEBUG_CFLAGS}" #those are CXXFLAGS
else
	DEFAULT_CFLAGS="-O2 ${DEFAULT_CFLAGS}"
fi

if test "x$PROF" = "xyes" ; then
	DEFAULT_CFLAGS="${DEFAULT_CFLAGS} -pg"
fi

if test "x$WARN" = "xyes" -o "x$WARN" = "xfull" ; then
	case x${CCNAM} in
		xicc)
		WARN_CFLAGS="${WARN_CFLAGS} -Wcheck"
		WARN_CFLAGS="${WARN_CFLAGS} -Wall -Wno-unused-parameter -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith  -Wwrite-strings -Wno-long-long"
		WARN_CFLAGS="${WARN_CFLAGS} -Wextra -ansi"
			;;
		xeko)
			WARN_CFLAGS="${WARN_CFLAGS} -Wno-unused-parameter"
			;;
		xgcc|xgcc44)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long  -Wno-variadic-macros -Wno-vla"
			fi
			if test "x${HAVE_CXX11}" = "x0" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -ansi"
			fi
			;;
		xgcc48)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long  -Wno-variadic-macros -Wno-vla"
				# WARN_CFLAGS="${WARN_CFLAGS} -fsanitize=address"
			fi
			if test "x${HAVE_CXX11}" = "x0" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -ansi"
			fi

			;;

		xclang)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long -Wno-vla-extension -Wno-variadic-macros"
				WARN_CFLAGS="${WARN_CFLAGS} -D__STRICT_ANSI__"
			fi
			;;
		xclang38)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long  -g -Wno-vla-extension -Wno-variadic-macros"
				WARN_CFLAGS="${WARN_CFLAGS} -D__STRICT_ANSI__"
				# WARN_CFLAGS="${WARN_CFLAGS} -fsanitize=address"
			fi
			;;

		*)
			echo
			echo "*******************************************************"
			echo "unsupported compiler ($CCNAM). Please file a bug."
			echo "*******************************************************"
			echo
			WARN_CFLAGS="${WARN_CFLAGS}"
	esac
fi


DEFAULT_CFLAGS="${DEFAULT_CFLAGS} ${WARN_CFLAGS} ${DEBUG_CFLAGS}"
AC_SUBST([DEFAULT_CFLAGS])

AC_HEADER_STDC
AC_PROG_LIBTOOL
AC_PROG_EGREP
AC_PROG_SED
# newer libtool...
LT_PREREQ([2.4.3])
LT_INIT


echo "DEFAULT_CFLAGS=${DEFAULT_CFLAGS}"
echo "DEBUG_CFLAGS=${DEBUG_CFLAGS}"
echo "TESTS_CFLAGS=${TESTS_CFLAGS}"
echo "-----------------------------------------------"
echo "     START  FFLAS-FFPACK CONFIG                "
echo "-----------------------------------------------"

FF_CHECK_OMP


# Checks which SIMD optimization flags to use and set SIMD_CFLAGS accordingly
INSTR_SET
AC_SUBST(SIMD_FLAGS)

dnl gcc-4.9.2 bug See https://trac.sagemath.org/ticket/17635#comment:178
AS_IF([ test  "x$CCNAM" = "xgcc492" ],[REQUIRED_FLAGS="${REQUIRED_FLAGS} -fpermissive"],[])

dnl With GCC 4.x, the default ABI version is 2. With this version, __m128 and
dnl __m256 are the same types and therefore we cannot have overloads for both
dnl types without linking error. It is fixed in ABI version 4.
dnl FIXME: Why do we set ABI version to 6 ?
AS_IF([test "x$CCNAM" = "xgccl5"],[REQUIRED_FLAGS+=" -fabi-version=6"])


AC_SUBST(REQUIRED_FLAGS)

PARFLAGS="${OMPFLAGS}"
PARLIBS="${OMPFLAGS}"

AC_SUBST(PARFLAGS)
AC_SUBST(PARLIBS)

echo "-----------------------------------------------"

# Machine characteristics

AC_CHECK_SIZEOF(char, 8)
AC_CHECK_SIZEOF(short, 16)
AC_CHECK_SIZEOF(int, 32)
AC_CHECK_SIZEOF(long, 32)
AC_CHECK_SIZEOF(long long, 64)
AC_CHECK_SIZEOF(__int64, 64)

AC_CHECK_TYPE([__int128_t], [AC_DEFINE(HAVE_INT128, 1, [Define that compiler allows int128_t types])])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stddef.h stdlib.h string.h sys/time.h stdint.h pthread.h])


# check endianness of the architecture
AC_C_BIGENDIAN(
  [AC_DEFINE(HAVE_BIG_ENDIAN, 1, [Define that architecture uses big endian storage])],
  [AC_DEFINE(HAVE_LITTLE_ENDIAN, 1, [Define that architecture uses little endian storage])],
  [])

# Create some useful data types of fixed, known lengths


# AC_DEFINE_UNQUOTED(INT8, $LINBOX_INT8, Canonical 8-bit data type)
# AC_DEFINE_UNQUOTED(INT16, $LINBOX_INT16, Canonical 16-bit data type)
# AC_DEFINE_UNQUOTED(INT32, $LINBOX_INT32, Canonical 32-bit data type)
# AC_DEFINE_UNQUOTED(INT64, $LINBOX_INT64, Canonical 64-bit data type)

echo "-----------------------------------------------"
# Feature checks
FF_MISC

AC_LANG([C++])


echo "-----------------------------------------------"

# Getting GMP from Givaro - AB 2014-12-12
#FF_CHECK_GMP

PKG_CHECK_MODULES([GIVARO],[givaro >= 4.1.1])
FF_CHECK_GIVARO_USABILITY()


dnl FF_CHECK_GIVARO(,,[
dnl echo '*******************************************************************************'
dnl echo ' WARNING: GIVARO not found!'
dnl echo
dnl echo ' GIVARO library is required for some tests in this library.'
dnl echo ' Please make sure GIVARO is installed and specify its location with the'
dnl echo ' option --with-givaro=<prefix> when running configure.'
dnl echo ' Do not forget to set/export LD_LIBRARY_PATH if necessary.'
dnl echo '*******************************************************************************'
dnl exit 1
dnl ])

BLAS_FOUND=false

FF_CHECK_BLAS_CFLAGS
FF_CHECK_BLAS_LIBS
FF_CHECK_MKL
FF_CHECK_USER_BLAS
FF_CHECK_USER_LAPACK

FF_OPENBLAS_NUM_THREADS

# FF_CHECK_BLAS

# FF_CHECK_GOTOBLAS

# FF_CHECK_GSL

# if test "$BLAS_FOUND" = "false" ; then
	# FF_CHECK_CBLAS
# fi

# if test "$BLAS_FOUND" = "false" ; then
	# FF_CHECK_OTHERBLAS
# fi

# FF_CHECK_LAPACK

# if test "$BLAS_FOUND" = "false" ; then
	# FF_CHECK_BLAS2
# fi



# BLAS_LIBS="${BLAS_LIBS}"
# BLAS_LIBS="-L/${BLAS_PATH} ${LAPACK_LIBS} ${BLAS_LIBS}"
# AC_SUBST(BLAS_LIBS)

# FF_CHECK_CUDA

# AM_CONDITIONAL(FFLASFFPACK_HAVE_BLAS,  test "x$BLAS_FOUND" != "xfalse")


# FF_BENCH


FF_DOC


# if test ! -d ./benchmarks/data  ; then
	# echo "Creating data dir in benchmark" ;
	# mkdir ./benchmarks/data  ;
# fi

CXXFLAGS="${CXXFLAGS} ${REQUIRED_FLAGS}"
CXXFLAGS="${CXXFLAGS} ${SIMD_CFLAGS}"
CXXFLAGS="${CXXFLAGS} ${GIVARO_CFLAGS}"
CXXFLAGS="${CXXFLAGS} ${BLAS_CFLAGS}"
CXXFLAGS="${CXXFLAGS} ${CUDA_CFLAGS}"
CXXFLAGS="${CXXFLAGS} ${PARFLAGS}"

AC_SUBST(CXXFLAGS)

FF_PRECOMPILE

echo "-----------------------------------------------"
echo "          END FFLAS-FFPACK CONFIG              "
echo "-----------------------------------------------"

AC_CONFIG_FILES([
Makefile
macros/Makefile
macros/CodeChunk/Makefile
fflas-ffpack-config
fflas-ffpack/Makefile
fflas-ffpack/fflas/Makefile
fflas-ffpack/fflas/fflas_fgemm/Makefile
fflas-ffpack/fflas/fflas_sparse/Makefile
fflas-ffpack/fflas/fflas_sparse/coo/Makefile
fflas-ffpack/fflas/fflas_sparse/csr/Makefile
fflas-ffpack/fflas/fflas_sparse/ell/Makefile
fflas-ffpack/fflas/fflas_sparse/ell_simd/Makefile
fflas-ffpack/fflas/fflas_sparse/csr_hyb/Makefile
fflas-ffpack/fflas/fflas_sparse/sell/Makefile
fflas-ffpack/fflas/fflas_sparse/hyb_zo/Makefile
fflas-ffpack/fflas/fflas_igemm/Makefile
fflas-ffpack/fflas/fflas_simd/Makefile
fflas-ffpack/ffpack/Makefile
fflas-ffpack/field/Makefile
fflas-ffpack/utils/Makefile
fflas-ffpack/paladin/Makefile
fflas-ffpack/interfaces/Makefile
fflas-ffpack/interfaces/libs/Makefile
fflas-ffpack/checkers/Makefile
autotune/Makefile
doc/Makefile
tests/Makefile
tests/data/Makefile
benchmarks/Makefile
examples/Makefile
tutorials/Makefile
fflas-ffpack.pc
],[chmod +x fflas-ffpack-config])
AC_OUTPUT

echo "" > fflas-ffpack/fflas-ffpack-thresholds.h
