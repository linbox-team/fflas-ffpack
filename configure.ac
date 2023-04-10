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


AC_INIT([FFLAS-FFPACK], [2.5.1],[ffpack-devel@googlegroups.com],[fflas-ffpack],
		[https://github.com/linbox-team/fflas-ffpack])


AC_CONFIG_MACRO_DIR([macros])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADERS([config.h])

AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([1.10 -Wall -Wno-portability foreign])
AX_PREFIX_CONFIG_H(fflas-ffpack/config.h, __FFLASFFPACK)

# Set CXXFLAGS to an empty string if not defined, and save it. We must do it
# before calling AC_PROG_CXX that sets it to '-O2 -g' if not defined
: ${CXXFLAGS=""}
ORIGINAL_CXXFLAGS="${CXXFLAGS}"

# We set the language to C++
AC_LANG([C++])
AC_PROG_CXX

#
AM_PROG_AR

# Libtool
AC_PROG_LIBTOOL
LT_PREREQ([2.4.3])
LT_INIT

# Look for headers
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stddef.h stdlib.h string.h sys/time.h stdint.h pthread.h])

#
AC_PROG_EGREP
AC_PROG_SED

#################################################
AS_BOX([       START FFLAS-FFPACK CONFIG       ])
#################################################

AC_COMPILER_NAME

# We need a C++11 compiler now - AB 2014-12-12
# clang-3.8 does not support __float128 without explicitly passing it -std=c++11
AS_IF([test "x${CCNAM}" = "xclang38"],
      [AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])],
      [AX_CXX_COMPILE_STDCXX_11([ext],[mandatory])]
      )
REQUIRED_FLAGS="${CXX11FLAGS}"
# Add the c++11 flags for the configure compilations as clang needs it to work properly with a glibc++ compiled with gcc
CXXFLAGS="${CXX11FLAGS} ${CXXFLAGS}"

AS_ECHO([---------------------------------------])

# Set OPTIM_FLAGS, DEBUG_FLAGS depending on compiler and command line arguments
SET_FLAGS

# Append -march=native or -mcpu=native (if recognized by the compiler) to
# OPTIM_FLAGS if not present in CXXFLAGS and not cross-compiling and
# --without-archnative is not set
ARCH_FLAGS

# Append -mfpmath=sse to OPTIM_FLAGS on i386 and i686 architecture with SSE
FPMATH_FLAGS

AS_ECHO([---------------------------------------])
# Machine characteristics

# Size of some types
AC_CHECK_SIZEOF(char)
AC_CHECK_SIZEOF(short)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(long long)
AC_CHECK_SIZEOF(__int64_t)

# Looking for int128
AC_CHECK_TYPE([__int128_t], [AC_DEFINE(HAVE_INT128, 1, [Define that compiler allows int128_t types])])


# check endianness of the architecture
AC_C_BIGENDIAN(
  [AC_DEFINE(HAVE_BIG_ENDIAN, 1, [Define that architecture uses big endian storage])],
  [AC_DEFINE(HAVE_LITTLE_ENDIAN, 1, [Define that architecture uses little endian storage])],
  [])

AS_ECHO([---------------------------------------])

# Looking for OpenMP
AC_OPENMP
FF_CHECK_OMP

PARFLAGS="${OPENMP_CXXFLAGS}"
PARLIBS="${OPENMP_CXXFLAGS}"
AC_SUBST(PARFLAGS)
AC_SUBST(PARLIBS)


# Feature checks
FF_MISC


# Looking for Givaro. We get the flags for GMP at the same time
PKG_CHECK_MODULES([GIVARO],[givaro >= 4.1.2])
FF_CHECK_GIVARO_USABILITY()

AS_ECHO([---------------------------------------])

BLAS_FOUND=false

FF_CHECK_BLAS_CFLAGS
FF_CHECK_BLAS_LIBS
FF_CHECK_MKL
FF_CHECK_USER_BLAS
FF_CHECK_USER_LAPACK

FF_OPENBLAS_NUM_THREADS

# FF_CHECK_CUDA

AS_ECHO([---------------------------------------])

FF_DOC
FF_PRECOMPILE

AS_ECHO([---------------------------------------])

CXXFLAGS="${ORIGINAL_CXXFLAGS}"
FFLASFFPACK_CXXFLAGS="${REQUIRED_FLAGS} ${OPTIM_FLAGS} ${DEBUG_FLAGS}"
AC_SUBST(FFLASFFPACK_CXXFLAGS)
AS_ECHO(["FFLASFFPACK_CXXFLAGS = ${FFLASFFPACK_CXXFLAGS}"])
AC_SUBST(REQUIRED_FLAGS)

#################################################
AS_BOX([        END FFLAS-FFPACK CONFIG        ])
#################################################

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
