/* Copyright (C) 2012 FFLAS-FFPACK
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

/*! @file fflas-ffpack/fflas-ffpack-config.h
 * @ingroup optimise
 * @brief Defaults for optimised values.
 * While \c fflas-ffpack-optimise.h is created by \c configure script,
 * (either left blank or filled by optimiser), this file produces the
 * defaults for the optimised values. If \c fflas-ffpack-optimise.h is not
 * empty, then its values preceeds the defaults here.
 */


#ifndef __FFLASFFPACK_fflas_ffpack_configuration_H
#define __FFLASFFPACK_fflas_ffpack_configuration_H

#ifndef GCC_VERSION
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif

#if defined __GNUC__ && __GNUC__>=6
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

#ifdef __CYGWIN__
# ifndef _GLIBCXX_USE_C99
#  define _GLIBCXX_USE_C99 true
#  ifndef _GLIBCXX_USE_C99_MATH_TR1
#    include <cstdlib>
#    include <string>
#    include <cmath>
#    undef fma
#    include <stdlib.h>
#    undef strtoull
#    undef strtoll
namespace std _GLIBCXX_VISIBILITY(default)
{
    _GLIBCXX_BEGIN_NAMESPACE_VERSION
    using ::fma;
    using ::strtoll;
    using ::strtoull;

    /*
       unsigned long      stoul( const std::string& str, std::size_t* pos = 0, int base = 10 ) {
       return std::strtoul(str.c_str(), NULL, base);
       }

       unsigned long long stoull( const std::string& str, std::size_t* pos = 0, int base = 10 ) {
       return std::strtoull(str.c_str(), NULL, base);
       }

       long      stol( const std::string& str, std::size_t* pos = 0, int base = 10 ) {
       return std::strtol(str.c_str(), NULL, base);
       }

       long long stoll( const std::string& str, std::size_t* pos = 0, int base = 10 ) {
       return std::strtoll(str.c_str(), NULL, base);
       }
       */

}
#  else
#    define _GLIBCXX_USE_C99 true
#    include <cstdlib>
#  endif
# endif
#endif

#include "fflas-ffpack/config.h"
#ifdef __FFLASFFPACK_USE_OPENMP
#  ifndef __GIVARO_USE_OPENMP
#    define __GIVARO_USE_OPENMP 1
#  endif
#endif

/* include definitions of thresholds set by the autotuning first */
#include "fflas-ffpack/fflas-ffpack-thresholds.h"

/* then include the default definitions */
#include "fflas-ffpack/fflas-ffpack-default-thresholds.h"

#if defined(_OPENMP) || defined(OMP_H) || defined(__OMP_H) || defined(__pmp_omp_h)
#ifndef __FFLASFFPACK_USE_OPENMP
#warning "openmp was not detected correctly at configure time, please report this bug"
#define __FFLASFFPACK_USE_OPENMP
#endif
#endif

#include "givaro/givconfig.h"

/* Define if sse instructions are supported */
#ifdef __SSE__
#define __FFLASFFPACK_HAVE_SSE_INSTRUCTIONS  1
#endif

/* Define if sse2 instructions are supported */
#ifdef __SSE2__
#define __FFLASFFPACK_HAVE_SSE2_INSTRUCTIONS  1
#endif

/* Define if sse3 instructions are supported */
#ifdef __SSE3__
#define __FFLASFFPACK_HAVE_SSE3_INSTRUCTIONS  1
#endif

/* Define if sse4.1 instructions are supported */
#ifdef __SSE4_1__
#define __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS  1
#endif

/* Define if sse4.2 instructions are supported */
#ifdef __SSE4_2__
#define __FFLASFFPACK_HAVE_SSE4_2_INSTRUCTIONS  1
#endif

/* 256 SIMD registers are not supported by gcc on CYGWIN
 * See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=54412
 */
#if not defined(__CYGWIN__) or not defined(__GNUC__)

/* Define if avx instructions are supported */
#ifdef __AVX__
#define __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS  1
#endif

/* Define if avx2 instructions are supported */
#ifdef __AVX2__
#define __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS  1
#endif

/* Define if avx512f instructions are supported */
#ifdef __AVX512F__
#define __FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS  1
#endif

/* Define if avx512dq instructions are supported */
#ifdef __AVX512DQ__
#define __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS 1
#endif

#endif // CYGWIN and GCC

/* Define if fma instructions are supported */
#ifdef __FMA__
#define __FFLASFFPACK_HAVE_FMA_INSTRUCTIONS  1
#endif

#ifdef __FFLASFFPACK_HAVE_INT128
#define int128_t __int128_t
#define uint128_t __uint128_t
#endif

#endif // __FFLASFFPACK_fflas_ffpack_configuration_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
