/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

#include "fflas-ffpack/fflas-ffpack-thresholds.h"

// winograd algorithm threshold (for double)
#ifndef __FFLASFFPACK_WINOTHRESHOLD
#define __FFLASFFPACK_WINOTHRESHOLD 1000
#endif

#ifndef __FFLASFFPACK_WINOTHRESHOLD_FLT
#define __FFLASFFPACK_WINOTHRESHOLD_FLT 2000
#endif

#ifndef __FFLASFFPACK_WINOTHRESHOLD_BAL
#define __FFLASFFPACK_WINOTHRESHOLD_BAL 1000
#endif

#ifndef __FFLASFFPACK_WINOTHRESHOLD_BAL_FLT
#define __FFLASFFPACK_WINOTHRESHOLD_BAL_FLT 2000
#endif


#if defined(_OPENMP) || defined(OMP_H) || defined(__OMP_H) || defined(__pmp_omp_h)
#ifndef __FFLASFFPACK_USE_OPENMP
#warning "openmp was not detected correctly at configure time, please report this bug"
#define __FFLASFFPACK_USE_OPENMP
#endif
#endif

#include "givaro/givconfig.h"

/* Define if mmx instructions are supported */
#ifdef __GIVARO_HAVE_MMX_INSTRUCTIONS
#define __FFLASFFPACK_HAVE_MMX_INSTRUCTIONS  1
#endif

/* Define if popcnt instructions are supported */
#ifdef __GIVARO_HAVE_POPCNT_INSTRUCTIONS
#define __FFLASFFPACK_HAVE_POPCNT_INSTRUCTIONS  1
#endif

/* Define if sse instructions are supported */
#ifdef __GIVARO_HAVE_SSE_INSTRUCTIONS
#define __FFLASFFPACK_HAVE_SSE_INSTRUCTIONS  1
#endif

/* Define if sse2 instructions are supported */
#ifdef __GIVARO_HAVE_SSE2_INSTRUCTIONS
#define __FFLASFFPACK_HAVE_SSE2_INSTRUCTIONS  1
#endif

/* Define if sse3 instructions are supported */
#ifdef __GIVARO_HAVE_SSE3_INSTRUCTIONS
#define __FFLASFFPACK_HAVE_SSE3_INSTRUCTIONS  1
#endif

/* Define if sse4a instructions are supported */
#ifdef __GIVARO_HAVE_SSE4A_INSTRUCTIONS
#define __FFLASFFPACK_HAVE_SSE4A_INSTRUCTIONS  1
#endif

/* Define if sse4.1 instructions are supported */
#ifdef __GIVARO_HAVE_SSE4_1_INSTRUCTIONS
#define __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS  1
#endif

/* Define if sse4.2 instructions are supported */
#ifdef __GIVARO_HAVE_SSE4_2_INSTRUCTIONS
#define __FFLASFFPACK_HAVE_SSE4_2_INSTRUCTIONS  1
#endif

/* Define if avx instructions are supported */
#ifdef __GIVARO_HAVE_AVX_INSTRUCTIONS
#define __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS  1
#endif

/* Define if avx2 instructions are supported */
#ifdef __GIVARO_HAVE_AVX2_INSTRUCTIONS
#define __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS  1
#endif

/* Define if avx512f instructions are supported */
#ifdef __GIVARO_HAVE_AVX512F_INSTRUCTIONS
#define __FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS  1
#endif

/* Define if bmi2 instructions are supported */
#ifdef __GIVARO_HAVE_BMI2_INSTRUCTIONS
#define __FFLASFFPACK_HAVE_BMI2_INSTRUCTIONS  1
#endif

/* Define if bmi instructions are supported */
#ifdef __GIVARO_HAVE_BMI_INSTRUCTIONS
#define __FFLASFFPACK_HAVE_BMI_INSTRUCTIONS  1
#endif

/* Define if fma4 instructions are supported */
#ifdef __GIVARO_HAVE_FMA4_INSTRUCTIONS
#define __FFLASFFPACK_HAVE_FMA4_INSTRUCTIONS  1
#endif

/* Define if fma instructions are supported */
#ifdef __GIVARO_HAVE_FMA_INSTRUCTIONS
#define __FFLASFFPACK_HAVE_FMA_INSTRUCTIONS  1
#endif

#ifdef __GIVARO_HAVE_INT128
#define __FFLASFFPACK_HAVE_INT128
#endif

#endif // __FFLASFFPACK_fflas_ffpack_configuration_H
