/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2012 FFLAS-FFPACK
 * Written by <brice.boyer@imag.fr>
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

#include "fflas-ffpack/config.h"
#ifdef __FFLASFFPACK_USE_OPENMP
#  ifndef __GIVARO_USE_OPENMP
#    define __GIVARO_USE_OPENMP 1
#  endif
#endif

#include "fflas-ffpack/fflas-ffpack-optimise.h"

#if defined(__FFLASFFPACK_USE_SSE) or defined(__FFLASFFPACK_USE_AVX) or defined(__FFLASFFPACK_USE_AVX2)
#define __FFLASFFPACK_USE_SIMD // see configure...
#endif

#ifdef __CYGWIN__
#  define _GLIBCXX_USE_C99
#  ifndef _GLIBCXX_USE_C99_MATH_TR1
#    include <math.h>
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

  unsigned long long stoull( const std::string& str) { return std::strtoull(str.c_str(), 0, 10); }
  long long stoll( const std::string& str) { return std::strtoll(str.c_str(), 0, 10); }

}
#  endif
#endif

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

#ifdef __x86_64__
#if defined(__GNUC__) || defined (__clang__) /* who supports __int128_t ? */
#define int128_t __int128_t
#define uint128_t unsigned __int128_t
#else /* hopefully this exists */
#define int128_t __int128
#define uint128_t unsigned __int128
#endif /* __int128_t */
#endif /* __x86_64__ */

#endif // __FFLASFFPACK_fflas_ffpack_configuration_H
