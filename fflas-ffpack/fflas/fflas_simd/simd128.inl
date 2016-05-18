/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   Bastien Vialla<bastien.vialla@lirmm.fr>
 * Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *
 *
 * ========LICENCE========
 * This file is part of the library FFLAS-FFPACK.
 *
 * FFLAS-FFPACK is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd128_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd128_INL

struct Simd128_base {

	/*
	* alias to 128 bit simd register
	*/
	using vect_t = __m128i;

	/*
	*  Return vector of type vect_t with all elements set to zero
	*  Return [0,0,0,0] int32_t
	*/
	static INLINE CONST vect_t zero() { return _mm_setzero_si128(); }

	/*
	* Compute the bitwise AND of packed 32-bits integer in a and b, and store the results in vect_t.
	* Args   : [a0, a1, a2, a3]	int32_t
	*	   [b0, b1, b2, b3]	int32_t
	* Return : [a0 AND b0, a1 AND b1, a2 AND b2, a3 AND b3]	int32_t
	*/
	static INLINE CONST vect_t vand(const vect_t a, const vect_t b) { return _mm_and_si128(b, a); }

	/*
	* Compute the bitwise OR of packed 32-bits integer in a and b, and store the results in vect_t.
	* Args   : [a0, a1, a2, a3]	int32_t
	*	   [b0, b1, b2, b3]	int32_t
	* Return : [a0 OR b0, a1 OR b1, a2 OR b2, a3 OR b3]	int32_t
	*/
	static INLINE CONST vect_t vor(const vect_t a, const vect_t b) { return _mm_or_si128(b, a); }

	/*
	* Compute the bitwise XOR of packed 32-bits integer in a and b, and store the results in vect_t.
	* Args   : [a0, a1, a2, a3]	int32_t
	*	   [b0, b1, b2, b3]	int32_t
	* Return : [a0 XOR b0, a1 XOR b1, a2 XOR b2, a3 XOR b3]	int32_t
	*/
	static INLINE CONST vect_t vxor(const vect_t a, const vect_t b) { return _mm_xor_si128(b, a); }

	/*
	* Compute the bitwise AND NOT of packed 32-bits integer in a and b, and store the results in vect_t.
	* Args   : [a0, a1, a2, a3]	int32_t
	*	   [b0, b1, b2, b3]	int32_t
	* Return : [a0 ANDNOT b0, a1 ANDNOT b1, a2 ANDNOT b2, a3 ANDNOT b3]	int32_t
	*/
	static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) { return _mm_andnot_si128(b, a); }

};

template <bool ArithType, bool Int, bool Signed, int Size> struct Simd128_impl;

#include "simd128_float.inl"
#include "simd128_double.inl"

#ifdef SIMD_INT
// Trop d'instructions SSE manquantes pour les int8_t

#include "simd128_int16.inl"
#include "simd128_int32.inl"
#include "simd128_int64.inl"

#endif //#ifdef SIMD_INT

template <class T>
using Simd128 =
    Simd128_impl<std::is_arithmetic<T>::value, std::is_integral<T>::value, std::is_signed<T>::value, sizeof(T)>;

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd128_INL
