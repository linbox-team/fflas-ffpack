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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd256_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd256_INL

struct Simd256_base {

	/*
	* alias to 256 bit simd register
	*/
	using vect_t = __m256i;

	/*
	* Compute the bitwise AND of packed 16-bits integer in a and b, and store the results in vect_t.
	* Args   :	[a0, ..., a15]		int16_t
			[b0, ..., b15]		int16_t
	* Return : [a0 AND b0, a1 AND b1, a2 AND b2, a3 AND b3, a4 AND b4, a5 AND b5, a6 AND b6, a7 AND b7,
	a8 AND b8, a9 AND b9, a10 AND b10, a11 AND b11, a12 AND b12, a13 AND b13, a14 AND b14, a15 AND b15]
	*/
	static INLINE CONST vect_t vand(const vect_t a, const vect_t b) { return _mm256_and_si256(b, a); }

	/*
	* Compute the bitwise OR of packed 16-bits integer in a and b, and store the results in vect_t.
	* Args   :	[a0, ..., a15]		int16_t
			[b0, ..., b15]		int16_t
	* Return : [a0 OR b0, a1 OR b1, a2 OR b2, a3 OR b3, a4 OR b4, a5 OR b5, a6 OR b6, a7 OR b7,
	a8 OR b8, a9 OR b9, a10 OR b10, a11 OR b11, a12 OR b12, a13 OR b13, a14 OR b14, a15 OR b15]
	*/
	static INLINE CONST vect_t vor(const vect_t a, const vect_t b) { return _mm256_or_si256(b, a); }

	/*
	* Compute the bitwise XOR of packed 16-bits integer in a and b, and store the results in vect_t.
	* Args   :	[a0, ..., a15]		int16_t
			[b0, ..., b15]		int16_t
	* Return : [a0 XOR b0, a1 XOR b1, a2 XOR b2, a3 XOR b3, a4 XOR b4, a5 XOR b5, a6 XOR b6, a7 XOR b7,
	a8 XOR b8, a9 XOR b9, a10 XOR b10, a11 XOR b11, a12 XOR b12, a13 XOR b13, a14 XOR b14, a15 XOR b15]
	*/
	static INLINE CONST vect_t vxor(const vect_t a, const vect_t b) { return _mm256_xor_si256(b, a); }

	/*
	* Compute the bitwise AND NOT of packed 16-bits integer in a and b, and store the results in vect_t.
	* Args   :	[a0, ..., a15]		int16_t
			[b0, ..., b15]		int16_t
	* Return : [a0 ANDNOT b0, a1 ANDNOT b1, a2 ANDNOT b2, a3 ANDNOT b3, a4 ANDNOT b4, a5 ANDNOT b5, a6 ANDNOT b6, a7
	ANDNOT b7,
	a8 ANDNOT b8, a9 ANDNOT b9, a10 ANDNOT b10, a11 ANDNOT b11, a12 ANDNOT b12, a13 ANDNOT b13, a14 ANDNOT b14, a15
	ANDNOT b15]
	*/
	static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) { return _mm256_andnot_si256(b, a); }

};

template <bool ArithType, bool Int, bool Signed, int Size> struct Simd256_impl;

#include "simd256_float.inl"
#include "simd256_double.inl"

#ifdef SIMD_INT
// Trop d'instructions SSE manquantes pour les int8_t

#if defined(__FFLASFFPACK_USE_AVX2)
#include "simd256_int16.inl"
#include "simd256_int32.inl"
#include "simd256_int64.inl"
#endif

#endif //#ifdef SIMD_INT

template <class T>
using Simd256 =
    Simd256_impl<std::is_arithmetic<T>::value, std::is_integral<T>::value, std::is_signed<T>::value, sizeof(T)>;

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd256_INL
