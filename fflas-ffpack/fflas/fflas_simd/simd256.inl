/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

struct Simd256i_base {

	/*
	* alias to 256 bit simd register
	*/
	using vect_t = __m256i;

	/*
	*  Return vector of type vect_t with all elements set to zero
	*  Return [0, ...,0]
	*/
	static INLINE CONST vect_t zero() { return _mm256_setzero_si256(); }

	/*
	* Compute the bitwise AND and store the results in vect_t.
	* Args   : [a0, ..., a255]
	*		   [b0, ..., b255]
	* Return : [a0 AND b0, ..., a255 AND b255]
	*/
	static INLINE CONST vect_t vand(const vect_t a, const vect_t b) { return _mm256_and_si256(b, a); }

	/*
	* Compute the bitwise OR and store the results in vect_t.
	* Args   : [a0, ..., a255]
	*		   [b0, ..., b255]
	* Return : [a0 OR b0, ..., a255 OR b255]
	*/
	static INLINE CONST vect_t vor(const vect_t a, const vect_t b) { return _mm256_or_si256(b, a); }

	/*
	* Compute the bitwise XOR and store the results in vect_t.
	* Args   : [a0, ..., a255]
	*		   [b0, ..., b255]
	* Return : [a0 XOR b0, ..., a255 XOR b255]
	*/
	static INLINE CONST vect_t vxor(const vect_t a, const vect_t b) { return _mm256_xor_si256(b, a); }

	/*
	* Compute the bitwise AND NOT and store the results in vect_t.
	* Args   : [a0, ..., a255]
	*		   [b0, ..., b255]
	* Return : [a0 AND (NOT b0), ..., a255 AND (NOT b255)]
	*/
	static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) { return _mm256_andnot_si256(b, a); }

	/*
	* Shuffle 128-bit integers in a and b using the control in imm8, and store the results in dst.
	* Args   :	[a0, a1] int128_t
	*			[b0, b1] int128_t
	* Return : [s[0..3]?a0:a1:b0:b1, s[4..7]?a0:a1:b0:b1] int128_t
	*/
	static INLINE CONST vect_t shuffle128(const vect_t a, const vect_t b, const int s) {
		static_assert(__builtin_constant_p(s),"Index s has to be a constant expression");
		return _mm256_permute2x128_si256(a, b, __builtin_constant_p(s)?s:0);
	}

	/*
	* Unpack and interleave 128-bit integers from the low half of a and b, and store the results in dst.
	* Args   : [a0, a1] int128_t
			   [b0, b1] int128_t
	* Return : [a0, b0] int128_t
	*/
	static INLINE CONST vect_t unpacklo128(const vect_t a, const vect_t b) { return shuffle128(a, b, 0x20); }

	/*
	* Unpack and interleave 128-bit integers from the high half of a and b, and store the results in dst.
	* Args   : [a0, a1] int128_t
			   [b0, b1] int128_t
	* Return : [a1, b1] int128_t
	*/
	static INLINE CONST vect_t unpackhi128(const vect_t a, const vect_t b) { return shuffle128(a, b, 0x31); }

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
