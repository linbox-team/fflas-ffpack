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

struct Simd128i_base {

    /*
     * alias to 128 bit simd register
     */
    using vect_t = __m128i;

    /*
     *  Return vector of type vect_t with all elements set to zero
     *  Return [0, ...,0]
     */
    static INLINE CONST vect_t zero() { return _mm_setzero_si128(); }

    /*
     * Shift packed 128-bit integers in a left by s bits while shifting in zeros, and store the results in vect_t.
     * Args   : [a0] int128_t
     * Return : [a0 << (s*8)] int128_t
     */
    template<uint8_t s>
    static INLINE CONST vect_t sll128(const vect_t a) { return _mm_slli_si128(a, s); }

    /*
     * Shift packed 128-bit integers in a right by s while shifting in zeros, and store the results in vect_t.
     * Args   : [a0] int128_t
     * Return : [a0 >> (s*8)] int128_t
     */
    template<uint8_t s>
    static INLINE CONST vect_t srl128(const vect_t a) { return _mm_srli_si128(a, s); }

    /*
     * Compute the bitwise AND and store the results in vect_t.
     * Args   : [a0, ..., a127]
     *		   [b0, ..., b127]
     * Return : [a0 AND b0, ..., a127 AND b127]
     */
    static INLINE CONST vect_t vand(const vect_t a, const vect_t b) { return _mm_and_si128(b, a); }

    /*
     * Compute the bitwise OR and store the results in vect_t.
     * Args   : [a0, ..., a127]
     *		   [b0, ..., b127]
     * Return : [a0 OR b0, ..., a127 OR b127]
     */
    static INLINE CONST vect_t vor(const vect_t a, const vect_t b) { return _mm_or_si128(b, a); }

    /*
     * Compute the bitwise XOR and store the results in vect_t.
     * Args   : [a0, ..., a127]
     *		   [b0, ..., b127]
     * Return : [a0 XOR b0, ..., a127 XOR b127]
     */
    static INLINE CONST vect_t vxor(const vect_t a, const vect_t b) { return _mm_xor_si128(b, a); }

    /*
     * Compute the bitwise NOT AND and store the results in vect_t.
     * Args   : [a0, ..., a127]
     *		   [b0, ..., b127]
     * Return : [NOT(a0) AND b0, ..., NOT(a127) AND b127]
     */
    static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) { return _mm_andnot_si128(a, b); }

};

template <bool ArithType, bool Int, bool Signed, int Size> struct Simd128_impl;

template <class T>
using Simd128 =
Simd128_impl<std::is_arithmetic<T>::value, std::is_integral<T>::value, std::is_signed<T>::value, sizeof(T)>;

#include "simd128_float.inl"
#include "simd128_double.inl"

#ifdef SIMD_INT
// Trop d'instructions SSE manquantes pour les int8_t

#include "simd128_int16.inl"
#include "simd128_int32.inl"
#ifdef __x86_64__
#include "simd128_int64.inl"
#endif
#endif //#ifdef SIMD_INT

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd128_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
