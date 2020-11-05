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

struct Simd256fp_base {
#if defined(__FFLASFFPACK_HAVE_AVX_INSTRUCTIONS)

    /*
     * Shuffle 128-bits selected by imm8 from a and b, and store the results in dst.
     * Args   :	[a0, a1]
     *			[b0, b1]
     * Return : [s[0..3]?a0:a1:b0:b1, s[4..7]?a0:a1:b0:b1]
     */
    template<int s>
    static INLINE CONST __m256d permute128(const __m256d a, const __m256d b) {
        return _mm256_permute2f128_pd(a, b, s);
    }

    template<int s>
    static INLINE CONST __m256 permute128(const __m256 a, const __m256 b) {
        return _mm256_permute2f128_ps(a, b, s);
    }
#endif
};

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

#if defined(__FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS)

    // CLANG < 3.8 does not implement m256_bslli_epi128 nor _mmm256_bsrli_epi128
#if defined(__clang__)
#if  __clang_major < 3 || (__clang_major__  == 3 && __clang_minor__ < 8)
#define _mm256_bslli_epi128(a, count) _mm256_slli_si256((a), (count))
#define _mm256_bsrli_epi128(a, count) _mm256_srli_si256((a), (count))
#endif
#endif

    /*
     * Shift packed 128-bit integers in a left by s bits while shifting in zeros, and store the results in vect_t.
     * Args   : [a0, a1] int128_t
     * Return : [a0 << (s*8), a1 << (s*8)] int128_t
     */
    template<uint8_t s>
    static INLINE CONST vect_t sll128(const vect_t a) { return _mm256_bslli_epi128(a, s); }

    /*
     * Shift packed 128-bit integers in a right by s while shifting in zeros, and store the results in vect_t.
     * Args   : [a0, a1] int128_t
     * Return : [a0 << (s*8), a1 << (s*8)] int128_t
     */
    template<uint8_t s>
    static INLINE CONST vect_t srl128(const vect_t a) { return _mm256_bsrli_epi128(a, s); }

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
     * Compute the bitwise NOT AND and store the results in vect_t.
     * Args   : [a0, ..., a255]
     *		   [b0, ..., b255]
     * Return : [(NOT a0) AND b0, ..., (NOT a255) AND b255]
     */
    static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) { return _mm256_andnot_si256(a, b); }

    /*
     * Shuffle 128-bit integers in a and b using the control in imm8, and store the results in dst.
     * Args   :	[a0, a1] int128_t
     *			[b0, b1] int128_t
     * Return : [s[0..3]?a0:a1:b0:b1, s[4..7]?a0:a1:b0:b1] int128_t
     */
    template<int s>
    static INLINE CONST vect_t permute128(const vect_t a, const vect_t b) {
        return _mm256_permute2x128_si256(a, b, s);
    }
#endif
};

template <bool ArithType, bool Int, bool Signed, int Size> struct Simd256_impl;

template <class T>
using Simd256 =
Simd256_impl<std::is_arithmetic<T>::value, std::is_integral<T>::value, std::is_signed<T>::value, sizeof(T)>;

#include "simd256_float.inl"
#include "simd256_double.inl"

#ifdef SIMD_INT
// To many missing insctructions on int8_t

#if defined(__FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS)
#ifdef __x86_64__
#include "simd256_int64.inl"
#endif
#include "simd256_int32.inl"
#include "simd256_int16.inl"
#endif

#endif //#ifdef SIMD_INT


#endif // __FFLASFFPACK_fflas_ffpack_utils_simd256_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
