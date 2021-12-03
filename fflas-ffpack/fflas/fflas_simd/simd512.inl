/*
 * Copyright (C) 2018 the FFLAS-FFPACK group
 *
 * Written by   Ozturk Ozden<ozden.ozturk@etu.univ-grenoble-alpes.fr>
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


#ifndef __FFLASFFPACK_simd512_INL
#define __FFLASFFPACK_simd512_INL

struct Simd512i_base {

    /*
     * alias to 512 bit simd register
     */
    using vect_t = __m512i;

    /*
     *  Return vector of type vect_t with all elements set to zero
     *  Return [0, ...,0]
     */
    static INLINE CONST vect_t zero() { return _mm512_setzero_si512(); }

    /*
     * Compute the bitwise OR and store the results in vect_t.
     * Args   : [a0, ..., a511]
     *		   [b0, ..., b511]
     * Return : [a0 OR b0, ..., a511 OR b511]
     */
    static INLINE CONST vect_t vor(const vect_t a, const vect_t b) { return _mm512_or_si512(b, a); }

    /*
     * Compute the bitwise XOR and store the results in vect_t.
     * Args   : [a0, ..., a511]
     *		   [b0, ..., b511]
     * Return : [a0 XOR b0, ..., a511 XOR b511]
     */
    static INLINE CONST vect_t vxor(const vect_t a, const vect_t b) { return _mm512_xor_si512(b, a); }

    /*
     * Compute the bitwise AND and store the results in vect_t.
     * Args   : [a0, ..., a511]
     *		   [b0, ..., b511]
     * Return : [a0 AND b0, ..., a511 AND b511]
     */
    static INLINE CONST vect_t vand(const vect_t a, const vect_t b) { return _mm512_and_si512(b, a); }

    /*
     * Compute the bitwise NOT AND and store the results in vect_t.
     * Args   : [a0, ..., a511]
     *		   [b0, ..., b511]
     * Return : [(NOT a0) AND b0, ..., (NOT a511) AND b511]
     */
    static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) { return _mm512_andnot_si512(a, b); }


};

template <bool ArithType, bool Int, bool Signed, int Size> struct Simd512_impl;

template <class T>
using Simd512 =
Simd512_impl<std::is_arithmetic<T>::value, std::is_integral<T>::value, std::is_signed<T>::value, sizeof(T)>;

#include "simd512_float.inl"
#include "simd512_double.inl"
#include "simd512_int64.inl"

#endif // __FFLASFFPACK_simd512_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
