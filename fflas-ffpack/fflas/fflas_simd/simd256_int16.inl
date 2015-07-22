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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd256_int16_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd256_int16_INL

/*
 * Simd256 specialized for int16_t
 */
template <> struct Simd256_impl<true, true, true, 2> {
#if defined(__FFLASFFPACK_USE_AVX2)
    /*
     * alias to 256 bit simd register
     */
    using vect_t = __m256i;

    /*
     * alias to 256 bit simd register
     */
    using half_t = __m128i;

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = int16_t;

    /*
     * Simd128 for scalar_t, to deal half_t
     */
    using simdHalf = Simd128<scalar_t>;

    /*
     *  number of scalar_t in a simd register
     */
    static const constexpr size_t vect_size = 16;

    /*
     *  alignement required by scalar_t pointer to be loaded in a vect_t
     */
    static const constexpr size_t alignment = 32;

    /*
     * Check if the pointer p is a multiple of alignemnt
     */
    template <class T> static constexpr bool valid(T *p) { return (int64_t)p % alignment == 0; }

    /*
     * Check if the number n is a multiple of vect_size
     */
    template <class T> static constexpr bool compliant(T n) { return n % vect_size == 0; }

    /*
     * Converter from vect_t to a tab.
     * exple:
     *      Converter conv;
     *      conv.v = a;
     *      scalart_t x = conv.t[1]
     */
    union Converter {
        vect_t v;
        scalar_t t[vect_size];
    };

    /*
     *  Return vector of type vect_t with all elements set to zero
     *  Return [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] int16_t
     */
    static INLINE CONST vect_t zero() { return _mm256_setzero_si256(); }

    /*
     *  Broadcast 16-bit integer a to all all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x] int16_t
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm256_set1_epi16(x); }

    /*
     *  Broadcast 16-bit integer a to all all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15] int16_t
     */
    static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1, const scalar_t x2, const scalar_t x3,
                                   const scalar_t x4, const scalar_t x5, const scalar_t x6, const scalar_t x7,
                                   const scalar_t x8, const scalar_t x9, const scalar_t x10, const scalar_t x11,
                                   const scalar_t x12, const scalar_t x13, const scalar_t x14, const scalar_t x15) {
        return _mm256_set_epi16(x15, x14, x13, x12, x11, x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0);
    }

    /*
     *  Gather 16-bit integer elements with indexes idx[0], ..., idx[15] from the address p in vect_t.
     *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]],
     p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]],
     p[idx[8]], p[idx[9]], p[idx[10]], p[idx[11]],
     p[idx[12]], p[idx[13]], p[idx[14]], p[idx[15]]] int16_t
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        return set(p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]], p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]], p[idx[8]],
                   p[idx[9]], p[idx[10]], p[idx[11]], p[idx[12]], p[idx[13]], p[idx[14]], p[idx[15]]);
    }

    /*
     * Load 256-bits of integer data from memory into dst.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     * Return [p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11]p[12],p[13],p[14],p[15]] int16_t
     */
    static INLINE PURE vect_t load(const scalar_t *const p) {
        return _mm256_load_si256(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Load 256-bits of integer data from memory into dst.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11]p[12],p[13],p[14],p[15]] int16_t
     */
    static INLINE PURE vect_t loadu(const scalar_t *const p) {
        return _mm256_loadu_si256(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Store 256-bits of integer data from a into memory.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void store(const scalar_t *p, vect_t v) {
        _mm256_store_si256(reinterpret_cast<vect_t *>(const_cast<scalar_t *>(p)), v);
    }

    /*
     * Store 256-bits of integer data from a into memory.
     * p does not need to be aligned on any particular boundary.
     */
    static INLINE void storeu(const scalar_t *p, vect_t v) {
        _mm256_storeu_si256(reinterpret_cast<vect_t *>(const_cast<scalar_t *>(p)), v);
    }

    /*
     * Store 256-bits of integer data from a into memory using a non-temporal memory hint.
     * p must be aligned on a 32-byte boundary or a general-protection exception may be generated.
     */
    static INLINE void stream(const scalar_t *p, const vect_t v) {
        _mm256_stream_si256(reinterpret_cast<vect_t *>(const_cast<scalar_t *>(p)), v);
    }

    /*
     * Shift packed 16-bit integers in a left by s while shifting in zeros, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15] int16_t
     * Return : [a0 << s, a1 << s, a2 << s, a3 << s, a4 << s, a5 << s, a6 << s, a7 << s,
     *           a8 << s, a9 << s, a10 << s, a11 << s, a12 << s, a13 << s, a14 << s, a15 << s] int16_t
     */
    static INLINE CONST vect_t sll(const vect_t a, const int s) { return _mm256_slli_epi16(a, s); }

    /*
     * Shift packed 16-bit integers in a right by s while shifting in zeros, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15] int16_t
     * Return : [a0 >> s, a1 >> s, a2 >> s, a3 >> s, a4 >> s, a5 >> s, a6 >> s, a7 >> s,
     *           a8 >> s, a9 >> s, a10 >> s, a11 >> s, a12 >> s, a13 >> s, a14 >> s, a15 >> s] int16_t
     */
    static INLINE CONST vect_t srl(const vect_t a, const int s) { return _mm256_srli_epi16(a, s); }


    static INLINE CONST vect_t sra(const vect_t a, const int s) { return _mm256_sra_epi16(a, Simd128<int>::set1(s)); }

    /*
     * Add packed 16-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15] int16_t
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15] int16_t
     * Return : [a0+b0, a1+b1, a2+b2, a3+b3, a4+b4, a5+b5, a6+b6, a7+b7,
     a8+b8, a9+b9, a10+b10, a11+b11, a12+b12, a13+b13, a14+b14, a15+b15]   int16_t
     */
    static INLINE CONST vect_t add(const vect_t a, const vect_t b) { return _mm256_add_epi16(a, b); }

    static INLINE vect_t addin(vect_t &a, const vect_t b) { return a = add(a, b); }

    /*
     * Subtract packed 16-bit integers in b from packed 16-bit integers in a, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15] int16_t
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15] int16_t
     * Return : [a0-b0, a1-b1, a2-b2, a3-b3, a4-b4, a5-b5, a6-b6, a7-b7,
     a8-b8, a9-b9, a10-b10, a11-b11, a12-b12, a13-b13, a14-b14, a15-b15]  int16_t
     */
    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) { return _mm256_sub_epi16(a, b); }

    static INLINE CONST vect_t subin(vect_t &a, const vect_t b) { return a = sub(a, b); }

    /*
     * Multiply the packed 16-bit integers in a and b, producing intermediate 32-bit integers, and store the low 16 bits
     of the intermediate integers in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]           int16_t
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]  		 int16_t
     * Return : [a0*b0 mod 2^16-1, a1*b1 mod 2^16-1, a2*b2 mod 2^16-1, a3*b3 mod 2^16-1,
     a4*b4 mod 2^16-1, a5*b5 mod 2^16-1, a6*b6 mod 2^16-1, a7*b7 mod 2^16-1,
     a8*b8 mod 2^16-1, a9*b9 mod 2^16-1, a10*b10 mod 2^16-1, a11*b11 mod 2^16-1,
     a12*b12 mod 2^16-1, a13*b13 mod 2^16-1, a14-b14 mod 2^16-1, a15*b15 mod 2^16-1] int16_t
     */
    static INLINE CONST vect_t mullo(const vect_t a, const vect_t b) { return _mm256_mullo_epi16(a, b); }

    /*
     * Multiply the packed 16-bit integers in a and b, producing intermediate 32-bit integers, and store the low 16 bits
     of the intermediate integers in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]           int16_t
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]           int16_t
     * Return : [a0*b0 mod 2^16-1, a1*b1 mod 2^16-1, a2*b2 mod 2^16-1, a3*b3 mod 2^16-1,
     a4*b4 mod 2^16-1, a5*b5 mod 2^16-1, a6*b6 mod 2^16-1, a7*b7 mod 2^16-1,
     a8*b8 mod 2^16-1, a9*b9 mod 2^16-1, a10*b10 mod 2^16-1, a11*b11 mod 2^16-1,
     a12*b12 mod 2^16-1, a13*b13 mod 2^16-1, a14-b14 mod 2^16-1, a15*b15 mod 2^16-1] int16_t
     */
    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) { return mullo(a, b); }

    /*
     * Multiply packed 16-bit integers in a and b, producing intermediate 32-bit integers, and add the low 16-bits of
     the intermediate with c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]           int16_t
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]           int16_t
     [c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15]           int16_t
     * Return : [(a0*b0 mod 2^16-1)+c0, (a1*b1 mod 2^16-1)+c1, (a2*b2 mod 2^16-1)+c2, (a3*b3 mod 2^16-1)+c3,
     (a4*b4 mod 2^16-1)+c4, (a5*b5 mod 2^16-1)+c5, (a6*b6 mod 2^16-1)+c6, (a7*b7 mod 2^16-1)+c7,
     (a8*b8 mod 2^16-1)+c8, (a9*b9 mod 2^16-1)+c9, (a10*b10 mod 2^16-1)+c10, (a11*b11 mod 2^16-1)+c11,
     (a12*b12 mod 2^16-1)+c12, (a13*b13 mod 2^16-1)+c13, (a14*b14 mod 2^16-1)+c14, (a15*b15 mod 2^16-1)+c15]
     */
    static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b) { return add(c, mul(a, b)); }

    static INLINE CONST vect_t fmaddin(vect_t c, const vect_t a, const vect_t b) { return c = fmadd(c, a, b); }

    /*
     * Multiply packed 16-bit integers in a and b, producing intermediate 32-bit integers, and substract elements of c
     to the low 16-bit of the intermiate result, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]           int16_t
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]           int16_t
     [c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15]           int16_t
     * Return : [-(a0*b0 mod 2^16-1)+c0, -(a1*b1 mod 2^16-1)+c1, -(a2*b2 mod 2^16-1)+c2, -(a3*b3 mod 2^16-1)+c3,
     -(a4*b4 mod 2^16-1)+c4, -(a5*b5 mod 2^16-1)+c5, -(a6*b6 mod 2^16-1)+c6, -(a7*b7 mod 2^16-1)+c7,
     -(a8*b8 mod 2^16-1)+c8, -(a9*b9 mod 2^16-1)+c9, -(a10*b10 mod 2^16-1)+c10, -(a11*b11 mod 2^16-1)+c11,
     -(a12*b12 mod 2^16-1)+c12, -(a13*b13 mod 2^16-1)+c13, -(a14*b14 mod 2^16-1)+c14, -(a15*b15 mod 2^16-1)+c15]
     */
    static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mul(a, b)); }

    static INLINE CONST vect_t fnmaddin(vect_t c, const vect_t a, const vect_t b) { return c = fnmadd(c, a, b); }

    /*
     * Multiply packed 16-bit integers in a and b, producing intermediate 32-bit integers, and substract the low 16-bits
     of the intermediate with c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]           int16_t
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]           int16_t
     [c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15]           int16_t
     * Return : [(a0*b0 mod 2^16-1)-c0, (a1*b1 mod 2^16-1)-c1, (a2*b2 mod 2^16-1)-c2, (a3*b3 mod 2^16-1)-c3,
     (a4*b4 mod 2^16-1)-c4, (a5*b5 mod 2^16-1)-c5, (a6*b6 mod 2^16-1)-c6, (a7*b7 mod 2^16-1)-c7,
     (a8*b8 mod 2^16-1)-c8, (a9*b9 mod 2^16-1)-c9, (a10*b10 mod 2^16-1)-c10, (a11*b11 mod 2^16-1)-c11,
     (a12*b12 mod 2^16-1)-c12, (a13*b13 mod 2^16-1)-c13, (a14*b14 mod 2^16-1)-c14, (a15*b15 mod 2^16-1)-c15]
     */
    static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mul(a, b)); }

    static INLINE CONST vect_t fsubin(vect_t c, const vect_t a, const vect_t b) { return c = fmsub(c, a, b); }

    /*
     * Multiply the packed 16-bit integers in a and b, producing intermediate 32-bit integers, and store the high 16
     bits of the intermediate integers in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15] int16_t
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15] int16_t
     * Return :
     */
    static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b) { return _mm256_mulhi_epi16(a, b); }

    /*
     * Multiply the low 8-bit integers from each packed 16-bit element in a and b, and store the signed 16-bit results
     in dst.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]    int16_t
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]    int16_t
     * Return : [a0*b0, a1*b1, a2*b2, a3*b3, a4*b4, a5*b5, a6*b6, a7*b7, a8*b8, a9*b9, a10*b10, a11*b11, a12*b12,
     a13*b13, a14*b14, a15*b15] int16_t
     */
    static INLINE CONST vect_t mulx(vect_t a, vect_t b) {
        vect_t mask = set1(0x00FF);
        a = vand(a, mask);
        b = vand(b, mask);
        return mullo(a, b);
    }

    /*
     * Compare packed 16-bits in a and b for equality, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15] int16_t
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15] int16_t
     * Return : [(a0==b0) ? 0xFFFF : 0, (a1==b1) ? 0xFFFF : 0,
     (a2==b2) ? 0xFFFF : 0, (a3==b3) ? 0xFFFF : 0,
     (a4==b4) ? 0xFFFF : 0, (a5==b5) ? 0xFFFF : 0,
     (a6==b6) ? 0xFFFF : 0, (a7==b7) ? 0xFFFF : 0,
     (a8==b8) ? 0xFFFF : 0, (a9==b9) ? 0xFFFF : 0,
     (a10==b10) ? 0xFFFF : 0, (a11==b11) ? 0xFFFF : 0,
     (a12==b12) ? 0xFFFF : 0, (a13==b13) ? 0xFFFF : 0,
     (a14==b14) ? 0xFFFF : 0, (a15==b15) ? 0xFFFF : 0]                     int16_t
     */
    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) { return _mm256_cmpeq_epi16(a, b); }

    /*
     * Compare packed 16-bits in a and b for greater-than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15] int16_t
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15] int16_t
     * Return : [(a0>b0) ? 0xFFFF : 0, (a1>b1) ? 0xFFFF : 0,
     (a2>b2) ? 0xFFFF : 0, (a3>b3) ? 0xFFFF : 0,
     (a4>b4) ? 0xFFFF : 0, (a5>b5) ? 0xFFFF : 0,
     (a6>b6) ? 0xFFFF : 0, (a7>b7) ? 0xFFFF : 0,
     (a8>b8) ? 0xFFFF : 0, (a9>b9) ? 0xFFFF : 0,
     (a10>b10) ? 0xFFFF : 0, (a11>b11) ? 0xFFFF : 0,
     (a12>b12) ? 0xFFFF : 0, (a13>b13) ? 0xFFFF : 0,
     (a14>b14) ? 0xFFFF : 0, (a15>b15) ? 0xFFFF : 0]					  int16_t
     */
    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) { return _mm256_cmpgt_epi16(a, b); }

    /*
     * Compare packed 16-bits in a and b for lesser-than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15] int16_t
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15] int16_t
     * Return : [(a0<b0) ? 0xFFFF : 0, (a1<b1) ? 0xFFFF : 0,
     (a2<b2) ? 0xFFFF : 0, (a3<b3) ? 0xFFFF : 0,
     (a4<b4) ? 0xFFFF : 0, (a5<b5) ? 0xFFFF : 0,
     (a6<b6) ? 0xFFFF : 0, (a7<b7) ? 0xFFFF : 0,
     (a8<b8) ? 0xFFFF : 0, (a9<b9) ? 0xFFFF : 0,
     (a10<b10) ? 0xFFFF : 0, (a11<b11) ? 0xFFFF : 0,
     (a12<b12) ? 0xFFFF : 0, (a13<b13) ? 0xFFFF : 0,
     (a14<b14) ? 0xFFFF : 0, (a15>b15) ? 0xFFFF : 0] 					  int16_t
     */
    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) { return _mm256_cmpgt_epi16(b, a); }

    /*
     * Compare packed 16-bits in a and b for greater or equal than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15] int16_t
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15] int16_t
     * Return : [(a0>=b0) ? 0xFFFF : 0, (a1>=b1) ? 0xFFFF : 0,
     (a2>=b2) ? 0xFFFF : 0, (a3>=b3) ? 0xFFFF : 0,
     (a4>=b4) ? 0xFFFF : 0, (a5>=b5) ? 0xFFFF : 0,
     (a6>=b6) ? 0xFFFF : 0, (a7>=b7) ? 0xFFFF : 0,
     (a8>=b8) ? 0xFFFF : 0, (a9>=b9) ? 0xFFFF : 0,
     (a10>=b10) ? 0xFFFF : 0, (a11>=b11) ? 0xFFFF : 0,
     (a12>=b12) ? 0xFFFF : 0, (a13>=b13) ? 0xFFFF : 0,
     (a14>=b14) ? 0xFFFF : 0, (a15>=b15) ? 0xFFFF : 0]					  int16_t
     */
    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

    /*
     * Compare packed 16-bits in a and b for lesser or equal than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15] int16_t
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15] int16_t
     * Return : [(a0<=b0) ? 0xFFFF : 0, (a1<=b1) ? 0xFFFF : 0,
     (a2<=b2) ? 0xFFFF : 0, (a3<=b3) ? 0xFFFF : 0,
     (a4<=b4) ? 0xFFFF : 0, (a5<=b5) ? 0xFFFF : 0,
     (a6<=b6) ? 0xFFFF : 0, (a7<=b7) ? 0xFFFF : 0,
     (a8<=b8) ? 0xFFFF : 0, (a9<=b9) ? 0xFFFF : 0,
     (a10<=b10) ? 0xFFFF : 0, (a11<=b11) ? 0xFFFF : 0,
     (a12<=b12) ? 0xFFFF : 0, (a13<=b13) ? 0xFFFF : 0,
     (a14<=b14) ? 0xFFFF : 0, (a15<=b15) ? 0xFFFF : 0] 					   int16_t
     */
    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

    /*
     * Compute the bitwise AND of packed 16-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [a0 AND b0, a1 AND b1, a2 AND b2, a3 AND b3, a4 AND b4, a5 AND b5, a6 AND b6, a7 AND b7,
     a8 AND b8, a9 AND b9, a10 AND b10, a11 AND b11, a12 AND b12, a13 AND b13, a14 AND b14, a15 AND b15]
     */
    static INLINE CONST vect_t vand(const vect_t a, const vect_t b) { return _mm256_and_si256(b, a); }

    /*
     * Compute the bitwise OR of packed 16-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [a0 OR b0, a1 OR b1, a2 OR b2, a3 OR b3, a4 OR b4, a5 OR b5, a6 OR b6, a7 OR b7,
     a8 OR b8, a9 OR b9, a10 OR b10, a11 OR b11, a12 OR b12, a13 OR b13, a14 OR b14, a15 OR b15]
     */
    static INLINE CONST vect_t vor(const vect_t a, const vect_t b) { return _mm256_or_si256(b, a); }

    /*
     * Compute the bitwise XOR of packed 16-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [a0 XOR b0, a1 XOR b1, a2 XOR b2, a3 XOR b3, a4 XOR b4, a5 XOR b5, a6 XOR b6, a7 XOR b7,
     a8 XOR b8, a9 XOR b9, a10 XOR b10, a11 XOR b11, a12 XOR b12, a13 XOR b13, a14 XOR b14, a15 XOR b15]
     */
    static INLINE CONST vect_t vxor(const vect_t a, const vect_t b) { return _mm256_xor_si256(b, a); }

    /*
     * Compute the bitwise AND NOT of packed 16-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [a0 ANDNOT b0, a1 ANDNOT b1, a2 ANDNOT b2, a3 ANDNOT b3, a4 ANDNOT b4, a5 ANDNOT b5, a6 ANDNOT b6, a7
     ANDNOT b7,
     a8 ANDNOT b8, a9 ANDNOT b9, a10 ANDNOT b10, a11 ANDNOT b11, a12 ANDNOT b12, a13 ANDNOT b13, a14 ANDNOT b14, a15
     ANDNOT b15]
     */
    static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) { return _mm256_andnot_si256(b, a); }

    /*
     * Horizontally add 16-bits elements of a.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]
     * Return : a0+a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14+a15
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        Converter ca;
        ca.v = a;
        return ca.t[0] + ca.t[1] + ca.t[2] + ca.t[3] + ca.t[4] + ca.t[5] + ca.t[6] + ca.t[7] + ca.t[8] + ca.t[9] +
               ca.t[10] + ca.t[11] + ca.t[12] + ca.t[13] + ca.t[14] + ca.t[15];
    }

    static INLINE PURE half_t load_half(const scalar_t *const p) {
        return _mm_load_si128(reinterpret_cast<const half_t *>(p));
    }

    static INLINE PURE half_t loadu_half(const scalar_t *const p) {
        return _mm_loadu_si128(reinterpret_cast<const half_t *>(p));
    }

    static INLINE void store_half(const scalar_t *p, half_t v) {
        _mm_store_si128(reinterpret_cast<half_t *>(const_cast<scalar_t *>(p)), v);
    }

    static INLINE void storeu_half(const scalar_t *p, half_t v) {
        _mm_storeu_si128(reinterpret_cast<half_t *>(const_cast<scalar_t *>(p)), v);
    }

    /*
     *
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]    int16_t
     [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]    int16_t
     [c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15]    int16_t
     * Return : [a0*b0+c0, a1*b1+c1, a2*b2+c2, a3*b3+c3, a4*b4+c4, a5*b5+c5, a6*b6+c6, a7*b7+c7, a8*b8+c8, a9*b9+c9,
     a10*b10+c10, a11*b11+c11, a12*b12+c12, a13*b13+c13, a14*b14+c14, a15*b15+c15] int16_t
     */
    static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b) { return add(c, mulx(a, b)); }

    static INLINE vect_t fmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmaddx(c, a, b); }

    static INLINE CONST vect_t fnmaddx(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mulx(a, b)); }

    static INLINE vect_t fnmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmaddx(c, a, b); }

    static INLINE CONST vect_t round(const vect_t a) { return a; }

    static INLINE CONST vect_t signbits(const vect_t x) {
        vect_t signBits = sub(zero(), srl(x, 4*sizeof(scalar_t)-1));
        return signBits;
    }

    static INLINE vect_t mod(vect_t &C, const vect_t &P, const vect_t &INVP, const vect_t &NEGP, const vect_t &MIN,
                             const vect_t &MAX, vect_t &Q, vect_t &T) {
#ifdef __INTEL_COMPILER
        C = _mm256_rem_epi16(C, P);
#else
        FFLASFFPACK_abort("pas implementé");
#endif
        NORML_MOD(C, P, NEGP, MIN, MAX, Q, T);
        return C;
    }

#else

#error "You need AVX2 instructions to perform 256bits operations on int16_t"

#endif // defined(__FFLASFFPACK_USE_AVX2)
};

// uint16_t
template <> struct Simd256_impl<true, true, false, 2> : public Simd256_impl<true, true, true, 2> {
    using scalar_t = uint16_t;

#if defined(__FFLASFFPACK_USE_AVX2)

    static INLINE CONST vect_t greater(vect_t a, vect_t b) {

        vect_t x;
        x = set1(-(static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1)));
        a = sub(x, a);
        b = sub(x, b);
        return _mm256_cmpgt_epi16(a, b);
    }

    static INLINE CONST vect_t lesser(vect_t a, vect_t b) {
        vect_t x;
        x = set1(-(static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1)));
        a = sub(x, a);
        b = sub(x, b);
        return _mm256_cmpgt_epi16(a, b);
    }

    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }
#else

#error "You need AVX2 instructions to perform 256bits operations on uint16_t"

#endif // defined(__FFLASFFPACK_USE_AVX2)
};

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd256_int16_INL
