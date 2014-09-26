/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   Bastien Vialla<bastien.vialla@lirmm.fr>
 * BB <bbboyer@ncsu.edu>
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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd128_int16_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd128_int16_INL

/*
 * Simd128 specialized for int16_t
 */
template<>
struct Simd128_impl<true, true, true, 2>{
#if defined(__FFLASFFPACK_USE_SIMD)
    /*
     * alias to 128 bit simd register
     */
    using vect_t = __m128i;

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = int16_t;

    /*
     *  number of scalar_t in a simd register
     */
    static const constexpr size_t vect_size = 8;

    /*
     *  alignement required by scalar_t pointer to be loaded in a vect_t
     */
    static const constexpr size_t alignment = 16;

    /*
     *  Return vector of type vect_t with all elements set to zero
     *  Return [0,0,0,0,0,0,0,0] int16_t
     */
    static INLINE CONST vect_t zero()
    {
        return _mm_setzero_si128();
    }

    /*
     *  Broadcast 16-bit integer a to all all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x,x,x,x,x,x,x,x] int16_t
     */
    static INLINE CONST vect_t set1(const scalar_t x)
    {
        return _mm_set1_epi16(x);
    }

    /*
     *  Broadcast 16-bit integer a to all all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x0,x1,x2,x3,x4,x5,x6,x7] int16_t
     */
    static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1, const scalar_t x2, const scalar_t x3,
                                   const scalar_t x4, const scalar_t x5, const scalar_t x6, const scalar_t x7)
    {
        return _mm_set_epi16(x7,x6,x5,x4,x3,x2,x1,x0);
    }

    /*
     *  Gather 16-bit integer elements with indexes idx[0], ..., idx[7] from the address p in vect_t.
     *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]],
                p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]]] int16_t
     */
    template<class T>
    static INLINE PURE vect_t gather(const scalar_t * const p, const T * const idx)
    {
        return set(p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]],
                   p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]]);
    }

    /*
     * Load 128-bits of integer data from memory into dst.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     * Return [p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]] int16_t
     */
    static INLINE PURE vect_t load(const scalar_t * const p)
    {
        return _mm_load_si128(reinterpret_cast<const vect_t*>(p));
    }

    /*
     * Load 128-bits of integer data from memory into dst.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]] int16_t
     */
    static INLINE PURE vect_t loadu(const scalar_t * const p)
    {
        return _mm_loadu_si128(reinterpret_cast<const vect_t*>(p));
    }

    /*
     * Store 128-bits of integer data from a into memory.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void store(const scalar_t * p, vect_t v)
    {
        _mm_store_si128(reinterpret_cast<vect_t*>(const_cast<scalar_t*>(p)), v);
    }

    /*
     * Store 128-bits of integer data from a into memory.
     * p does not need to be aligned on any particular boundary.
     */
    static INLINE void storeu(const scalar_t * p, vect_t v)
    {
        _mm_storeu_si128(reinterpret_cast<vect_t*>(const_cast<scalar_t*>(p)), v);
    }

    /*
     * Add packed 16-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int16_t
                [b0, b1, b2, b3, b4, b5, b6, b7] int16_t
     * Return : [a0+b0, a1+b1, a2+b2, a3+b3, a4+b4, a5+b5, a6+b6, a7+b7]   int16_t
     */
    static INLINE CONST vect_t add(const vect_t a, const vect_t b)
    {
        return _mm_add_epi16(a, b);
    }

    static INLINE vect_t addin(vect_t &a, const vect_t b)
    {
        return a = add(a,b);
    }

    /*
     * Subtract packed 16-bit integers in b from packed 16-bit integers in a, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int16_t
                [b0, b1, b2, b3, b4, b5, b6, b7] int16_t
     * Return : [a0-b0, a1-b1, a2-b2, a3-b3, a4-b4, a5-b5, a6-b6, a7-b7]  int16_t
     */
    static INLINE CONST vect_t sub(const vect_t a, const vect_t b)
    {
        return _mm_sub_epi16(a, b);
    }

    static INLINE CONST vect_t subin(vect_t &a, const vect_t b)
    {
        return a = sub(a,b);
    }

    /*
     * Multiply the packed 16-bit integers in a and b, producing intermediate 32-bit integers, and store the low 16 bits of the intermediate integers in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]           int16_t
                [b0, b1, b2, b3, b4, b5, b6, b7]           int16_t
     * Return : [a0*b0 mod 2^16-1, a1*b1 mod 2^16-1, a2*b2 mod 2^16-1, a3*b3 mod 2^16-1,
                 a4*b4 mod 2^16-1, a5*b5 mod 2^16-1, a6*b6 mod 2^16-1, a7*b7 mod 2^16-1] int16_t
     */
    static INLINE CONST vect_t mullo(const vect_t a, const vect_t b)
    {
        return _mm_mullo_epi16(a, b);
    }

    /*
     * Multiply the packed 16-bit integers in a and b, producing intermediate 32-bit integers, and store the high 16 bits of the intermediate integers in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int16_t
                [b0, b1, b2, b3, b4, b5, b6, b7] int16_t
     * Return :
     */
    static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b)
    {
        return _mm_mulhi_epi16(a, b);
    }

    /*
     * Multiply the low 16-bit integers from each packed 32-bit element in a and b, and store the signed 32-bit results in dst.
     * Args   : [0, a1, 0, a3, 0, a5, 0, a7]    int16_t
                [0, b1, 0, b3, 0, b5, 0, b7]    int16_t
     * Return : [a1*b1, a3*b2, a5*b5, a7*b7]    int32_t
     */
    static INLINE CONST vect_t mulx(const vect_t a, const vect_t b)
    {
        vect_t ah, al;
        ah = mulhi(a, b);
        al = mullo(a, b);
        return set(_mm_extract_epi16(ah,1),_mm_extract_epi16(al,1),_mm_extract_epi16(ah,3),_mm_extract_epi16(al,3),_mm_extract_epi16(ah,5),_mm_extract_epi16(al,5),_mm_extract_epi16(ah,7),_mm_extract_epi16(al,7));
    }

    /*
     *
     * Args   : [0, a1, 0, a3, 0, a5, 0, a7] int16_t
                [0, b1, 0, b3, 0, b5, 0, b7] int16_t
                [c0, c1, c2, c3]                            int32_t
     * Return : [c0+a1*b1, c1+a3*b2, c2+a5*b5, c3+a7*b7] int32_t
     */
    static INLINE CONST vect_t maddx(vect_t c, const vect_t a, const vect_t b)
    {
        return simd128<int32_t>::add(c, mulx(a, b));
    }

    /*
     * Compare packed 16-bits in a and b for equality, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int16_t
                [b0, b1, b2, b3, b4, b5, b6, b7] int16_t
     * Return : [(a0==b0) ? 0xFFFF : 0, (a1==b1) ? 0xFFFF : 0,
                 (a2==b2) ? 0xFFFF : 0, (a3==b3) ? 0xFFFF : 0,
                 (a4==b4) ? 0xFFFF : 0, (a5==b5) ? 0xFFFF : 0,
                 (a6==b6) ? 0xFFFF : 0, (a7==b7) ? 0xFFFF : 0]                     int16_t
     */
    static INLINE CONST vect_t eq(const vect_t a, const vect_t b)
    {
        return _mm_cmpeq_epi16(a, b);
    }

    /*
     * Compare packed 16-bits in a and b for greater-than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int16_t
                [b0, b1, b2, b3, b4, b5, b6, b7] int16_t
     * Return : [(a0>b0) ? 0xFFFF : 0, (a1>b1) ? 0xFFFF : 0,
                 (a2>b2) ? 0xFFFF : 0, (a3>b3) ? 0xFFFF : 0,
                 (a4>b4) ? 0xFFFF : 0, (a5>b5) ? 0xFFFF : 0,
                 (a6>b6) ? 0xFFFF : 0, (a7>b7) ? 0xFFFF : 0]                      int16_t
     */
    static INLINE CONST vect_t greater(const vect_t a, const vect_t b)
    {
        return _mm_cmpgt_epi16(a, b);
    }

    /*
     * Compare packed 16-bits in a and b for lesser-than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int16_t
                [b0, b1, b2, b3, b4, b5, b6, b7] int16_t
     * Return : [(a0<b0) ? 0xFFFF : 0, (a1<b1) ? 0xFFFF : 0,
                 (a2<b2) ? 0xFFFF : 0, (a3<b3) ? 0xFFFF : 0,
                 (a4<b4) ? 0xFFFF : 0, (a5<b5) ? 0xFFFF : 0,
                 (a6<b6) ? 0xFFFF : 0, (a7<b7) ? 0xFFFF : 0]                      int16_t
     */
    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b)
    {
        return _mm_cmpgt_epi16(b, a);
    }

    /*
     * Compare packed 16-bits in a and b for greater or equal than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int16_t
                [b0, b1, b2, b3, b4, b5, b6, b7] int16_t
     * Return : [(a0>=b0) ? 0xFFFF : 0, (a1>=b1) ? 0xFFFF : 0,
                 (a2>=b2) ? 0xFFFF : 0, (a3>=b3) ? 0xFFFF : 0,
                 (a4>=b4) ? 0xFFFF : 0, (a5>=b5) ? 0xFFFF : 0,
                 (a6>=b6) ? 0xFFFF : 0, (a7>=b7) ? 0xFFFF : 0]                    int16_t
     */
    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b)
    {
        return vor(greater(a, b), eq(a, b));
    }

    /*
     * Compare packed 16-bits in a and b for lesser or equal than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int16_t
                [b0, b1, b2, b3, b4, b5, b6, b7] int16_t
     * Return : [(a0<=b0) ? 0xFFFF : 0, (a1<=b1) ? 0xFFFF : 0,
                 (a2<=b2) ? 0xFFFF : 0, (a3<=b3) ? 0xFFFF : 0,
                 (a4<=b4) ? 0xFFFF : 0, (a5<=b5) ? 0xFFFF : 0,
                 (a6<=b6) ? 0xFFFF : 0, (a7<=b7) ? 0xFFFF : 0]                     int16_t
     */
    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b)
    {
        return vor(lesser(a, b), eq(a, b));
    }

    /*
     * Compute the bitwise AND of packed 16-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
                [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0 AND b0, a1 AND b1, a2 AND b2, a3 AND b3, a4 AND b4, a5 AND b5, a6 AND b6, a7 AND b7]
     */
    static INLINE CONST vect_t vand(const vect_t a, const vect_t b)
    {
        return _mm_and_si128(b, a);
    }

    /*
     * Compute the bitwise OR of packed 16-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
                [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0 OR b0, a1 OR b1, a2 OR b2, a3 OR b3, a4 OR b4, a5 OR b5, a6 OR b6, a7 OR b7]
     */
    static INLINE CONST vect_t vor(const vect_t a, const vect_t b)
    {
        return _mm_or_si128(b, a);
    }

    /*
     * Compute the bitwise XOR of packed 16-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
                [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0 XOR b0, a1 XOR b1, a2 XOR b2, a3 XOR b3, a4 XOR b4, a5 XOR b5, a6 XOR b6, a7 XOR b7]
     */
    static INLINE CONST vect_t vxor(const vect_t a, const vect_t b)
    {
        return _mm_xor_si128(b, a);
    }

    /*
     * Compute the bitwise AND NOT of packed 16-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
                [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0 ANDNOT b0, a1 ANDNOT b1, a2 ANDNOT b2, a3 ANDNOT b3, a4 ANDNOT b4, a5 ANDNOT b5, a6 ANDNOT b6, a7 ANDNOT b7]
     */
    static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b)
    {
        return _mm_andnot_si128(b, a);
    }

#else
#error "You need SSE instructions to perform 128 bits operations on int16"
#endif // defined(__FFLASFFPACK_USE_AVX2)
};

/*
 * Simd128 specialized for uint16_t
 */
template<>
struct Simd128_impl<true, true, false, 2>{
};

#endif
