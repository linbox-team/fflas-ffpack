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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd256_float_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd256_float_INL

/*
 * Simd256 specialized for float
 */
template <> struct Simd256_impl<true, false, true, 4> {
#if defined(__FFLASFFPACK_USE_AVX) or defined(__FFLASFFPACK_USE_AVX2)
    /*
     * alias to 256 bit simd register
     */
    using vect_t = __m256;

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = float;

    /*
     *	number of scalar_t in a simd register
     */
    static const constexpr size_t vect_size = 8;

    /*
     *	alignement required by scalar_t pointer to be loaded in a vect_t
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
     *	Return vector of type vect_t with all elements set to zero
     *  Return [0,0,0,0,0,0,0,0]
     */
    static INLINE CONST vect_t zero() { return _mm256_setzero_ps(); }

    /*
     *	Broadcast single-precision (32-bit) floating-point value x to all elements of vect_t.
     *  Return [x,x,x,x,x,x,x,x]
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm256_set1_ps(x); }

    /*
     *	Set packed single-precision (32-bit) floating-point elements in vect_t with the supplied values.
     *  Return [x1,x2,x3,x4,x5,x6,x7,x8]
     */
    static INLINE CONST vect_t set(const scalar_t x1, const scalar_t x2, const scalar_t x3, const scalar_t x4,
                                   const scalar_t x5, const scalar_t x6, const scalar_t x7, const scalar_t x8) {
        return _mm256_set_ps(x8, x7, x6, x5, x4, x3, x2, x1);
    }

    /*
     *	Gather single-precision (32-bit) floating-point elements with indexes idx[0], ..., idx[3] from the address p in
     *vect_t.
     *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]], p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]]]
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        // TODO AVX2 Gather
        return _mm256_set_ps(p[idx[7]], p[idx[6]], p[idx[5]], p[idx[4]], p[idx[3]], p[idx[2]], p[idx[1]], p[idx[0]]);
    }

    /*
     * Load 256-bits (composed of 8 packed single-precision (32-bit) floating-point elements) from memory into vect_t.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     * Return [p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]]
     */
    static INLINE PURE vect_t load(const scalar_t *const p) { return _mm256_load_ps(p); }

    /*
     * Load 256-bits (composed of 8 packed single-precision (32-bit) floating-point elements) from memory into vect_t.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]]
     */
    static INLINE PURE vect_t loadu(const scalar_t *const p) { return _mm256_loadu_ps(p); }

    /*
     * Store 256-bits (composed of 8 packed single-precision (32-bit) floating-point elements) from a into memory.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void store(const scalar_t *p, const vect_t v) { _mm256_store_ps(const_cast<scalar_t *>(p), v); }

    /*
     * Store 256-bits (composed of 8 packed single-precision (32-bit) floating-point elements) from a into memory.
     * p does not need to be aligned on any particular boundary.
     */
    static INLINE void storeu(const scalar_t *p, const vect_t v) { _mm256_storeu_ps(const_cast<scalar_t *>(p), v); }

    /*
     * Store 256-bits (composed of 8 packed double-precision (32-bit) floating-point elements) from a into memory using
     * a non-temporal memory hint.
     * p must be aligned on a 32-byte boundary or a general-protection exception may be generated.
     */
    static INLINE void stream(const scalar_t *p, const vect_t v) { _mm256_stream_ps(const_cast<scalar_t *>(p), v); }

    /*
     * Add packed single-precision (32-bit) floating-point elements in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0+b0, a1+b1, a2+b2, a3+b3, a4+b4, a5+b5, a6+b6, a7+b7]
     */
    static INLINE CONST vect_t add(const vect_t a, const vect_t b) { return _mm256_add_ps(a, b); }

    static INLINE vect_t addin(vect_t &a, const vect_t b) { return a = add(a, b); }

    /*
     * Subtract packed single-precision (32-bit) floating-point elements in b from packed single-precision (32-bit)
     * floating-point elements in a, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0-b0, a1-b1, a2-b2, a3-b3, a4-b4, a5-b5, a6-b6, a7-b7]
     */
    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) { return _mm256_sub_ps(a, b); }

    static INLINE CONST vect_t subin(vect_t &a, const vect_t b) { return a = sub(a, b); }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0*b0, a1*b1, a2*b2, a3*b3, a4*b4, a5*b5, a6*b6, a7*b7]
     */
    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) { return _mm256_mul_ps(a, b); }

    static INLINE CONST vect_t mulin(vect_t &a, const vect_t b) { return a = mul(a, b); }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, add the intermediate result to
     * packed elements in c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
     * Return : [a0*b0+c0, a1*b1+c1, a2*b2+c2, a3*b3+c3, a4*b4+c4, a5*b5+c5, a6*b6+c6, a7*b7+c7]
     */
    static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __FMA__
        return _mm256_fmadd_ps(a, b, c);
#else
        return add(c, mul(a, b));
#endif
    }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, add the intermediate result to
     * packed elements in c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
     * Return : [a0*b0+c0, a1*b1+c1, a2*b2+c2, a3*b3+c3, a4*b4+c4, a5*b5+c5, a6*b6+c6, a7*b7+c7]
     */
    static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b) { return fmadd(c, a, b); }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, add the intermediate result to
     * packed elements in c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
     * Return : [a0*b0+c0, a1*b1+c1, a2*b2+c2, a3*b3+c3, a4*b4+c4, a5*b5+c5, a6*b6+c6, a7*b7+c7]
     */
    static INLINE CONST vect_t maddx(const vect_t c, const vect_t a, const vect_t b) { return fmadd(c, a, b); }

    static INLINE CONST vect_t fmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fmadd(c, a, b); }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, add the negated intermediate result
     * to packed elements in c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
     * Return : [-(a0*b0)+c0, -(a1*b1)+c1, -(a2*b2)+c2, -(a3*b3)+c3, -(a4*b4)+c4, -(a5*b5)+c5, -(a6*b6)+c6, -(a7*b7)+c7]
     */
    static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __FMA__
        return _mm256_fnmadd_ps(a, b, c);
#else
        return sub(c, mul(a, b));
#endif
    }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, add the negated intermediate result
     * to packed elements in c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
     * Return : [-(a0*b0)+c0, -(a1*b1)+c1, -(a2*b2)+c2, -(a3*b3)+c3, -(a4*b4)+c4, -(a5*b5)+c5, -(a6*b6)+c6, -(a7*b7)+c7]
     */
    static INLINE CONST vect_t nmadd(const vect_t c, const vect_t a, const vect_t b) { return fnmadd(c, a, b); }

    static INLINE CONST vect_t fnmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmadd(c, a, b); }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, subtract packed elements in c from
     * the intermediate result, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
     * Return : [a0*b0-c0, a1*b1-c1, a2*b2-c2, a3*b3-c3, a4*b4-c4, a5*b5-c5, a6*b6-c6, a7*b7-c7]
     */
    static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __FMA__
        return _mm256_fmsub_ps(a, b, c);
#else
        return sub(mul(a, b), c);
#endif
    }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, subtract packed elements in c from
     * the intermediate result, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
     * Return : [a0*b0-c0, a1*b1-c1, a2*b2-c2, a3*b3-c3, a4*b4-c4, a5*b5-c5, a6*b6-c6, a7*b7-c7]
     */
    static INLINE CONST vect_t msub(const vect_t c, const vect_t a, const vect_t b) { return fmsub(c, a, b); }

    static INLINE CONST vect_t fmsubin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsub(c, a, b); }

    /*
     * Compare packed single-precision (32-bit) floating-point elements in a and b for equality, and store the results
     in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [(a0==b0) ? 0xFFFFFFFF : 0,
     (a1==b1) ? 0xFFFFFFFF : 0,
     (a2==b2) ? 0xFFFFFFFF : 0,
     (a3==b3) ? 0xFFFFFFFF : 0,
     (a4==b4) ? 0xFFFFFFFF : 0,
     (a5==b5) ? 0xFFFFFFFF : 0,
     (a6==b6) ? 0xFFFFFFFF : 0,
     (a7==b7) ? 0xFFFFFFFF : 0]
     */
    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) { return _mm256_cmp_ps(a, b, _CMP_EQ_OQ); }

    /*
     * Compare packed single-precision (32-bit) floating-point elements in a and b for lesser-than, and store the
     results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [(a0<b0) ? 0xFFFFFFFF : 0,
     (a1<b1) ? 0xFFFFFFFF : 0,
     (a2<b2) ? 0xFFFFFFFF : 0,
     (a3<b3) ? 0xFFFFFFFF : 0,
     (a4<b4) ? 0xFFFFFFFF : 0,
     (a5<b5) ? 0xFFFFFFFF : 0,
     (a6<b6) ? 0xFFFFFFFF : 0,
     (a7<b7) ? 0xFFFFFFFF : 0]
     */
    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) { return _mm256_cmp_ps(a, b, _CMP_LT_OS); }

    /*
     * Compare packed single-precision (32-bit) floating-point elements in a and b for lesser or equal than, and store
     the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [(a0<=b0) ? 0xFFFFFFFF : 0,
     (a1<=b1) ? 0xFFFFFFFF : 0,
     (a2<=b2) ? 0xFFFFFFFF : 0,
     (a3<=b3) ? 0xFFFFFFFF : 0,
     (a4<=b4) ? 0xFFFFFFFF : 0,
     (a5<=b5) ? 0xFFFFFFFF : 0,
     (a6<=b6) ? 0xFFFFFFFF : 0,
     (a7<=b7) ? 0xFFFFFFFF : 0]
     */
    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return _mm256_cmp_ps(a, b, _CMP_LE_OS); }

    /*
     * Compare packed single-precision (32-bit) floating-point elements in a and b for greater-than, and store the
     results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [(a0>b0) ? 0xFFFFFFFF : 0,
     (a1>b1) ? 0xFFFFFFFF : 0,
     (a2>b2) ? 0xFFFFFFFF : 0,
     (a3>b3) ? 0xFFFFFFFF : 0,
     (a4>b4) ? 0xFFFFFFFF : 0,
     (a5>b5) ? 0xFFFFFFFF : 0,
     (a6>b6) ? 0xFFFFFFFF : 0,
     (a7>b7) ? 0xFFFFFFFF : 0]
     */
    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) { return _mm256_cmp_ps(a, b, _CMP_GT_OS); }

    /*
     * Compare packed single-precision (32-bit) floating-point elements in a and b for greater or equal than, and store
     the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [(a0>=b0) ? 0xFFFFFFFF : 0,
     (a1>=b1) ? 0xFFFFFFFF : 0,
     (a2>=b2) ? 0xFFFFFFFF : 0,
     (a3>=b3) ? 0xFFFFFFFF : 0,
     (a4>=b4) ? 0xFFFFFFFF : 0,
     (a5>=b5) ? 0xFFFFFFFF : 0,
     (a6>=b6) ? 0xFFFFFFFF : 0,
     (a7>=b7) ? 0xFFFFFFFF : 0]
     */
    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return _mm256_cmp_ps(a, b, _CMP_GE_OS); }

    /*
     * Compute the bitwise AND of packed single-precision (32-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0 AND b0, a1 AND b1, a2 AND b2, a3 AND b3, a4 AND b4, a5 AND b5, a6 AND b6, a7 AND b7]
     */
    static INLINE CONST vect_t vand(const vect_t a, const vect_t b) { return _mm256_and_ps(a, b); }

    /*
     * Compute the bitwise OR of packed single-precision (32-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0 OR b0, a1 OR b1, a2 OR b2, a3 OR b3, a4 OR b4, a5 OR b5, a6 OR b6, a7 OR b7]
     */
    static INLINE CONST vect_t vor(const vect_t a, const vect_t b) { return _mm256_or_ps(a, b); }

    /*
     * Compute the bitwise XOR of packed single-precision (32-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0 XOR b0, a1 XOR b1, a2 XOR b2, a3 XOR b3, a4 XOR b4, a5 XOR b5, a6 XOR b6, a7 XOR b7]
     */
    static INLINE CONST vect_t vxor(const vect_t a, const vect_t b) { return _mm256_xor_ps(a, b); }

    /*
     * Compute the bitwise AND NOT of packed single-precision (32-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0 ANDNOT b0, a1 ANDNOT b1, a2 ANDNOT b2, a3 ANDNOT b3, a4 ANDNOT b4, a5 ANDNOT b5, a6 ANDNOT b6, a7
     * ANDNOT b7]
     */
    static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) { return _mm256_andnot_ps(a, b); }

    /*
     * Round the packed single-precision (32-bit) floating-point elements in a down to an integer value, and store the
     * results as packed double-precision floating-point elements in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
     * Return : [floor(a0), floor(a1), floor(a2), floor(a3), floor(a4), floor(a5), floor(a6), floor(a7)]
     */
    static INLINE CONST vect_t floor(const vect_t a) { return _mm256_floor_ps(a); }

    /*
     * Round the packed single-precision (32-bit) floating-point elements in a up to an integer value, and store the
     * results as packed single-precision floating-point elements in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
     * Return : [ceil(a0), ceil(a1), ceil(a2), ceil(a3), ceil(a4), ceil(a5), ceil(a6), ceil(a7)]
     */
    static INLINE CONST vect_t ceil(const vect_t a) { return _mm256_ceil_ps(a); }

    /*
     * Round the packed single-precision (32-bit) floating-point elements in a, and store the results as packed
     * single-precision floating-point elements in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
     * Return : [round(a0), round(a1), round(a2), round(a3), round(a4), round(a5), round(a6), round(a7)]
     */
    static INLINE CONST vect_t round(const vect_t a) {
        return _mm256_round_ps(a, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    }

    /*
     * Horizontally add adjacent pairs of single-precision (32-bit) floating-point elements in a and b, and pack the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0+a1, b0+b1, a2+a3, b2+b3, a4+a5, b4+b5, a6+a7, b6+b7]
     */
    static INLINE CONST vect_t hadd(const vect_t a, const vect_t b) { return _mm256_hadd_ps(a, b); }

    /*
     * Horizontally add single-precision (32-bit) floating-point elements in a.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
     * Return : a0+a1+a2+a3+a4+a5+a6+a7
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        return ((const scalar_t *)&a)[0] + ((const scalar_t *)&a)[1] + ((const scalar_t *)&a)[2] +
               ((const scalar_t *)&a)[3] + ((const scalar_t *)&a)[4] + ((const scalar_t *)&a)[5] +
               ((const scalar_t *)&a)[6] + ((const scalar_t *)&a)[7];
    }

    static INLINE vect_t mod(vect_t &C, const vect_t &P, const vect_t &INVP, const vect_t &NEGP, const vect_t &MIN,
                             const vect_t &MAX, vect_t &Q, vect_t &T) {
        FLOAT_MOD(C, P, INVP, Q);
        NORML_MOD(C, P, NEGP, MIN, MAX, Q, T);

        return C;
    }

#else // __AVX__
#error "You need AVX instructions to perform 256bits operations on float"
#endif
};

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd256_float_INL
