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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd128_double_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd128_double_INL

#include "givaro/givtypestring.h"
#include "fflas-ffpack/utils/align-allocator.h"
#include <vector>
#include <type_traits>

/*
 * Simd128 specialized for double
 */
template <> struct Simd128_impl<true, false, true, 8> {
#if defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS)

    /*
     * alias to 128 bit simd register
     */
    using vect_t = __m128d;

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = double;

    /*
     *  number of scalar_t in a simd register
     */
    static const constexpr size_t vect_size = 2;

    /*
     *  string describing the Simd struct
     */
    static const std::string type_string () {
        return "Simd" + std::to_string(8*vect_size*sizeof(scalar_t)) + "<"
                      + Givaro::TypeString<scalar_t>::get() + ">";
    }

    /*
     *  alignement required by scalar_t pointer to be loaded in a vect_t
     */
    static const constexpr size_t alignment = 16;
    using aligned_allocator = AlignedAllocator<scalar_t, Alignment(alignment)>;
    using aligned_vector = std::vector<scalar_t, aligned_allocator>;

    /* To check compatibility with Modular struct */
    template <class Field>
    using is_same_element = std::is_same<typename Field::Element, scalar_t>;

    union Converter {
        vect_t v;
        scalar_t t[vect_size];
    };

    /*
     * Check if the pointer p is a multiple of alignemnt
     */
    template <class T> static constexpr bool valid(T *p) { return (int64_t)p % alignment == 0; }

    /*
     * Check if the number n is a multiple of vect_size
     */
    template <class T> static constexpr bool compliant(T n) { return n % vect_size == 0; }

    /*
     * Return vector of type vect_t with all elements set to zero.
     * Return [0,0]
     */
    static INLINE CONST vect_t zero() { return _mm_setzero_pd(); }

    /*
     * Broadcast double-precision (64-bit) floating-point value a to all elements of vect_t.
     * Return [x,x]
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm_set1_pd(x); }

    /*
     *  Set packed double-precision (64-bit) floating-point elements in vect_t with the supplied values.
     *  Return [x1,x2]
     */
    static INLINE CONST vect_t set(const scalar_t x1, const scalar_t x2) { return _mm_set_pd(x2, x1); }

    /*
     *  Gather double-precision (64-bit) floating-point elements with indexes idx[0], ..., idx[3] from the address p in
     * vect_t.
     *  Return [p[idx[0]], p[idx[1]]]
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        return _mm_set_pd(p[idx[1]], p[idx[0]]);
    }

    /*
     * Load 128-bits (composed of 2 packed double-precision (64-bit) floating-point elements) from memory into vect_t.
     * p must be aligned on a 16-byte boundary or a general-protection exception will be generated.
     * Return [p[0], p[1]]
     */
    static INLINE PURE vect_t load(const scalar_t *const p) { return _mm_load_pd(p); }

    /*
     * Load 128-bits (composed of 2 packed double-precision (64-bit) floating-point elements) from memory into vect_t.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0], p[1]]
     */
    static INLINE PURE vect_t loadu(const scalar_t *const p) { return _mm_loadu_pd(p); }

    /*
     * Store 128-bits (composed of 2 packed double-precision (64-bit) floating-point elements) from p into memory.
     * p must be aligned on a 16-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void store(const scalar_t *p, const vect_t v) { _mm_store_pd(const_cast<scalar_t *>(p), v); }

    /*
     * Store 128-bits (composed of 2 packed double-precision (64-bit) floating-point elements) from p into memory.
     * p must be aligned on a 16-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void storeu(const scalar_t *p, const vect_t v) { _mm_storeu_pd(const_cast<scalar_t *>(p), v); }

    /*
     * Store 128-bits (composed of 2 packed double-precision (64-bit) floating-point elements) from a into memory using
     * a non-temporal memory hint.
     * p must be aligned on a 16-byte boundary or a general-protection exception may be generated.
     */
    static INLINE void stream(const scalar_t *p, const vect_t v) { _mm_stream_pd(const_cast<scalar_t *>(p), v); }

    /*
     * Shuffle double-precision (64-bit) floating-point elements using the control in s,
     * and store the results in dst.
     * Args   : [a0, a1] double
     * Return : [a[s[0]], a[s[1]]] double
     */
#if defined(__FFLASFFPACK_HAVE_AVX_INSTRUCTIONS)
    template<uint8_t s>
    static INLINE CONST vect_t shuffle(const vect_t a) {
        return _mm_permute_pd(a, s);
    }
#endif

    /*
     * Unpack and interleave double-precision (64-bit) floating-point elements
     * from the low half of a and b.
     * Args: a = [ a0, a1 ]
     *       b = [ b0, b1 ]
     * Return:   [ a0, b0 ]
     */
    static INLINE CONST vect_t
    unpacklo_intrinsic (const vect_t a, const vect_t b) {
        return _mm_unpacklo_pd(a, b);
    }

    /*
     * Unpack and interleave double-precision (64-bit) floating-point elements
     * from the high half of a and b.
     * Args: a = [ a0, a1 ]
     *       b = [ b0, b1 ]
     * Return:   [ a1, b1 ]
     */
    static INLINE CONST vect_t
    unpackhi_intrinsic (const vect_t a, const vect_t b) {
        return _mm_unpackhi_pd(a, b);
    }

    /*
     * Unpack and interleave double-precision (64-bit) floating-point elements
     * from the low half of a and b.
     * Args: a = [ a0, a1 ]
     *       b = [ b0, b1 ]
     * Return:   [ a0, b0 ]
     */
    static INLINE CONST vect_t unpacklo(const vect_t a, const vect_t b) {
        return unpacklo_intrinsic(a, b);
    }

    /*
     * Unpack and interleave double-precision (64-bit) floating-point elements
     * from the high half of a and b.
     * Args: a = [ a0, a1 ]
     *       b = [ b0, b1 ]
     * Return:   [ a1, b1 ]
     */
    static INLINE CONST vect_t unpackhi(const vect_t a, const vect_t b) {
        return unpackhi_intrinsic(a, b);
    }

    /*
     * Perform unpacklo and unpackhi with a and b and store the results in lo
     * and hi.
     * Args: a = [ a0, a1  ]
     *       b = [ b0, b1  ]
     * Return: lo = [ a0, b0 ]
     *         hi = [ a1, b1 ]
     */
    static INLINE void
    unpacklohi (vect_t& lo, vect_t& hi, const vect_t a, const vect_t b) {
        lo = unpacklo (a, b);
        hi = unpackhi (a, b);
    }

    /*
     * Pack double-precision (64-bit) floating-point elements from the even
     * positions of a and b.
     * Args: a = [ a0, a1 ]
     *       b = [ b0, b1 ]
     * Return:   [ a0, b0 ]
     */
    static INLINE CONST vect_t pack_even (const vect_t a, const vect_t b) {
        return unpacklo (a, b); /* same as unpacklo for vect_size = 2 */
    }

    /*
     * Pack double-precision (64-bit) floating-point elements from the odd
     * positions of a and b.
     * Args: a = [ a0, a1 ]
     *       b = [ b0, b1 ]
     * Return:   [ a1, b1 ]
     */
    static INLINE CONST vect_t pack_odd (const vect_t a, const vect_t b) {
        return unpackhi (a, b); /* same as unpackhi for vect_size = 2 */
    }

    /*
     * Perform pack_even and pack_odd with a and b and store the results in even
     * and odd.
     * Args: a = [ a0, a1 ]
     *       b = [ b0, b1 ]
     * Return: even = [ a0, b0 ]
     *         odd  = [ a1, b1 ]
     */
    static INLINE void
    pack (vect_t& even, vect_t& odd, const vect_t a, const vect_t b) {
        even = pack_even (a, b);
        odd = pack_odd (a, b);
    }

    /*
     * Transpose the 2x2 matrix formed by the 2 rows of double-precision
     * (64-bit) floating-point elements in r0 and r1, and store the transposed
     * matrix in these vectors.
     * Args: r0 = [ r00, r01 ]
     *       r1 = [ r10, r11 ]
     * Return: r0 = [ r00, r10 ]
     *         r1 = [ r01, r11 ]
     */
    static INLINE void
    transpose (vect_t& r0, vect_t& r1) {
        unpacklohi (r0, r1, r0, r1);
    }

    /*
     * Blend packed double-precision (64-bit) floating-point elements from a and
     * b using control mask s.
     * Args: a = [ a0, a1 ]
     *       b = [ b0, a1 ]
     *       s = a 2-bit immediate integer
     * Return: [ s[0] ? a0 : b0, s[1] ? a1 : b1 ]
     */
    template<uint8_t s>
    static INLINE CONST vect_t blend(const vect_t a, const vect_t b) {
        return _mm_blend_pd(a, b, s);
    }

    /*
     * Blend packed double-precision (64-bit) floating-point elements from a and
     * b using the vector mask as control.
     * Args: a = [ a0, a1 ]
     *       b = [ b0, b1 ]
     *       mask
     * Return: [ mask[63] ? a0 : b0, mask[127] ? a1 : b1 ]
     */
    static INLINE CONST vect_t blendv(const vect_t a, const vect_t b, const vect_t mask) {
        return _mm_blendv_pd(a, b, mask);
    }

    /*
     * Add packed double-precision (64-bit) floating-point elements in a and b, and store the results in vect_t.
     * Args   : [a0, a1], [b0, b1]
     * Return : [a0+b0, a1+b1]
     */
    static INLINE CONST vect_t add(const vect_t a, const vect_t b) { return _mm_add_pd(a, b); }

    static INLINE vect_t addin(vect_t &a, const vect_t b) { return a = add(a, b); }

    /*
     * Subtract packed double-precision (64-bit) floating-point elements in b from packed double-precision (64-bit)
     * floating-point elements in a, and store the results in vect_t.
     * Args   : [a0, a1], [b0, b1]
     * Return : [a0-b0, a1-b1]
     */
    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) { return _mm_sub_pd(a, b); }

    static INLINE CONST vect_t subin(vect_t &a, const vect_t b) { return a = sub(a, b); }

    /*
     * Multiply packed double-precision (64-bit) floating-point elements in a and b, and store the results in vect_t.
     * Args   : [a0, a1], [b0, b1]
     * Return : [a0*b0, a1*b1]
     */
    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) { return _mm_mul_pd(a, b); }

    static INLINE CONST vect_t mulin(vect_t &a, const vect_t b) { return a = mul(a, b); }

    /*
     * Divide packed double-precision (64-bit) floating-point elements in a by packed elements in b,
     * and store the results in dst.
     * Args   : [a0, a1], [b0, b1]
     * Return : [a0/b0, a1/b1]
     */
    static INLINE CONST vect_t div(const vect_t a, const vect_t b) { return _mm_div_pd(a, b); }

    /*
     * Multiply packed double-precision (64-bit) floating-point elements in a and b, add the intermediate result to
     * packed elements in c, and store the results in vect_t.
     * Args   : [a0, a1], [b0, b1], [c0, c1]
     * Return : [a0*b0+c0, a1*b1+c1]
     */
    static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __FMA__
        return _mm_fmadd_pd(a, b, c);
#else
	Converter ca, cb, cc;
        ca.v = a;
        cb.v = b;
        cc.v = c;
        return set(std::fma (ca.t[0], cb.t[0], cc.t[0]),
                   std::fma (ca.t[1], cb.t[1], cc.t[1]) );
#endif
    }

    static INLINE CONST vect_t fmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fmadd(c, a, b); }

    /*
     * Multiply packed double-precision (64-bit) floating-point elements in a and b, add the negated intermediate result
     * to packed elements in c, and store the results in vect_t.
     * Args   : [a0, a1], [b0, b1], [c0, c1]
     * Return : [-(a0*b0)+c0, -(a1*b1)+c1]
     */
    static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __FMA__
        return _mm_fnmadd_pd(a, b, c);
#else
	return fmadd (c, sub (zero(), a), b);
#endif
    }

    static INLINE CONST vect_t fnmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmadd(c, a, b); }

    /*
     * Multiply packed double-precision (64-bit) floating-point elements in a and b, subtract packed elements in c from
     * the intermediate result, and store the results in vect_t.
     * Args   : [a0, a1], [b0, b1], [c0, c1]
     * Return : [a0*b0-c0, a1*b1-c1]
     */
    static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __FMA__
        return _mm_fmsub_pd(a, b, c);
#else
	return fmadd (sub (zero(), c), a, b);
#endif
    }

    static INLINE CONST vect_t fmsubin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsub(c, a, b); }

    /*
     * Compare packed double-precision (64-bit) floating-point elements in a and b for equality, and store the results
     in vect_t.
     * Args   : [a0, a1], [b0, b1]
     * Return : [(a0==b0) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a1==b1) ? 0xFFFFFFFFFFFFFFFF : 0]
     */
    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) { return _mm_cmpeq_pd(a, b); }

    /*
     * Compare packed double-precision (64-bit) floating-point elements in a and b for lesser-than, and store the
     results in vect_t.
     * Args   : [a0, a1], [b0, b1]
     * Return : [(a0<b0) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a1<b1) ? 0xFFFFFFFFFFFFFFFF : 0]
     */
    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) { return _mm_cmplt_pd(a, b); }

    /*
     * Compare packed double-precision (64-bit) floating-point elements in a and b for lesser or equal than, and store
     the results in vect_t.
     * Args   : [a0, a1], [b0, b1]
     * Return : [(a0<=b0) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a1<=b1) ? 0xFFFFFFFFFFFFFFFF : 0]
     */
    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return _mm_cmple_pd(a, b); }

    /*
     * Compare packed double-precision (64-bit) floating-point elements in a and b for greater-than, and store the
     results in vect_t.
     * Args   : [a0, a1], [b0, b1]
     * Return : [(a0>b0) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a1>b1) ? 0xFFFFFFFFFFFFFFFF : 0]
     */
    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) { return _mm_cmpgt_pd(a, b); }

    /*
     * Compare packed double-precision (64-bit) floating-point elements in a and b for greater or equal than, and store
     the results in vect_t.
     * Args   : [a0, a1], [b0, b1]
     * Return : [(a0>=b0) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a1>=b1) ? 0xFFFFFFFFFFFFFFFF : 0]
     */
    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return _mm_cmpge_pd(a, b); }

    /*
     * Compute the bitwise AND of packed double-precision (64-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1], [b0, b1]
     * Return : [a0 AND b0, a1 AND b1]
     */
    static INLINE CONST vect_t vand(const vect_t a, const vect_t b) { return _mm_and_pd(a, b); }

    /*
     * Compute the bitwise OR of packed double-precision (64-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1], [b0, b1]
     * Return : [a0 OR b0, a1 OR b1]
     */
    static INLINE CONST vect_t vor(const vect_t a, const vect_t b) { return _mm_or_pd(a, b); }

    /*
     * Compute the bitwise XOR of packed double-precision (64-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1], [b0, b1]
     * Return : [a0 XOR b0, a1 XOR b1]
     */
    static INLINE CONST vect_t vxor(const vect_t a, const vect_t b) { return _mm_xor_pd(a, b); }

    /*
     * Compute the bitwise NOT AND of packed double-precision (64-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1], [b0, b1]
     * Return : [NOT(a0) AND b0, NOT(a1) AND b1]
     */
    static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) { return _mm_andnot_pd(a, b); }

    /*
     * Round the packed double-precision (64-bit) floating-point elements in a down to an integer value, and store the
     * results as packed double-precision floating-point elements in vect_t.
     * Args   : [a0, a1]
     * Return : [floor(a0), floor(a1)]
     */
    static INLINE CONST vect_t floor(const vect_t a) { return _mm_floor_pd(a); }

    /*
     * Round the packed double-precision (64-bit) floating-point elements in a up to an integer value, and store the
     * results as packed double-precision floating-point elements in vect_t.
     * Args   : [a0, a1]
     * Return : [ceil(a0), ceil(a1)]
     */
    static INLINE CONST vect_t ceil(const vect_t a) { return _mm_ceil_pd(a); }

    /*
     * Round the packed double-precision (64-bit) floating-point elements in a, and store the results as packed
     * double-precision floating-point elements in vect_t.
     * Args   : [a0, a1]
     * Return : [round(a0), round(a1)]
     */
    static INLINE CONST vect_t round(const vect_t a) {
        return _mm_round_pd(a, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    }

    /*
     * Horizontally add adjacent pairs of double-precision (64-bit) floating-point elements in a and b, and pack the
     * results in vect_t.
     * Args   : [a0, a1], [b0, b1]
     * Return : [a0+a1, b0+b1]
     */
    static INLINE CONST vect_t hadd(const vect_t a, const vect_t b) { return _mm_hadd_pd(a, b); }

    /*
     * Horizontally add double-precision (64-bit) floating-point elements in a.
     * Args   : [a0, a1]
     * Return : a0+a1
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        return ((const scalar_t *)&a)[0] + ((const scalar_t *)&a)[1];
    }

    static INLINE vect_t mod(vect_t &C, const vect_t &P, const vect_t &INVP, const vect_t &NEGP, const vect_t &MIN,
                             const vect_t &MAX, vect_t &Q, vect_t &T) {
        FLOAT_MOD(C, P, INVP, Q);
        NORML_MOD(C, P, NEGP, MIN, MAX, Q, T);

        return C;
    }

#else // __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS
#error "You need SSE instructions to perform 128bits operations on double"
#endif
};

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd128_double_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
