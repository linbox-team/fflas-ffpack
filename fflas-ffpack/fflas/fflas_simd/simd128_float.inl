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
 * License along with this library; if not, see
 * <https://www.gnu.org/licenses/>.
 * ========LICENCE========
 *.
 */

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd128_float_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd128_float_INL

#include "givaro/givtypestring.h"
#include "fflas-ffpack/utils/align-allocator.h"
#include <vector>
#include <type_traits>

/*
 * Simd128 specialized for float
 */
template <> struct Simd128_impl<true, false, true, 4> {
#if defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS)

    /*
     * alias to 128 bit simd register
     */
    using vect_t = __m128;

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = float;

    /*
     *  number of scalar_t in a simd register
     */
    static const constexpr size_t vect_size = 4;

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
     *  Return vector of type vect_t with all elements set to zero
     *  Return [0,0,0,0]
     */
    static INLINE CONST vect_t zero() { return _mm_setzero_ps(); }

    /*
     *  Broadcast single-precision (32-bit) floating-point value x to all elements of vect_t.
     *  Return [x,x,x,x]
     */
    static INLINE CONST vect_t set1(const scalar_t x) {
#ifdef __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS
        // return _mm_broadcast_ss(&x);
        return _mm_set1_ps(x);
#else
        return _mm_set1_ps(x);
#endif
    }

    /*
     *  Set packed single-precision (32-bit) floating-point elements in vect_t with the supplied values.
     *  Return [x1,x2,x3,x4]
     */
    static INLINE CONST vect_t set(const scalar_t x1, const scalar_t x2, const scalar_t x3, const scalar_t x4) {
        return _mm_set_ps(x4, x3, x2, x1);
    }

    /*
     *  Gather single-precision (32-bit) floating-point elements with indexes idx[0], ..., idx[3] from the address p in
     * vect_t.
     *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]]]
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        return _mm_set_ps(p[idx[3]], p[idx[2]], p[idx[1]], p[idx[0]]);
    }

    /*
     * Load 128-bits (composed of 4 packed single-precision (32-bit) floating-point elements) from memory into vect_t.
     * p must be aligned on a 16-byte boundary or a general-protection exception will be generated.
     * Return [p[0], p[1], p[2], p[3]]
     */
    static INLINE PURE vect_t load(const scalar_t *const p) { return _mm_load_ps(p); }

    /*
     * Load 128-bits (composed of 4 packed single-precision (32-bit) floating-point elements) from memory into vect_t.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0], p[1], p[2], p[3]]
     */
    static INLINE PURE vect_t loadu(const scalar_t *const p) { return _mm_loadu_ps(p); }

    /*
     * Store 128-bits (composed of 4 packed single-precision (32-bit) floating-point elements) from a into memory.
     * p must be aligned on a 16-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void store(const scalar_t *p, const vect_t v) { _mm_store_ps(const_cast<scalar_t *>(p), v); }

    /*
     * Store 128-bits (composed of 4 packed single-precision (32-bit) floating-point elements) from a into memory.
     * p does not need to be aligned on any particular boundary.
     */
    static INLINE void storeu(const scalar_t *p, const vect_t v) { _mm_storeu_ps(const_cast<scalar_t *>(p), v); }

    /*
     * Store 128-bits (composed of 4 packed double-precision (32-bit) floating-point elements) from a into memory using
     * a non-temporal memory hint.
     * p must be aligned on a 16-byte boundary or a general-protection exception may be generated.
     */
    static INLINE void stream(const scalar_t *p, const vect_t v) { _mm_stream_ps(const_cast<scalar_t *>(p), v); }

    /*
     * Shuffle single-precision (32-bit) floating-point elements in a using the control in s,
     * and store the results in dst.
     * Args   :	[a0, a1, a2, a3] float
     * Return :	[a[s[0..1]], ..., a[s[6..7]] float
     */
    template<uint8_t s>
    static INLINE CONST vect_t shuffle(const vect_t a) {
        return _mm_castsi128_ps (_mm_shuffle_epi32 (_mm_castps_si128 (a), s));
    }

    /*
     * Unpack and interleave single-precision (32-bit) floating-point elements
     * from the low half of a and b.
     * Args: a = [ a0, a1, a2, a3 ]
     *       b = [ b0, b1, b2, b3 ]
     * Return:   [ a0, b0, a1, b1 ]
     */
    static INLINE CONST vect_t
    unpacklo_intrinsic (const vect_t a, const vect_t b) {
        return _mm_unpacklo_ps(a, b);
    }

    /*
     * Unpack and interleave single-precision (32-bit) floating-point elements
     * from the high half of a and b.
     * Args: a = [ a0, a1, a2, a3 ]
     *       b = [ b0, b1, b2, b3 ]
     * Return:   [ a2, b2, a3, b3 ]
     */
    static INLINE CONST vect_t
    unpackhi_intrinsic (const vect_t a, const vect_t b) {
        return _mm_unpackhi_ps(a, b);
    }

    /*
     * Unpack and interleave single-precision (32-bit) floating-point elements
     * from the low half of a and b.
     * Args: a = [ a0, a1, a2, a3 ]
     *       b = [ b0, b1, b2, b3 ]
     * Return:   [ a0, b0, a1, b1 ]
     */
    static INLINE CONST vect_t unpacklo(const vect_t a, const vect_t b) {
        return unpacklo_intrinsic(a, b);
    }

    /*
     * Unpack and interleave single-precision (32-bit) floating-point elements
     * from the high half of a and b.
     * Args: a = [ a0, a1, a2, a3 ]
     *       b = [ b0, b1, b2, b3 ]
     * Return:   [ a2, b2, a3, b3 ]
     */
    static INLINE CONST vect_t unpackhi(const vect_t a, const vect_t b) {
        return unpackhi_intrinsic(a, b);
    }

    /*
     * Perform unpacklo and unpackhi with a and b and store the results in lo
     * and hi.
     * Args: a = [ a0, a1, a2, a3 ]
     *       b = [ b0, b1, b2, b3 ]
     * Return: lo = [ a0, b0, a1, b1 ]
     *         hi = [ a2, b2, a3, b3 ]
     */
    static INLINE void
    unpacklohi (vect_t& lo, vect_t& hi, const vect_t a, const vect_t b) {
        lo = unpacklo (a, b);
        hi = unpackhi (a, b);
    }

    /*
     * Pack single-precision (32-bit) floating-point elements from the even
     * positions of a and b.
     * Args: a = [ a0, a1, a2, a3 ]
     *       b = [ b0, b1, b2, b3 ]
     * Return:   [ a0, a2, b0, b2 ]
     */
    static INLINE CONST vect_t pack_even (const vect_t a, const vect_t b) {
        /* 0xd8 = 3120 base_4 */
        __m128i t1 = _mm_shuffle_epi32 (_mm_castps_si128 (a), 0xd8);
        __m128i t2 = _mm_shuffle_epi32 (_mm_castps_si128 (b), 0xd8);
        return _mm_castsi128_ps (_mm_unpacklo_epi64 (t1, t2));
    }

    /*
     * Pack single-precision (32-bit) floating-point elements from the odd
     * positions of a and b.
     * Args: a = [ a0, a1, a2, a3 ]
     *       b = [ b0, b1, b2, b3 ]
     * Return:   [ a1, a3, b1, b3 ]
     */
    static INLINE CONST vect_t pack_odd (const vect_t a, const vect_t b) {
        /* 0xd8 = 3120 base_4 */
        __m128i t1 = _mm_shuffle_epi32 (_mm_castps_si128 (a), 0xd8);
        __m128i t2 = _mm_shuffle_epi32 (_mm_castps_si128 (b), 0xd8);
        return _mm_castsi128_ps (_mm_unpackhi_epi64 (t1, t2));
    }

    /*
     * Perform pack_even and pack_odd with a and b and store the results in even
     * and odd.
     * Args: a = [ a0, a1, a2, a3 ]
     *       b = [ b0, b1, b2, b3 ]
     * Return: even = [ a0, a2, b0, b2 ]
     *         odd  = [ a1, a3, b1, b3 ]
     */
    static INLINE void
    pack (vect_t& even, vect_t& odd, const vect_t a, const vect_t b) {
        /* 0xd8 = 3120 base_4 */
        __m128i t1 = _mm_shuffle_epi32 (_mm_castps_si128 (a), 0xd8);
        __m128i t2 = _mm_shuffle_epi32 (_mm_castps_si128 (b), 0xd8);
        even = _mm_castsi128_ps (_mm_unpacklo_epi64 (t1, t2));
        odd = _mm_castsi128_ps (_mm_unpackhi_epi64 (t1, t2));
    }

    /*
     * Transpose the 4x4 matrix formed by the 4 rows of single-precision
     * (32-bit) floating-point elements in r0, r1, r2 and r3, and store the
     * transposed matrix in these vectors.
     * Args: r0 = [ r00, r01, r02, r03 ]
     *       r1 = [ r10, r11, r12, r13 ]
     *       r2 = [ r20, r21, r22, r23 ]
     *       r3 = [ r30, r31, r32, r33 ]
     * Return: r0 = [ r00, r10, r20, r30 ]
     *         r1 = [ r01, r11, r21, r31 ]
     *         r2 = [ r02, r12, r22, r32 ]
     *         r3 = [ r03, r13, r23, r33 ]
     * Note: taken from https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=transpose&expand=5961
     */
    static INLINE void
    transpose (vect_t& r0, vect_t& r1, vect_t& r2, vect_t& r3) {
        vect_t t0, t1, t2, t3;
        t0 = unpacklo_intrinsic (r0, r1);
        t2 = unpacklo_intrinsic (r2, r3);
        t1 = unpackhi_intrinsic (r0, r1);
        t3 = unpackhi_intrinsic (r2, r3);
        r0 = _mm_movelh_ps (t0, t2);
        r1 = _mm_movehl_ps (t2, t0);
        r2 = _mm_movelh_ps (t1, t3);
        r3 = _mm_movehl_ps (t3, t1);
    }

    /*
     * Blend packed single-precision (32-bit) floating-point elements from a and
     * b using control mask s.
     * Args: a = [ a0, ..., a3 ]
     *       b = [ b0, ..., b3 ]
     *       s = a 4-bit immediate integer
     * Return: [ s[0] ? a0 : b0, ..., s[3] ? a3 : b3 ]
     */
    template<uint8_t s>
    static INLINE CONST vect_t blend(const vect_t a, const vect_t b) {
        return _mm_blend_ps(a, b, s);
    }

    /*
     * Blend packed single-precision (32-bit) floating-point elements from a and
     * b using the vector mask as control.
     * Args: a = [ a0, ..., a3 ]
     *       b = [ b0, ..., b3 ]
     *       mask
     * Return: [ mask[31] ? a0 : b0, ..., mask[127] ? a3 : b3 ]
     */
    static INLINE CONST vect_t blendv(const vect_t a, const vect_t b, const vect_t mask) {
        return _mm_blendv_ps(a, b, mask);
    }

    /*
     * Add packed single-precision (32-bit) floating-point elements in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
     * Return : [a0+b0, a1+b1, a2+b2, a3+b3]
     */
    static INLINE CONST vect_t add(const vect_t a, const vect_t b) { return _mm_add_ps(a, b); }

    static INLINE vect_t addin(vect_t &a, const vect_t b) { return a = add(a, b); }

    /*
     * Subtract packed single-precision (32-bit) floating-point elements in b from packed single-precision (32-bit)
     * floating-point elements in a, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
     * Return : [a0-b0, a1-b1, a2-b2, a3-b3]
     */
    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) { return _mm_sub_ps(a, b); }

    static INLINE CONST vect_t subin(vect_t &a, const vect_t b) { return a = sub(a, b); }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
     * Return : [a0*b0, a1*b1, a2*b2, a3*b3]
     */
    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) { return _mm_mul_ps(a, b); }

    static INLINE CONST vect_t mulin(vect_t &a, const vect_t b) { return a = mul(a, b); }

    /*
     * Divide packed single-precision (32-bit) floating-point elements in a by packed elements in b,
     * and store the results in dst.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3] float
     * Return : [a0/b0, a1/b1, a2/b2, a3/b3] float
     */
    static INLINE CONST vect_t div(const vect_t a, const vect_t b) { return _mm_div_ps(a, b); }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, add the intermediate result to
     * packed elements in c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3], [c0, c1, c2, c3]
     * Return : [a0*b0+c0, a1*b1+c1, a2*b2+c2, a3*b3+c3]
     */
    static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __FMA__
        return _mm_fmadd_ps(a, b, c);
#else
	Converter ca, cb, cc;
        ca.v = a;
        cb.v = b;
        cc.v = c;
        return set(std::fma (ca.t[0], cb.t[0], cc.t[0]),
                   std::fma (ca.t[1], cb.t[1], cc.t[1]),
                   std::fma (ca.t[2], cb.t[2], cc.t[2]),
                   std::fma (ca.t[3], cb.t[3], cc.t[3]) );

#endif
    }

    static INLINE CONST vect_t fmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fmadd(c, a, b); }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, add the intermediate result to
     * packed elements in c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3], [c0, c1, c2, c3]
     * Return : [a0*b0+c0, a1*b1+c1, a2*b2+c2, a3*b3+c3, a4*b4+c4]
     */
    static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __FMA__
        return _mm_fnmadd_ps(a, b, c);
#else
	return fmadd (c, sub (zero(), a), b);
#endif
    }

    static INLINE CONST vect_t fnmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmadd(c, a, b); }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, subtract packed elements in c from
     * the intermediate result, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3], [c0, c1, c2, c3]
     * Return : [a0*b0-c0, a1*b1-c1, a2*b2-c2, a3*b3-c3]
     */
    static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __FMA__
        return _mm_fmsub_ps(a, b, c);
#else
	return fmadd (sub (zero(), c), a, b);
#endif
    }

    static INLINE CONST vect_t fmsubin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsub(c, a, b); }

    /*
     * Compare packed single-precision (32-bit) floating-point elements in a and b for equality, and store the results
     in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
     * Return : [(a0==b0) ? 0xFFFFFFFF : 0,
     (a1==b1) ? 0xFFFFFFFF : 0,
     (a2==b2) ? 0xFFFFFFFF : 0,
     (a3==b3) ? 0xFFFFFFFF : 0]
     */
    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) { return _mm_cmpeq_ps(a, b); }

    /*
     * Compare packed single-precision (32-bit) floating-point elements in a and b for lesser-than, and store the
     results in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
     * Return : [(a0<b0) ? 0xFFFFFFFF : 0,
     (a1<b1) ? 0xFFFFFFFF : 0,
     (a2<b2) ? 0xFFFFFFFF : 0,
     (a3<b3) ? 0xFFFFFFFF : 0]
     */
    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) { return _mm_cmplt_ps(a, b); }

    /*
     * Compare packed single-precision (32-bit) floating-point elements in a and b for lesser or equal than, and store
     the results in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
     * Return : [(a0<=b0) ? 0xFFFFFFFF : 0,
     (a1<=b1) ? 0xFFFFFFFF : 0,
     (a2<=b2) ? 0xFFFFFFFF : 0,
     (a3<=b3) ? 0xFFFFFFFF : 0]
     */
    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return _mm_cmple_ps(a, b); }

    /*
     * Compare packed single-precision (32-bit) floating-point elements in a and b for greater-than, and store the
     results in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
     * Return : [(a0>b0) ? 0xFFFFFFFF : 0,
     (a1>b1) ? 0xFFFFFFFF : 0,
     (a2>b2) ? 0xFFFFFFFF : 0,
     (a3>b3) ? 0xFFFFFFFF : 0]
     */
    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) { return _mm_cmpgt_ps(a, b); }

    /*
     * Compare packed single-precision (32-bit) floating-point elements in a and b for greater or equal than, and store
     the results in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
     * Return : [(a0>=b0) ? 0xFFFFFFFF : 0,
     (a1>=b1) ? 0xFFFFFFFF : 0,
     (a2>=b2) ? 0xFFFFFFFF : 0,
     (a3>=b3) ? 0xFFFFFFFF : 0]
     */
    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return _mm_cmpge_ps(a, b); }

    /*
     * Compute the bitwise AND of packed single-precision (32-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
     * Return : [a0 AND b0, a1 AND b1, a2 AND b2, a3 AND b3]
     */
    static INLINE CONST vect_t vand(const vect_t a, const vect_t b) { return _mm_and_ps(a, b); }

    /*
     * Compute the bitwise OR of packed single-precision (32-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
     * Return : [a0 OR b0, a1 OR b1, a2 OR b2, a3 OR b3]
     */
    static INLINE CONST vect_t vor(const vect_t a, const vect_t b) { return _mm_or_ps(a, b); }

    /*
     * Compute the bitwise XOR of packed single-precision (32-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
     * Return : [a0 XOR b0, a1 XOR b1, a2 XOR b2, a3 XOR b3]
     */
    static INLINE CONST vect_t vxor(const vect_t a, const vect_t b) { return _mm_xor_ps(a, b); }

    /*
     * Compute the bitwise NOT AND of packed single-precision (32-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
     * Return : [NOT(a0) AND b0, NOT(a1) AND b1, NOT(a2) AND b2, NOT(a3) AND b3]
     */
    static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) { return _mm_andnot_ps(a, b); }

    /*
     * Round the packed single-precision (32-bit) floating-point elements in a down to an integer value, and store the
     * results as packed double-precision floating-point elements in vect_t.
     * Args   : [a0, a1, a2, a3]
     * Return : [floor(a0), floor(a1), floor(a2), floor(a3)]
     */
    static INLINE CONST vect_t floor(const vect_t a) { return _mm_floor_ps(a); }

    /*
     * Round the packed single-precision (32-bit) floating-point elements in a up to an integer value, and store the
     * results as packed single-precision floating-point elements in vect_t.
     * Args   : [a0, a1, a2, a3]
     * Return : [ceil(a0), ceil(a1), ceil(a2), ceil(a3)]
     */
    static INLINE CONST vect_t ceil(const vect_t a) { return _mm_ceil_ps(a); }

    /*
     * Round the packed single-precision (32-bit) floating-point elements in a, and store the results as packed
     * single-precision floating-point elements in vect_t.
     * Args   : [a0, a1, a2, a3]
     * Return : [round(a0), round(a1), round(a2), round(a3)]
     */
    static INLINE CONST vect_t round(const vect_t a) {
        return _mm_round_ps(a, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    }

    /*
     * Horizontally add adjacent pairs of single-precision (32-bit) floating-point elements in a and b, and pack the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
     * Return : [a0+a1, b0+b1, a2+a3, b2+b3]
     */
    static INLINE CONST vect_t hadd(const vect_t a, const vect_t b) { return _mm_hadd_ps(a, b); }

    /*
     * Horizontally add single-precision (32-bit) floating-point elements in a.
     * Args   : [a0, a1, a2, a3]
     * Return : a0+a1+a2+a3
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        return ((const scalar_t *)&a)[0] + ((const scalar_t *)&a)[1] + ((const scalar_t *)&a)[2] +
        ((const scalar_t *)&a)[3];
    }

    static INLINE vect_t mod(vect_t &C, const vect_t &P, const vect_t &INVP, const vect_t &NEGP, const vect_t &MIN,
                             const vect_t &MAX, vect_t &Q, vect_t &T) {
        FLOAT_MOD(C, P, INVP, Q);
        NORML_MOD(C, P, NEGP, MIN, MAX, Q, T);

        return C;
    }
#else // __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS
#error "You need SSE instructions to perform 128bits operations on double"
#endif // __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
};

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd128_float_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
