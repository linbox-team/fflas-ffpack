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

#ifndef __FFLASFFPACK_simd512_float_INL
#define __FFLASFFPACK_simd512_float_INL

#include "givaro/givtypestring.h"
#include "fflas-ffpack/utils/align-allocator.h"
#include <vector>
#include <type_traits>

/*
 * Simd512 specialized for float
 */
template <> struct Simd512_impl<true, false, true, 4> {
#if defined(__FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS)
    /*
     * alias to 512 bit simd register
     */
    using vect_t = __m512;

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = float;

    /*
     *	number of scalar_t in a simd register
     */
    static const constexpr size_t vect_size = 16;

    /*
     *  string describing the Simd struct
     */
    static const std::string type_string () {
        return "Simd" + std::to_string(8*vect_size*sizeof(scalar_t)) + "<"
                      + Givaro::TypeString<scalar_t>::get() + ">";
    }

    /*
     *	alignement required by scalar_t pointer to be loaded in a vect_t
     */
    static const constexpr size_t alignment = 64;
    using aligned_allocator = AlignedAllocator<scalar_t, Alignment(alignment)>;
    using aligned_vector = std::vector<scalar_t, aligned_allocator>;

    /* To check compatibility with Modular struct */
    template <class Field>
    using is_same_element = std::is_same<typename Field::Element, scalar_t>;

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
     *  Return [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
     */
    static INLINE CONST vect_t zero() { return _mm512_setzero_ps(); }

    /*
     *	Broadcast single-precision (32-bit) floating-point value x to all elements of vect_t.
     *  Return [x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x]
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm512_set1_ps(x); }

    /*
     *	Set packed single-precision (32-bit) floating-point elements in vect_t with the supplied values.
     *  Return [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16]
     */
    static INLINE CONST vect_t set(const scalar_t x1, const scalar_t x2, const scalar_t x3, const scalar_t x4,
                                   const scalar_t x5, const scalar_t x6, const scalar_t x7, const scalar_t x8,
                                   const scalar_t x9, const scalar_t x10, const scalar_t x11, const scalar_t x12,
                                   const scalar_t x13, const scalar_t x14, const scalar_t x15, const scalar_t x16) {
        return _mm512_set_ps(x16, x15, x14, x13, x12, x11, x10, x9, x8, x7, x6, x5, x4, x3, x2, x1);
    }

    /*
     *	Gather single-precision (32-bit) floating-point elements with indexes idx[0], ..., idx[15] from the address p in
     *vect_t.
     *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]], p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]], p[idx[8]],
     p[idx[9]], p[idx[10]], p[idx[11]], p[idx[12]], p[idx[13]], p[idx[14]], p[idx[15]]]
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        // TODO AVX2 Gather
        return _mm512_set_ps(p[idx[15]], p[idx[14]], p[idx[13]], p[idx[12]], p[idx[11]], p[idx[10]], p[idx[9]], p[idx[8]], p[idx[7]], p[idx[6]], p[idx[5]], p[idx[4]], p[idx[3]], p[idx[2]], p[idx[1]], p[idx[0]]);
    }

    /*
     * Load 512-bits (composed of 16 packed single-precision (32-bit) floating-point elements) from memory into vect_t.
     * p must be aligned on a 64-byte boundary or a general-protection exception will be generated.
     * Return [p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15]]
     */
    static INLINE PURE vect_t load(const scalar_t *const p) { return _mm512_load_ps(p); }

    /*
     * Load 512-bits (composed of 16 packed single-precision (32-bit) floating-point elements) from memory into vect_t.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15]]
     */
    static INLINE PURE vect_t loadu(const scalar_t *const p) { return _mm512_loadu_ps(p); }

    /*
     * Store 512-bits (composed of 16 packed single-precision (32-bit) floating-point elements) from a into memory.
     * p must be aligned on a 64-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void store(const scalar_t *p, const vect_t v) { _mm512_store_ps(const_cast<scalar_t *>(p), v); }

    /*
     * Store 512-bits (composed of 16 packed single-precision (32-bit) floating-point elements) from a into memory.
     * p does not need to be aligned on any particular boundary.
     */
    static INLINE void storeu(const scalar_t *p, const vect_t v) { _mm512_storeu_ps(const_cast<scalar_t *>(p), v); }

    /*
     * Store 512-bits (composed of 16 packed double-precision (32-bit) floating-point elements) from a into memory using
     * a non-temporal memory hint.
     * p must be aligned on a 64-byte boundary or a general-protection exception may be generated.
     */
    static INLINE void stream(const scalar_t *p, const vect_t v) { _mm512_stream_ps(const_cast<scalar_t *>(p), v); }

    /*
     * Shuffle single-precision (32-bit) floating-point elements in a within 128-bit lanes using the control in s,
     * and store the results in dst.
     * Args   :	[a0, ..., a7] float
     [b0, ..., b7] float
     * Return :	[a[s[0..3]], ..., a[s[28..31]]] float
     */
    template<uint8_t s>
    static INLINE CONST vect_t shuffle_twice(const vect_t a) {
        return _mm512_permute_ps(a, s);
    }

    /*
     * Unpack and interleave single-precision (32-bit) floating-point elements
     * from the low half of each 128-bit lane in a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, ..., a13, a14, a15 ]
     *       b = [ b0, b1, b2, b3, b4, b5, ..., b13, b14, b15 ]
     * Return:	 [ a0, b0, a1, b1, a4, b4, a5, b5, ..., a12, b12, a13, b13 ]
     */
    static INLINE CONST vect_t
    unpacklo_intrinsic (const vect_t a, const vect_t b) {
        return _mm512_unpacklo_ps(a,b);
    }

    /*
     * Unpack and interleave single-precision (32-bit) floating-point elements
     * from the high half of each 128-bit lane in a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, ..., a13, a14, a15 ]
     *       b = [ b0, b1, b2, b3, b4, b5, ..., b13, b14, b15 ]
     * Return:	 [ a2, b2, a3, b3, a6, b6, a7, b7, ..., a14, b14, a15, b15 ]
     */
    static INLINE CONST vect_t
    unpackhi_intrinsic (const vect_t a, const vect_t b) {
        return _mm512_unpackhi_ps(a,b);
    }

    /*
     * Unpack and interleave single-precision (32-bit) floating-point elements
     * from the low half of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, ..., a13, a14, a15 ]
     *       b = [ b0, b1, b2, b3, b4, b5, ..., b13, b14, b15 ]
     * Return:   [ a0, b0, a1, b1, ..., a6, b7, a7, b7 ]
     */
    static INLINE CONST vect_t unpacklo(const vect_t a, const vect_t b) {
        int32_t permute_idx[16] = { 0, 1, 8, 9, 2, 3, 10, 11, 4, 5, 12, 13, 6, 7, 14, 15  };
        __m512i s = _mm512_loadu_si512 (permute_idx);
        vect_t ta = _mm512_permutexvar_ps (s, a);
        vect_t tb = _mm512_permutexvar_ps (s, b);
        return _mm512_unpacklo_ps (ta, tb);
    }

    /*
     * Unpack and interleave single-precision (32-bit) floating-point elements
     * from the high half of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, ..., a13, a14, a15 ]
     *       b = [ b0, b1, b2, b3, b4, b5, ..., b13, b14, b15 ]
     * Return:   [ a8, b8, a9, b9, ..., a14, b14, a15, b15 ]
     */
    static INLINE CONST vect_t unpackhi(const vect_t a, const vect_t b) {
        int32_t permute_idx[16] = { 0, 1, 8, 9, 2, 3, 10, 11, 4, 5, 12, 13, 6, 7, 14, 15  };
        __m512i s = _mm512_loadu_si512 (permute_idx);
        vect_t ta = _mm512_permutexvar_ps (s, a);
        vect_t tb = _mm512_permutexvar_ps (s, b);
        return _mm512_unpackhi_ps (ta, tb);
    }

    /*
     * Perform unpacklo and unpackhi with a and b and store the results in lo
     * and hi.
     * Args: a = [ a0, a1, a2, a3, a4, a5, ..., a13, a14, a15 ]
     *       b = [ b0, b1, b2, b3, b4, b5, ..., b13, b14, b15 ]
     * Return: lo = [ a0, b0, a1, b1, ..., a6, b7, a7, b7 ]
     *         hi = [ a8, b8, a9, b9, ..., a14, b14, a15, b15 ]
     */
    static INLINE void
    unpacklohi (vect_t& lo, vect_t& hi, const vect_t a, const vect_t b) {
        int32_t permute_idx[16] = { 0, 1, 8, 9, 2, 3, 10, 11, 4, 5, 12, 13, 6, 7, 14, 15  };
        __m512i s = _mm512_loadu_si512 (permute_idx);
        vect_t ta = _mm512_permutexvar_ps (s, a);
        vect_t tb = _mm512_permutexvar_ps (s, b);
        lo = _mm512_unpacklo_ps (ta, tb);
        hi = _mm512_unpackhi_ps (ta, tb);
    }

    /*
     * Pack single-precision (32-bit) floating-point elements from the even
     * positions of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, ..., a13, a14, a15 ]
     *       b = [ b0, b1, b2, b3, b4, b5, ..., b13, b14, b15 ]
     * Return:   [ a0, a2, ..., a12, a14, b0, b2, ..., b12, b14 ]
     */
    static INLINE CONST vect_t pack_even (const vect_t a, const vect_t b) {
        int32_t permute_idx[16] = { 0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15 };
        __m512i s = _mm512_loadu_si512 (permute_idx);
        vect_t t1 = _mm512_shuffle_ps (a, b, 0x88); /* 0x88 = 2020 base_4 */
        return _mm512_permutexvar_ps (s, t1);
    }

    /*
     * Pack single-precision (32-bit) floating-point elements from the odd
     * positions of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, ..., a13, a14, a15 ]
     *       b = [ b0, b1, b2, b3, b4, b5, ..., b13, b14, b15 ]
     * Return:   [ a1, a3, ..., a13, a15, b1, b3, ..., b13, b15 ]
     */
    static INLINE CONST vect_t pack_odd (const vect_t a, const vect_t b) {
        int32_t permute_idx[16] = { 0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15 };
        __m512i s = _mm512_loadu_si512 (permute_idx);
        vect_t t2 = _mm512_shuffle_ps (a, b, 0xdd); /* 0xdd = 3131 base_4 */
        return _mm512_permutexvar_ps (s, t2);
    }

    /*
     * Perform pack_even and pack_odd with a and b and store the results in even
     * and odd.
     * Args: a = [ a0, a1, a2, a3, a4, a5, ..., a13, a14, a15 ]
     *       b = [ b0, b1, b2, b3, b4, b5, ..., b13, b14, b15 ]
     * Return: even = [ a0, a2, ..., a12, a14, b0, b2, ..., b12, b14 ]
     *         odd = [ a1, a3, ..., a13, a15, b1, b3, ..., b13, b15 ]
     */
    static INLINE void
    pack (vect_t& even, vect_t& odd, const vect_t a, const vect_t b) {
        int32_t permute_idx[16] = { 0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15 };
        __m512i s = _mm512_loadu_si512 (permute_idx);
        vect_t t1 = _mm512_shuffle_ps (a, b, 0x88); /* 0x88 = 2020 base_4 */
        vect_t t2 = _mm512_shuffle_ps (a, b, 0xdd); /* 0xdd = 3131 base_4 */
        even = _mm512_permutexvar_ps (s, t1);
        odd = _mm512_permutexvar_ps (s, t2);
    }

    /*
     * Transpose the 16x16 matrix formed by the 16 rows of single-precision
     * (32-bit) floating-point elements in r0, r1, r2, r3, ..., r14 and r15, and
     * store the transposed matrix in these vectors.
     * Args: r0 = [ r00, r01, r02, r03, r04, r05, r06, r07, ... ]
     *       r1 = [ r10, r11, r12, r13, r14, r15, r16, r17, ... ]
     *       ...                                            ...
     * Return: r0 = [ r00, r10, r20, r30, r40, r50, r60, r70, ... ]
     *         r1 = [ r01, r11, r21, r31, r41, r51, r61, r71, ... ]
     *         ...                                            ...
     */
    static INLINE void
    transpose (vect_t& r0, vect_t& r1, vect_t& r2, vect_t& r3, vect_t& r4,
               vect_t& r5, vect_t& r6, vect_t& r7, vect_t& r8, vect_t& r9,
               vect_t& r10, vect_t& r11, vect_t& r12, vect_t& r13, vect_t& r14,
               vect_t& r15) {
        vect_t t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11,t12,t13,t14,t15;
        vect_t v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11,v12,v13,v14,v15;

        int32_t permute_idx1[16] = { 0x00, 0x01, 0x02, 0x03,
                                     0x10, 0x11, 0x12, 0x13,
                                     0x08, 0x09, 0x0a, 0x0b,
                                     0x18, 0x19, 0x1a, 0x1b };
        int32_t permute_idx2[16] = { 0x04, 0x05, 0x06, 0x07,
                                     0x14, 0x15, 0x16, 0x17,
                                     0x0c, 0x0d, 0x0e, 0x0f,
                                     0x1c, 0x1d, 0x1e, 0x1f };
        int32_t permute_idx3[16] = { 0x00, 0x01, 0x02, 0x03,
                                     0x04, 0x05, 0x06, 0x07,
                                     0x10, 0x11, 0x12, 0x13,
                                     0x14, 0x15, 0x16, 0x17 };
        int32_t permute_idx4[16] = { 0x08, 0x09, 0x0a, 0x0b,
                                     0x0c, 0x0d, 0x0e, 0x0f,
                                     0x18, 0x19, 0x1a, 0x1b,
                                     0x1c, 0x1d, 0x1e, 0x1f };

        __m512i i1 = _mm512_loadu_si512 (permute_idx1);

        v0 = unpacklo_intrinsic (r0, r2);
        v1 = unpacklo_intrinsic (r1, r3);
        v4 = unpacklo_intrinsic (r4, r6);
        v5 = unpacklo_intrinsic (r5, r7);
        v8 = unpacklo_intrinsic (r8, r10);
        v9 = unpacklo_intrinsic (r9, r11);
        v12 = unpacklo_intrinsic (r12, r14);
        v13 = unpacklo_intrinsic (r13, r15);
        v2 = unpackhi_intrinsic (r0, r2);
        v3 = unpackhi_intrinsic (r1, r3);
        v6 = unpackhi_intrinsic (r4, r6);
        v7 = unpackhi_intrinsic (r5, r7);
        v10 = unpackhi_intrinsic (r8, r10);
        v11 = unpackhi_intrinsic (r9, r11);
        v14 = unpackhi_intrinsic (r12, r14);
        v15 = unpackhi_intrinsic (r13, r15);

        t0 = unpacklo_intrinsic (v0, v1);
        t2 = unpacklo_intrinsic (v2, v3);
        t4 = unpacklo_intrinsic (v4, v5);
        t6 = unpacklo_intrinsic (v6, v7);
        t8 = unpacklo_intrinsic (v8, v9);
        t10 = unpacklo_intrinsic (v10, v11);
        t12 = unpacklo_intrinsic (v12, v13);
        t14 = unpacklo_intrinsic (v14, v15);
        t1 = unpackhi_intrinsic (v0, v1);
        t3 = unpackhi_intrinsic (v2, v3);
        t5 = unpackhi_intrinsic (v4, v5);
        t7 = unpackhi_intrinsic (v6, v7);
        t9 = unpackhi_intrinsic (v8, v9);
        t11 = unpackhi_intrinsic (v10, v11);
        t13 = unpackhi_intrinsic (v12, v13);
        t15 = unpackhi_intrinsic (v14, v15);


        __m512i i2 = _mm512_loadu_si512 (permute_idx2);
        v0 = _mm512_permutex2var_ps (t0, i1, t4);
        v1 = _mm512_permutex2var_ps (t1, i1, t5);
        v2 = _mm512_permutex2var_ps (t2, i1, t6);
        v3 = _mm512_permutex2var_ps (t3, i1, t7);
        v8 = _mm512_permutex2var_ps (t8, i1, t12);
        v9 = _mm512_permutex2var_ps (t9, i1, t13);
        v10 = _mm512_permutex2var_ps (t10, i1, t14);
        v11 = _mm512_permutex2var_ps (t11, i1, t15);
        __m512i i3 = _mm512_loadu_si512 (permute_idx3);
        v4 = _mm512_permutex2var_ps (t0, i2, t4);
        v5 = _mm512_permutex2var_ps (t1, i2, t5);
        v6 = _mm512_permutex2var_ps (t2, i2, t6);
        v7 = _mm512_permutex2var_ps (t3, i2, t7);
        v12 = _mm512_permutex2var_ps (t8, i2, t12);
        v13 = _mm512_permutex2var_ps (t9, i2, t13);
        v14 = _mm512_permutex2var_ps (t10, i2, t14);
        v15 = _mm512_permutex2var_ps (t11, i2, t15);

        __m512i i4 = _mm512_loadu_si512 (permute_idx4);
        r0 = _mm512_permutex2var_ps (v0, i3, v8);
        r1 = _mm512_permutex2var_ps (v1, i3, v9);
        r2 = _mm512_permutex2var_ps (v2, i3, v10);
        r3 = _mm512_permutex2var_ps (v3, i3, v11);
        r4 = _mm512_permutex2var_ps (v4, i3, v12);
        r5 = _mm512_permutex2var_ps (v5, i3, v13);
        r6 = _mm512_permutex2var_ps (v6, i3, v14);
        r7 = _mm512_permutex2var_ps (v7, i3, v15);
        r8 = _mm512_permutex2var_ps (v0, i4, v8);
        r9 = _mm512_permutex2var_ps (v1, i4, v9);
        r10 = _mm512_permutex2var_ps (v2, i4, v10);
        r11 = _mm512_permutex2var_ps (v3, i4, v11);
        r12 = _mm512_permutex2var_ps (v4, i4, v12);
        r13 = _mm512_permutex2var_ps (v5, i4, v13);
        r14 = _mm512_permutex2var_ps (v6, i4, v14);
        r15 = _mm512_permutex2var_ps (v7, i4, v15);
    }

    /*
     * Blend packed single-precision (32-bit) floating-point elements from a and
     * b using control mask s.
     * Args: a = [ a0, ..., a15 ]
     *       b = [ b0, ..., b15 ]
     *       s = a 16-bit mask
     * Return: [ s[0] ? a0 : b0, ..., s[15] ? a15 : b15 ]
     */
    template<uint16_t s>
    static INLINE CONST vect_t blend(const vect_t a, const vect_t b) {
        return _mm512_mask_blend_ps(s, a, b);
    }

    /*
     * Blend packed single-precision (32-bit) floating-point elements from a and
     * b using the vector mask as control.
     * Args: a = [ a0, ..., a15 ]
     *       b = [ b0, ..., b15 ]
     *       mask
     * Return: [ mask[31] ? a0 : b0, ..., mask[511] ? a15 : b15 ]
     */
    static INLINE CONST vect_t blendv(const vect_t a, const vect_t b, const vect_t mask) {
#ifdef __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS
        __mmask16 k = _mm512_movepi32_mask (_mm512_castps_si512 (mask));
#else
        __mmask16 k = _mm512_cmplt_epi32_mask (_mm512_castps_si512 (mask), _mm512_setzero_si512());
#endif
        return _mm512_mask_blend_ps (k, a, b);
    }

    /*
     * Add packed single-precision (32-bit) floating-point elements in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [a0+b0, a1+b1, a2+b2, a3+b3, a4+b4, a5+b5, a6+b6, a7+b7, a8+b8, a9+b9, a10+b10, a11+b11, a12+b12, a13+b13, a14+b14, a15+b15]
     */
    static INLINE CONST vect_t add(const vect_t a, const vect_t b) { return _mm512_add_ps(a, b); }

    static INLINE vect_t addin(vect_t &a, const vect_t b) { return a = add(a, b); }

    /*
     * Subtract packed single-precision (32-bit) floating-point elements in b from packed single-precision (32-bit)
     * floating-point elements in a, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [a0-b0, a1-b1, a2-b2, a3-b3, a4-b4, a5-b5, a6-b6, a7-b7, a8-b8, a9-b9, a10-b10, a11-b11, a12-b12, a13-b13, a14-b14, a15-b15]
     */
    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) { return _mm512_sub_ps(a, b); }

    static INLINE CONST vect_t subin(vect_t &a, const vect_t b) { return a = sub(a, b); }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [a0*b0, a1*b1, a2*b2, a3*b3, a4*b4, a5*b5, a6*b6, a7*b7, a8*b8, a9*b9, a10*b10, a11*b11, a12*b12, a13*b13, a14*b14, a15*b15]
     */
    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) { return _mm512_mul_ps(a, b); }

    static INLINE CONST vect_t mulin(vect_t &a, const vect_t b) { return a = mul(a, b); }

    /*
     * Divide packed single-precision (32-bit) floating-point elements in a by packed elements in b,
     * and store the results in dst.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [a0/b0, a1/b1, a2/b2, a3/b3, a4/b4, a5/b5, a6/b6, a7/b7, a8/b8, a9/b9, a10/b10, a11/b11, a12/b12, a13/b13, a14/b14, a15/b15]
     */
    static INLINE CONST vect_t div(const vect_t a, const vect_t b) { return _mm512_div_ps(a, b); }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, add the intermediate result to
     * packed elements in c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15],
     *			[c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15]
     * Return : [a0*b0+c0, a1*b1+c1, a2*b2+c2, a3*b3+c3, a4*b4+c4, a5*b5+c5, a6*b6+c6, a7*b7+c7,
     *			a8*b8+c8, a9*b9+c9, a10*b10+c10, a11*b11+c11, a12*b12+c12, a13*b13+c13, a14*b14+c14, a15*b15+c15]
     */
    static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b) {
        return _mm512_fmadd_ps(a, b, c);
    }

    static INLINE CONST vect_t fmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fmadd(c, a, b); }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, add the negated intermediate result
     * to packed elements in c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15],
     *			[c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15]
     * Return : [-(a0*b0)+c0, -(a1*b1)+c1, -(a2*b2)+c2, -(a3*b3)+c3, -(a4*b4)+c4, -(a5*b5)+c5, -(a6*b6)+c6, -(a7*b7)+c7,
     *			-(a8*b8)+c8, -(a9*b9)+c9, -(a10*b10)+c10, -(a11*b11)+c11, -(a12*b12)+c12, -(a13*b13)+c13, -(a14*b14)+c14, -(a15*b15)+c15]
     */
    static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b) {
        return _mm512_fnmadd_ps(a, b, c);
    }

    static INLINE CONST vect_t fnmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmadd(c, a, b); }

    /*
     * Multiply packed single-precision (32-bit) floating-point elements in a and b, subtract packed elements in c from
     * the intermediate result, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15],
     *			[c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15]
     * Return : [a0*b0-c0, a1*b1-c1, a2*b2-c2, a3*b3-c3, a4*b4-c4, a5*b5-c5, a6*b6-c6, a7*b7-c7,
     *			a8*b8-c8, a9*b9-c9, a10*b10-c10, a11*b11-c11, a12*b12-c12, a13*b13-c13, a14*b14-c14, a15*b15-c15]
     */
    static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b) {
        return _mm512_fmsub_ps(a, b, c);
    }

    static INLINE CONST vect_t fmsubin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsub(c, a, b); }


    /*
     * Compare packed single-precision (32-bit) floating-point elements in a and b for equality, and store the results
     in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [(a0==b0) ? 0xFFFFFFFF : 0,
     (a1==b1) ? 0xFFFFFFFF : 0,
     (a2==b2) ? 0xFFFFFFFF : 0,
     (a3==b3) ? 0xFFFFFFFF : 0,
     (a4==b4) ? 0xFFFFFFFF : 0,
     (a5==b5) ? 0xFFFFFFFF : 0,
     (a6==b6) ? 0xFFFFFFFF : 0,
     (a7==b7) ? 0xFFFFFFFF : 0,
     (a8==b8) ? 0xFFFFFFFF : 0,
     (a9==b9) ? 0xFFFFFFFF : 0,
     (a10==b10) ? 0xFFFFFFFF : 0,
     (a11==b11) ? 0xFFFFFFFF : 0,
     (a12==b12) ? 0xFFFFFFFF : 0,
     (a13==b13) ? 0xFFFFFFFF : 0,
     (a14==b14) ? 0xFFFFFFFF : 0,
     (a15==b15) ? 0xFFFFFFFF : 0]


*/

    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) {
        int32_t i = 0xFFFFFFFF;
        __m512i c = _mm512_set1_epi32(i);
        return _mm512_maskz_expand_ps(_mm512_cmp_ps_mask(a, b, _CMP_EQ_OQ), _mm512_castsi512_ps(c));
    }

    /*
     * Compare packed single-precision (32-bit) floating-point elements in a and b for lesser-than, and store the
     results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [(a0<b0) ? 0xFFFFFFFF : 0,
     (a1<b1) ? 0xFFFFFFFF : 0,
     (a2<b2) ? 0xFFFFFFFF : 0,
     (a3<b3) ? 0xFFFFFFFF : 0,
     (a4<b4) ? 0xFFFFFFFF : 0,
     (a5<b5) ? 0xFFFFFFFF : 0,
     (a6<b6) ? 0xFFFFFFFF : 0,
     (a7<b7) ? 0xFFFFFFFF : 0,
     (a8<b8) ? 0xFFFFFFFF : 0,
     (a9<b9) ? 0xFFFFFFFF : 0,
     (a10<b10) ? 0xFFFFFFFF : 0,
     (a11<b11) ? 0xFFFFFFFF : 0,
     (a12<b12) ? 0xFFFFFFFF : 0,
     (a13<b13) ? 0xFFFFFFFF : 0,
     (a14<b14) ? 0xFFFFFFFF : 0,
     (a15<b15) ? 0xFFFFFFFF : 0]
     */
    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {
        int32_t i = 0xFFFFFFFF;
        __m512i c = _mm512_set1_epi32(i);
        return _mm512_maskz_expand_ps(_mm512_cmp_ps_mask(a, b, _CMP_LT_OS), _mm512_castsi512_ps(c));
    }

    /*
     * Compare packed single-precision (32-bit) floating-point elements in a and b for lesser or equal than, and store
     the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [(a0<=b0) ? 0xFFFFFFFF : 0,
     (a1<=b1) ? 0xFFFFFFFF : 0,
     (a2<=b2) ? 0xFFFFFFFF : 0,
     (a3<=b3) ? 0xFFFFFFFF : 0,
     (a4<=b4) ? 0xFFFFFFFF : 0,
     (a5<=b5) ? 0xFFFFFFFF : 0,
     (a6<=b6) ? 0xFFFFFFFF : 0,
     (a7<=b7) ? 0xFFFFFFFF : 0,
     (a8<=b8) ? 0xFFFFFFFF : 0,
     (a9<=b9) ? 0xFFFFFFFF : 0,
     (a10<=b10) ? 0xFFFFFFFF : 0,
     (a11<=b11) ? 0xFFFFFFFF : 0,
     (a12<=b12) ? 0xFFFFFFFF : 0,
     (a13<=b13) ? 0xFFFFFFFF : 0,
     (a14<=b14) ? 0xFFFFFFFF : 0,
     (a15<=b15) ? 0xFFFFFFFF : 0]
     */
    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) {
        int32_t i = 0xFFFFFFFF;
        __m512i c = _mm512_set1_epi32(i);
        return _mm512_maskz_expand_ps(_mm512_cmp_ps_mask(a, b, _CMP_LE_OS), _mm512_castsi512_ps(c));
    }

    /*
     * Compare packed single-precision (32-bit) floating-point elements in a and b for greater-than, and store the
     results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [(a0>b0) ? 0xFFFFFFFF : 0,
     (a1>b1) ? 0xFFFFFFFF : 0,
     (a2>b2) ? 0xFFFFFFFF : 0,
     (a3>b3) ? 0xFFFFFFFF : 0,
     (a4>b4) ? 0xFFFFFFFF : 0,
     (a5>b5) ? 0xFFFFFFFF : 0,
     (a6>b6) ? 0xFFFFFFFF : 0,
     (a7>b7) ? 0xFFFFFFFF : 0,
     (a8>b8) ? 0xFFFFFFFF : 0,
     (a9>b9) ? 0xFFFFFFFF : 0,
     (a10>b10) ? 0xFFFFFFFF : 0,
     (a11>b11) ? 0xFFFFFFFF : 0,
     (a12>b12) ? 0xFFFFFFFF : 0,
     (a13>b13) ? 0xFFFFFFFF : 0,
     (a14>b14) ? 0xFFFFFFFF : 0,
     (a15>b15) ? 0xFFFFFFFF : 0]
     */
    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {
        int32_t i = 0xFFFFFFFF;
        __m512i c = _mm512_set1_epi32(i);
        return _mm512_maskz_expand_ps(_mm512_cmp_ps_mask(a, b, _CMP_GT_OS), _mm512_castsi512_ps(c));
    }

    /*
     * Compare packed single-precision (32-bit) floating-point elements in a and b for greater or equal than, and store
     the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [(a0>=b0) ? 0xFFFFFFFF : 0,
     (a1>=b1) ? 0xFFFFFFFF : 0,
     (a2>=b2) ? 0xFFFFFFFF : 0,
     (a3>=b3) ? 0xFFFFFFFF : 0,
     (a4>=b4) ? 0xFFFFFFFF : 0,
     (a5>=b5) ? 0xFFFFFFFF : 0,
     (a6>=b6) ? 0xFFFFFFFF : 0,
     (a7>=b7) ? 0xFFFFFFFF : 0,
     (a8>=b8) ? 0xFFFFFFFF : 0,
     (a9>=b9) ? 0xFFFFFFFF : 0,
     (a10>=b10) ? 0xFFFFFFFF : 0,
     (a11>=b11) ? 0xFFFFFFFF : 0,
     (a12>=b12) ? 0xFFFFFFFF : 0,
     (a13>=b13) ? 0xFFFFFFFF : 0,
     (a14>=b14) ? 0xFFFFFFFF : 0,
     (a15>=b15) ? 0xFFFFFFFF : 0]
     */
    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) {
        int32_t i = 0xFFFFFFFF;
        __m512i c = _mm512_set1_epi32(i);
        return _mm512_maskz_expand_ps(_mm512_cmp_ps_mask(a, b, _CMP_GE_OS), _mm512_castsi512_ps(c));
    }

#ifdef __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS
    /*
     * Compute the bitwise AND of packed single-precision (32-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Without AVX512_DQ need to cast m512 into m512i and use _mm512_and_si512
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [a0 AND b0, a1 AND b1, a2 AND b2, a3 AND b3, a4 AND b4, a5 AND b5, a6 AND b6, a7 AND b7,
     *			a8 AND b8, a9 AND b9, a10 AND b10, a11 AND b11, a12 AND b12, a13 AND b13, a14 AND b14, a15 AND b15]
     */
    static INLINE CONST vect_t vand(const vect_t a, const vect_t b) { return _mm512_and_ps(a, b); }

    /*
     * Compute the bitwise OR of packed single-precision (32-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [a0 OR b0, a1 OR b1, a2 OR b2, a3 OR b3, a4 OR b4, a5 OR b5, a6 OR b6, a7 OR b7,
     *			a8 OR b8, a9 OR b9, a10 OR b10, a11 OR b11, a12 OR b12, a13 OR b13, a14 OR b14, a15 OR b15]
     */
    static INLINE CONST vect_t vor(const vect_t a, const vect_t b) { return _mm512_or_ps(a, b); }

    /*
     * Compute the bitwise XOR of packed single-precision (32-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [a0 XOR b0, a1 XOR b1, a2 XOR b2, a3 XOR b3, a4 XOR b4, a5 XOR b5, a6 XOR b6, a7 XOR b7,
     *			a8 XOR b8, a9 XOR b9, a10 XOR b10, a11 XOR b11, a12 XOR b12, a13 XOR b13, a14 XOR b14, a15 XOR b15]
     */
    static INLINE CONST vect_t vxor(const vect_t a, const vect_t b) { return _mm512_xor_ps(a, b); }

    /*
     * Compute the bitwise NOT AND of packed single-precision (32-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [NOT(a0) AND b0, NOT(a1) AND b1, NOT(a2) AND b2, NOT(a3) AND b3, NOT(a4) AND b4, NOT(a5) AND b5, NOT(a6) AND b6, NOT(a7) AND b7,
     *			NOT(a8) AND b8, NOT(a9) AND b9, NOT(a10) AND b10, NOT(a11) AND b11, NOT(a12) AND b12, NOT(a13) AND b13, NOT(a14) AND b14, NOT(a15) AND b15]
     */
    static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) { return _mm512_andnot_ps(a, b); }
#endif /* __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS */
    /*
     * Round the packed single-precision (32-bit) floating-point elements in a down to an integer value, and store the
     * results as packed double-precision floating-point elements in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]
     * Return : [floor(a0), floor(a1), floor(a2), floor(a3), floor(a4), floor(a5), floor(a6), floor(a7),
     *			floor(a8), floor(a9), floor(a10), floor(a11), floor(a12), floor(a13), floor(a14), floor(a15)]
     */
    static INLINE CONST vect_t floor(const vect_t a) { return _mm512_floor_ps(a); }

    /*
     * Round the packed single-precision (32-bit) floating-point elements in a up to an integer value, and store the
     * results as packed single-precision floating-point elements in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]
     * Return : [ceil(a0), ceil(a1), ceil(a2), ceil(a3), ceil(a4), ceil(a5), ceil(a6), ceil(a7),
     *			ceil(a8), ceil(a9), ceil(a10), ceil(a11), ceil(a12), ceil(a13), ceil(a14), ceil(a15)]
     */
    static INLINE CONST vect_t ceil(const vect_t a) { return _mm512_ceil_ps(a); }

    /*
     * Round the packed single-precision (32-bit) floating-point elements in a, and store the results as packed
     * single-precision floating-point elements in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]
     * Return : [round(a0), round(a1), round(a2), round(a3), round(a4), round(a5), round(a6), round(a7),
     *			round(a8), round(a9), round(a10), round(a11), round(a12), round(a13), round(a14), round(a15)]
     */
    static INLINE CONST vect_t round(const vect_t a) {
        return _mm512_roundscale_ps(a, _MM_FROUND_TO_NEAREST_INT);
    }

    /*
     * Horizontally add adjacent pairs of single-precision (32-bit) floating-point elements in a and b, and pack the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15],
     *			[b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]
     * Return : [a0+a1, b0+b1, a2+a3, b2+b3, a4+a5, b4+b5, a6+a7, b6+b7,
     *			a8+a9, b8+b9, a10+a11, b10+b11, a12+a13, b12+b13, a14+a15, b14+b15]

     A TESTER
     */
    static INLINE CONST vect_t hadd(const vect_t a, const vect_t b) {
        __m256 lowa  = _mm512_castps512_ps256(a);
        __m256 higha = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(a),1));

        __m256 lowb  = _mm512_castps512_ps256(b);
        __m256 highb = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(b),1));

        __m256 reslow = _mm256_hadd_ps(lowa, lowb);
        __m256 reshigh = _mm256_hadd_ps(higha, highb);

        __m512 res = _mm512_castps256_ps512(reslow);
        res = _mm512_insertf32x8(res, reshigh, 1);

        return res;
    }

    /*
     * Horizontally add single-precision (32-bit) floating-point elements in a.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]
     * Return : a0+a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14+a15
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        return ((const scalar_t *)&a)[0] + ((const scalar_t *)&a)[1] + ((const scalar_t *)&a)[2] +
        ((const scalar_t *)&a)[3] + ((const scalar_t *)&a)[4] + ((const scalar_t *)&a)[5] +
        ((const scalar_t *)&a)[6] + ((const scalar_t *)&a)[7] +  ((const scalar_t *)&a)[8] +
        ((const scalar_t *)&a)[9] + ((const scalar_t *)&a)[10] + ((const scalar_t *)&a)[11] +
        ((const scalar_t *)&a)[12] + ((const scalar_t *)&a)[13] + ((const scalar_t *)&a)[14] +
        ((const scalar_t *)&a)[15];
    }

#ifdef __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS
    /* Call NORML_MOD which needs vand which is not defined without AVX512DQ */
    static INLINE vect_t mod(vect_t &C, const vect_t &P, const vect_t &INVP, const vect_t &NEGP, const vect_t &MIN,
                             const vect_t &MAX, vect_t &Q, vect_t &T) {
        FLOAT_MOD(C, P, INVP, Q);
        NORML_MOD(C, P, NEGP, MIN, MAX, Q, T);

        return C;
    }
#endif /* __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS */

#else // __FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS
#error "You need AVX512 instructions to perform 512bits operations on float"
#endif
};

#endif // __FFLASFFPACK_simd512_float_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
