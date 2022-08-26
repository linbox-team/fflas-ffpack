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

#ifndef __FFLASFFPACK_simd512_double_INL
#define __FFLASFFPACK_simd512_double_INL

#if not (defined(__FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS))
#error "You need AVX512 instructions to perform 512bits operations on double"
#endif

#include "givaro/givtypestring.h"
#include "fflas-ffpack/utils/align-allocator.h"
#include <vector>
#include <type_traits>

/*
 * Simd512 specialized for double
 */
template <> struct Simd512_impl<true, false, true, 8> {
    /*
     * alias to 512 bit simd register
     */
    using vect_t = __m512d;

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = double;

    /*
     *	number of scalar_t in a simd register
     */
    static const constexpr size_t vect_size = 8;

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
     *  Return [0,0,0,0,0,0,0,0]
     */
    static INLINE CONST vect_t zero() {
        return _mm512_setzero_pd();
    }

    /*
     *	Broadcast double-precision (64-bit) floating-point value x to all elements of vect_t.
     *  Return [x,x,x,x,x,x,x,x]
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm512_set1_pd(x); }

    /*
     *	Set packed double-precision (64-bit) floating-point elements in vect_t with the supplied values.
     *  Return [x1,x2,x3,x4,x5,x6,x7,x8]
     */
    static INLINE CONST vect_t set(const scalar_t x1, const scalar_t x2, const scalar_t x3, const scalar_t x4, const scalar_t x5, const scalar_t x6, const scalar_t x7, const scalar_t x8) {
        return _mm512_set_pd(x8, x7, x6, x5, x4, x3, x2, x1);
    }

    /*
     *	Gather double-precision (64-bit) floating-point elements with indexes idx[0], ..., idx[7] from the address p in
     *vect_t.
     *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]], p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]]]
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        return _mm512_set_pd(p[idx[7]], p[idx[6]], p[idx[5]], p[idx[4]] ,p[idx[3]], p[idx[2]], p[idx[1]], p[idx[0]]);
    }

    /*
     * Load 512-bits (composed of 8 packed double-precision (64-bit) floating-point elements) from memory into vect_t.
     * p must be aligned on a 64-byte boundary or a general-protection exception will be generated.
     * Return [p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]]
     */
    static INLINE PURE vect_t load(const scalar_t *const p) { return _mm512_load_pd(p); }

    /*
     * Load 512-bits (composed of 8 packed double-precision (64-bit) floating-point elements) from memory into vect_t.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]]
     */
    static INLINE PURE vect_t loadu(const scalar_t *const p) { return _mm512_loadu_pd(p); }

    /*
     * Store 512-bits (composed of 8 packed double-precision (64-bit) floating-point elements) from p into memory.
     * p must be aligned on a 64-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void store(const scalar_t *p, const vect_t v) { _mm512_store_pd(const_cast<scalar_t *>(p), v); }

    /*
     * Store 512-bits (composed of 8 packed double-precision (64-bit) floating-point elements) from p into memory.
     * p does not need to be aligned on any particular boundary.
     */
    static INLINE void storeu(const scalar_t *p, const vect_t v) { _mm512_storeu_pd(const_cast<scalar_t *>(p), v); }

    /*
     * Store 512-bits (composed of 8 packed double-precision (64-bit) floating-point elements) from a into memory using
     * a non-temporal memory hint.
     * p must be aligned on a 64-byte boundary or a general-protection exception may be generated.
     */
    static INLINE void stream(const scalar_t *p, const vect_t v) { _mm512_stream_pd(const_cast<scalar_t *>(p), v); }

    /*
       Shuffle double-precision (64-bit) floating-point elements within 128-bit lanes using the control in imm8,
     * and store the results in dst.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] double
     *		   [b0, b1, b2, b3, b4, b5, b6, b7] double
     * Return : [a[s[0..1]], ..., a[s[14..15]]] double

     PROBLEME DANS L'IMPLEMENTATION, DES FONCTIONS EXISTENT POUR PERMUTER LE m512d MAIS CHAQUE PARTIE DU m512d NE PEUT PRENDRE QUE 2 VALEURS
     ALORS QUE POUR LE m256d IL POUVAIT EN PRENDRE 4 (DONC IL DEVRAIT POUVOIR PRENDRE LES 8 VALEURS DANS CETTE VERSION)*/


    template<uint8_t s>
    static INLINE CONST vect_t shuffle(const vect_t a) {
        return _mm512_permute_pd(a, s);
    }

    /*
     * Unpack and interleave double-precision (64-bit) floating-point elements
     * from the low half of each 128-bit lane in a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a0, b0, a2, b2, a4, b4, a6, b6 ]
     */
    static INLINE CONST vect_t
    unpacklo_intrinsic (const vect_t a, const vect_t b) {
        return _mm512_unpacklo_pd(a,b);
    }

    /*
     * Unpack and interleave double-precision (64-bit) floating-point elements
     * from the high half of each 128-bit lane in a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a1, b1, a3, b3, a5, b5, a7, b7 ]
     */
    static INLINE CONST vect_t
    unpackhi_intrinsic (const vect_t a, const vect_t b) {
        return _mm512_unpackhi_pd(a,b);
    }

    /*
     * Unpack and interleave double-precision (64-bit) floating-point elements
     * from the low half of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a0, b0, a1, b1, a2, b2, a3, b3 ]
     */
    static INLINE CONST vect_t unpacklo(const vect_t a, const vect_t b) {
        int64_t permute_idx[8] = { 0, 4, 1, 5, 2, 6, 3, 7 };
        __m512i s = _mm512_loadu_si512 (permute_idx);
        vect_t ta = _mm512_permutexvar_pd (s, a);
        vect_t tb = _mm512_permutexvar_pd (s, b);
        return _mm512_unpacklo_pd (ta, tb);
    }

    /*
     * Unpack and interleave double-precision (64-bit) floating-point elements
     * from the high half of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a4, b4, a5, b5, a6, b6, a7, b7 ]
     */
    static INLINE CONST vect_t unpackhi(const vect_t a, const vect_t b) {
        int64_t permute_idx[8] = { 0, 4, 1, 5, 2, 6, 3, 7 };
        __m512i s = _mm512_loadu_si512 (permute_idx);
        vect_t ta = _mm512_permutexvar_pd (s, a);
        vect_t tb = _mm512_permutexvar_pd (s, b);
        return _mm512_unpackhi_pd (ta, tb);
    }

    /*
     * Perform unpacklo and unpackhi with a and b and store the results in lo
     * and hi.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return: lo = [ a0, b0, a1, b1, a2, b2, a3, b3 ]
     *         hi = [ a4, b4, a5, b5, a6, b6, a7, b7 ]
     */
    static INLINE void
    unpacklohi (vect_t& lo, vect_t& hi, const vect_t a, const vect_t b) {
        int64_t permute_idx[8] = { 0, 4, 1, 5, 2, 6, 3, 7 };
        __m512i s = _mm512_loadu_si512 (permute_idx);
        vect_t ta = _mm512_permutexvar_pd (s, a);
        vect_t tb = _mm512_permutexvar_pd (s, b);
        lo = _mm512_unpacklo_pd (ta, tb);
        hi = _mm512_unpackhi_pd (ta, tb);
    }

    /*
     * Pack double-precision (64-bit) floating-point elements from the even
     * positions of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a0, a2, a4, a6, b0, b2, b4, b6 ]
     */
    static INLINE CONST vect_t pack_even (const vect_t a, const vect_t b) {
        int64_t permute_idx[8] = { 0, 2, 4, 6, 1, 3, 5, 7 };
        __m512i s = _mm512_loadu_si512 (permute_idx);
        vect_t lo = _mm512_unpacklo_pd (a, b);
        return _mm512_permutexvar_pd (s, lo);
    }

    /*
     * Pack double-precision (64-bit) floating-point elements from the odd
     * positions of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a1, a3, a5, a7, b1, b3, b5, b7 ]
     */
    static INLINE CONST vect_t pack_odd (const vect_t a, const vect_t b) {
        int64_t permute_idx[8] = { 0, 2, 4, 6, 1, 3, 5, 7 };
        __m512i s = _mm512_loadu_si512 (permute_idx);
        vect_t hi = _mm512_unpackhi_pd (a, b);
        return _mm512_permutexvar_pd (s, hi);
    }

    /*
     * Perform pack_even and pack_odd with a and b and store the results in even
     * and odd.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return: even = [ a0, a2, a4, a6, b0, b2, b4, b6 ]
     *         odd = [ a1, a3, a5, a7, b1, b3, b5, b7 ]
     */
    static INLINE void
    pack (vect_t& even, vect_t& odd, const vect_t a, const vect_t b) {
        int64_t permute_idx[8] = { 0, 2, 4, 6, 1, 3, 5, 7 };
        __m512i s = _mm512_loadu_si512 (permute_idx);
        vect_t lo = _mm512_unpacklo_pd (a, b);
        vect_t hi = _mm512_unpackhi_pd (a, b);
        even = _mm512_permutexvar_pd (s, lo);
        odd  = _mm512_permutexvar_pd (s, hi);
    }

    /*
     * Transpose the 8x8 matrix formed by the 8 rows of double-precision
     * (64-bit) floating-point elements in r0, r1, r2, r3, r4, r5, r6 and r7,
     * and store the transposed matrix in these vectors.
     * Args: r0 = [ r00, r01, r02, r03, r04, r05, r06, r07 ]
     *       r1 = [ r10, r11, r12, r13, r14, r15, r16, r17 ]
     *       ...                   ...                   ...
     *       r6 = [ r60, r61, r62, r63, r64, r65, r66, r67 ]
     *       r7 = [ r70, r71, r72, r73, r74, r75, r76, r77 ]
     * Return: r0 = [ r00, r10, r20, r30, r40, r50, r60, r70 ]
     *         r1 = [ r01, r11, r21, r31, r41, r51, r61, r71 ]
     *         ...                   ...                   ...
     *         r6 = [ r06, r16, r26, r36, r46, r56, r66, r76 ]
     *         r7 = [ r07, r17, r27, r37, r47, r57, r67, r77 ]
     */
    static INLINE void
    transpose (vect_t& r0, vect_t& r1, vect_t& r2, vect_t& r3, vect_t& r4,
               vect_t& r5, vect_t& r6, vect_t& r7) {
        vect_t t0, t1, t2, t3, t4, t5, t6, t7;
        vect_t v0, v1, v2, v3, v4, v5, v6, v7;

        int64_t permute_idx1[8] = { 0x0, 0x1, 0x8, 0x9, 0x4, 0x5, 0xc, 0xd };
        int64_t permute_idx2[8] = { 0x2, 0x3, 0xa, 0xb, 0x6, 0x7, 0xe, 0xf };
        int64_t permute_idx3[8] = { 0x0, 0x1, 0x2, 0x3, 0x8, 0x9, 0xa, 0xb };
        int64_t permute_idx4[8] = { 0x4, 0x5, 0x6, 0x7, 0xc, 0xd, 0xe, 0xf };

        __m512i i1 = _mm512_loadu_si512 (permute_idx1);
        t0 = unpacklo_intrinsic (r0, r1);
        t2 = unpacklo_intrinsic (r2, r3);
        t4 = unpacklo_intrinsic (r4, r5);
        t6 = unpacklo_intrinsic (r6, r7);
        t1 = unpackhi_intrinsic (r0, r1);
        t3 = unpackhi_intrinsic (r2, r3);
        t5 = unpackhi_intrinsic (r4, r5);
        t7 = unpackhi_intrinsic (r6, r7);

        __m512i i2 = _mm512_loadu_si512 (permute_idx2);
        v0 = _mm512_permutex2var_pd (t0, i1, t2);
        v1 = _mm512_permutex2var_pd (t1, i1, t3);
        v4 = _mm512_permutex2var_pd (t4, i1, t6);
        v5 = _mm512_permutex2var_pd (t5, i1, t7);
        __m512i i3 = _mm512_loadu_si512 (permute_idx3);
        v2 = _mm512_permutex2var_pd (t0, i2, t2);
        v3 = _mm512_permutex2var_pd (t1, i2, t3);
        v6 = _mm512_permutex2var_pd (t4, i2, t6);
        v7 = _mm512_permutex2var_pd (t5, i2, t7);

        __m512i i4 = _mm512_loadu_si512 (permute_idx4);
        r0 = _mm512_permutex2var_pd (v0, i3, v4);
        r1 = _mm512_permutex2var_pd (v1, i3, v5);
        r2 = _mm512_permutex2var_pd (v2, i3, v6);
        r3 = _mm512_permutex2var_pd (v3, i3, v7);
        r4 = _mm512_permutex2var_pd (v0, i4, v4);
        r5 = _mm512_permutex2var_pd (v1, i4, v5);
        r6 = _mm512_permutex2var_pd (v2, i4, v6);
        r7 = _mm512_permutex2var_pd (v3, i4, v7);
    }

    /*
     * Blend packed double-precision (64-bit) floating-point elements from a and
     * b using control mask s.
     * Args: a = [ a0, ..., a7 ]
     *       b = [ b0, ..., b7 ]
     *       s = a 8-bit mask
     * Return: [ s[0] ? a0 : b0, ..., s[7] ? a7 : b7 ]
     */
    template<uint8_t s>
    static INLINE CONST vect_t blend(const vect_t a, const vect_t b) {
        return _mm512_mask_blend_pd(s, a, b);
    }

    /*
     * Blend packed double-precision (64-bit) floating-point elements from a and
     * b using the vector mask as control.
     * Args: a = [ a0, ..., a7 ]
     *       b = [ b0, ..., b7 ]
     *       mask
     * Return: [ mask[31] ? a0 : b0, ..., mask[511] ? a7 : b7 ]
     */
    static INLINE CONST vect_t blendv(const vect_t a, const vect_t b, const vect_t mask) {
#ifdef __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS
        __mmask8 k = _mm512_movepi64_mask (_mm512_castpd_si512 (mask));
#else
        __mmask8 k = _mm512_cmplt_epi64_mask (_mm512_castpd_si512 (mask), _mm512_setzero_si512());
#endif
        return _mm512_mask_blend_pd (k, a, b);
    }

    /*
     * Add packed double-precision (64-bit) floating-point elements in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0+b0, a1+b1, a2+b2, a3+b3, a4+b4, a5+b5, a6+b6, a7+b7]
     */
    static INLINE CONST vect_t add(const vect_t a, const vect_t b) { return _mm512_add_pd(a, b); }

    static INLINE vect_t addin(vect_t &a, const vect_t b) { return a = add(a, b); }

    /*
     * Subtract packed double-precision (64-bit) floating-point elements in b from packed double-precision (64-bit)
     * floating-point elements in a, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0-b0, a1-b1, a2-b2, a3-b3, a4-b4, a5-b5, a6-b6, a7-b7]
     */
    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) {
        return _mm512_sub_pd(a, b); }

    static INLINE CONST vect_t subin(vect_t &a, const vect_t b) { return a = sub(a, b); }

    /*
     * Multiply packed double-precision (64-bit) floating-point elements in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0*b0, a1*b1, a2*b2, a3*b3, a4*b4, a5*b5, a6*b6, a7*b7]
     */
    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) { return _mm512_mul_pd(a, b); }

    static INLINE CONST vect_t mulin(vect_t &a, const vect_t b) { return a = mul(a, b); }

    /*
     * Divide packed double-precision (64-bit) floating-point elements in a by packed elements in b,
     * and store the results in dst.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0/b0, a1/b1, a2/b2, a3/b3, a4/b4, a5/b5, a6/b6, a7/b7]
     */
    static INLINE CONST vect_t div(const vect_t a, const vect_t b) { return _mm512_div_pd(a, b); }

    /*
     * Multiply packed double-precision (64-bit) floating-point elements in a and b, add the intermediate result to
     * packed elements in c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
     * Return : [a0*b0+c0, a1*b1+c1, a2*b2+c2, a3*b3+c3, a4*b4+c4, a5*b5+c5, a6*b6+c6, a7*b7+c7]
     */
    static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b) {
        return _mm512_fmadd_pd(a, b, c);
    }

    static INLINE CONST vect_t fmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fmadd(c, a, b); }

    /*
     * Multiply packed double-precision (64-bit) floating-point elements in a and b, add the negated intermediate result
     * to packed elements in c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
     * Return : [-(a0*b0)+c0, -(a1*b1)+c1, -(a2*b2)+c2, -(a3*b3)+c3, -(a4*b4)+c4, -(a5*b5)+c5, -(a6*b6)+c6, -(a7*b7)+c7]
     */
    static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b) {
        return _mm512_fnmadd_pd(a, b, c);
    }

    static INLINE CONST vect_t fnmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmadd(c, a, b); }

    /*
     * Multiply packed double-precision (64-bit) floating-point elements in a and b, subtract packed elements in c from
     * the intermediate result, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
     * Return : [a0*b0-c0, a1*b1-c1, a2*b2-c2, a3*b3-c3, a4*b4-c4, a5*b5-c5, a6*b6-c6, a7*b7-c7]
     */
    static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b) {
        return _mm512_fmsub_pd(a, b, c);
    }

    static INLINE CONST vect_t fmsubin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsub(c, a, b); }

    /*
     * Compare packed double-precision (64-bit) floating-point elements in a and b for equality, and store the results
     in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7],
     * 			[b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [(a0==b0) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a1==b1) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a2==b2) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a3==b3) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a4==b4) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a5==b5) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a6==b6) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a7==b7) ? 0xFFFFFFFFFFFFFFFF : 0]
     */
    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) {
        int64_t i = 0xFFFFFFFFFFFFFFFF;
        __m512i c = _mm512_set1_epi64(i);
        return _mm512_maskz_expand_pd(_mm512_cmp_pd_mask(a, b, _CMP_EQ_OQ), _mm512_castsi512_pd(c));
    }

    /*
     * Compare packed double-precision (64-bit) floating-point elements in a and b for lesser-than, and store the
     results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7],
     * 			[b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [(a0<b0) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a1<b1) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a2<b2) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a3<b3) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a4<b4) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a5<b5) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a6<b6) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a7<b7) ? 0xFFFFFFFFFFFFFFFF : 0]
     */
    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {
        int64_t i = 0xFFFFFFFFFFFFFFFF;
        __m512i c = _mm512_set1_epi64(i);
        return _mm512_maskz_expand_pd(_mm512_cmp_pd_mask(a, b, _CMP_LT_OS), _mm512_castsi512_pd(c));
    }

    /*
     * Compare packed double-precision (64-bit) floating-point elements in a and b for lesser or equal than, and store
     the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7],
     * 			[b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [(a0<=b0) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a1<=b1) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a2<=b2) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a3<=b3) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a4<=b4) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a5<=b5) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a6<=b6) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a7<=b7) ? 0xFFFFFFFFFFFFFFFF : 0]
     */
    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) {
        int64_t i = 0xFFFFFFFFFFFFFFFF;
        __m512i c = _mm512_set1_epi64(i);
        return _mm512_maskz_expand_pd(_mm512_cmp_pd_mask(a, b, _CMP_LE_OS), _mm512_castsi512_pd(c));
    }

    /*
     * Compare packed double-precision (64-bit) floating-point elements in a and b for greater-than, and store the
     results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7],
     * 			[b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [(a0>b0) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a1>b1) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a2>b2) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a3>b3) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a4>b4) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a5>b5) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a6>b6) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a7>b7) ? 0xFFFFFFFFFFFFFFFF : 0]
     */
    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {
        int64_t i = 0xFFFFFFFFFFFFFFFF;
        __m512i c = _mm512_set1_epi64(i);
        return _mm512_maskz_expand_pd(_mm512_cmp_pd_mask(a, b, _CMP_GT_OS), _mm512_castsi512_pd(c));
    }

    /*
     * Compare packed double-precision (64-bit) floating-point elements in a and b for greater or equal than, and store
     the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7],
     * 			[b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [(a0>=b0) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a1>=b1) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a2>=b2) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a3>=b3) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a4>=b4) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a5>=b5) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a6>=b6) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a7>=b7) ? 0xFFFFFFFFFFFFFFFF : 0]
     */
    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) {
        int64_t i = 0xFFFFFFFFFFFFFFFF;
        __m512i c = _mm512_set1_epi64(i);
        return _mm512_maskz_expand_pd(_mm512_cmp_pd_mask(a, b, _CMP_GE_OS), _mm512_castsi512_pd(c));
    }

#ifdef __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS
    /*
     * Compute the bitwise AND of packed double-precision (64-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0 AND b0, a1 AND b1, a2 AND b2, a3 AND b3, a4 AND b4, a5 AND b5, a6 AND b6, a7 AND b7]
     */
    static INLINE CONST vect_t vand(const vect_t a, const vect_t b) { return _mm512_and_pd(a, b); }

    /*
     * Compute the bitwise OR of packed double-precision (64-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0 OR b0, a1 OR b1, a2 OR b2, a3 OR b3, a4 OR b4, a5 OR b5, a6 OR b6, a7 OR b7]
     */
    static INLINE CONST vect_t vor(const vect_t a, const vect_t b) { return _mm512_or_pd(a, b); }

    /*
     * Compute the bitwise XOR of packed double-precision (64-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0 XOR b0, a1 XOR b1, a2 XOR b2, a3 XOR b3, a4 XOR b4, a5 XOR b5, a6 XOR b6, a7 XOR b7]
     */
    static INLINE CONST vect_t vxor(const vect_t a, const vect_t b) { return _mm512_xor_pd(a, b); }

    /*
     * Compute the bitwise NOT AND of packed double-precision (64-bit) floating-point elements in a and b, and store the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [NOT(a0) AND b0, NOT(a1) AND b1, NOT(a2) AND b2, NOT(a3) AND b3, NOT(a4) AND b4,
     * NOT(a5) AND b5, NOT(a6) AND b6, NOT(a7) AND b7]
     */
    static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) { return _mm512_andnot_pd(a, b); }
#endif /* __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS */

    /*
     * Round the packed double-precision (64-bit) floating-point elements in a down to an integer value, and store the
     * results as packed double-precision floating-point elements in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
     * Return : [floor(a0), floor(a1), floor(a2), floor(a3), floor(a4), floor(a5), floor(a6), floor(a7)]
     */
    static INLINE CONST vect_t floor(const vect_t a) { return _mm512_floor_pd(a); }

    /*
     * Round the packed double-precision (64-bit) floating-point elements in a up to an integer value, and store the
     * results as packed double-precision floating-point elements in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
     * Return : [ceil(a0), ceil(a1), ceil(a2), ceil(a3), ceil(a4), ceil(a5), ceil(a6), ceil(a7)]
     */
    static INLINE CONST vect_t ceil(const vect_t a) { return _mm512_ceil_pd(a); }

    /*
     * Round the packed double-precision (64-bit) floating-point elements in a, and store the results as packed
     * double-precision floating-point elements in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
     * Return : [round(a0), round(a1), round(a2), round(a3), round(a4), round(a5), round(a6), round(a7)]
     */
    static INLINE CONST vect_t round(const vect_t a) {
        return _mm512_roundscale_pd(a, _MM_FROUND_TO_NEAREST_INT);
    }

    /*
     * Horizontally add adjacent pairs of double-precision (64-bit) floating-point elements in a and b, and pack the
     * results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0+a1, b0+b1, a2+a3, b2+b3, a4+a5, b4+b5, a6+a7, b6+b7]

     A TESTER
     */
    static INLINE CONST vect_t hadd(const vect_t a, const vect_t b) {

        __m256d lowa  = _mm512_castpd512_pd256(a); //sépare le m512d en 2 m256d
        __m256d higha = _mm512_extractf64x4_pd(a,1);

        __m256d lowb  = _mm512_castpd512_pd256(b);
        __m256d highb = _mm512_extractf64x4_pd(b,1);

        __m256d reslow = _mm256_hadd_pd(lowa, lowb); //fait le hadd sur les deux m256d
        __m256d reshigh = _mm256_hadd_pd(higha, highb);

        __m512d res = _mm512_castpd256_pd512(reslow); //met les 2 m256d dans un m512d
        res = _mm512_insertf64x4(res, reshigh, 1);
        return res;
    }

    /*
     * Horizontally add double-precision (64-bit) floating-point elements in a.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
     * Return : a0+a1+a2+a3+a4+a5+a6+a7
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        return ((const scalar_t *)&a)[0] + ((const scalar_t *)&a)[1] + ((const scalar_t *)&a)[2] +
        ((const scalar_t *)&a)[3] + ((const scalar_t *)&a)[4] + ((const scalar_t *)&a)[5] +
        ((const scalar_t *)&a)[6] + ((const scalar_t *)&a)[7];
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

};

#endif // __FFLASFFPACK_simd256_double_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
