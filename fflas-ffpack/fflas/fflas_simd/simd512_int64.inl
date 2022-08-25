/*
 * Copyright (C) 2018 the FFLAS-FFPACK group
 *
 * Written by   Ozden Ozturk <ozden.ozturk@etu.univ-grenoble-alpes.fr>
 * Pierre Karpman <pierre.karpman@univ-grenoble-alpes.fr>
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

#ifndef _simd512_int64_INL
#define _simd512_int64_INL

#ifndef __FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS
#error "You need AVX512 instructions to perform 512bits operations on int64_t"
#endif

#include "givaro/givtypestring.h"
#include "fflas-ffpack/utils/align-allocator.h"
#include "fflas-ffpack/utils/bit_manipulation.h"
#include <vector>
#include <type_traits>

/*
 * Simd512 specialized for int64_t
 */
template <> struct Simd512_impl<true, true, true, 8> : public Simd512i_base {

    /*
     * alias to 512 bit simd register
     */
    using vect_t = __m512i;

    /*
     * alias to 512 bit simd register
     */
    using half_t = __m256i;

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = int64_t;

    /*
     * Simd256 for scalar_t, to deal half_t
     */
    using simdHalf = Simd256<scalar_t>;

    /*
     *  number of scalar_t in a simd register
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
     *  alignement required by scalar_t pointer to be loaded in a vect_t
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
     * Converter from vect_t to a tab.
     * exple:
     *	Converter conv;
     *	conv.v = a;
     *	scalart_t x = conv.t[1]
     */
    union Converter {
        vect_t v;
        scalar_t t[vect_size];
    };

    /*
     *  Broadcast 64-bit integer a to all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x,x,x,x,x,x,x,x] int64_t
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm512_set1_epi64(x); }

    /*
     *  Set packed 64-bit integers in dst with the supplied values.
     *  Return [x0,x1,x2,x3,x4,x5,x6,x7] int64_t
     */
    static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1, const scalar_t x2, const scalar_t x3, const scalar_t x4, const scalar_t x5, const scalar_t x6, const scalar_t x7) {
        return _mm512_set_epi64(x7, x6, x5, x4, x3, x2, x1, x0);
    }
    /*
     *  Set packed 64-bit integers in dst with the supplied values, and padd with 0s
     *  Return [x0,x1,x2,x3,0,0,0,0] int64_t
     */
    static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1, const scalar_t x2, const scalar_t x3) {
        return _mm512_set_epi64(scalar_t(0), scalar_t(0), scalar_t(0), scalar_t(0), x3, x2, x1, x0);
    }

    /*
     *  Gather 64-bit integer elements with indexes idx[0], ..., idx[7] from the address p in vect_t.
     *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]], p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]]] int64_t
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        return set(p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]], p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]]);
    }

    /*
     * Load 512-bits of integer data from memory into dst.
     * p must be aligned on a 64-byte boundary or a general-protection exception will be generated.
     * Return [p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]] int32_t
     */
    static INLINE PURE vect_t load(const scalar_t *const p) {
        //std::cerr<<"load simd512_int64"<<std::endl;
        return _mm512_load_si512(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Load 512-bits of integer data from memory into dst.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]] int64_t
     */
    static INLINE PURE vect_t loadu(const scalar_t *const p) {
        return _mm512_loadu_si512(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Store 512-bits of integer data from a into memory.
     * p must be aligned on a 64-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void store(scalar_t *p, vect_t v) {
        _mm512_store_si512(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Store 512-bits of integer data from a into memory following mask.
     * p must be aligned on a 64-byte boundary or a general-protection exception will be generated.
     */
    template<uint8_t k>
    static INLINE void maskstore(scalar_t *p, vect_t v) {
        _mm512_mask_store_epi64(p, k, v);
    }


    /*
     * Store 512-bits of integer data from a into memory.
     * p does not need to be aligned on any particular boundary.
     */
    static INLINE void storeu(scalar_t *p, vect_t v) {
        _mm512_storeu_si512(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Store 512-bits of integer data from a into memory using a non-temporal memory hint.
     * p must be aligned on a 64-byte boundary or a general-protection exception may be generated.
     */
    static INLINE void stream(scalar_t *p, const vect_t v) {
        _mm512_stream_si512(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Shift packed 64-bit integers in a left by s while shifting in zeros, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     * Return : [a0 << s, a1 << s, a2 << s, a3 << s, a4 << s, a5 << s, a6 << s, a7 << s] int64_t
     */
    template <int s>
    static INLINE CONST vect_t sll(const vect_t a) { return _mm512_slli_epi64(a, s); }

    /*
     * Shift packed 64-bit integers in a right by s while shifting in zeros, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     * Return : [a0 >> s, a1 >> s, a2 >> s, a3 >> s, a4 >> s, a5 >> s, a6 >> s, a7 >> s] int64_t
     */
    template <int s>
    static INLINE CONST vect_t srl(const vect_t a) { return _mm512_srli_epi64(a, s); }

    /*
     * Shift packed 64-bit integers in a right by s while shifting in sign bits, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     * Return : [a0 >> s, a1 >> s, a2 >> s, a3 >> s, a4 >> s, a5 >> s, a6 >> s, a7 >> s] int64_t
     */
    template <int s>
    static INLINE CONST vect_t sra(const vect_t a) { return _mm512_srai_epi64(a, s); }

    /*
     * Shuffle 64-bit integers in a using the control in imm8, and store the results in dst.
     * Args   : [a0, ..., a7] int32_t
     * Return : [a[s[0..1]], ..., a[s[14..15]]] int32_t
     */

    template<uint8_t s>
    static INLINE CONST vect_t shuffle(const vect_t a) {
        return _mm512_permutex_epi64(a, s);
    }

    /*
     * Unpack and interleave 64-bit integers from the low half of each 128-bit
     * lane in a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a0, b0, a2, b2, a4, b4, a6, b6 ]
     */
    static INLINE CONST vect_t
    unpacklo_intrinsic (const vect_t a, const vect_t b) {
        return _mm512_unpacklo_epi64(a, b);
    }

    /*
     * Unpack and interleave 64-bit integers from the high half of each 128-bit
     * lane in a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a1, b1, a3, b3, a5, b5, a7, b7 ]
     */
    static INLINE CONST vect_t
    unpackhi_intrinsic (const vect_t a, const vect_t b) {
        return _mm512_unpackhi_epi64(a, b);
    }

    /*
     * Unpack and interleave 64-bit integers from the low half of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a0, b0, a1, b1, a2, b2, a3, b3 ]
     */
    static INLINE CONST vect_t unpacklo(const vect_t a, const vect_t b) {
        int64_t permute_idx[8] = { 0, 4, 1, 5, 2, 6, 3, 7 };
        vect_t s = loadu (permute_idx);
        vect_t ta = _mm512_permutexvar_epi64 (s, a);
        vect_t tb = _mm512_permutexvar_epi64 (s, b);
        return _mm512_unpacklo_epi64 (ta, tb);
    }

    /*
     * Unpack and interleave 64-bit integers from the high half of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a4, b4, a5, b5, a6, b6, a7, b7 ]
     */
    static INLINE CONST vect_t unpackhi(const vect_t a, const vect_t b) {
        int64_t permute_idx[8] = { 0, 4, 1, 5, 2, 6, 3, 7 };
        vect_t s = loadu (permute_idx);
        vect_t ta = _mm512_permutexvar_epi64 (s, a);
        vect_t tb = _mm512_permutexvar_epi64 (s, b);
        return _mm512_unpackhi_epi64 (ta, tb);
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
        vect_t s = loadu (permute_idx);
        vect_t ta = _mm512_permutexvar_epi64 (s, a);
        vect_t tb = _mm512_permutexvar_epi64 (s, b);
        lo = _mm512_unpacklo_epi64 (ta, tb);
        hi = _mm512_unpackhi_epi64 (ta, tb);
    }

    /*
     * Pack 64-bit integers from the even positions of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a0, a2, a4, a6, b0, b2, b4, b6 ]
     */
    static INLINE CONST vect_t pack_even (const vect_t a, const vect_t b) {
        int64_t permute_idx[8] = { 0, 2, 4, 6, 1, 3, 5, 7 };
        vect_t s = loadu (permute_idx);
        vect_t lo = _mm512_unpacklo_epi64 (a, b);
        return _mm512_permutexvar_epi64 (s, lo);
    }

    /*
     * Pack 64-bit integers from the odd positions of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a1, a3, a5, a7, b1, b3, b5, b7 ]
     */
    static INLINE CONST vect_t pack_odd (const vect_t a, const vect_t b) {
        int64_t permute_idx[8] = { 0, 2, 4, 6, 1, 3, 5, 7 };
        vect_t s = loadu (permute_idx);
        vect_t hi = _mm512_unpackhi_epi64 (a, b);
        return _mm512_permutexvar_epi64 (s, hi);
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
        vect_t s = loadu (permute_idx);
        vect_t lo = _mm512_unpacklo_epi64 (a, b);
        vect_t hi = _mm512_unpackhi_epi64 (a, b);
        even = _mm512_permutexvar_epi64 (s, lo);
        odd  = _mm512_permutexvar_epi64 (s, hi);
    }

    /*
     * Transpose the 8x8 matrix formed by the 8 rows of 64-bit integers in r0,
     * r1, r2, r3, r4, r5, r6 and r7, and store the transposed matrix in these
     * vectors.
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

        vect_t i1 = loadu (permute_idx1);
        t0 = unpacklo_intrinsic (r0, r1);
        t2 = unpacklo_intrinsic (r2, r3);
        t4 = unpacklo_intrinsic (r4, r5);
        t6 = unpacklo_intrinsic (r6, r7);
        t1 = unpackhi_intrinsic (r0, r1);
        t3 = unpackhi_intrinsic (r2, r3);
        t5 = unpackhi_intrinsic (r4, r5);
        t7 = unpackhi_intrinsic (r6, r7);

        vect_t i2 = loadu (permute_idx2);
        v0 = _mm512_permutex2var_epi64 (t0, i1, t2);
        v1 = _mm512_permutex2var_epi64 (t1, i1, t3);
        v4 = _mm512_permutex2var_epi64 (t4, i1, t6);
        v5 = _mm512_permutex2var_epi64 (t5, i1, t7);
        vect_t i3 = loadu (permute_idx3);
        v2 = _mm512_permutex2var_epi64 (t0, i2, t2);
        v3 = _mm512_permutex2var_epi64 (t1, i2, t3);
        v6 = _mm512_permutex2var_epi64 (t4, i2, t6);
        v7 = _mm512_permutex2var_epi64 (t5, i2, t7);

        vect_t i4 = loadu (permute_idx4);
        r0 = _mm512_permutex2var_epi64 (v0, i3, v4);
        r1 = _mm512_permutex2var_epi64 (v1, i3, v5);
        r2 = _mm512_permutex2var_epi64 (v2, i3, v6);
        r3 = _mm512_permutex2var_epi64 (v3, i3, v7);
        r4 = _mm512_permutex2var_epi64 (v0, i4, v4);
        r5 = _mm512_permutex2var_epi64 (v1, i4, v5);
        r6 = _mm512_permutex2var_epi64 (v2, i4, v6);
        r7 = _mm512_permutex2var_epi64 (v3, i4, v7);
    }

    /*
     * Blend 64-bit integers from a and b using control mask s.
     * Args: a = [ a0, ..., a7 ]
     *       b = [ b0, ..., b7 ]
     *       s = a 8-bit mask
     * Return: [ s[0] ? a0 : b0, ..., s[7] ? a7 : b7 ]
     */
    template<uint8_t s>
    static INLINE CONST vect_t blend(const vect_t a, const vect_t b) {
        return _mm512_mask_blend_epi64(s, a, b);
    }

    /*
     * Add packed 64-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, b4, b5, b6, b7]	int64_t
     * Return : [a0+b0, a1+b1, a2+b2, a3+b3, a4+b4, a5+b5, a6+b6, a7+b7]   int64_t
     */
    static INLINE CONST vect_t add(const vect_t a, const vect_t b) {
        return _mm512_add_epi64(a, b);
    }

    static INLINE vect_t addin(vect_t &a, const vect_t b) { return a = add(a, b); }

    /*
     * Subtract packed 64-bits integers in b from packed 64-bits integers in a, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, b4, b5, b6, b7]	int64_t
     * Return : [a0-b0, a1-b1, a2-b2, a3-b3, a4-b4, a5-b5, a6-b6, a7-b7]  int64_t
     */
    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) {
        //std::cerr<<"sub in simd512_int64"<<std::endl;
        return _mm512_sub_epi64(a, b); }

    static INLINE vect_t subin(vect_t &a, const vect_t b) { return a = sub(a, b); }

    /*
     * Multiply the packed 64-bits integers in a and b, producing intermediate 128-bit integers, and store the low 64
     bits of the intermediate integers in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, a4, a5, a6, a7]	int64_t
     * Return : [a0*b0 smod 2^32, ..., a7*b7 smod 2^32]	int64_t
     *	   where (a smod p) is the signed representant of a modulo p, that is -p/2 <= (a smod p) < p/2
     */
    static INLINE CONST vect_t mullo(vect_t a, vect_t b) {
        //std::cerr<<"mullo simd512_int64"<<std::endl;
        return _mm512_mullo_epi64(a, b);
    }

    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) { return mullo(a, b); }

    /*
     * Multiply the packed 64-bits integers in a and b, producing intermediate 128-bit integers, and store the high 64
     bits of the intermediate integers in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, a4, a5, a6, a7]	int64_t
     * Return : [Floor(a0*b0/2^64), ..., Floor(a7*b7/2^64)] int64_t
     */

    static INLINE CONST vect_t mulhi(vect_t a, vect_t b) {
        //#pragma warning "The simd mulhi function is emulate, it may impact the performances."
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set(mulhi_64(ca.t[0], cb.t[0]), mulhi_64 (ca.t[1], cb.t[1]),
                   mulhi_64(ca.t[2], cb.t[2]), mulhi_64 (ca.t[3], cb.t[3]),
                   mulhi_64(ca.t[4], cb.t[4]), mulhi_64 (ca.t[5], cb.t[5]),
                   mulhi_64(ca.t[6], cb.t[6]), mulhi_64 (ca.t[7], cb.t[7]));
    }

    /*
     * Multiply the low 32-bits integers from each packed 64-bit element in a and b, and store the signed 64-bit results
     in dst.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, a4, a5, a6, a7]	int64_t
     * Return : [(a0 smod 2^32)*(b0 smod 2^32), ..., (a7 smod 2^32)*(b7 smod 2^32)]	int64_t
     *	   where (a smod p) is the signed representant of a modulo p, that is -p/2 <= (a smod p) < p/2
     */
    static INLINE CONST vect_t mulx(const vect_t a, const vect_t b) {
        //std::cerr<<"mulx simd512_int64"<<std::endl;
        return _mm512_mul_epi32(a, b); }

    /*
     * Multiply packed 64-bit integers in a and b, producing intermediate 128-bit integers, and add the low 64-bits of
     the intermediate with c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, a4, a5, a6, a7]	int64_t
     *	    	[c0, c1, c2, c3, c4, c5, c6, c7]	int64_t
     * Return : [(a0*b0+c0) smod 2^64, ..., (a7*b7+c7) smod 2^64]	int64_t
     */
    static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b) {
        return add(c, mul(a, b)); }

    static INLINE vect_t fmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fmadd(c, a, b); }

    /*
     * Multiply the low 32-bit integers from each packed 64-bit element in a and b,
     * keep the signed 64-bit results and add the low 64-bits of c.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, a4, a5, a6, a7]	int64_t
     *	    	[c0, c1, c2, c3, c4, c5, c6, c7]	int64_t
     * Return :	[((a0 smod 2^32)*(b0 smod 2^32)+c0) smod 2^64, ...,
     *		 ((a7 smod 2^32)*(b7 smod 2^32)+c7) smod 2^64]	int64_t
     */
    static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b) { return add(c, mulx(a, b)); }

    static INLINE vect_t fmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmaddx(c, a, b); }

    /*
     * Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers,
     * and substract the low 64 bits of the intermediate from elements of c.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, a4, a5, a6, a7]	int64_t
     *	    	[c0, c1, c2, c3, c4, c5, c6, c7]	int64_t
     * Return :	[(-a0*b0+c0) smod 2^64, ..., (-a7*b7+c7) smod 2^64]	int64_t
     */
    static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mul(a, b)); }

    static INLINE vect_t fnmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmadd(c, a, b); }

    /*
     * Multiply the low 32-bit integers from each packed 64-bit element in a and b,
     * keep the signed 64-bit results and substract them from elements of c.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, a4, a5, a6, a7]	int64_t
     *	    	[c0, c1, c2, c3, c4, c5, c6, c7]	int64_t
     * Return :	[(-(a0 smod 2^32)*(b0 smod 2^32)+c0) smod 2^64, ...,
     *		 (-(a7 smod 2^32)*(b7 smod 2^32)+c7) smod 2^64]	int64_t
     */
    static INLINE CONST vect_t fnmaddx(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mulx(a, b)); }

    static INLINE vect_t fnmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmaddx(c, a, b); }

    /*
     * Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers,
     * and substract elements of c to the low 64-bits of the intermediate.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, a4, a5, a6, a7]	int64_t
     *	    	[c0, c1, c2, c3, c4, c5, c6, c7]	int64_t
     * Return :	[(a0*b0-c0) smod 2^64, ..., (a7*b7-c7) smod 2^64]	int64_t
     */
    static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b) {
        //std::cerr<<"fmsub in simd512_int64"<<std::endl;
        return sub(mul(a, b), c); }

    static INLINE vect_t fmsubin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsub(c, a, b); }

    /*
     * Multiply the low 32-bit integers from each packed 64-bit element in a and b,
     * keep the signed 64-bit results and substract elements of c from them.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, a4, a5, a6, a7]	int64_t
     *	    	[c0, c1, c2, c3, c4, c5, c6, c7]	int64_t
     * Return :	[((a0 smod 2^32)*(b0 smod 2^32)-c0) smod 2^64, ...,
     *		 ((a7 smod 2^32)*(b7 smod 2^32)-c7) smod 2^64]	int64_t
     */
    static INLINE CONST vect_t fmsubx(const vect_t c, const vect_t a, const vect_t b) {

        //std::cerr<<"fmsubx in simd512_int64"<<std::endl;
        return sub(mulx(a, b), c); }

    static INLINE vect_t fmsubxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsubx(c, a, b); }

    /*
     * Compare packed 64-bits in a and b for equality, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, b4, b5, b6, b7]	int64_t
     * Return : [(a0==b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1==b1) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a2==b2) ? 0xFFFFFFFFFFFFFFFF : 0, (a3==b3) ? 0xFFFFFFFFFFFFFFFF : 0, ..., (a7==b7) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
     */
    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) {
        int64_t i = 0xFFFFFFFFFFFFFFFF;
        __m512i c = _mm512_set1_epi64(i);
        __m512i d = _mm512_maskz_expand_epi64(_mm512_cmpeq_epi64_mask(a, b), c);
        return d;
    }

    /*
     * Compare packed 64-bits in a and b for greater-than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, b4, b5, b6, b7]	int64_t
     * Return : [(a0>b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1>b1) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a2>b2) ? 0xFFFFFFFFFFFFFFFF : 0, (a3>b3) ? 0xFFFFFFFFFFFFFFFF : 0, ..., (a7>b7) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
     */
    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {
        //		std::cerr<<"COUCOU greater simd512_int64.inl"<<std::endl;
        int64_t i = 0xFFFFFFFFFFFFFFFF;
        __m512i c = _mm512_set1_epi64(i);
        __m512i d = _mm512_maskz_expand_epi64(_mm512_cmpgt_epi64_mask(a, b), c);

        return d;
    }

    /*
     * Compare packed 64-bits in a and b for lesser-than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, b4, b5, b6, b7]	int64_t
     * Return : [(a0<b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1<b1) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a2<b2) ? 0xFFFFFFFFFFFFFFFF : 0, (a3<b3) ? 0xFFFFFFFFFFFFFFFF : 0, ..., (a7<b7) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
     */
    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {
        int64_t i = 0xFFFFFFFFFFFFFFFF;
        __m512i c = _mm512_set1_epi64(i);
        __m512i d = _mm512_maskz_expand_epi64(_mm512_cmpgt_epi64_mask(b, a), c);
        return d;
    }

    /*
     * Compare packed 64-bits in a and b for greater or equal than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, b4, b5, b6, b7]	int64_t
     * Return : [(a0>=b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1>=b1) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a2>=b2) ? 0xFFFFFFFFFFFFFFFF : 0, (a3>=b3) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a4>=b4) ? 0xFFFFFFFFFFFFFFFF : 0, (a5>=b5) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a6>=b6) ? 0xFFFFFFFFFFFFFFFF : 0, (a7>=b7) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
     */
    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

    /*
     * Compare packed 64-bits in a and b for lesser or equal than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     *	    	[b0, b1, b2, b3, b4, b5, b6, b7]	int64_t
     * Return : [(a0<=b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1<=b1) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a2<=b2) ? 0xFFFFFFFFFFFFFFFF : 0, (a3<=b3) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a4<=b4) ? 0xFFFFFFFFFFFFFFFF : 0, (a5<=b5) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a6<=b6) ? 0xFFFFFFFFFFFFFFFF : 0, (a7<=b7) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
     */
    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

    /*
     * Horizontally add 64-bits elements of a.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int64_t
     * Return : a0+a1+a2+a3+a4+a5+a6+a7	int64_t
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        Converter ca;
        ca.v = a;
        return scalar_t(ca.t[0] + ca.t[1] + ca.t[2] + ca.t[3] + ca.t[4] + ca.t[5] + ca.t[6] + ca.t[7]);
    }

    static INLINE CONST vect_t round(const vect_t a) { return a; }

    // mask the high 32 bits of a 64 bits, that is 00000000FFFFFFFF
    static INLINE CONST vect_t mask_high() { return srl<32>(_mm512_set1_epi8(-1)); }

    static INLINE CONST vect_t mulhi_fast(vect_t x, vect_t y);

    static INLINE vect_t mod(vect_t &C, const __m512d &P, const __m512d &INVP, const __m512d &NEGP, const vect_t &POW50REM,
                             const __m512d &MIN, const __m512d &MAX, __m512d &Q, __m512d &T);

protected:
    /* return the sign where vect_t is seen as sixteen int32_t */
    static INLINE CONST vect_t signbits(const vect_t x) {
        vect_t signBits = sub(zero(), srl<4*sizeof(scalar_t)-1>(x));
        return signBits;
    }
}; // Simd512_impl<true, true, true, 8>

/*
 * Simd512 specialized for uint64_t
 */
template <> struct Simd512_impl<true, true, false, 8> : public Simd512_impl<true, true, true, 8> {

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = uint64_t;

    /*
     *  string describing the Simd struct
     */
    static const std::string type_string () {
        return "Simd" + std::to_string(8*vect_size*sizeof(scalar_t)) + "<"
                      + Givaro::TypeString<scalar_t>::get() + ">";
    }

    using aligned_allocator = AlignedAllocator<scalar_t, Alignment(alignment)>;
    using aligned_vector = std::vector<scalar_t, aligned_allocator>;

    /* To check compatibility with Modular struct */
    template <class Field>
    using is_same_element = std::is_same<typename Field::Element, scalar_t>;

    /*
     * Simd128 for scalar_t, to deal half_t
     */
    using simdHalf = Simd256<scalar_t>;

    /*
     * Converter from vect_t to a tab.
     * exple:
     *	Converter conv;
     *	conv.v = a;
     *	scalart_t x = conv.t[1]
     */
    union Converter {
        vect_t v;
        scalar_t t[vect_size];
    };

    /*
     *  Broadcast 64-bit unsigned integer a to all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x,x,x,x,x,x,x,x] uint64_t
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm512_set1_epi64(x); }

    /*
     *  Set packed 64-bit unsigned integers in dst with the supplied values.
     *  Return [x0,x1,x2,x3,x4,x5,x6,x7] uint64_t
     */
    static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1, const scalar_t x2, const scalar_t x3, const scalar_t x4, const scalar_t x5, const scalar_t x6, const scalar_t x7) {
        return _mm512_set_epi64(x7, x6, x5, x4, x3, x2, x1, x0);
    }

    /*
     *  Gather 64-bit unsigned integer elements with indexes idx[0], ..., idx[7] from the address p in vect_t.
     *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]], ..., p[idx[7]]] uint64_t
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        return set(p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]], p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]]);
    }

    /*
     * Load 512-bits of unsigned integer data from memory into dst.
     * p must be aligned on a 64-byte boundary or a general-protection exception will be generated.
     * Return [p[0],p[1],p[2],p[3], ..., p[7]] uint64_t
     */
    static INLINE PURE vect_t load(const scalar_t *const p) {
        return _mm512_load_si512(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Load 512-bits of unsigned integer data from memory into dst.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0],p[1],p[2],p[3], ..., p[7]] uint64_t
     */
    static INLINE PURE vect_t loadu(const scalar_t *const p) {
        return _mm512_loadu_si512(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Store 512-bits of unsigned integer data from a into memory.
     * p must be aligned on a 64-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void store(scalar_t *p, vect_t v) {
        _mm512_store_si512(reinterpret_cast<vect_t *>(p), v);
    }

    /* Store 512-bits of integer data from a into memory following mask.
     * p must be aligned on a 64-byte boundary or a general-protection exception will be generated.
     */
    template<uint8_t k>
    static INLINE void maskstore(scalar_t *p, vect_t v) {
        _mm512_mask_store_epi64(p, k, v);
    }

    /*
     * Store 512-bits of unsigned integer data from a into memory.
     * p does not need to be aligned on any particular boundary.
     */
    static INLINE void storeu(scalar_t *p, vect_t v) {
        _mm512_storeu_si512(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Store 512-bits of unsigned integer data from a into memory using a non-temporal memory hint.
     * p must be aligned on a 64-byte boundary or a general-protection exception may be generated.
     */
    static INLINE void stream(scalar_t *p, const vect_t v) {
        _mm512_stream_si512(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Shift packed 64-bit unsigned integers in a right by s while shifting in sign bits, and store the results in vect_t.
     * Args   : [a0, ..., a7]			uint64_t
     * Return : [Floor(a0/2^s), ..., Floor(a7/2^s)]	uint64_t
     */
    template<int s>
    static INLINE CONST vect_t sra(const vect_t a) { return _mm512_srli_epi64(a, s); }

    static INLINE CONST vect_t greater(vect_t a, vect_t b) {
        __m512i c = _mm512_set1_epi64(0xFFFFFFFFFFFFFFFFLL);
        return _mm512_maskz_expand_epi64(_mm512_cmpgt_epu64_mask(a, b), c);
    }

    static INLINE CONST vect_t lesser(vect_t a, vect_t b) {
        __m512i c = _mm512_set1_epi64(0xFFFFFFFFFFFFFFFFLL);
        return _mm512_maskz_expand_epi64(_mm512_cmplt_epu64_mask(a, b), c);
    }

    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

    /*
     * Multiply the packed 64-bits integers in a and b, producing intermediate 128-bit integers, and store the low 64
     bits of the intermediate integers in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	uint64_t
     [b0, b1, b2, b3, b4, b5, b6, b7]	uint64_t
     * Return : [a0*b0 mod 2^64, a1*b1 mod 2^64, a2*b2 mod 2^64, a3*b3 mod 2^64, ..., a7*b7 mod 2^64]		uint64_t
     */
    static INLINE CONST vect_t mullo(vect_t a, vect_t b) {
        //#pragma warning "The simd mullo function is emulate, it may impact the performances."
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set(ca.t[0] * cb.t[0], ca.t[1] * cb.t[1], ca.t[2] * cb.t[2], ca.t[3] * cb.t[3], ca.t[4] * cb.t[4], ca.t[5] * cb.t[5], ca.t[6] * cb.t[6], ca.t[7] * cb.t[7]);
    }

    /*
     * Multiply the packed 64-bits integers in a and b, producing intermediate 128-bit integers, and store the high 64
     bits of the intermediate integers in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	uint64_t
     [b0, b1, b2, b3, b4, b5, b6, b7]  	uint64_t
     * Return :
     */
    static INLINE CONST vect_t mulhi(vect_t a, vect_t b) {
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set(mulhi_u64(ca.t[0], cb.t[0]), mulhi_u64 (ca.t[1], cb.t[1]),
                   mulhi_u64(ca.t[2], cb.t[2]), mulhi_u64 (ca.t[3], cb.t[3]),
                   mulhi_u64(ca.t[4], cb.t[4]), mulhi_u64 (ca.t[5], cb.t[5]),
                   mulhi_u64(ca.t[6], cb.t[6]), mulhi_u64 (ca.t[7], cb.t[7]));
    }

    /*
     * Multiply the low 32-bits integers from each packed 64-bit element in a and b, and store the unsigned 64-bit
     results in dst.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	uint64_t
     [b0, b1, b2, b3, b4, b5, b6, b7]	uint64_t
     * Return : [a0*b0, a1*b1, a2*b2, a3*b3, ..., a7*b7] uint64_t
     */
    static INLINE CONST vect_t mulx(const vect_t a, const vect_t b) { return _mm512_mul_epu32(a, b); }

    static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b) { return add(c, mulx(a, b)); }

    static INLINE vect_t fmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmaddx(c, a, b); }

    static INLINE CONST vect_t fnmaddx(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mulx(a, b)); }

    static INLINE vect_t fnmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmaddx(c, a, b); }

    static INLINE CONST vect_t fmsubx(const vect_t c, const vect_t a, const vect_t b) {
        //std::cerr<<"fmsubx in simd512_int64"<<std::endl;
        return sub(mulx(a, b), c); }

    static INLINE vect_t fmsubxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsubx(c, a, b); }

    /*
     * Horizontally add 64-bits elements of a.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
     * Return : a0+a1+a2+a3+a4+a5+a6+a7
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        Converter ca;
        ca.v = a;
        return ca.t[0] + ca.t[1] + ca.t[2] + ca.t[3] + ca.t[4] + ca.t[5] + ca.t[6] + ca.t[7];
    }
}; // Simd512_impl<true, true, false, 8>

#define vect_t Simd512_impl<true, true, true, 8>::vect_t

// warning : may be off by 1 multiple, but we save a mul...
INLINE CONST vect_t Simd512_impl<true, true, true, 8>::mulhi_fast(vect_t x, vect_t y) {
    // unsigned mulhi starts:
    // x1 = xy_high = mulhiu_fast(x,y)
    const vect_t mask = mask_high();

    vect_t x0 = vand(x, mask), x1 = srl<32>(x);
    vect_t y0 = vand(y, mask), y1 = srl<32>(y);

    x0 = Simd512_impl<true, true, false, 8>::mulx(x0, y1); // x0y1
    y0 = Simd512_impl<true, true, false, 8>::mulx(x1, y0); // x1y0
    y1 = Simd512_impl<true, true, false, 8>::mulx(x1, y1); // x1y1

    x1 = vand(y0, mask);
    y0 = srl<32>(y0); // x1y0_lo = x1 // y1yo_hi = y0
    x1 = srl<32>(add(x1, x0));
    y0 = add(y1, y0);

    x1 = add(x1, y0);
    // unsigned mulhi ends

    // fixing signs
    x0 = vand(signbits(x), y);
    x1 = sub(x1, x0);
    x0 = vand(signbits(y), x);
    x1 = sub(x1, x0);
    // end fixing
    return x1;
}

// FIXME why cannot use Simd512<double>::vect_t in the declaration instead of __m512d?
// --**~~~~ Only suitable for use with Modular<int64_t> or ModularBalanced<int64_t>, so p <= max_cardinality < 2**33 ~~~~~**--
INLINE vect_t Simd512_impl<true, true, true, 8>::mod(vect_t &C, const __m512d &P, const __m512d &INVP, const __m512d &NEGP, const vect_t &POW50REM,
                                                     const __m512d &MIN, const __m512d &MAX, __m512d &Q, __m512d &T) {
    vect_t Cq50, Cr50, Ceq;
    __m512d nCmod;

    // nothing so special with 50; could be something else

    Cq50 = sra<50>(C);                      // Cq50[i] < 2**14
    Cr50 = set1(0x3FFFFFFFFFFFFLL);
    Cr50 = vand(C, Cr50);                   // Cr50[i] < 2**50

    Ceq = fmadd(Cr50, Cq50, POW50REM);      // Ceq[i] < 2**47 + 2**50 < 2**51; Ceq[i] ~ C[i] mod p

    nCmod = _mm512_cvtepi64_pd(Ceq);
    nCmod = Simd512<double>::mod(nCmod, P, INVP, NEGP, MIN, MAX, Q, T);

    C = _mm512_cvtpd_epi64(nCmod);

    return C;
}

#undef vect_t

#endif // __FFLASFFPACK_simd512_int64_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
