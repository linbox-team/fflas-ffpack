/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   Bastien Vialla<bastien.vialla@lirmm.fr>
 * Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 * Romain Lebreton <romain.lebreton@lirmm.fr>
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

#ifndef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
#error "You need SSE instructions to perform 128 bits operations on int16"
#endif

#include "givaro/givtypestring.h"
#include "fflas-ffpack/utils/align-allocator.h"
#include <vector>
#include <type_traits>

/*
 * Simd128 specialized for int16_t
 */
template <> struct Simd128_impl<true, true, true, 2> : public Simd128i_base {

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
     *  Broadcast 16-bit integer a to all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x,x,x,x,x,x,x,x] int16_t
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm_set1_epi16(x); }

    /*
     *  Set packed 16-bit integers in dst with the supplied values.
     *  Return [x0, ..., x7] int16_t
     */
    static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1, const scalar_t x2, const scalar_t x3,
                                   const scalar_t x4, const scalar_t x5, const scalar_t x6, const scalar_t x7) {
        return _mm_set_epi16(x7, x6, x5, x4, x3, x2, x1, x0);
    }

    /*
     *  Gather 16-bit integer elements with indexes idx[0], ..., idx[7] from the address p in vect_t.
     *  Return [p[idx[0]], ..., p[idx[7]]] int16_t
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        return set(p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]], p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]]);
    }

    /*
     * Load 128-bits of integer data from memory into dst.
     * p must be aligned on a 16-byte boundary or a general-protection exception will be generated.
     * Return [p[0], ..., p[7]] int16_t
     */
    static INLINE PURE vect_t load(const scalar_t *const p) {
        return _mm_load_si128(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Load 128-bits of integer data from memory into dst.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0], ..., p[7]] int16_t
     */
    static INLINE PURE vect_t loadu(const scalar_t *const p) {
        return _mm_loadu_si128(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Store 128-bits of integer data from a into memory.
     * p must be aligned on a 16-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void store(scalar_t *p, vect_t v) {
        _mm_store_si128(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Store 128-bits of integer data from a into memory.
     * p does not need to be aligned on any particular boundary.
     */
    static INLINE void storeu(scalar_t *p, vect_t v) {
        _mm_storeu_si128(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Store 128-bits of integer data from a into memory using a non-temporal memory hint.
     * p must be aligned on a 16-byte boundary or a general-protection exception may be generated.
     */
    static INLINE void stream(scalar_t *p, const vect_t v) {
        _mm_stream_si128(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Shift packed 16-bit integers in a left by s while shifting in zeros, and store the results in vect_t.
     * Args   :	[a0, ..., a7] int16_t
     * Return :	[a0 << s, a1 << s, a2 << s, a3 << s, a4 << s, a5 << s, a6 << s, a7 << s] int16_t
     */
    template<int s>
    static INLINE CONST vect_t sll(const vect_t a) { return _mm_slli_epi16(a, s); }

    /*
     * Shift packed 16-bit integers in a right by s while shifting in zeros, and store the results in vect_t.
     * Args   :	[a0, ..., a7] int16_t
     * Return :	[a0 >> s, a1 >> s, a2 >> s, a3 >> s, a4 >> s, a5 >> s, a6 >> s, a7 >> s] int16_t
     */
    template<int s>
    static INLINE CONST vect_t srl(const vect_t a) { return _mm_srli_epi16(a, s); }

    /*
     * Shift packed 16-bit integers in a right by s while shifting in sign bits, and store the results in vect_t.
     * Args   :	[a0, ..., a7] int16_t
     * Return :	[a0 >> s, a1 >> s, a2 >> s, a3 >> s, a4 >> s, a5 >> s, a6 >> s, a7 >> s] int16_t
     */
    template<int s>
    static INLINE CONST vect_t sra(const vect_t a) { return _mm_srai_epi16(a, s); }

    /*
     * Shuffle 16-bit integers in a using the control in imm8, and store the results in dst.
     * Args   :	[a0, ..., a7] int16_t
     * Return :	[a[s[0..3]], ..., a[s[28..31]] int16_t
     */
    template<uint32_t s>
    static INLINE CONST vect_t shuffle(const vect_t a) {
        //#pragma warning "The simd shuffle function is emulated, it may impact the performances.";
        Converter conv;
        conv.v = a;
        return set (conv.t[( s      & 0x0000000F)], conv.t[( s      & 0x000000F0)],
                    conv.t[((s>> 8) & 0x0000000F)], conv.t[((s>> 8) & 0x000000F0)],
                    conv.t[((s>>16) & 0x0000000F)], conv.t[((s>>16) & 0x000000F0)],
                    conv.t[((s>>24) & 0x0000000F)], conv.t[((s>>24) & 0x000000F0)]);
    }

    /*
     * Unpack and interleave 16-bit integers from the low half of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a0, b0, a1, b1, a2, b2, a3, b3 ]
     */
    static INLINE CONST vect_t
    unpacklo_intrinsic (const vect_t a, const vect_t b) {
        return _mm_unpacklo_epi16(a, b);
    }

    /*
     * Unpack and interleave 16-bit integers from the high half of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a4, b4, a5, b5, a6, b6, a7, b7 ]
     */
    static INLINE CONST vect_t
    unpackhi_intrinsic (const vect_t a, const vect_t b) {
        return _mm_unpackhi_epi16(a, b);
    }

    /*
     * Unpack and interleave 16-bit integers from the low half of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a0, b0, a1, b1, a2, b2, a3, b3 ]
     */
    static INLINE CONST vect_t unpacklo(const vect_t a, const vect_t b) {
        return unpacklo_intrinsic(a, b);
    }

    /*
     * Unpack and interleave 16-bit integers from the high half of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a4, b4, a5, b5, a6, b6, a7, b7 ]
     */
    static INLINE CONST vect_t unpackhi(const vect_t a, const vect_t b) {
        return unpackhi_intrinsic(a, b);
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
        lo = unpacklo (a, b);
        hi = unpackhi (a, b);
    }
    /*
     * Pack 16-bit integers from the even positions of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a0, a2, a4, a6, b0, b2, b4, b6 ]
     */
    static INLINE CONST vect_t pack_even (const vect_t a, const vect_t b) {
        vect_t idx = _mm_set_epi8 (15,14,11,10,7,6,3,2,13,12,9,8,5,4,1,0);
        vect_t t1 = _mm_shuffle_epi8 (a, idx);
        vect_t t2 = _mm_shuffle_epi8 (b, idx);
        return _mm_unpacklo_epi64 (t1, t2);
    }

    /*
     * Pack 16-bit integers from the odd positions of a and b.
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return:   [ a1, a3, a5, a7, b1, b3, b5, b7 ]
     */
    static INLINE CONST vect_t pack_odd (const vect_t a, const vect_t b) {
        vect_t idx = _mm_set_epi8 (15,14,11,10,7,6,3,2,13,12,9,8,5,4,1,0);
        vect_t t1 = _mm_shuffle_epi8 (a, idx);
        vect_t t2 = _mm_shuffle_epi8 (b, idx);
        return _mm_unpackhi_epi64 (t1, t2);
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
        vect_t idx = _mm_set_epi8 (15,14,11,10,7,6,3,2,13,12,9,8,5,4,1,0);
        vect_t t1 = _mm_shuffle_epi8 (a, idx);
        vect_t t2 = _mm_shuffle_epi8 (b, idx);
        even = _mm_unpacklo_epi64 (t1, t2);
        odd = _mm_unpackhi_epi64 (t1, t2);
    }

    /*
     * Transpose the 8x8 matrix formed by the 8 rows of 16-bit integers in r0,
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
        t0 = unpacklo_intrinsic (r0, r4);
        t1 = unpacklo_intrinsic (r1, r5);
        t2 = unpacklo_intrinsic (r2, r6);
        t3 = unpacklo_intrinsic (r3, r7);
        t4 = unpackhi_intrinsic (r0, r4);
        t5 = unpackhi_intrinsic (r1, r5);
        t6 = unpackhi_intrinsic (r2, r6);
        t7 = unpackhi_intrinsic (r3, r7);

        v0 = unpacklo_intrinsic (t0, t2);
        v1 = unpacklo_intrinsic (t1, t3);
        v4 = unpacklo_intrinsic (t4, t6);
        v5 = unpacklo_intrinsic (t5, t7);
        v2 = unpackhi_intrinsic (t0, t2);
        v3 = unpackhi_intrinsic (t1, t3);
        v6 = unpackhi_intrinsic (t4, t6);
        v7 = unpackhi_intrinsic (t5, t7);

        r0 = unpacklo_intrinsic (v0, v1);
        r2 = unpacklo_intrinsic (v2, v3);
        r4 = unpacklo_intrinsic (v4, v5);
        r6 = unpacklo_intrinsic (v6, v7);
        r1 = unpackhi_intrinsic (v0, v1);
        r3 = unpackhi_intrinsic (v2, v3);
        r5 = unpackhi_intrinsic (v4, v5);
        r7 = unpackhi_intrinsic (v6, v7);
    }

    /*
     * Blend 16-bit integers from a and b using control mask s.
     * Args: a = [ a0, ..., a7 ]
     *       b = [ b0, ..., b7 ]
     *       s = a 8-bit immediate integer
     * Return: [ s[0] ? a0 : b0, ..., s[7] ? a7 : b7 ]
     */
    template<uint8_t s>
    static INLINE CONST vect_t blend(const vect_t a, const vect_t b) {
        return _mm_blend_epi16(a, b, s);
    }

    /*
     * Add packed 16-bits integer in a and b, and store the results in vect_t.
     * Args   :	[a0, ..., a7] int16_t
     [b0, ..., b7] int16_t
     * Return :	[a0+b0, ..., a7+b7]   int16_t
     */
    static INLINE CONST vect_t add(const vect_t a, const vect_t b) { return _mm_add_epi16(a, b); }

    static INLINE vect_t addin(vect_t &a, const vect_t b) { return a = add(a, b); }

    /*
     * Subtract packed 16-bit integers in b from packed 16-bit integers in a, and store the results in vect_t.
     * Args   :	[a0, ..., a7] int16_t
     [b0, ..., b7] int16_t
     * Return :	[a0-b0, ..., a7-b7]  int16_t
     */
    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) { return _mm_sub_epi16(a, b); }

    static INLINE vect_t subin(vect_t &a, const vect_t b) { return a = sub(a, b); }

    /*
     * Multiply the packed 16-bit integers in a and b, producing intermediate 32-bit integers, and store the low 16 bits
     of the intermediate integers in vect_t.
     * Args   :	[a0, ..., a7] int16_t
     [b0, ..., b7] int16_t
     * Return :	[a0*b0 smod 2^16, ..., a7*b7 smod 2^16]	int16_t
     *	   where (a smod p) is the signed representant of a modulo p, that is -p/2 <= (a smod p) < p/2
     */
    static INLINE CONST vect_t mullo(const vect_t a, const vect_t b) { return _mm_mullo_epi16(a, b); }

    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) { return mullo(a, b); }

    /*
     * Multiply the packed 16-bit integers in a and b, producing intermediate 32-bit integers, and store the high 16
     bits of the intermediate integers in vect_t.
     * Args   :	[a0, ..., a7] int16_t
     [b0, ..., b7] int16_t
     * Return :	[Floor(a0*b0/2^16), ..., Floor(a7*b7/2^16)] int16_t
     */
    static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b) { return _mm_mulhi_epi16(a, b); }

    /*
     * Multiply the low 8-bit integers from each packed 16-bit element in a and b, and store the signed 16-bit results
     in vect_t.
     * Args   :	[a0, ..., a7] int16_t
     [b0, ..., b7] int16_t
     * Return :	[(a0 smod 2^8)*(b0 smod 2^8), ..., (a7 smod 2^8)*(b7 smod 2^8)]	int16_t
     *	   where (a smod p) is the signed representant of a modulo p, that is -p/2 <= (a smod p) < p/2
     */
    static INLINE CONST vect_t mulx(const vect_t a, const vect_t b) {
        //#pragma warning "The simd mulx function is emulated, it may impact the performances."
        vect_t a1, b1, mask1, mask2;
        mask1 = set1(0x00FF);
        mask2 = set1(0x0080);
        a1 = add(a,mask2);
        a1 = vand(a1,mask1);
        a1 = sub(a1,mask2);
        b1 = add(b,mask2);
        b1 = vand(b1,mask1);
        b1 = sub(b1,mask2);
        return mul(a1,b1);
    }

    /*
     * Multiply the packed 16-bit integers in a and b, producing intermediate 32-bit integers,
     * keep the low 16 bits of the intermediate and add the low 16-bits of c.
     * Args   :	[a0, ..., a7]		int16_t
     [b0, ..., b7]		int16_t
     [c0, ..., c7]		int16_t
     * Return :	[(a0*b0+c0) smod 2^16, ..., (a7*b7+c7) smod 2^16]	int16_t
     */
    static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b) { return add(c, mul(a, b)); }

    static INLINE vect_t fmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fmadd(c, a, b); }

    /*
     * Multiply the low 8-bit integers from each packed 16-bit element in a and b,
     * keep the signed 16-bit results and add the low 16-bits of c.
     * Args   :	[a0, ..., a7]		int16_t
     [b0, ..., b7]		int16_t
     [c0, ..., c7]		int16_t
     * Return :	[((a0 smod 2^8)*(b0 smod 2^8)+c0) smod 2^16, ...,
     *		 ((a7 smod 2^8)*(b7 smod 2^8)+c7) smod 2^16]	int16_t
     */
    static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b) { return add(c, mulx(a, b)); }

    static INLINE vect_t fmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmaddx(c, a, b); }

    /*
     * Multiply the packed 16-bit integers in a and b, producing intermediate 32-bit integers,
     * and substract the low 16 bits of the intermediate from elements of c.
     * Args   :	[a0, ..., a7]		int16_t
     [b0, ..., b7]		int16_t
     [c0, ..., c7]		int16_t
     * Return :	[(-a0*b0+c0) smod 2^16, ..., (-a7*b7+c7) smod 2^16]	int16_t
     */
    static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mul(a, b)); }

    static INLINE vect_t fnmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmadd(c, a, b); }

    /*
     * Multiply the low 8-bit integers from each packed 16-bit element in a and b,
     * keep the signed 16-bit results and substract them from elements of c.
     * Args   :	[a0, ..., a7]		int16_t
     [b0, ..., b7]		int16_t
     [c0, ..., c7]		int16_t
     * Return :	[(-(a0 smod 2^8)*(b0 smod 2^8)+c0) smod 2^16, ...,
     *		 (-(a7 smod 2^8)*(b7 smod 2^8)+c7) smod 2^16]		int16_t
     */
    static INLINE CONST vect_t fnmaddx(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mulx(a, b)); }

    static INLINE vect_t fnmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmaddx(c, a, b); }

    /*
     * Multiply packed 16-bit integers in a and b, producing intermediate 32-bit integers,
     * and substract elements of c to the low 16-bits of the intermediate.
     * Args   :	[a0, ..., a7]		int16_t
     [b0, ..., b7]		int16_t
     [c0, ..., c7]		int16_t
     * Return :	[(a0*b0-c0) smod 2^16, ..., (a7*b7-c7) smod 2^16]	int16_t
     */
    static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b) { return sub(mul(a, b), c); }

    static INLINE vect_t fmsubin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsub(c, a, b); }

    /*
     * Multiply the low 8-bit integers from each packed 16-bit element in a and b,
     * keep the signed 16-bit results and substract elements of c from them.
     * Args   :	[a0, ..., a7]		int16_t
     [b0, ..., b7]		int16_t
     [c0, ..., c7]		int16_t
     * Return :	[((a0 smod 2^8)*(b0 smod 2^8)-c0) smod 2^16, ...,
     *		 ((a7 smod 2^8)*(b7 smod 2^8)-c7) smod 2^16]		int16_t
     */
    static INLINE CONST vect_t fmsubx(const vect_t c, const vect_t a, const vect_t b) { return sub(mulx(a, b), c); }

    static INLINE vect_t fmsubxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsubx(c, a, b); }

    /*
     * Compare packed 16-bits in a and b for equality, and store the results in vect_t.
     * Args   :	[a0, ..., a7]		int16_t
     [b0, ..., b7]		int16_t
     * Return :	[(a0==b0) ? 0xFFFF : 0, ..., (a7==b7) ? 0xFFFF : 0]			int16_t
     */
    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) { return _mm_cmpeq_epi16(a, b); }

    /*
     * Compare packed 16-bits in a and b for greater-than, and store the results in vect_t.
     * Args   :	[a0, ..., a7]		int16_t
     [b0, ..., b7]		int16_t
     * Return :	[(a0>b0) ? 0xFFFF : 0, ..., (a7>b7) ? 0xFFFF : 0]			int16_t
     */
    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) { return _mm_cmpgt_epi16(a, b); }

    /*
     * Compare packed 16-bits in a and b for lesser-than, and store the results in vect_t.
     * Args   :	[a0, ..., a7]		int16_t
     [b0, ..., b7]		int16_t
     * Return :	[(a0<b0) ? 0xFFFF : 0, ..., (a7<b7) ? 0xFFFF : 0]			int16_t
     */
    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) { return _mm_cmplt_epi16(a, b); }

    /*
     * Compare packed 16-bits in a and b for greater or equal than, and store the results in vect_t.
     * Args   :	[a0, ..., a7]		int16_t
     [b0, ..., b7]		int16_t
     * Return :	[(a0>=b0) ? 0xFFFF : 0, ..., (a7>=b7) ? 0xFFFF : 0]			int16_t
     */
    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

    /*
     * Compare packed 16-bits in a and b for lesser or equal than, and store the results in vect_t.
     * Args   :	[a0, ..., a7]		int16_t
     [b0, ..., b7]		int16_t
     * Return :	[(a0<=b0) ? 0xFFFF : 0, ..., (a7<=b7) ? 0xFFFF : 0]			int16_t
     */
    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

    /*
     * Horizontally add 16-bits elements of a.
     * Args   :	[a0, a1, a2, a3, a4, a5, a6, a7]
     * Return :	a0+a1+a2+a3
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        Converter conv;
        conv.v = a;
        return scalar_t(conv.t[0] + conv.t[1] + conv.t[2] + conv.t[3] + conv.t[4] + conv.t[5] + conv.t[6] + conv.t[7]);
    }

    static INLINE CONST vect_t round(const vect_t a) { return a; }

    static INLINE vect_t mod(vect_t &C, const vect_t &P, const __m64 &INVP, const vect_t &NEGP, const vect_t &MIN,
                             const vect_t &MAX, vect_t &Q, vect_t &T) {
#ifdef __INTEL_COMPILER
        C = _mm_rem_epi16(C, P);
#else
        FFLASFFPACK_abort("not implemented");
#endif
        NORML_MOD(C, P, NEGP, MIN, MAX, Q, T);
        return C;
    }
};

/*
 * Simd128 specialized for uint16_t
 */
template <> struct Simd128_impl<true, true, false, 2> : public Simd128_impl<true, true, true, 2> {

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = uint16_t;

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
     *  Broadcast 16-bit unsigned integer a to all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x,x,x,x,x,x,x,x] uint16_t
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm_set1_epi16(x); }

    /*
     *  Broadcast 16-bit unsigned integer a to all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x0, ..., x7] uint16_t
     */
    static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1, const scalar_t x2, const scalar_t x3,
                                   const scalar_t x4, const scalar_t x5, const scalar_t x6, const scalar_t x7) {
        return _mm_set_epi16(x7, x6, x5, x4, x3, x2, x1, x0);
    }

    /*
     *  Gather 16-bit unsigned integer elements with indexes idx[0], ..., idx[7] from the address p in vect_t.
     *  Return [p[idx[0]],..., p[idx[7]]] uint16_t
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        return set(p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]], p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]]);
    }

    /*
     * Load 128-bits of unsigned integer data from memory into dst.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     * Return [p[idx[0]],..., p[idx[7]]] uint16_t
     */
    static INLINE PURE vect_t load(const scalar_t *const p) {
        return _mm_load_si128(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Load 128-bits of unsigned integer data from memory into dst.
     * p does not need to be aligned on any particular boundary.
     * Return [p[idx[0]],..., p[idx[7]]] uint16_t
     */
    static INLINE PURE vect_t loadu(const scalar_t *const p) {
        return _mm_loadu_si128(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Store 128-bits of unsigned integer data from a into memory.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void store(scalar_t *p, vect_t v) {
        _mm_store_si128(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Store 128-bits of unsigned integer data from a into memory.
     * p does not need to be aligned on any particular boundary.
     */
    static INLINE void storeu(scalar_t *p, vect_t v) {
        _mm_storeu_si128(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Store 128-bits of unsigned integer data from a into memory using a non-temporal memory hint.
     * p must be aligned on a 16-byte boundary or a general-protection exception may be generated.
     */
    static INLINE void stream(scalar_t *p, const vect_t v) {
        _mm_stream_si128(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Shift packed 16-bit unsigned integers in a right by s while shifting in sign bits, and store the results in vect_t.
     * Args   :	[a0, ..., a7]			uint16_t
     * Return :	[Floor(a0/2^s), ..., Floor(a7/2^s)]	int16_t
     */
    template<int s>
    static INLINE CONST vect_t sra(const vect_t a) { return _mm_srli_epi16(a, s); }

    static INLINE CONST vect_t greater(vect_t a, vect_t b) {
        vect_t x;
        x = set1(static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1));
        a = vxor(x, a);
        b = vxor(x, b);
        return _mm_cmpgt_epi16(a, b);
    }

    static INLINE CONST vect_t lesser(vect_t a, vect_t b) {
        vect_t x;
        x = set1(static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1));
        a = vxor(x, a);
        b = vxor(x, b);
        return _mm_cmplt_epi16(a, b);
    }

    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

    /*
     * Multiply the packed unsigned 16-bit integers in a and b, producing intermediate 32-bit integers,
     * and store the high 16 bits of the intermediate integers in vect_t.
     * Args   :	[a0, ..., a7] uint16_t
     *			[b0, ..., b7] uint16_t
     * Return :	[Floor(a0*b0/2^16), ..., Floor(a7*b7/2^16)] uint16_t
     */
    static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b) { return _mm_mulhi_epu16(a, b); }

    /*
     * Multiply the low unsigned 8-bit integers from each packed 16-bit element in a and b,
     * and store the signed 16-bit results in vect_t.
     * Args   :	[a0, ..., a7] uint16_t
     *			[b0, ..., b7] uint16_t
     * Return :	[(a0 mod 2^8)*(b0 mod 2^8), ..., (a7 mod 2^8)*(b7 mod 2^8)] uint16_t
     */
    static INLINE CONST vect_t mulx(const vect_t a, const vect_t b) {
        //#pragma warning "The simd mulx function is emulated, it may impact the performances."
        vect_t a1, b1, mask1;
        mask1 = set1(0x00FF);
        a1 = vand(a,mask1);
        b1 = vand(b,mask1);
        return mul(a1,b1);
    }

    static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b) { return add(c, mulx(a, b)); }

    static INLINE vect_t fmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmaddx(c, a, b); }

    static INLINE CONST vect_t fnmaddx(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mulx(a, b)); }

    static INLINE vect_t fnmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmaddx(c, a, b); }

    static INLINE CONST vect_t fmsubx(const vect_t c, const vect_t a, const vect_t b) { return sub(mulx(a, b), c); }

    static INLINE vect_t fmsubxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsubx(c, a, b); }

    /*
     * Horizontally add 16-bits elements of a.
     * Args   :	[a0, a1, a2, a3, a4, a5, a6, a7]
     * Return :	a0+a1+a2+a3
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        Converter conv;
        conv.v = a;
        return scalar_t(conv.t[0] + conv.t[1] + conv.t[2] + conv.t[3] + conv.t[4] + conv.t[5] + conv.t[6] + conv.t[7]);
    }
}; //Simd128_impl<true,true,false,2>

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd128_int16_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
