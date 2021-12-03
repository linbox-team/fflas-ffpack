/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   Bastien Vialla<bastien.vialla@lirmm.fr>
 * Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 * Romain Lebreton <romain.lebreton@lirmm.fr>
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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd128_int64_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd128_int64_INL

#ifndef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
#error "You need SSE instructions to perform 128 bits operations on int64"
#endif

#include "givaro/givtypestring.h"
#include "fflas-ffpack/utils/align-allocator.h"
#include "fflas-ffpack/utils/bit_manipulation.h"
#include <vector>
#include <type_traits>

/*
 * Simd128 specialized for int64_t
 */
template <> struct Simd128_impl<true, true, true, 8> : public Simd128i_base {

    /*
     * alias to 128 bit simd register
     */
    using vect_t = __m128i;

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = int64_t;

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
     *  Return [x,x] int64_t
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm_set1_epi64x(x); }

    /*
     *  Set packed 64-bit integers in dst with the supplied values.
     *  Return [x0,x1] int64_t
     */
    static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1) { return _mm_set_epi64x(x1, x0); }

    /*
     *  Gather 64-bit integer elements with indexes idx[0], idx[1] from the address p in vect_t.
     *  Return [p[idx[0]], p[idx[1]]] int64_t
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        return set(p[idx[0]], p[idx[1]]);
    }

    /*
     *  Extract one 64-bit integer from src at *_immediate_* index idx
     *  Return v[idx] int64_t
     */
    template<int idx>
    static INLINE CONST scalar_t get(vect_t v) {
        return _mm_extract_epi64(v, idx);
    }

    /*
     * Load 128-bits of integer data from memory into dst.
     * p must be aligned on a 16-byte boundary or a general-protection exception will be generated.
     * Return [p[0],p[1]] int64_t
     */
    static INLINE PURE vect_t load(const scalar_t *const p) {
        return _mm_load_si128(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Load 128-bits of integer data from memory into dst.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0],p[1]] int64_t
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
     * Shift packed 64-bit integers in a left by s while shifting in zeros, and store the results in vect_t.
     * Args   : [a0, a1] int64_t
     * Return : [a0 << s, a1 << s] int64_t
     */
    template<int s>
    static INLINE CONST vect_t sll(const vect_t a) { return _mm_slli_epi64(a, s); }

    /*
     * Shift packed 64-bit integers in a right by s while shifting in zeros, and store the results in vect_t.
     * Args   : [a0, a1] int64_t
     * Return : [a0 >> s, a1 >> s] int64_t
     */
    template<int s>
    static INLINE CONST vect_t srl(const vect_t a) { return _mm_srli_epi64(a, s); }

    /*
     * Shift packed 64-bit integers in a right by s while shifting in sign bits, and store the results in vect_t.
     * Args   : [a0, a1] int64_t
     * Return : [a0 >> s, a1 >> s] int64_t
     */
    template<int s>
    static INLINE CONST vect_t sra(const vect_t a) {
#if defined(__FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS) and defined(__FFLASFFPACK_HAVE_AVX512VL_INSTRUCTIONS)
        return _mm_srai_epi64(a, s);
#else
        vect_t m = sll<63-s>(set1(1));
        vect_t x = srl<s>(a);
        vect_t result = sub(vxor(x, m), m); // result = x^m - m
        return result;
#endif // __FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS and __FFLASFFPACK_HAVE_AVX512VL_INSTRUCTIONS
    }

    /*
     * Shuffle 64-bit integers in a using the control in imm8, and store the results in dst.
     * Args   : [a0, a1] int64_t
     * Return : [a[s[0]], a[s[1]]] int64_t
     */
    template<uint8_t s>
    static INLINE CONST vect_t shuffle(const vect_t a) {
        // Transform s = [d1 d0]_base2 to s1 = [2*d1+1 2*d1 2*d0+1 2*d0]_base4
        constexpr uint8_t s1 = ((s & 1)?(3*4+2):(1*4+0))+16*((s & 2)?(3*4+2):(1*4+0));
        return _mm_shuffle_epi32(a, s1);
    }

    /*
     * Unpack and interleave 64-bit integers from the low half of a and b.
     * Args: a = [ a0, a1 ]
     *       b = [ b0, b1 ]
     * Return:   [ a0, b0 ]
     */
    static INLINE CONST vect_t
    unpacklo_intrinsic (const vect_t a, const vect_t b) {
        return _mm_unpacklo_epi64(a, b);
    }

    /*
     * Unpack and interleave 64-bit integers from the high half of a and b.
     * Args: a = [ a0, a1 ]
     *       b = [ b0, b1 ]
     * Return:   [ a1, b1 ]
     */
    static INLINE CONST vect_t
    unpackhi_intrinsic (const vect_t a, const vect_t b) {
        return _mm_unpackhi_epi64(a, b);
    }

    /*
     * Unpack and interleave 64-bit integers from the low half of a and b.
     * Args: a = [ a0, a1 ]
     *       b = [ b0, b1 ]
     * Return:   [ a0, b0 ]
     */
    static INLINE CONST vect_t unpacklo(const vect_t a, const vect_t b) {
        return unpacklo_intrinsic(a, b);
    }

    /*
     * Unpack and interleave 64-bit integers from the high half of a and b.
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
     * Pack 64-bit integers from the even positions of a and b.
     * Args: a = [ a0, a1 ]
     *       b = [ b0, b1 ]
     * Return:   [ a0, b0 ]
     */
    static INLINE CONST vect_t pack_even (const vect_t a, const vect_t b) {
        return unpacklo (a, b); /* same as unpacklo for vect_size = 2 */
    }

    /*
     * Pack 64-bit integers from the odd positions of a and b.
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
     * Transpose the 2x2 matrix formed by the 2 rows of 64-bit integers in r0
     * and r1, and store the transposed matrix in these vectors.
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
     * Blend 64-bit integers from a and b using control mask s.
     * Args: a = [ a0, a1 ]
     *       b = [ b0, b1 ]
     *       s = a 2-bit immediate integer
     * Return: [ s[0] ? a0 : b0, s[1] ? a1 : b1 ]
     */
    template<uint8_t s>
    static INLINE CONST vect_t blend(const vect_t a, const vect_t b) {
#ifdef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
        /* If we have AVX2, we can use _mm_blend_epi32, which is faster.
         * We need to transform s = [d1 d0]_base2
         * into s1 = [d1 d1 d0 d0]_base2
         */
        constexpr uint8_t s1 = (s & 0x1) * 3 + (((s & 0x2) << 1)*3);
        return _mm_blend_epi32(a, b, s1);
#else
        /* If we only have SSE4.1, we need to use _mm_blend_epi16.
         * We need to transform s = [d1 d0]_base2
         * into s1 = [d1 d1 d1 d1 d0 d0 d0 d0]_base2
         */
        constexpr uint8_t s1 = (s & 0x1) * 15 + ((s & 0x2) << 3) * 15;
        return _mm_blend_epi16(a, b, s1);
#endif
    }

    /*
     * Add packed 64-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     * Return : [a0+b0, a1+b1]   int64_t
     */
    static INLINE CONST vect_t add(const vect_t a, const vect_t b) { return _mm_add_epi64(a, b); }

    static INLINE vect_t addin(vect_t &a, const vect_t b) { return a = add(a, b); }

    /*
     * Subtract packed 64-bit integers in b from packed 64-bit integers in a, and store the results in vect_t.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     * Return : [a0-b0, a1-b1]  int64_t
     */
    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) { return _mm_sub_epi64(a, b); }

    static INLINE vect_t subin(vect_t &a, const vect_t b) { return a = sub(a, b); }

    /*
     * Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers, and store the low 64
     bits of the intermediate integers in vect_t.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     * Return : [a0*b0 smod 2^64, a1*b1 smod 2^64] int64_t
     *	   where (a smod p) is the signed representant of a modulo p, that is -p/2 <= (a smod p) < p/2
     */
    static INLINE CONST vect_t mullo(const vect_t x0, const vect_t x1) {
#if defined(__FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS) and defined(__FFLASFFPACK_HAVE_AVX512VL_INSTRUCTIONS)
        return _mm_mullo_epi64(x0, x1);
#else
        // _mm_mullo_epi64 emul
        //#pragma warning "The simd mullo function is emulate, it may impact the performances."
        Converter c0, c1;
        c0.v = x0;
        c1.v = x1;
        return set((scalar_t)(c0.t[0] * c1.t[0]), (scalar_t)(c0.t[1] * c1.t[1]));
#endif //  __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS and __FFLASFFPACK_HAVE_AVX512VL_INSTRUCTIONS
    }

    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) { return mullo(a, b); }

    /*
     * Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers, and store the high 64
     bits of the intermediate integers in vect_t.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     * Return : [Floor(a0*b0/2^64), Floor(a1*b1/2^64)] int64_t
     */
    static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b) {
        //#pragma warning "The simd mulhi function is emulated, it may impact the performances."
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set(mulhi_64(ca.t[0], cb.t[0]), mulhi_64 (ca.t[1], cb.t[1]));
    }

    /*
     * Multiply the low 32-bit integers from each packed 64-bit element in a and b, and store the signed 64-bit results
     in vect_t.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     * Return : [(a0 smod 2^32)*(b0 smod 2^32), (a1 smod 2^32)*(b1 smod 2^32)]	int64_t
     *	   where (a smod p) is the signed representant of a modulo p, that is -p/2 <= (a smod p) < p/2
     */
    static INLINE CONST vect_t mulx(const vect_t a, const vect_t b) { return _mm_mul_epi32(a, b); }

    /*
     * Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers,
     * keep the low 64 bits of the intermediate and add the low 64-bits of c.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     [c0, c1] int64_t
     * Return : [(a0*b0+c0) smod 2^64, (a1*b1+c1) smod 2^64]	int64_t
     */
    static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b) { return add(c, mul(a, b)); }

    static INLINE vect_t fmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fmadd(c, a, b); }

    /*
     * Multiply the low 32-bit integers from each packed 64-bit element in a and b,
     * keep the signed 64-bit results and add the low 64-bits of c.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     [c0, c1] int64_t
     * Return : [((a0 smod 2^32)*(b0 smod 2^32)+c0) smod 2^64,
     *		 ((a1 smod 2^32)*(b1 smod 2^32)+c1) smod 2^64]	int64_t
     */
    static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b) { return add(c, mulx(a, b)); }

    static INLINE vect_t fmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmaddx(c, a, b); }

    /*
     * Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers,
     * and substract the low 64 bits of the intermediate from elements of c.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     [c0, c1] int64_t
     * Return : [(-a0*b0+c0) smod 2^64, (-a1*b1+c1) smod 2^64]	int64_t
     */
    static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mul(a, b)); }

    static INLINE vect_t fnmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmadd(c, a, b); }

    /*
     * Multiply the low 32-bit integers from each packed 64-bit element in a and b,
     * keep the signed 64-bit results and substract them from elements of c.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     [c0, c1] int64_t
     * Return : [(-(a0 smod 2^32)*(b0 smod 2^32)+c0) smod 2^64,
     *		 (-(a1 smod 2^32)*(b1 smod 2^32)+c1) smod 2^64]	int64_t
     */
    static INLINE CONST vect_t fnmaddx(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mulx(a, b)); }

    static INLINE vect_t fnmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmaddx(c, a, b); }

    /*
     * Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers,
     * and substract elements of c to the low 64-bits of the intermediate.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     [c0, c1] int64_t
     * Return : [(a0*b0-c0) smod 2^64, (a1*b1-c1) smod 2^64]	int64_t
     */
    static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b) { return sub(mul(a, b), c); }

    static INLINE vect_t fmsubin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsub(c, a, b); }

    /*
     * Multiply the low 32-bit integers from each packed 64-bit element in a and b,
     * keep the signed 64-bit results and substract elements of c from them.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     [c0, c1] int64_t
     * Return : [(-(a0 smod 2^32)*(b0 smod 2^32)+c0) smod 2^64,
     *		 (-(a1 smod 2^32)*(b1 smod 2^32)+c1) smod 2^64]	int64_t
     */
    static INLINE CONST vect_t fmsubx(const vect_t c, const vect_t a, const vect_t b) { return sub(mulx(a, b), c); }

    static INLINE vect_t fmsubxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsubx(c, a, b); }

    /*
     * Compare packed 64-bits in a and b for equality, and store the results in vect_t.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     * Return : [(a0==b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1==b1) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
     */
    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) { return _mm_cmpeq_epi64(a, b); }

    /*
     * Compare packed 64-bits in a and b for greater-than, and store the results in vect_t.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     * Return : [(a0>b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1>b1) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
     */
    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {
#ifdef __FFLASFFPACK_HAVE_SSE4_2_INSTRUCTIONS
        return _mm_cmpgt_epi64(a, b);
#else
        //#warning "The simd greater function is emulate, it may impact the performances."
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set((ca.t[0] > cb.t[0]) ? 0xFFFFFFFFFFFFFFFF : 0, (ca.t[1] > cb.t[1]) ? 0xFFFFFFFFFFFFFFFF : 0);
#endif // __FFLASFFPACK_HAVE_SSE4_2_INSTRUCTIONS
    }

    /*
     * Compare packed 64-bits in a and b for lesser-than, and store the results in vect_t.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     * Return : [(a0<b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1<b1) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
     */
    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {
#ifdef __FFLASFFPACK_HAVE_SSE4_2_INSTRUCTIONS
        return _mm_cmpgt_epi64(b, a);
#else
        //#warning "The simd lesser function is emulate, it may impact the performances."
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set((ca.t[0] < cb.t[0]) ? 0xFFFFFFFFFFFFFFFF : 0, (ca.t[1] < cb.t[1]) ? 0xFFFFFFFFFFFFFFFF : 0);
#endif // __FFLASFFPACK_HAVE_SSE4_2_INSTRUCTIONS
    }

    /*
     * Compare packed 64-bits in a and b for greater or equal than, and store the results in vect_t.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     * Return : [(a0>=b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1>=b1) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
     */
    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

    /*
     * Compare packed 64-bits in a and b for lesser or equal than, and store the results in vect_t.
     * Args   : [a0, a1] int64_t
     [b0, b1] int64_t
     * Return : [(a0<=b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1<=b1) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
     */
    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

    /*
     * Horizontally add 64-bits elements of a.
     * Args   : [a0, a1]	int64_t
     * Return : a0+a1	int64_t
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        Converter conv;
        conv.v = a;
        return scalar_t(conv.t[0] + conv.t[1]);
    }

    static INLINE CONST vect_t round(const vect_t a) { return a; }

    // mask the high 32 bits of a 64 bits, that is 00000000FFFFFFFF
    static INLINE CONST vect_t mask_high() { return srl<32>(_mm_set1_epi8(-1)); }

    static INLINE CONST vect_t mulhi_fast(vect_t x, vect_t y);

    static INLINE vect_t mod(vect_t &C, const __m128d &P, const __m128d &INVP, const __m128d &NEGP, const vect_t &POW50REM,
                             const __m128d &MIN, const __m128d &MAX, __m128d &Q, __m128d &T);

protected:
    /* return the sign where vect_t is seen as four int32_t */
    static INLINE CONST vect_t signbits(const vect_t x) {
        vect_t signBits = sub(zero(), srl< 4*sizeof(scalar_t)-1>(x));
        return signBits;
    }
}; // Simd128_impl<true, true, true, 8>

/*
 * Simd128 specialized for uint64_t
 */
template <> struct Simd128_impl<true, true, false, 8> : public Simd128_impl<true, true, true, 8> {

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
     *  Return [x,x] uint64_t
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm_set1_epi64x(x); }

    /*
     *  Set packed 64-bit integers in dst with the supplied values.
     *  Return [x0,x1] uint64_t
     */
    static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1) { return _mm_set_epi64x(x1, x0); }

    /*
     *  Gather 64-bit unsigned integer elements with indexes idx[0], ..., idx[1] from the address p in vect_t.
     *  Return [p[idx[0]], p[idx[1]]] uint64_t
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        return set(p[idx[0]], p[idx[1]]);
    }

    /*
     * Load 128-bits of unsigned integer data from memory into dst.
     * p must be aligned on a 16-byte boundary or a general-protection exception will be generated.
     * Return [p[0],p[1]] uint64_t
     */
    static INLINE PURE vect_t load(const scalar_t *const p) {
        return _mm_load_si128(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Load 128-bits of unsigned integer data from memory into dst.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0],p[1]] uint64_t
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
     * Shift packed 64-bit unsigned integers in a right by s while shifting in sign bits, and store the results in vect_t.
     * Args   : [a0, a1]				uint64_t
     * Return : [Floor(a0/2^s), Floor(a1/2^s)]	uint64_t
     */
    template<int s>
    static INLINE CONST vect_t sra(const vect_t a) { return _mm_srli_epi64(a, s); }

    static INLINE CONST vect_t greater(vect_t a, vect_t b) {
#ifdef __FFLASFFPACK_HAVE_SSE4_2_INSTRUCTIONS
        vect_t x;
        x = set1(static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1));
        a = vxor(x, a);
        b = vxor(x, b);
        return _mm_cmpgt_epi64 (a, b);
#else
        //#pragma warning "The simd greater function is emulated, it may impact the performances."
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set((ca.t[0] > cb.t[0]) ? 0xFFFFFFFFFFFFFFFF : 0, (ca.t[1] > cb.t[1]) ? 0xFFFFFFFFFFFFFFFF : 0);
#endif // __FFLASFFPACK_HAVE_SSE4_2_INSTRUCTIONS
    }

    static INLINE CONST vect_t lesser(vect_t a, vect_t b) {
#ifdef __FFLASFFPACK_HAVE_SSE4_2_INSTRUCTIONS
        vect_t x;
        x = set1(static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1));
        a = vxor(x, a);
        b = vxor(x, b);
        return _mm_cmpgt_epi64 (b, a);
#else
        //#pragma warning "The simd greater function is emulated, it may impact the performances."
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set((ca.t[0] < cb.t[0]) ? 0xFFFFFFFFFFFFFFFF : 0, (ca.t[1] < cb.t[1]) ? 0xFFFFFFFFFFFFFFFF : 0);
#endif // __FFLASFFPACK_HAVE_SSE4_2_INSTRUCTIONS
    }

    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

    /*
     * Multiply the packed 64-bit unsigned integers in a and b, producing intermediate 128-bit integers, and store the low 64
     bits of the intermediate integers in vect_t.
     * Args   : [a0, a1] uint64_t
     [b0, b1] uint64_t
     * Return : [a0*b0 mod 2^64, a1*b1 mod 2^64] uint64_t
     */
    static INLINE CONST vect_t mullo(const vect_t x0, const vect_t x1) {
        // _mm_mullo_epi32 emul
        //#pragma warning "The simd mullo function is emulated, it may impact the performances."
        Converter c0, c1;
        c0.v = x0;
        c1.v = x1;
        return set((scalar_t)(c0.t[0] * c1.t[0]), (scalar_t)(c0.t[1] * c1.t[1]));
    }

    /*
     * Multiply the packed unsigned 64-bit integers in a and b, producing intermediate 128-bit integers,
     * and store the high 64 bits of the intermediate integers in vect_t.
     * Args   : [a0, a1] uint64_t
     [b0, b1] uint64_t
     * Return : [Floor(a0*b0/2^16), Floor(a1*b1/2^16)] uint64_t
     */
    static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b) {
        //#pragma warning "The simd mulhi function is emulate, it may impact the performances."
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set(mulhi_u64(ca.t[0], cb.t[0]), mulhi_u64 (ca.t[1], cb.t[1]));
    }

    /*
     * Multiply the low unsigned 32-bit integers from each packed 64-bit element in a and b, and store the signed 64-bit results
     in vect_t.
     * Args   : [a0, a1] uint64_t
     [b0, b1] uint64_t
     * Return : [(a0 mod 2^32)*(b0 mod 2^32), (a1 mod 2^32)*(b1 mod 2^32)]	uint64_t
     */
    static INLINE CONST vect_t mulx(const vect_t a, const vect_t b) { return _mm_mul_epu32(a, b); }

    static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b) { return add(c, mulx(a, b)); }

    static INLINE vect_t fmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmaddx(c, a, b); }

    static INLINE CONST vect_t fnmaddx(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mulx(a, b)); }

    static INLINE vect_t fnmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmaddx(c, a, b); }

    static INLINE CONST vect_t fmsubx(const vect_t c, const vect_t a, const vect_t b) { return sub(mulx(a, b), c); }

    static INLINE vect_t fmsubxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsubx(c, a, b); }

    /*
     * Horizontally add 64-bits elements of a.
     * Args   : [a0, a1, a2, a3]
     * Return : a0+a1+a2+a3
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        Converter c;
        c.v = a;
        return c.t[0] + c.t[1];
    }
}; //Simd128_impl<true,true,false,8>

#define vect_t Simd128_impl<true,true,true,8>::vect_t

// warning : may be off by 1 multiple, but we save a mul...
INLINE CONST vect_t Simd128_impl<true,true,true,8>::mulhi_fast(vect_t x, vect_t y) {
    // unsigned mulhi starts:
    // x1 = xy_high = mulhiu_fast(x,y)
    const vect_t mask = mask_high();

    vect_t x0 = vand(x, mask), x1 = srl<32>(x);
    vect_t y0 = vand(y, mask), y1 = srl<32>(y);

    x0 = Simd128_impl<true, true, false, 8>::mulx(x0, y1); // x0y1
    y0 = Simd128_impl<true, true, false, 8>::mulx(x1, y0); // x1y0
    y1 = Simd128_impl<true, true, false, 8>::mulx(x1, y1); // x1y1

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

// FIXME why cannot use Simd128<double>::vect_t in the declaration instead of __m128d?
// --**~~~~ Only suitable for use with Modular<int64_t> or ModularBalanced<int64_t>, so p <= max_cardinality < 2**33 ~~~~~**--
INLINE vect_t Simd128_impl<true, true, true, 8>::mod(vect_t &C, const __m128d &P, const __m128d &INVP, const __m128d &NEGP, const vect_t &POW50REM,
                                                     const __m128d &MIN, const __m128d &MAX, __m128d &Q, __m128d &T) {
    vect_t Cq50, Cr50, Ceq;
    __m128d nCmod;

    // nothing so special with 50; could be something else

    Cq50 = sra<50>(C);                      // Cq50[i] < 2**14
    Cr50 = set1(0x3FFFFFFFFFFFFLL);
    Cr50 = vand(C, Cr50);                   // Cr50[i] < 2**50

    Ceq = fmadd(Cr50, Cq50, POW50REM);      // Ceq[i] < 2**47 + 2**50 < 2**51; Ceq[i] ~ C[i] mod p

#if defined(__FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS) and defined(__FFLASFFPACK_HAVE_AVX512VL_INSTRUCTIONS)
    nCmod = _mm_cvtepi64_pd(Ceq);
#else
    // ><
    Converter cC;
    cC.v = Ceq;
    nCmod = _mm_set_pd(static_cast<double>(cC.t[1]),static_cast<double>(cC.t[0]));
#endif

    nCmod = Simd128<double>::mod(nCmod, P, INVP, NEGP, MIN, MAX, Q, T);

#if defined(__FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS) and defined(__FFLASFFPACK_HAVE_AVX512VL_INSTRUCTIONS)
    C = _mm_cvtpd_epi64(nCmod);
#else
    double r[2];
    _mm_storeu_pd(r, nCmod); // could be changed to store if guaranteed to be aligned
    C = set(static_cast<int64_t>(r[0]),static_cast<int64_t>(r[1]));
#endif

    return C;
}

#undef vect_t

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd128_int64_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
