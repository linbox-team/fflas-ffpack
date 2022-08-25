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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd128_int32_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd128_int32_INL

#ifndef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
#error "You need SSE instructions to perform 128 bits operations on int32"
#endif

#ifdef __x86_64__
#include "fflas-ffpack/fflas/fflas_simd/simd128_int64.inl"
#endif

#include "givaro/givtypestring.h"
#include "fflas-ffpack/utils/align-allocator.h"
#include <vector>
#include <type_traits>

/*
 * Simd128 specialized for int32_t
 */
template <> struct Simd128_impl<true, true, true, 4> : public Simd128i_base {

    /*
     * alias to 128 bit simd register
     */
    using vect_t = __m128i;

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = int32_t;

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
     *  Broadcast 32-bit integer a to all elements of dst. This intrinsic may generate vpbroadcastd.
     *  Return [x,x,x,x] int32_t
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm_set1_epi32(x); }

    /*
     *  Set packed 32-bit integers in dst with the supplied values.
     *  Return [x0,x1,x2,x3] int32_t
     */
    static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1, const scalar_t x2, const scalar_t x3) {
        return _mm_set_epi32(x3, x2, x1, x0);
    }

    /*
     *  Gather 32-bit integer elements with indexes idx[0], ..., idx[3] from the address p in vect_t.
     *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]]] int32_t
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        return set(p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]]);
    }

    /*
     * Load 128-bits of integer data from memory into dst.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     * Return [p[0],p[1],p[2],p[3]] int32_t
     */
    static INLINE PURE vect_t load(const scalar_t *const p) {
        return _mm_load_si128(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Load 128-bits of integer data from memory into dst.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]] int32_t
     */
    static INLINE PURE vect_t loadu(const scalar_t *const p) {
        return _mm_loadu_si128(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Store 128-bits of integer data from a into memory.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
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
     * Args   : [a0, a1, a2, a3] int32_t
     * Return : [a0 << s, a1 << s, a2 << s, a3 << s] int32_t
     */
    template<int s>
    static INLINE CONST vect_t sll(const vect_t a) { return _mm_slli_epi32(a, s); }

    /*
     * Shift packed 64-bit integers in a right by s while shifting in zeros, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3] int32_t
     * Return : [a0 >> s, a1 >> s, a2 >> s, a3 >> s] int32_t
     */
    template<int s>
    static INLINE CONST vect_t srl(const vect_t a) { return _mm_srli_epi32(a, s); }

    /*
     * Shift packed 32-bit integers in a right by s while shifting in sign bits, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3] int32_t
     * Return : [a0 >> s, a1 >> s, a2 >> s, a3 >> s] int32_t
     */
    template<int s>
    static INLINE CONST vect_t sra(const vect_t a) { return _mm_srai_epi32(a, s); }

    /*
     * Shuffle 32-bit integers in a using the control in imm8, and store the results in dst.
     * Args   : [a0, a1, a2, a3] int32_t
     * Return : [a[s[0..1]], ..., a[s[6..7]] int32_t
     */
    template<uint8_t s>
    static INLINE CONST vect_t shuffle(const vect_t a) {
        return _mm_shuffle_epi32(a, s);
    }

    /*
     * Unpack and interleave 32-bit integers from the low half of a and b.
     * Args: a = [ a0, a1, a2, a3 ]
     *       b = [ b0, b1, b2, b3 ]
     * Return:   [ a0, b0, a1, b1 ]
     */
    static INLINE CONST vect_t
    unpacklo_intrinsic (const vect_t a, const vect_t b) {
        return _mm_unpacklo_epi32(a, b);
    }

    /*
     * Unpack and interleave 32-bit integers from the high half of a and b.
     * Args: a = [ a0, a1, a2, a3 ]
     *       b = [ b0, b1, b2, b3 ]
     * Return:   [ a2, b2, a3, b3 ]
     */
    static INLINE CONST vect_t
    unpackhi_intrinsic (const vect_t a, const vect_t b) {
        return _mm_unpackhi_epi32(a, b);
    }

    /*
     * Unpack and interleave 32-bit integers from the low half of a and b.
     * Args: a = [ a0, a1, a2, a3 ]
     *       b = [ b0, b1, b2, b3 ]
     * Return:   [ a0, b0, a1, b1 ]
     */
    static INLINE CONST vect_t unpacklo(const vect_t a, const vect_t b) {
        return unpacklo_intrinsic(a, b);
    }

    /*
     * Unpack and interleave 32-bit integers from the high half of a and b.
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
     * Pack 32-bit integers from the even positions of a and b.
     * Args: a = [ a0, a1, a2, a3 ]
     *       b = [ b0, b1, b2, b3 ]
     * Return:   [ a0, a2, b0, b2 ]
     */
    static INLINE CONST vect_t pack_even (const vect_t a, const vect_t b) {
        /* 0xd8 = 3120 base_4 */
        vect_t t1 = _mm_shuffle_epi32 (a, 0xd8);
        vect_t t2 = _mm_shuffle_epi32 (b, 0xd8);
        return _mm_unpacklo_epi64 (t1, t2);
    }

    /*
     * Pack 32-bit integers from the odd positions of a and b.
     * Args: a = [ a0, a1, a2, a3 ]
     *       b = [ b0, b1, b2, b3 ]
     * Return:   [ a1, a3, b1, b3 ]
     */
    static INLINE CONST vect_t pack_odd (const vect_t a, const vect_t b) {
        /* 0xd8 = 3120 base_4 */
        vect_t t1 = _mm_shuffle_epi32 (a, 0xd8);
        vect_t t2 = _mm_shuffle_epi32 (b, 0xd8);
        return _mm_unpackhi_epi64 (t1, t2);
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
        vect_t t1 = _mm_shuffle_epi32 (a, 0xd8);
        vect_t t2 = _mm_shuffle_epi32 (b, 0xd8);
        even = _mm_unpacklo_epi64 (t1, t2);
        odd = _mm_unpackhi_epi64 (t1, t2);
    }

    /*
     * Transpose the 4x4 matrix formed by the 4 rows of 32-bit integers in r0,
     * r1, r2 and r3, and store the transposed matrix in these vectors.
     * Args: r0 = [ r00, r01, r02, r03 ]
     *       r1 = [ r10, r11, r12, r13 ]
     *       r2 = [ r20, r21, r22, r23 ]
     *       r3 = [ r30, r31, r32, r33 ]
     * Return: r0 = [ r00, r10, r20, r30 ]
     *         r1 = [ r01, r11, r21, r31 ]
     *         r2 = [ r02, r12, r22, r32 ]
     *         r3 = [ r03, r13, r23, r33 ]
     */
    static INLINE void
    transpose (vect_t& r0, vect_t& r1, vect_t& r2, vect_t& r3) {
        vect_t t0, t1, t2, t3;
        t0 = unpacklo_intrinsic (r0, r1);
        t2 = unpacklo_intrinsic (r2, r3);
        t1 = unpackhi_intrinsic (r0, r1);
        t3 = unpackhi_intrinsic (r2, r3);
        r0 = _mm_castps_si128 (_mm_movelh_ps (_mm_castsi128_ps (t0),
                                              _mm_castsi128_ps (t2)));
        r1 = _mm_castps_si128 (_mm_movehl_ps (_mm_castsi128_ps (t2),
                                              _mm_castsi128_ps (t0)));
        r2 = _mm_castps_si128 (_mm_movelh_ps (_mm_castsi128_ps (t1),
                                              _mm_castsi128_ps (t3)));
        r3 = _mm_castps_si128 (_mm_movehl_ps (_mm_castsi128_ps (t3),
                                              _mm_castsi128_ps (t1)));
    }

    /*
     * Blend 32-bit integers from a and b using control mask s.
     * Args: a = [ a0, ..., a3 ]
     *       b = [ b0, ..., b3 ]
     *       s = a 4-bit immediate integer
     * Return: [ s[0] ? a0 : b0, ..., s[3] ? a3 : b3 ]
     */
    template<uint8_t s>
    static INLINE CONST vect_t blend(const vect_t a, const vect_t b) {
#ifdef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
        /* If we have AVX2, we can use _mm_blend_epi32 */
        return _mm_blend_epi32(a, b, s);
#else
        /* If we only have SSE4.1, we need to use _mm_blend_epi16.
         * We need to transform s = [d3 d2 d1 d0]_base2
         * into s1 = [d3 d3 d2 d2 d1 d1 d0 d0]_base2
         */
        constexpr uint8_t s1 = (s & 0x1) * 3 + (((s & 0x2) << 1)*3)
                                + (((s & 0x4) << 2)*3) + (((s & 0x8) << 3)*3);
        return _mm_blend_epi16(a, b, s1);
#endif
    }

    /*
     * Add packed 32-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3] int32_t
     [b0, b1, b2, b3] int32_t
     * Return : [a0+b0, a1+b1, a2+b2, a3+b3]   int32_t
     */
    static INLINE CONST vect_t add(const vect_t a, const vect_t b) { return _mm_add_epi32(a, b); }

    static INLINE vect_t addin(vect_t &a, const vect_t b) { return a = add(a, b); }

    /*
     * Subtract packed 32-bit integers in b from packed 32-bit integers in a, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3] int32_t
     [b0, b1, b2, b3] int32_t
     * Return : [a0-b0, a1-b1, a2-b2, a3-b3]  int32_t
     */
    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) { return _mm_sub_epi32(a, b); }

    static INLINE vect_t subin(vect_t &a, const vect_t b) { return a = sub(a, b); }

    /*
     * Multiply the packed 32-bit integers in a and b, producing intermediate 64-bit integers, and store the low 32 bits
     of the intermediate integers in vect_t.
     * Args   : [a0, a1, a2, a3] int32_t
     [b0, b1, b2, b3] int32_t
     * Return : [a0*b0 smod 2^32, ..., a3*b3 smod 2^32]	int32_t
     *	   where (a smod p) is the signed representant of a modulo p, that is -p/2 <= (a smod p) < p/2
     */
    static INLINE CONST vect_t mullo(const vect_t a, const vect_t b) { return _mm_mullo_epi32(a, b); }

    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) { return mullo(a, b); }

    /*
     * Multiply the packed 32-bit integers in a and b, producing intermediate 64-bit integers, and store the high 32
     bits of the intermediate integers in vect_t.
     * Args   : [a0, a1, a2, a3] int32_t
     [b0, b1, b2, b3] int32_t
     * Return : [Floor(a0*b0/2^32), ..., Floor(a3*b3/2^32)] int32_t
     */

    static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b) {
        // _mm_mulhi_epi32 emul
        //#pragma warning "The simd mulhi function is emulated, it may impact the performances."
#ifndef __x86_64__
        vect_t a1, a2, b1, b2;
        a1 = set(_mm_extract_epi32(a, 0), 0, _mm_extract_epi32(a, 2), 0);
        a2 = set(_mm_extract_epi32(a, 1), 0, _mm_extract_epi32(a, 3), 0);
        b1 = set(_mm_extract_epi32(b, 0), 0, _mm_extract_epi32(b, 2), 0);
        b2 = set(_mm_extract_epi32(b, 1), 0, _mm_extract_epi32(b, 3), 0);
        a1 = _mm_mul_epi32(a1, b1);
        a2 = _mm_mul_epi32(a2, b2);
        return set(_mm_extract_epi32(a1, 1), _mm_extract_epi32(a2, 1), _mm_extract_epi32(a1, 3),
                   _mm_extract_epi32(a2, 3));
#else
        typedef Simd128_impl<true, true, true, 8> Simd128_64;
        vect_t C,A1,B1;
        C  = Simd128_64::mulx(a,b);
        A1 = Simd128_64::srl<32>(a);
        B1 = Simd128_64::srl<32>(b);
        A1 = Simd128_64::mulx(A1,B1);
        C  = Simd128_64::srl<32>(C);
        A1 = Simd128_64::srl<32>(A1);
        A1 = Simd128_64::sll<32>(A1);
        return Simd128_64::vor(C,A1);
#endif
    }

    /*
     * Multiply the low 16-bit integers from each packed 32-bit element in a and b, and store the signed 32-bit results
     in vect_t.
     * Args   : [a0, a1, a2, a3] int32_t
     [b0, b1, b2, b3] int32_t
     * Return : [(a0 smod 2^16)*(b0 smod 2^16), (a1 smod 2^16)*(b1 smod 2^16),
     *	    (a2 smod 2^16)*(b2 smod 2^16), (a3 smod 2^16)*(b3 smod 2^16)]	int32_t
     *	   where (a smod p) is the signed representant of a modulo p, that is -p/2 <= (a smod p) < p/2
     */
    static INLINE CONST vect_t mulx(const vect_t a, const vect_t b) {
        //#pragma warning "The simd mulx function is emulated, it may impact the performances."
        vect_t a1, b1, mask1, mask2;
        mask1 = set1(0x0000FFFF);
        mask2 = set1(0x00008000);
        a1 = add(a,mask2);
        a1 = vand(a1,mask1);
        a1 = sub(a1,mask2);
        b1 = add(b,mask2);
        b1 = vand(b1,mask1);
        b1 = sub(b1,mask2);
        return mul(a1,b1);
    }

    /*
     * Multiply the packed 32-bit integers in a and b, producing intermediate 64-bit integers,
     * keep the low 32 bits of the intermediate and add the low 32-bits of c.
     * Args   :	[a0, a1, a2, a3]		int32_t
     [b0, b1, b2, b3]		int32_t
     [c0, c1, c2, c3]		int32_t
     * Return :	[(a0*b0+c0) smod 2^32, ..., (a3*b3+c3) smod 2^32]	int32_t
     */
    static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b) { return add(c, mul(a, b)); }

    static INLINE vect_t fmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fmadd(c, a, b); }

    /*
     * Multiply the low 16-bit integers from each packed 32-bit element in a and b,
     * keep the signed 32-bit results and add the low 32-bits of c.
     * Args   :	[a0, a1, a2, a3]		int32_t
     [b0, b1, b2, b3]		int32_t
     [c0, c1, c2, c3]		int32_t
     * Return :	[((a0 smod 2^16)*(b0 smod 2^16)+c0) smod 2^32, ...,
     *		 ((a3 smod 2^16)*(b3 smod 2^16)+c3) smod 2^32]	int32_t
     */
    static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b) { return add(c, mulx(a, b)); }

    static INLINE vect_t fmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmaddx(c, a, b); }

    /*
     * Multiply the packed 32-bit integers in a and b, producing intermediate 64-bit integers,
     * and substract the low 32 bits of the intermediate from elements of c.
     * Args   :	[a0, a1, a2, a3]		int32_t
     [b0, b1, b2, b3]		int32_t
     [c0, c1, c2, c3]		int32_t
     * Return :	[(-a0*b0+c0) smod 2^32, ..., (-a3*b3+c3) smod 2^32]	int32_t
     */
    static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mul(a, b)); }

    static INLINE vect_t fnmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmadd(c, a, b); }

    /*
     * Multiply the low 16-bit integers from each packed 32-bit element in a and b,
     * keep the signed 32-bit results and add the low 32-bits of c and substract them from elements of c.
     * Args   :	[a0, a1, a2, a3]		int32_t
     [b0, b1, b2, b3]		int32_t
     [c0, c1, c2, c3]		int32_t
     * Return :	[(-(a0 smod 2^16)*(b0 smod 2^16)+c0) smod 2^32, ...,
     *		 (-(a3 smod 2^16)*(b3 smod 2^16)+c3) smod 2^32]	int32_t
     */
    static INLINE CONST vect_t fnmaddx(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mulx(a, b)); }

    static INLINE vect_t fnmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmaddx(c, a, b); }

    /*
     * Multiply the packed 32-bit integers in a and b, producing intermediate 64-bit integers,
     * and substract elements of c to the low 32-bits of the intermediate.
     * Args   :	[a0, a1, a2, a3]		int32_t
     [b0, b1, b2, b3]		int32_t
     [c0, c1, c2, c3]		int32_t
     * Return : [(a0*b0-c0) smod 2^32, ..., (a3*b3-c3) smod 2^32]	int32_t
     */
    static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b) { return sub(mul(a, b), c); }

    static INLINE vect_t fmsubin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsub(c, a, b); }

    /*
     * Multiply the low 16-bit integers from each packed 32-bit element in a and b,
     * keep the signed 32-bit results and substract elements of c from them.
     * Args   :	[a0, a1, a2, a3]		int32_t
     [b0, b1, b2, b3]		int32_t
     [c0, c1, c2, c3]		int32_t
     * Return :	[((a0 smod 2^16)*(b0 smod 2^16)-c0) smod 2^32, ...,
     *		 ((a3 smod 2^16)*(b3 smod 2^16)-c3) smod 2^32]	int32_t
     */
    static INLINE CONST vect_t fmsubx(const vect_t c, const vect_t a, const vect_t b) { return sub(mulx(a, b), c); }

    static INLINE vect_t fmsubxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsubx(c, a, b); }

    /*
     * Compare packed 32-bits in a and b for equality, and store the results in vect_t.
     * Args   :	[a0, a1, a2, a3]		int32_t
     [b0, b1, b2, b3]		int32_t
     * Return : [(a0==b0) ? 0xFFFFFFFF : 0, (a1==b1) ? 0xFFFFFFFF : 0,
     (a2==b2) ? 0xFFFFFFFF : 0, (a3==b3) ? 0xFFFFFFFF : 0]			int32_t
     */
    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) { return _mm_cmpeq_epi32(a, b); }

    /*
     * Compare packed 32-bits in a and b for greater-than, and store the results in vect_t.
     * Args   :	[a0, a1, a2, a3]		int32_t
     [b0, b1, b2, b3]		int32_t
     * Return : [(a0>b0) ? 0xFFFFFFFF : 0, (a1>b1) ? 0xFFFFFFFF : 0,
     (a2>b2) ? 0xFFFFFFFF : 0, (a3>b3) ? 0xFFFFFFFF : 0]			int32_t
     */
    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) { return _mm_cmpgt_epi32(a, b); }

    /*
     * Compare packed 32-bits in a and b for lesser-than, and store the results in vect_t.
     * Args   :	[a0, a1, a2, a3]		int32_t
     [b0, b1, b2, b3]		int32_t
     * Return : [(a0<b0) ? 0xFFFFFFFF : 0, (a1<b1) ? 0xFFFFFFFF : 0,
     (a2<b2) ? 0xFFFFFFFF : 0, (a3<b3) ? 0xFFFFFFFF : 0]			int32_t
     */
    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) { return _mm_cmplt_epi32(a, b); }

    /*
     * Compare packed 32-bits in a and b for greater or equal than, and store the results in vect_t.
     * Args   :	[a0, a1, a2, a3]		int32_t
     [b0, b1, b2, b3]		int32_t
     * Return : [(a0>=b0) ? 0xFFFFFFFF : 0, (a1>=b1) ? 0xFFFFFFFF : 0,
     (a2>=b2) ? 0xFFFFFFFF : 0, (a3>=b3) ? 0xFFFFFFFF : 0]			int32_t
     */
    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

    /*
     * Compare packed 32-bits in a and b for lesser or equal than, and store the results in vect_t.
     * Args   :	[a0, a1, a2, a3]		int32_t
     [b0, b1, b2, b3]		int32_t
     * Return : [(a0<=b0) ? 0xFFFFFFFF : 0, (a1<=b1) ? 0xFFFFFFFF : 0,
     (a2<=b2) ? 0xFFFFFFFF : 0, (a3<=b3) ? 0xFFFFFFFF : 0]			int32_t
     */
    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

    /*
     * Horizontally add 32-bits elements of a.
     * Args   : [a0, a1, a2, a3]
     * Return : a0+a1+a2+a3
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        Converter conv;
        conv.v = a;
        return scalar_t(conv.t[0] + conv.t[1] + conv.t[2] + conv.t[3]);
    }

    static INLINE CONST vect_t round(const vect_t a) { return a; }

    static INLINE vect_t mod(vect_t &C, const vect_t &P, const vect_t &INVP, const vect_t &NEGP, const vect_t &MIN,
                             const vect_t &MAX, vect_t &Q, vect_t &T) {
#ifdef __INTEL_COMPILER
        C = _mm_rem_epi32(C, P);
#else
        FFLASFFPACK_abort("pas implementé");
        // C = fnmadd(C,_mm_castps_si128(_mm_floor_ps(_mm_mul_ps(INVP,_mm_castsi128_ps(C)))),P);
#endif
        NORML_MOD(C, P, NEGP, MIN, MAX, Q, T);
        return C;
    }

};

/*
 * Simd128 specialized for uint32_t
 */
template <> struct Simd128_impl<true, true, false, 4> : public Simd128_impl<true, true, true, 4> {

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = uint32_t;

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
     *  Broadcast 32-bit unsigned integer a to all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x,x,x,x] uint32_t
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm_set1_epi32(x); }

    /*
     *  Set packed 32-bit unsigned integers in dst with the supplied values.
     *  Return [x0,x1,x2,x3] uint32_t
     */
    static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1, const scalar_t x2, const scalar_t x3) {
        return _mm_set_epi32(x3, x2, x1, x0);
    }

    /*
     *  Gather 32-bit unsigned integer elements with indexes idx[0], ..., idx[3] from the address p in vect_t.
     *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]]] uint32_t
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        return set(p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]]);
    }

    /*
     * Load 128-bits of unsigned integer data from memory into dst.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     * Return [p[0],p[1],p[2],p[3]] uint32_t
     */
    static INLINE PURE vect_t load(const scalar_t *const p) {
        return _mm_load_si128(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Load 128-bits of unsigned integer data from memory into dst.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]] uint32_t
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
     * Shift packed 32-bit unsigned integers in a right by s while shifting in sign bits, and store the results in vect_t.
     * Args   : [a0, ..., a3]			int32_t
     * Return : [Floor(a0/2^s), ..., Floor(a3/2^s)]	int32_t
     */
    template<int s>
    static INLINE CONST vect_t sra(const vect_t a) { return _mm_srli_epi32(a, s); }

    static INLINE CONST vect_t greater(vect_t a, vect_t b) {
        vect_t x;
        x = set1(static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1));
        a = vxor(x, a);
        b = vxor(x, b);
        return _mm_cmpgt_epi32(a, b);
    }

    static INLINE CONST vect_t lesser(vect_t a, vect_t b) {
        vect_t x;
        x = set1(static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1));
        a = vxor(x, a);
        b = vxor(x, b);
        return _mm_cmplt_epi32(a, b);
    }

    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

    /*
     * Multiply the packed unsigned 32-bit integers in a and b, producing intermediate 64-bit integers,
     * and store the high 32	bits of the intermediate integers in vect_t.
     * Args   : [a0, a1, a2, a3]		 uint32_t
     *		   [b0, b1, b2, b3]		 uint32_t
     * Return : [Floor(a0*b0/2^32), ..., Floor(a3*b3/2^32)] uint32_t
     */
    static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b) {
        // _mm_mulhi_epi32 emul
        //#pragma warning "The simd mulhi function is emulated, it may impact the performances."
#ifndef __x86_64__
        vect_t a1, a2, b1, b2;
        a1 = set(_mm_extract_epi32(a, 0), 0, _mm_extract_epi32(a, 2), 0);
        a2 = set(_mm_extract_epi32(a, 1), 0, _mm_extract_epi32(a, 3), 0);
        b1 = set(_mm_extract_epi32(b, 0), 0, _mm_extract_epi32(b, 2), 0);
        b2 = set(_mm_extract_epi32(b, 1), 0, _mm_extract_epi32(b, 3), 0);
        a1 = _mm_mul_epu32(a1, b1);
        a2 = _mm_mul_epu32(a2, b2);
        return set(_mm_extract_epi32(a1, 1), _mm_extract_epi32(a2, 1), _mm_extract_epi32(a1, 3),
                   _mm_extract_epi32(a2, 3));
#else
        typedef Simd128_impl<true, true, false, 8> Simd128_64;
        vect_t C,A1,B1;
        C  = Simd128_64::mulx(a,b);
        A1 = Simd128_64::srl<32>(a);
        B1 = Simd128_64::srl<32>(b);
        A1 = Simd128_64::mulx(A1,B1);
        C  = Simd128_64::srl<32>(C);
        A1 = Simd128_64::srl<32>(A1);
        A1 = Simd128_64::sll<32>(A1);
        return Simd128_64::vor(C,A1);
#endif
    }

    /*
     * Multiply the low unsigned 16-bit integers from each packed 32-bit element in a and b,
     * and store the signed 32-bit results in vect_t.
     * Args   : [a0, a1, a2, a3]		 uint32_t
     *		   [b0, b1, b2, b3]		 uint32_t
     * Return : [(a0 mod 2^16)*(b0 mod 2^16), (a1 mod 2^16)*(b1 mod 2^16),
     *	    (a2 mod 2^16)*(b2 mod 2^16), (a3 mod 2^16)*(b3 mod 2^16)]	uint32_t
     */
    static INLINE CONST vect_t mulx(const vect_t a, const vect_t b) {
        //#pragma warning "The simd mulx function is emulated, it may impact the performances."
        vect_t a1, b1, mask1;
        mask1 = set1(0x0000FFFF);
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
     * Horizontally add 32-bits elements of a.
     * Args   : [a0, a1, a2, a3]
     * Return : a0+a1+a2+a3
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        Converter conv;
        conv.v = a;
        return conv.t[0] + conv.t[1] + conv.t[2] + conv.t[3];
    }
}; //Simd128_impl<true,true,false,4>

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd128_int32_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
