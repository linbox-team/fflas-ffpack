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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd256_int32_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd256_int32_INL

#ifndef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
#error "You need AVX2 instructions to perform 256bits operations on int32_t"
#endif

#ifdef __x86_64__
#include "fflas-ffpack/fflas/fflas_simd/simd256_int64.inl"
#endif

#include "fflas-ffpack/utils/align-allocator.h"
#include <vector>
#include <type_traits>

/*
 * Simd256 specialized for int32_t
 */
template <> struct Simd256_impl<true, true, true, 4> : public Simd256i_base {

    /*
     * alias to 256 bit simd register
     */
    using vect_t = __m256i;

    /*
     * alias to 256 bit simd register
     */
    using half_t = __m128i;

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = int32_t;

    /*
     * Simd128 for scalar_t, to deal half_t
     */
    using simdHalf = Simd128<scalar_t>;

    /*
     *  number of scalar_t in a simd register
     */
    static const constexpr size_t vect_size = 8;

    /*
     *  alignement required by scalar_t pointer to be loaded in a vect_t
     */
    static const constexpr size_t alignment = 32;
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
     *  Broadcast 32-bit integer a to all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x,x,x,x,x,x,x,x] int32_t
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm256_set1_epi32(x); }

    /*
     *  Set packed 32-bit integers in dst with the supplied values.
     *  Return [x0,x1,x2,x3,x4,x5,x6,x7] int32_t
     */
    static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1, const scalar_t x2, const scalar_t x3,
                                   const scalar_t x4, const scalar_t x5, const scalar_t x6, const scalar_t x7) {
        return _mm256_set_epi32(x7, x6, x5, x4, x3, x2, x1, x0);
    }

    /*
     *  Gather 32-bit integer elements with indexes idx[0], ..., idx[7] from the address p in vect_t.
     *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]], p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]]] int32_t
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        return set(p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]], p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]]);
    }

    /*
     * Load 256-bits of integer data from memory into dst.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     * Return [p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]] int32_t
     */
    static INLINE PURE vect_t load(const scalar_t *const p) {
        return _mm256_load_si256(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Load 256-bits of integer data from memory into dst.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]] int32_t
     */
    static INLINE PURE vect_t loadu(const scalar_t *const p) {
        return _mm256_loadu_si256(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Store 256-bits of integer data from a into memory.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void store(scalar_t *p, vect_t v) {
        _mm256_store_si256(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Store 256-bits of integer data from a into memory.
     * p does not need to be aligned on any particular boundary.
     */
    static INLINE void storeu(scalar_t *p, vect_t v) {
        _mm256_storeu_si256(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Store 256-bits of integer data from a into memory using a non-temporal memory hint.
     * p must be aligned on a 32-byte boundary or a general-protection exception may be generated.
     */
    static INLINE void stream(scalar_t *p, const vect_t v) {
        _mm256_stream_si256(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Shift packed 32-bit integers in a left by s while shifting in zeros, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int32_t
     * Return : [a0 << s, a1 << s, a2 << s, a3 << s, a4 << s, a5 << s, a6 << s, a7 << s] int32_t
     */
    template<int s>
    static INLINE CONST vect_t sll(const vect_t a) { return _mm256_slli_epi32(a, s); }

    /*
     * Shift packed 32-bit integers in a right by s while shifting in zeros, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int32_t
     * Return : [a0 >> s, a1 >> s, a2 >> s, a3 >> s, a4 >> s, a5 >> s, a6 >> s, a7 >> s] int32_t
     */
    template<int s>
    static INLINE CONST vect_t srl(const vect_t a) { return _mm256_srli_epi32(a, s); }

    /*
     * Shift packed 32-bit integers in a right by s while shifting in sign bits, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int32_t
     * Return : [a0 >> s, a1 >> s, a2 >> s, a3 >> s, a4 >> s, a5 >> s, a6 >> s, a7 >> s] int32_t
     */
    template<int s>
    static INLINE CONST vect_t sra(const vect_t a) { return _mm256_srai_epi32(a, s); }

    /*
     * Shuffle 32-bit integers in a within 128-bit lanes using the control in imm8, and store the results in dst.
     * Args   : [a0, ..., a7] int32_t
     * Return : [a[s[0..1]], ..., a[s[6..7]],a[4+s[0..1]], ..., a[4+s[6..7]]] int32_t
     */
    template<uint8_t s>
    static INLINE CONST vect_t shuffle_twice(const vect_t a) {
        return _mm256_shuffle_epi32(a, s);
    }

    /*
     * Shuffle 32-bit integers in a using the control in imm8, and store the results in dst.
     * Args   : [a0, ..., a7] int32_t
     * Return : [a[s[0..3]], ..., a[s[28..31]]] int32_t
     */
    template<uint32_t s>
    static INLINE CONST vect_t shuffle(const vect_t a) {
        //#pragma warning "The simd shuffle function is emulated, it may impact the performances."
        Converter conv;
        conv.v = a;
        return set (conv.t[( s      & 0x0000000F)], conv.t[( s      & 0x000000F0)],
                    conv.t[((s>> 8) & 0x0000000F)], conv.t[((s>> 8) & 0x000000F0)],
                    conv.t[((s>>16) & 0x0000000F)], conv.t[((s>>16) & 0x000000F0)],
                    conv.t[((s>>24) & 0x0000000F)], conv.t[((s>>24) & 0x000000F0)]);
    }

    /*
     * Unpack and interleave 32-bit integers from the low half of a and b within 128-bit lanes, and store the results in dst.
     * Args   :	[a0, ..., a7] int32_t
     [b0, ..., b7] int32_t
     * Return :	[a0, b0, a1, b1, a4, b4, a5, b5] int32_t
     */
    static INLINE CONST vect_t unpacklo_twice(const vect_t a, const vect_t b) { return _mm256_unpacklo_epi32(a, b); }

    /*
     * Unpack and interleave 32-bit integers from the high half of a and b within 128-bit lanes, and store the results in dst.
     * Args   :	[a0, ..., a7] int32_t
     [b0, ..., b7] int32_t
     * Return :	[a2, b2, a3, b3, a6, b6, a7, b7] int32_t
     */
    static INLINE CONST vect_t unpackhi_twice(const vect_t a, const vect_t b) { return _mm256_unpackhi_epi32(a, b); }

    /*
     * Unpack and interleave 32-bit integers from the low half of a and b, and store the results in dst.
     * Args   :	[a0, ..., a7] int32_t
     [b0, ..., b7] int32_t
     * Return :	[a0, b0, ..., a3, b3] int32_t
     */
    static INLINE CONST vect_t unpacklo(const vect_t a, const vect_t b) {
        //using Simd256_64 = Simd256<uint64_t>;
        //Simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4 so a -> [a0,a2,a1,a3] uint64
        //Simd256_64::template shuffle<0xD8>(b); // 0xD8 = 3120 base_4
        vect_t a1 = _mm256_permute4x64_epi64(a, 0xD8);
        vect_t b1 = _mm256_permute4x64_epi64(a, 0xD8);
        return unpacklo_twice(a1, b1);
    }

    /*
     * Unpack and interleave 32-bit integers from the high half of a and b, and store the results in dst.
     * Args   :	[a0, ..., a7] int32_t
     [b0, ..., b7] int32_t
     * Return :	[a4, b4, ..., a7, b7] int32_t
     */
    static INLINE CONST vect_t unpackhi(const vect_t a, const vect_t b) {
        // using Simd256_64 = Simd256<uint64_t>;
        // vect_t a1 = Simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
        // vect_t b1 = Simd256_64::template shuffle<0xD8>(b); // 0xD8 = 3120 base_4
        vect_t a1 = _mm256_permute4x64_epi64(a, 0xD8);
        vect_t b1 = _mm256_permute4x64_epi64(a, 0xD8);
        return unpackhi_twice(a1, b1);
    }

    /* unpacklohi:
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return: r1 = [ a0, b0, a1, b1, a2, b2, a3, b3 ]
     *         r2 = [ a4, b4, a5, b5, a6, b6, a7, b7 ]
     */
    static INLINE void
    unpacklohi (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
        vect_t t1, t2;
        /* 0xd8 = 3120 base_4 */
        t1 = _mm256_permute4x64_epi64 (a, 0xd8);
        t2 = _mm256_permute4x64_epi64 (b, 0xd8);
        r1 = _mm256_unpacklo_epi32 (t1, t2);
        r2 = _mm256_unpackhi_epi32 (t1, t2);
    }

    /* pack:
     * Args: a = [ a0, a1, a2, a3, a4, a5, a6, a7 ]
     *       b = [ b0, b1, b2, b3, b4, b5, b6, b7 ]
     * Return: r1 = [ a0, a2, a4, a6, b0, b2, b4, b6 ]
     *         r2 = [ a1, a3, a5, a7, b1, b3, b5, b7 ]
     */
    static INLINE void
    pack (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
        vect_t t1, t2;
        /* 0xd8 = 3120 base_4 */
        r1 = _mm256_shuffle_epi32 (a, 0xd8);
        r2 = _mm256_shuffle_epi32 (b, 0xd8);
        t1 = _mm256_unpacklo_epi64 (r1, r2);
        t2 = _mm256_unpackhi_epi64 (r1, r2);
        /* 0xd8 = 3120 base_4 */
        r1 = _mm256_permute4x64_epi64 (t1, 0xd8);
        r2 = _mm256_permute4x64_epi64 (t2, 0xd8);
    }

    /*
     * Blend packed 32-bit integers from a and b using control mask imm8, and store the results in dst.
     * Args   :	[a0, ..., a7] int32_t
     [b0, ..., b7] int32_t
     * Return :	[s[0]?a0:b0,   , s[7]?a7:b7] int32_t
     */
    template<uint8_t s>
    static INLINE CONST vect_t blend(const vect_t a, const vect_t b) {
        return _mm256_blend_epi32(a, b, s);
    }

    /*
     * Add packed 32-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] 						int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7] 						int32_t
     * Return : [a0+b0, a1+b1, a2+b2, a3+b3, a4+b4, a5+b5, a6+b6, a7+b7]   int32_t
     */
    static INLINE CONST vect_t add(const vect_t a, const vect_t b) { return _mm256_add_epi32(a, b); }

    static INLINE vect_t addin(vect_t &a, const vect_t b) { return a = add(a, b); }

    /*
     * Subtract packed 32-bits integers in b from packed 32-bits integers in a, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] 						int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7] 						int32_t
     * Return : [a0-b0, a1-b1, a2-b2, a3-b3, a4-b4, a5-b5, a6-b6, a7-b7]  int32_t
     */
    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) { return _mm256_sub_epi32(a, b); }

    static INLINE vect_t subin(vect_t &a, const vect_t b) { return a = sub(a, b); }

    /*
     * Multiply the packed 32-bits integers in a and b, producing intermediate 64-bit integers, and store the low 32
     bits of the intermediate integers in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]						int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7]	 					int32_t
     * Return : [a0*b0 smod 2^32, ..., a7*b7 smod 2^32]	int32_t
     *	   where (a smod p) is the signed representant of a modulo p, that is -p/2 <= (a smod p) < p/2
     */
    static INLINE CONST vect_t mullo(const vect_t a, const vect_t b) { return _mm256_mullo_epi32(a, b); }

    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) { return mullo(a, b); }

    /*
     * Multiply the packed 32-bits integers in a and b, producing intermediate 64-bit integers, and store the high 32
     bits of the intermediate integers in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7] int32_t
     * Return : [Floor(a0*b0/2^32), ..., Floor(a7*b7/2^32)] int32_t
     */
    static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b) {
        //#pragma warning "The simd mulhi function is emulated, it may impact the performances."
#if 0
        typedef Simd256_impl<true, true, true, 8> Simd256_64;
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        vect_t a1, a2, b1, b2, c1, c2;
        a1 = set(ca.t[0], 0, ca.t[1], 0, ca.t[2], 0, ca.t[3], 0);
        a2 = set(ca.t[4], 0, ca.t[5], 0, ca.t[6], 0, ca.t[7], 0);
        b1 = set(cb.t[0], 0, cb.t[1], 0, cb.t[2], 0, cb.t[3], 0);
        b2 = set(cb.t[4], 0, cb.t[5], 0, cb.t[6], 0, cb.t[7], 0);
        c1 = _mm256_mul_epi32(a1, b1); //Simd256_64::mulx(a1, b1);
        c2 = _mm256_mul_epi32(a1, b2); //Simd256_64::mulx(a2, b2);
        ca.v = c1;
        cb.v = c2;
        return set(ca.t[1], ca.t[3], ca.t[5], ca.t[7], cb.t[1], cb.t[3], cb.t[5], cb.t[7]);
#else
        //typedef Simd256_impl<true, true, true, 8> Simd256_64;
        vect_t C,A1,B1;
        //C  = Simd256_64::mulx(a,b);
        C  = _mm256_mul_epi32(a,b);
        //A1 = Simd256_64::srl(a,32);
        A1 = _mm256_srli_epi64(a, 32);
        //B1 = Simd256_64::srl(b,32);
        B1 = _mm256_srli_epi64(b, 32);
        //A1 = Simd256_64::mulx(A1,B1);
        A1 =  _mm256_mul_epi32(A1,B1);
        //C  = Simd256_64::srl(C,32);
        C  = _mm256_srli_epi64(C, 32);
        //A1 = Simd256_64::srl(A1,32);
        A1 = _mm256_srli_epi64(A1, 32);
        //A1 = Simd256_64::sll(A1,32);
        A1 = _mm256_slli_epi64(A1, 32);
        //return Simd256_64::vor(C,A1);
        return _mm256_or_si256(A1, C);
#endif
    }

    /*
     * Multiply the low 16-bit integers from each packed 32-bit element in a and b, and store the signed 32-bit results
     in dst.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7]	int32_t
     * Return : [(a0 smod 2^16)*(b0 smod 2^16), ..., (a7 smod 2^16)*(b7 smod 2^16)]	int32_t
     *	where (a smod p) is the signed representant of a modulo p, that is -p/2 <= (a smod p) < p/2
     */
    static INLINE CONST vect_t mulx(vect_t a, vect_t b) {
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
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7]	int32_t
     *	   [c0, c1, c2, c3, c4, c5, c6, c7]	int32_t
     * Return :	[(a0*b0+c0) smod 2^32, ..., (a7*b7+c7) smod 2^32]	int32_t
     */
    static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b) { return add(c, mul(a, b)); }

    static INLINE vect_t fmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fmadd(c, a, b); }

    /*
     * Multiply the low 16-bit integers from each packed 32-bit element in a and b,
     * keep the signed 32-bit results and add the low 32-bits of c.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7]	int32_t
     *	   [c0, c1, c2, c3, c4, c5, c6, c7]	int32_t
     * Return :	[((a0 smod 2^16)*(b0 smod 2^16)+c0) smod 2^32, ...,
     *		 ((a7 smod 2^16)*(b7 smod 2^16)+c7) smod 2^32]	int32_t
     */
    static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b) { return add(c, mulx(a, b)); }

    static INLINE vect_t fmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmaddx(c, a, b); }

    /*
     * Multiply the packed 32-bit integers in a and b, producing intermediate 64-bit integers,
     * and substract the low 32 bits of the intermediate from elements of c.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7]	int32_t
     *	   [c0, c1, c2, c3, c4, c5, c6, c7]	int32_t
     * Return :	[(-a0*b0+c0) smod 2^32, ..., (-a7*b7+c7) smod 2^32]	int32_t
     */
    static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mul(a, b)); }

    static INLINE vect_t fnmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmadd(c, a, b); }

    /*
     * Multiply the low 16-bit integers from each packed 32-bit element in a and b,
     * keep the signed 32-bit results and substract them from elements of c.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7]	int32_t
     *	   [c0, c1, c2, c3, c4, c5, c6, c7]	int32_t
     * Return :	[(-(a0 smod 2^16)*(b0 smod 2^16)+c0) smod 2^32, ...,
     *		 (-(a7 smod 2^16)*(b7 smod 2^16)+c7) smod 2^32]	int32_t
     */
    static INLINE CONST vect_t fnmaddx(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mulx(a, b)); }

    static INLINE vect_t fnmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmaddx(c, a, b); }

    /*
     * Multiply the packed 32-bit integers in a and b, producing intermediate 64-bit integers,
     * and substract elements of c to the low 32-bits of the intermediate.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7]	int32_t
     *	   [c0, c1, c2, c3, c4, c5, c6, c7]	int32_t
     * Return : [(a0*b0-c0) smod 2^32, ..., (a7*b7-c7) smod 2^32]	int32_t
     */
    static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b) { return sub(mul(a, b), c); }

    static INLINE vect_t fmsubin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsub(c, a, b); }

    /*
     * Multiply the low 16-bit integers from each packed 32-bit element in a and b,
     * keep the signed 32-bit results and substract elements of c from them.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7]	int32_t
     *	   [c0, c1, c2, c3, c4, c5, c6, c7]	int32_t
     * Return :	[((a0 smod 2^16)*(b0 smod 2^16)-c0) smod 2^32, ...,
     *		 ((a7 smod 2^16)*(b7 smod 2^16)-c7) smod 2^32]	int32_t
     */
    static INLINE CONST vect_t fmsubx(const vect_t c, const vect_t a, const vect_t b) { return sub(mulx(a, b), c); }

    static INLINE vect_t fmsubxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsubx(c, a, b); }

    /*
     * Compare packed 32-bits in a and b for equality, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7]	int32_t
     * Return : [(a0==b0) ? 0xFFFFFFFF : 0, (a1==b1) ? 0xFFFFFFFF : 0,
     *	    (a2==b2) ? 0xFFFFFFFF : 0, (a3==b3) ? 0xFFFFFFFF : 0,
     *	    (a4==b4) ? 0xFFFFFFFF : 0, (a5==b5) ? 0xFFFFFFFF : 0,
     *	    (a6==b6) ? 0xFFFFFFFF : 0, (a7==b7) ? 0xFFFFFFFF : 0]	int32_t
     */
    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) { return _mm256_cmpeq_epi32(a, b); }

    /*
     * Compare packed 32-bits in a and b for greater-than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7]	int32_t
     * Return : [(a0>b0) ? 0xFFFFFFFF : 0, (a1>b1) ? 0xFFFFFFFF : 0,
     *	    (a2>b2) ? 0xFFFFFFFF : 0, (a3>b3) ? 0xFFFFFFFF : 0,
     *	    (a4>b4) ? 0xFFFFFFFF : 0, (a5>b5) ? 0xFFFFFFFF : 0,
     *	    (a6>b6) ? 0xFFFFFFFF : 0, (a7>b7) ? 0xFFFFFFFF : 0]		int32_t
     */
    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) { return _mm256_cmpgt_epi32(a, b); }

    /*
     * Compare packed 32-bits in a and b for lesser-than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7]	int32_t
     * Return : [(a0<b0) ? 0xFFFFFFFF : 0, (a1<b1) ? 0xFFFFFFFF : 0,
     *	    (a2<b2) ? 0xFFFFFFFF : 0, (a3<b3) ? 0xFFFFFFFF : 0,
     *	    (a4<b4) ? 0xFFFFFFFF : 0, (a5<b5) ? 0xFFFFFFFF : 0,
     *	    (a6<b6) ? 0xFFFFFFFF : 0, (a7<b7) ? 0xFFFFFFFF : 0]		int32_t
     */
    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) { return _mm256_cmpgt_epi32(b, a); }

    /*
     * Compare packed 32-bits in a and b for greater or equal than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7]	int32_t
     * Return : [(a0>=b0) ? 0xFFFFFFFF : 0, (a1>=b1) ? 0xFFFFFFFF : 0,
     *	    (a2>=b2) ? 0xFFFFFFFF : 0, (a3>=b3) ? 0xFFFFFFFF : 0,
     *	    (a4>=b4) ? 0xFFFFFFFF : 0, (a5>=b5) ? 0xFFFFFFFF : 0,
     *	    (a6>=b6) ? 0xFFFFFFFF : 0, (a7>=b7) ? 0xFFFFFFFF : 0]	int32_t
     */
    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

    /*
     * Compare packed 32-bits in a and b for lesser or equal than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	int32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7]	int32_t
     * Return : [(a0<=b0) ? 0xFFFFFFFF : 0, (a1<=b1) ? 0xFFFFFFFF : 0,
     *	    (a2<=b2) ? 0xFFFFFFFF : 0, (a3<=b3) ? 0xFFFFFFFF : 0,
     *	    (a4<=b4) ? 0xFFFFFFFF : 0, (a5<=b5) ? 0xFFFFFFFF : 0,
     *	    (a6<=b6) ? 0xFFFFFFFF : 0, (a7<=b7) ? 0xFFFFFFFF : 0]	int32_t
     */
    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

    /*
     * Horizontally add 32-bits elements of a.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
     * Return : a0+a1+a2+a3+a4+a5+a6+a7
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        Converter ca;
        ca.v = a;
        return scalar_t(ca.t[0] + ca.t[1] + ca.t[2] + ca.t[3] + ca.t[4] + ca.t[5] + ca.t[6] + ca.t[7]);
    }

    static INLINE CONST vect_t round(const vect_t a) { return a; }

    static INLINE vect_t mod(vect_t &C, const vect_t &P, const vect_t &INVP, const vect_t &NEGP, const vect_t &MIN,
                             const vect_t &MAX, vect_t &Q, vect_t &T) {
#ifdef __INTEL_COMPILER
        C = _mm256_rem_epi32(C, P);
#else
        FFLASFFPACK_abort("pas implementé");
        // C = fnmadd(C,_mm256_castps_si128(_mm256_floor_ps(_mm256_mul_ps(INVP,_mm256_castsi128_ps(C)))),P);
#endif
        NORML_MOD(C, P, NEGP, MIN, MAX, Q, T);
        return C;
    }
};

/*
 * Simd256 specialized for uint32_t
 */
template <> struct Simd256_impl<true, true, false, 4> : public Simd256_impl<true, true, true, 4> {

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = uint32_t;

    using aligned_allocator = AlignedAllocator<scalar_t, Alignment(alignment)>;
    using aligned_vector = std::vector<scalar_t, aligned_allocator>;

    /* To check compatibility with Modular struct */
    template <class Field>
    using is_same_element = std::is_same<typename Field::Element, scalar_t>;

    /*
     * Simd128 for scalar_t, to deal half_t
     */
    using simdHalf = Simd128<scalar_t>;

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
     *  Return [x,x,x,x,x,x,x,x] uint32_t
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm256_set1_epi32(x); }

    /*
     *  Set packed 32-bit unsigned integers in dst with the supplied values.
     *  Return [x0,x1,x2,x3,x4,x5,x6,x7] uint32_t
     */
    static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1, const scalar_t x2, const scalar_t x3,
                                   const scalar_t x4, const scalar_t x5, const scalar_t x6, const scalar_t x7) {
        return _mm256_set_epi32(x7, x6, x5, x4, x3, x2, x1, x0);
    }

    /*
     *  Gather 32-bit unsigned integer elements with indexes idx[0], ..., idx[7] from the address p in vect_t.
     *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]], p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]]] uint32_t
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        return set(p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]], p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]]);
    }

    /*
     * Load 256-bits of unsigned integer data from memory into dst.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     * Return [p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]] uint32_t
     */
    static INLINE PURE vect_t load(const scalar_t *const p) {
        return _mm256_load_si256(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Load 256-bits of unsigned integer data from memory into dst.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]] uint32_t
     */
    static INLINE PURE vect_t loadu(const scalar_t *const p) {
        return _mm256_loadu_si256(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Store 256-bits of unsigned integer data from a into memory.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void store(scalar_t *p, vect_t v) {
        _mm256_store_si256(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Store 256-bits of unsigned integer data from a into memory.
     * p does not need to be aligned on any particular boundary.
     */
    static INLINE void storeu(scalar_t *p, vect_t v) {
        _mm256_storeu_si256(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Store 256-bits of unsigned integer data from a into memory using a non-temporal memory hint.
     * p must be aligned on a 32-byte boundary or a general-protection exception may be generated.
     */
    static INLINE void stream(scalar_t *p, const vect_t v) {
        _mm256_stream_si256(reinterpret_cast<vect_t *>(p), v);
    }

    /*
     * Shift packed 32-bit unsigned integers in a right by s while shifting in sign bits, and store the results in vect_t.
     * Args   : [a0, ..., a7]			int32_t
     * Return : [Floor(a0/2^s), ..., Floor(a7/2^s)]	int32_t
     */
    template<int s>
    static INLINE CONST vect_t sra(const vect_t a) { return _mm256_srli_epi32(a, s); }

    static INLINE CONST vect_t greater(vect_t a, vect_t b) {
        vect_t x;
        x = set1((static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1)));
        a = sub(a,x);
        b = sub(b,x);
        return _mm256_cmpgt_epi32(a, b);
    }

    static INLINE CONST vect_t lesser(vect_t a, vect_t b) {
        vect_t x;
        x = set1((static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1)));
        a = sub(a,x);
        b = sub(b,x);
        return _mm256_cmpgt_epi32(b, a);
    }

    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

    /*
     * Multiply the packed unsigned 32-bit integers in a and b, producing intermediate 64-bit integers,
     * and store the high 32	bits of the intermediate integers in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] uint32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7] uint32_t
     * Return : [Floor(a0*b0/2^32), ..., Floor(a7*b7/2^32)] uint32_t
     */
    static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b) {
        // //#pragma warning "The simd mulhi function is emulated, it may impact the performances."
        // typedef Simd256_impl<true, true, false, 8> Simd256_64;
        // vect_t C,A1,B1;
        // C  = Simd256_64::mulx(a,b);
        // A1 = Simd256_64::srl(a,32);
        // B1 = Simd256_64::srl(b,32);
        // A1 = Simd256_64::mulx(A1,B1);
        // C  = Simd256_64::srl(C,32);
        // A1 = Simd256_64::srl(A1,32);
        // A1 = Simd256_64::sll(A1,32);
        // return Simd256_64::vor(C,A1);
        vect_t C,A1,B1;
        //C  = Simd256_64::mulx(a,b);
        C  = _mm256_mul_epu32(a,b);
        //A1 = Simd256_64::srl(a,32);
        A1 = _mm256_srli_epi64(a, 32);
        //B1 = Simd256_64::srl(b,32);
        B1 = _mm256_srli_epi64(b, 32);
        //A1 = Simd256_64::mulx(A1,B1);
        A1 =  _mm256_mul_epu32(A1,B1);
        //C  = Simd256_64::srl(C,32);
        C  = _mm256_srli_epi64(C, 32);
        //A1 = Simd256_64::srl(A1,32);
        A1 = _mm256_srli_epi64(A1, 32);
        //A1 = Simd256_64::sll(A1,32);
        A1 = _mm256_slli_epi64(A1, 32);
        //return Simd256_64::vor(C,A1);
        return _mm256_or_si256(A1, C);
    }

    /*
     * Multiply the low unsigned 16-bit integers from each packed 32-bit element in a and b,
     * and store the signed 32-bit results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]	uint32_t
     *	   [b0, b1, b2, b3, b4, b5, b6, b7]	uint32_t
     * Return : [(a0 mod 2^16)*(b0 mod 2^16), ..., (a7 mod 2^16)*(b7 mod 2^16)]	uint32_t
     */
    static INLINE CONST vect_t mulx(vect_t a, vect_t b) {
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
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
     * Return : a0+a1+a2+a3+a4+a5+a6+a7
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        Converter ca;
        ca.v = a;
        return scalar_t(ca.t[0] + ca.t[1] + ca.t[2] + ca.t[3] + ca.t[4] + ca.t[5] + ca.t[6] + ca.t[7]);
    }
}; //Simd256_impl<true,true,false,4>

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd256_int32_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
