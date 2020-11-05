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
     * Unpack and interleave 64-bit integers from the low half of a and b within 128-bit lanes, and store the results in dst.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int64_t
     [b0, b1, b2, b3, a4, a5, a6, a7] int64_t
     * Return : [a0, b0, a2, b2, a4, b4, a6, b6] int64_t
     */
    static INLINE CONST vect_t unpacklo_twice(const vect_t a, const vect_t b) {
        //std::cerr<<"unpacklo_twice simd512_int64"<<std::endl;
        return _mm512_unpacklo_epi64(a, b); }

    /*
     * Unpack and interleave 64-bit integers from the high half of a and b within 128-bit lanes, and store the results in dst.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int64_t
     [b0, b1, b2, b3, a4, a5, a6, a7] int64_t
     * Return : [a0, b0, a2, b2, a4, b4, a6, b6] int64_t
     */
    static INLINE CONST vect_t unpackhi_twice(const vect_t a, const vect_t b) {
        //std::cerr<<"unpackhi_twice simd512_int64"<<std::endl;
        return _mm512_unpackhi_epi64(a, b); }

    /*
     * Unpack and interleave 64-bit integers from the low half of a and b, and store the results in dst.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int64_t
     [b0, b1, b2, b3, b4, b5, b6, b7] int64_t
     * Return : [a0, b0, a1, b1, a2, b2, a3, b3] int64_t
     */

    static INLINE CONST vect_t unpacklo(const vect_t a, const vect_t b) {
        //std::cerr<<"unpacklo simd512_int64"<<std::endl;
        __m256i a1 = _mm512_castsi512_si256(a); // a1 = [a0, a1, a2, a3]
        __m256i b1 = _mm512_castsi512_si256(b); // b1 = [b0, b1, b2, b3]
        __m256i a2 = _mm256_permute4x64_epi64(a1, 0xD8); // a2 = [a0, a2, a1, a3] (0xD8 a1 <-> a2)
        __m256i b2 = _mm256_permute4x64_epi64(b1, 0xD8); // b2 = [b0, b2, b1, b3] (0xD8 b1 <-> b2)
        __m256i low = _mm256_unpacklo_epi64(a2, b2); // low = [a0, bo, a1, b1]
        __m256i high = _mm256_unpackhi_epi64(a2, b2); // high = [a2, b2, a3, b3]
        __m512i res = _mm512_castsi256_si512(low);
        res = _mm512_inserti64x4(res, high, 1);
        return res;
    }

    /*
     * Unpack and interleave 64-bit integers from the high half of a and b, and store the results in dst.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int64_t
     [b0, b1, b2, b3, b4, b5, b6, b7] int64_t
     * Return : [a4, b4, a5, b5, a6, b6, a7, b7] int64_t
     */

    static INLINE CONST vect_t unpackhi(const vect_t a, const vect_t b) {
        //std::cerr<<"unpackhi simd512_int64"<<std::endl;
        __m256i a1 = _mm512_extracti64x4_epi64(a,1); // a1 = [a4, a5, a6, a7]
        __m256i b1 = _mm512_extracti64x4_epi64(b,1); // b1 = [b4, b5, b6, b7]
        __m256i a2 = _mm256_permute4x64_epi64(a1, 0xD8); // a2 = [a4, a6, a5, a7] (0xD8 a5 <-> a6)
        __m256i b2 = _mm256_permute4x64_epi64(b1, 0xD8); // b2 = [b4, b6, b5, b7] (0xD8 b5 <-> b6)
        __m256i low = _mm256_unpacklo_epi64(a2, b2); // low = [a0, bo, a1, b1]
        __m256i high = _mm256_unpackhi_epi64(a2, b2); // high = [a2, b2, a3, b3]
        __m512i res = _mm512_castsi256_si512(low);
        res = _mm512_inserti64x4(res, high, 1);
        return res;
    }

    /*
     * Unpack and interleave 64-bit integers from the low then high half of a and b, and store the results in dst.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int64_t
     [b0, b1, b2, b3, b4, b5, b6, b7] int64_t
     * Return : [a0, b0, a1, b1, a2, b2, a3, b3] int64_t
     *		   [a4, b4, a5, b5, a6, b6, a7, b7] int64_t
     */

    static INLINE void unpacklohi(vect_t& l, vect_t& h, const vect_t a, const vect_t b) {
        l = unpacklo(a, b);
        h = unpackhi(a, b);
    }

    /*
     * Blend packed 64-bit integers from a and b using control mask imm8, and store the results in dst.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int64_t
     [b0, b1, b2, b3, b4, b5, b6, b7] int64_t
     * Return : [s[0]?a0:b0,   , s[7]?a7:b7] int64_t
     */
    template<uint8_t s>
    static INLINE CONST vect_t blend(const vect_t a, const vect_t b) {
        // _mm_blend_epi16 is faster than _mm_blend_epi32 and require SSE4.1 instead of AVX2
        // We have to transform s = [d3 d2 d1 d0]_base2 to s1 = [d3 d3 d2 d2 d1 d1 d0 d0]_base2
        //constexpr uint8_t s1 = (s & 0x1) * 3 + (((s & 0x2) << 1)*3)  + (((s & 0x4) << 2)*3) + (((s & 0x8) << 3)*3);
        //std::cerr<<"blend simd512_int64"<<std::endl;
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

#ifdef __FFLASFFPACK_HAVE_INT128
    static INLINE CONST vect_t mulhi(vect_t a, vect_t b) {
        //#pragma warning "The simd mulhi function is emulate, it may impact the performances."
        // ugly solution, but it works.
        // tested with gcc, clang, icc
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set((scalar_t)((int128_t(ca.t[0]) * cb.t[0]) >> 64), (scalar_t)((int128_t(ca.t[1]) * cb.t[1]) >> 64),
                   (scalar_t)((int128_t(ca.t[2]) * cb.t[2]) >> 64), (scalar_t)((int128_t(ca.t[3]) * cb.t[3]) >> 64),
                   (scalar_t)((int128_t(ca.t[4]) * cb.t[4]) >> 64), (scalar_t)((int128_t(ca.t[5]) * cb.t[5]) >> 64),
                   (scalar_t)((int128_t(ca.t[6]) * cb.t[6]) >> 64), (scalar_t)((int128_t(ca.t[7]) * cb.t[7]) >> 64));
    }
#endif

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
#ifdef __FFLASFFPACK_HAVE_INT128
    static INLINE CONST vect_t mulhi(vect_t a, vect_t b) {
        //#pragma warning "The simd mulhi function is emulate, it may impact the performances."
        // ugly solution, but it works.
        // tested with gcc, clang, icc
        Converter c0, c1;
        c0.v = a;
        c1.v = b;
        return set((scalar_t)(((uint128_t)(c0.t[0]) * c1.t[0]) >> 64), (scalar_t)(((uint128_t)(c0.t[1]) * c1.t[1]) >> 64),
                   (scalar_t)(((uint128_t)(c0.t[2]) * c1.t[2]) >> 64), (scalar_t)(((uint128_t)(c0.t[3]) * c1.t[3]) >> 64),
                   (scalar_t)(((uint128_t)(c0.t[4]) * c1.t[4]) >> 64), (scalar_t)(((uint128_t)(c0.t[5]) * c1.t[5]) >> 64),
                   (scalar_t)(((uint128_t)(c0.t[6]) * c1.t[6]) >> 64), (scalar_t)(((uint128_t)(c0.t[7]) * c1.t[7]) >> 64));
    }
#endif

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
