/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd128_int64_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd128_int64_INL

/*
 * Simd128 specialized for int64_t
 */
template <> struct Simd128_impl<true, true, true, 8> {

#if defined(__FFLASFFPACK_USE_SIMD)
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
     *  alignement required by scalar_t pointer to be loaded in a vect_t
     */
    static const constexpr size_t alignment = 16;

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
     *      Converter conv;
     *      conv.v = a;
     *      scalart_t x = conv.t[1]
     */
    union Converter {
        vect_t v;
        scalar_t t[vect_size];
    };

    /*
     *  Return vector of type vect_t with all elements set to zero
     *  Return [0,0] int64_t
     */
    static INLINE CONST vect_t zero() { return _mm_setzero_si128(); }

    /*
     *  Broadcast 64-bit integer a to all all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x,x] int64_t
     */
    static INLINE CONST vect_t set1(const scalar_t x) { return _mm_set1_epi64x(x); }

    /*
     *  Broadcast 64-bit integer a to all all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x0,x1] int64_t
     */
    static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1) { return _mm_set_epi64x(x1, x0); }

    /*
     *  Gather 64-bit integer elements with indexes idx[0], ..., idx[1] from the address p in vect_t.
     *  Return [p[idx[0]], p[idx[1]]] int64_t
     */
    template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
        return set(p[idx[0]], p[idx[1]]);
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
    static INLINE void store(const scalar_t *p, vect_t v) {
        _mm_store_si128(reinterpret_cast<vect_t *>(const_cast<scalar_t *>(p)), v);
    }

    /*
     * Store 128-bits of integer data from a into memory.
     * p does not need to be aligned on any particular boundary.
     */
    static INLINE void storeu(const scalar_t *p, vect_t v) {
        _mm_storeu_si128(reinterpret_cast<vect_t *>(const_cast<scalar_t *>(p)), v);
    }

    /*
     * Store 128-bits of integer data from a into memory using a non-temporal memory hint.
     * p must be aligned on a 16-byte boundary or a general-protection exception may be generated.
     */
    // static INLINE void stream(scalar_t *p, const vect_t v) { _mm_stream_si128(static_cast<vect_t *>(p), v); }

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
     * Shift packed 64-bit integers in a left by s while shifting in zeros, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3] int64_t
     * Return : [a0 << s, a1 << s, a2 << s, a3 << s] int64_t
     */
    static INLINE CONST vect_t sll(const vect_t a, const int s) { return _mm_slli_epi64(a, s); }

    /*
     * Shift packed 64-bit integers in a right by s while shifting in zeros, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3] int64_t
     * Return : [a0 >> s, a1 >> s, a2 >> s, a3 >> s] int64_t
     */
    static INLINE CONST vect_t srl(const vect_t a, const int s) { return _mm_srli_epi64(a, s); }

    static INLINE CONST vect_t sra(const vect_t a, const int s) {
#ifdef __AVX512__
        return _mm_sra_epi64(a, set1(s));
#else
        const int b = 63 - s;
        vect_t m = sll(set1(1), b);
        vect_t x = srl(a, s);
        vect_t result = sub(vxor(x, m), m); // result = x^m - m
        return result;
#endif // 512
    }

    /*
     * Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers, and store the low 64
     bits of the intermediate integers in vect_t.
     * Args   : [a0, a1]           int64_t
     [b0, b1]           int64_t
     * Return : [a0*b0 mod 2^16-1, a1*b1 mod 2^16-1] int64_t
     */
    static INLINE CONST vect_t mullo(const vect_t x0, const vect_t x1) {
        // _mm_mullo_epi32 emul
        // #pragma warning "The simd mullo function is emulate, it may impact the performances."

        Converter c0, c1;
        c0.v = x0;
        c1.v = x1;
        return set((scalar_t)(c0.t[0] * c1.t[0]), (scalar_t)(c0.t[1] * c1.t[1]));
    }

    static INLINE CONST vect_t mullox(const vect_t x0, const vect_t x1) { return _mm_mullo_epi32(x0, x1); }

    /*
     * Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers, and store the low 64
     bits of the intermediate integers in vect_t.
     * Args   : [a0, a1]           int64_t
     [b0, b1]           int64_t
     * Return : [a0*b0 mod 2^16-1, a1*b1 mod 2^16-1] int64_t
     */
    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) { return mullo(a, b); }

    static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b) {
// #pragma warning "The simd mulhi function is emulate, it may impact the performances."
#ifdef __X86_64__
        Converter c0, c1;
        c0.v = a;
        c1.v = b;
        return set((scalar_t)((int128_t(c0.t[0]) * c1.t[0]) >> 64), (scalar_t)((int128_t(c0.t[1]) * c1.t[1]) >> 64));
#else
        return zero();
#endif
    }

    static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b) { return add(c, mul(a, b)); }

    static INLINE vect_t fmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fmadd(c, a, b); }

    static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mul(a, b)); }

    static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b) { return sub(mul(a, b), c); }

    static INLINE CONST vect_t mulx(const vect_t a, const vect_t b) { return _mm_mul_epi32(a, b); }

    static INLINE CONST vect_t mulux(const vect_t a, const vect_t b) { return _mm_mul_epu32(a, b); }

    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) { return _mm_cmpeq_epi64(a, b); }

    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {
#ifdef __SSE4_2__
        return _mm_cmpgt_epi64(a, b);
#else
#warning "The simd greater function is emulate, it may impact the performances."
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set((ca.t[0] > cb.t[0]) ? 0xFFFFFFFFFFFFFFFF : 0, (ca.t[1] > cb.t[1]) ? 0xFFFFFFFFFFFFFFFF : 0);
#endif // __SSE4_2__
    }

    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {
#ifdef __SSE4_2__
        return _mm_cmpgt_epi64(b, a);
#else
#warning "The simd lesser function is emulate, it may impact the performances."
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set((ca.t[0] < cb.t[0]) ? 0xFFFFFFFFFFFFFFFF : 0, (ca.t[1] < cb.t[1]) ? 0xFFFFFFFFFFFFFFFF : 0);
#endif // __SSE4_2__
    }

    /*
     * Compare packed 64-bits in a and b for greater or equal than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int64_t
     [b0, b1, b2, b3, b4, b5, b6, b7] int64_t
     * Return : [(a0>=b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1>=b1) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a2>=b2) ? 0xFFFFFFFFFFFFFFFF : 0, (a3>=b3) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a4>=b4) ? 0xFFFFFFFFFFFFFFFF : 0, (a5>=b5) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a6>=b6) ? 0xFFFFFFFFFFFFFFFF : 0, (a7>=b7) ? 0xFFFFFFFFFFFFFFFF : 0]                    int64_t
     */
    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

    /*
     * Compare packed 64-bits in a and b for lesser or equal than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] int64_t
     [b0, b1, b2, b3, b4, b5, b6, b7] int64_t
     * Return : [(a0<=b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1<=b1) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a2<=b2) ? 0xFFFFFFFFFFFFFFFF : 0, (a3<=b3) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a4<=b4) ? 0xFFFFFFFFFFFFFFFF : 0, (a5<=b5) ? 0xFFFFFFFFFFFFFFFF : 0,
     (a6<=b6) ? 0xFFFFFFFFFFFFFFFF : 0, (a7<=b7) ? 0xFFFFFFFFFFFFFFFF : 0]                     int64_t
     */
    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

    /*
     * Compute the bitwise AND of packed 64-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
     [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0 AND b0, a1 AND b1, a2 AND b2, a3 AND b3, a4 AND b4, a5 AND b5, a6 AND b6, a7 AND b7]
     */
    static INLINE CONST vect_t vand(const vect_t a, const vect_t b) { return _mm_and_si128(a, b); }

    /*
     * Compute the bitwise OR of packed 64-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
     [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0 OR b0, a1 OR b1, a2 OR b2, a3 OR b3, a4 OR b4, a5 OR b5, a6 OR b6, a7 OR b7]
     */
    static INLINE CONST vect_t vor(const vect_t a, const vect_t b) { return _mm_or_si128(a, b); }

    /*
     * Compute the bitwise XOR of packed 64-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
     [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0 XOR b0, a1 XOR b1, a2 XOR b2, a3 XOR b3, a4 XOR b4, a5 XOR b5, a6 XOR b6, a7 XOR b7]
     */
    static INLINE CONST vect_t vxor(const vect_t a, const vect_t b) { return _mm_xor_si128(b, a); }

    /*
     * Compute the bitwise AND NOT of packed 64-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
     [b0, b1, b2, b3, b4, b5, b6, b7]
     * Return : [a0 ANDNOT b0, a1 ANDNOT b1, a2 ANDNOT b2, a3 ANDNOT b3, a4 ANDNOT b4, a5 ANDNOT b5, a6 ANDNOT b6, a7
     ANDNOT b7]
     */
    static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) { return _mm_andnot_si128(b, a); }

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

    static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b) { return add(c, mulx(a, b)); }

    static INLINE vect_t fmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmaddx(c, a, b); }

    static INLINE CONST vect_t fnmaddx(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mulx(a, b)); }

    static INLINE vect_t fnmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmaddx(c, a, b); }

    static INLINE CONST vect_t round(const vect_t a) { return a; }

    // mask the high 32 bits of a 64 bits, that is 00000000FFFFFFFF
    static INLINE CONST vect_t mask_high() { return srl(_mm_set1_epi8(-1), 32); }

    static INLINE CONST vect_t signbits(const vect_t x) {
        vect_t signBits = sub(zero(), srl(x, 4*sizeof(scalar_t)-1));
        return signBits;
    }

    // warning : may be off by 1 multiple, but we save a mul...
    static INLINE CONST vect_t mulhi_fast(vect_t x, vect_t y) {
        // unsigned mulhi starts:
        // x1 = xy_high = mulhiu_fast(x,y)
        const vect_t mask = mask_high();

        vect_t x0 = vand(x, mask), x1 = srl(x, 32);
        vect_t y0 = vand(y, mask), y1 = srl(y, 32);

        x0 = mulux(x0, y1); // x0y1
        y0 = mulux(x1, y0); // x1y0
        y1 = mulux(x1, y1); // x1y1

        x1 = vand(y0, mask);
        y0 = srl(y0, 32); // x1y0_lo = x1 // y1yo_hi = y0
        x1 = srl(add(x1, x0), 32);
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

    template <bool overflow, bool poweroftwo>
    static INLINE vect_t mod(vect_t &C, const vect_t &P, const int8_t &shifter, const vect_t &magic, const vect_t &NEGP,
                             const vect_t &MIN, const vect_t &MAX, vect_t &Q, vect_t &T) {
#ifdef __INTEL_COMPILER
        // Works fine with ICC 15.0.1 - A.B.
        // #warning "not tested"
        C = _mm_rem_epi64(C, P);
#else
        if (poweroftwo) {
            Q = srl(C, 63);
            vect_t un = set1(1);
            T = sub(sll(un, shifter), un);
            Q = add(C, vand(Q, T));
            Q = sll(srl(Q, shifter), shifter);
            C = sub(C, Q);
            Q = vand(greater(zero(), Q), P);
            C = add(C, Q);
        } else {
            Q = mulhi_fast(C, magic);
            if (overflow) {
                Q = add(Q, C);
            }
            Q = sra(Q, shifter);
            vect_t q1 = mulux(Q, P);
            vect_t q2 = sll(mulux(srl(Q, 32), P), 32);
            C = sub(C, add(q1, q2));
            T = greater_eq(C, P);
            C = sub(C, vand(T, P));
        }
#endif
        NORML_MOD(C, P, NEGP, MIN, MAX, Q, T);
        return C;
    }

#else

#error "You need SSE instructions to perform 128 bits operations on int64"

#endif // __FFLASFFPACK_USE_SIMD
};

// uint64_t
template <> struct Simd128_impl<true, true, false, 8> : public Simd128_impl<true, true, true, 8> {
    using scalar_t = uint64_t;

    /*
    * Load 128-bits of unsigned integer data from memory into dst.
    * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
    * Return [p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]] int16_t
    */
    static INLINE PURE vect_t load(const scalar_t *const p) {
        return _mm_load_si128(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Load 128-bits of unsigned integer data from memory into dst.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]] int16_t
     */
    static INLINE PURE vect_t loadu(const scalar_t *const p) {
        return _mm_loadu_si128(reinterpret_cast<const vect_t *>(p));
    }

    /*
     * Store 128-bits of unsigned integer data from a into memory.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void store(const scalar_t *p, vect_t v) {
        _mm_store_si128(reinterpret_cast<vect_t *>(const_cast<scalar_t *>(p)), v);
    }

    static INLINE CONST vect_t greater(vect_t a, vect_t b) {
#ifdef __SSE4_2__
        vect_t x;
        x = set1(-(static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1)));
        a = sub(x, a);
        b = sub(x, b);
        return _mm_cmpgt_epi64(a, b);
#else
#warning "The simd greater function is emulate, it may impact the performances."
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set((ca.t[0] > cb.t[0]) ? 0xFFFFFFFFFFFFFFFF : 0, (ca.t[1] > cb.t[1]) ? 0xFFFFFFFFFFFFFFFF : 0);
#endif
    }

    static INLINE CONST vect_t lesser(vect_t a, vect_t b) {
#ifdef __SSE4_2__
        vect_t x;
        x = set1(-(static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1)));
        a = sub(x, a);
        b = sub(x, b);
        return _mm_cmpgt_epi64(a, b);
#else
#warning "The simd greater function is emulate, it may impact the performances."
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set((ca.t[0] < cb.t[0]) ? 0xFFFFFFFFFFFFFFFF : 0, (ca.t[1] < cb.t[1]) ? 0xFFFFFFFFFFFFFFFF : 0);
#endif
    }

    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }
};

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd128_int64_INL
