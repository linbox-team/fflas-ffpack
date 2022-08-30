/*
 * Copyright (C) 2014 FFLAS-FFPACK group
 *
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *
 * Part of this code is taken from http://libdivide.com/
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


#ifndef __FFLASFFPACK_utils_bit_manipulation_H
#define __FFLASFFPACK_utils_bit_manipulation_H

#ifndef __has_builtin
#define __has_builtin(x) 0  // Compatibility with non-clang compilers.
#endif

#include <givaro/udl.h>
#include "fflas-ffpack/fflas-ffpack-config.h"

// count leading zeros
inline int32_t clz(uint64_t val) {
#if __GNUC__ || __has_builtin(__builtin_clzll)
    return __builtin_clzll(val);
#else
    if (! val) return 64 ;
    int32_t result = 0;
    while (! (val & (1_ui64 << 63))) {
        val <<= 1;
        result++;
    }
    return result;
#endif
}

inline int32_t clz(uint32_t val) {
#if __GNUC__ || __has_builtin(__builtin_clzll)
    return __builtin_clz(val);
#else
    if (! val) return 32 ;
    int32_t result = 0;
    while (! (val & (1 << 31))) {
        val <<= 1;
        result++;
    }
    return result;
#endif
}

// count trailing zeros
inline int32_t ctz(uint32_t val) {
#if __GNUC__ || __has_builtin(__builtin_ctz)
    return __builtin_ctz(val);
#else
    if (!val) return 32;
    int32_t result = 0;
    val = (val ^ (val - 1)) >> 1;  // Set v's trailing 0s to 1s and zero rest
    while (val) {
        val >>= 1;
        result++;
    }
    return result;
#endif
}

// count trailing zeros
inline int32_t ctz(uint64_t val) {
#if __GNUC__ || __has_builtin(__builtin_ctzll)
    return __builtin_ctzll(val);
#else
    if (!val) return 64;
    uint32_t lo = val & 0xFFFFFFFF;
    if (lo != 0) return ctz(lo);
    return 32 + ctz(val >> 32);
#endif
}



#if defined (__FFLASFFPACK_HAVE_INT128) && defined(__x86_64__)
// division 128bits by 64 bits
// int128_t(u1,u0) = u1*2^64+u0, div v, rem r
// return quo
static uint64_t divide_128(uint64_t u1, uint64_t u0, uint64_t v, uint64_t *r)
{
    // u0 -> rax
    // u1 -> rdx
    // divq
    uint64_t result;
    __asm__("divq %[v]"
            : "=a"(result), "=d"(*r)
            : [v] "r"(v), "a"(u0), "d"(u1)
           );
    return result;
}
#endif

static inline uint64_t getpoweroftwoden_128(uint32_t d, uint64_t q, uint64_t *r) {
#if defined (__FFLASFFPACK_HAVE_INT128) && defined(__x86_64__)
    return divide_128(1_ui64 << (d - 1), 0, q, r);
#else
    lldiv_t ta;
    ta = lldiv(1ULL<<63,q);
    lldiv_t br;
    br = lldiv(ta.rem<<d,q);
    *r = br.rem;
    return (ta.quot<<d)+br.quot;
#endif
}



static inline uint32_t mullhi_u32(uint32_t x, uint32_t y) {
    uint64_t xl = x, yl = y;
    uint64_t rl = xl * yl;
    return (uint32_t)(rl >> 32);
}

static inline int64_t mulhi_64(int64_t x, int64_t y) {
#ifdef __FFLASFFPACK_HAVE_INT128
    int128_t xl = x, yl = y;
    int128_t rl = xl * yl;
    return (int64_t)(rl >> 64);
#else
    const uint32_t mask = 0xFFFFFFFF;
    const uint32_t x0 = (uint32_t)(x & mask), y0 = (uint32_t)(y & mask);
    const int32_t x1 = (int32_t)(x >> 32), y1 = (int32_t)(y >> 32);
    const uint32_t x0y0_hi = mullhi_u32(x0, y0);
    const int64_t t = x1*(int64_t)y0 + x0y0_hi;
    const int64_t w1 = x0*(int64_t)y1 + (t & mask);
    return x1*(int64_t)y1 + (t >> 32) + (w1 >> 32);
#endif
}

static inline uint64_t mulhi_u64(uint64_t x, uint64_t y) {
#ifdef __FFLASFFPACK_HAVE_INT128
    uint128_t xl = x, yl = y;
    uint128_t rl = xl * yl;
    return (uint64_t)(rl >> 64);
#else
    const uint64_t mask_lo = 0xFFFFFFFFULL;

    uint64_t x0 = x & mask_lo, x1 = x >> 32;
    uint64_t y0 = y & mask_lo, y1 = y >> 32;

    uint64_t xy_hi  = x1 * y1;
    uint64_t xy_mid = x1 * y0;
    uint64_t yx_mid = x0 * y1;
    uint64_t xy_lo  = x0 * y0;

    uint64_t carry_bit = ((xy_mid & mask_lo) +
                          (yx_mid & mask_lo) +
                          (xy_lo >> 32)) >> 32;

    return xy_hi + (xy_mid >> 32) + (yx_mid >> 32) + carry_bit;
#endif
}

static inline int64_t mulhi_fast_64(int64_t x, int64_t y) {
#ifdef __FFLASFFPACK_HAVE_INT128
    int128_t xl = x, yl = y;
    int128_t rl = xl * yl;
    return (int64_t)(rl >> 64);
#else
    const uint32_t mask = 0xFFFFFFFF;
    const uint32_t x0 = (uint32_t)(x & mask), y0 = (uint32_t)(y & mask);
    const int32_t x1 = (int32_t)(x >> 32), y1 = (int32_t)(y >> 32);
    // const uint32_t x0y0_hi = libdivide__mullhi_u32(x0, y0);
    const int64_t t = x1*(int64_t)y0  ; // + x0y0_hi;
    const int64_t w1 = x0*(int64_t)y1 ; // + (t & mask);
    return x1*(int64_t)y1 + (t >> 32) + (w1 >> 32);
#endif
}



#endif // __FFLASFFPACK_utils_bit_manipulation_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
