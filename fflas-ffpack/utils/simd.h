/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   Bastien Vialla<bastien.vialla@lirmm.fr>
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

#ifndef _FFLASFFPAC_simd_h
#define _FFLASFFPAC_simd_h

#include <immintrin.h>

#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define INLINE __attribute__((always_inline)) inline
#else
#define INLINE inline
#endif

#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define CONST __attribute__((const))
#else
#define CONST
#endif

#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define PURE __attribute__((pure))
#else
#define PURE
#endif

template<class T>
 struct simdToType;

 template<>
 struct simdToType<__m256d>
 {
    using type = double;
 };

template<>
 struct simdToType<__m128d>
 {
    using type = double;
 };

 template<>
 struct simdToType<__m256>
 {
    using type = float;
 };

 template<>
 struct simdToType<__m128>
 {
    using type = float;
 };


/*
 * is_simd trait
 */
template<class T>
 struct is_simd
 {
    static const constexpr bool value = false;
    using type = std::integral_constant<bool, false>;
 };

 template<>
 struct is_simd<__m256d>
 {
     static const constexpr bool value = true;
     using type = std::integral_constant<bool, true>;
 };

template<>
 struct is_simd<__m256>
 {
     static const constexpr bool value = true;
     using type = std::integral_constant<bool, true>;
 };

template<>
 struct is_simd<__m128d>
 {
     static const constexpr bool value = true;
     using type = std::integral_constant<bool, true>;
 };

template<>
 struct is_simd<__m128>
 {
     static const constexpr bool value = true;
     using type = std::integral_constant<bool, true>;
 };

/*
 * Simd functor
 */

template<class T>
struct Simd;

template<>
struct Simd<double>
{
#ifdef __AVX__
	using vect_t = __m256d;

	 static const constexpr size_t vect_size = 4;

     static const constexpr size_t alignment = 32;

     static INLINE CONST vect_t zero() {return _mm256_setzero_pd();}

     static INLINE CONST vect_t set1(const double x) {return _mm256_set1_pd(x);}

     static INLINE CONST vect_t set(const double x1, const double x2, const double x3, const double x4) {return _mm256_set_pd(x4, x3, x2, x1);}

     template<class T>
     static INLINE PURE vect_t gather(const double * const p, const T * const idx) {
        // TODO AVX2 Gather
        return _mm256_set_pd(p[idx[3]], p[idx[2]], p[idx[1]], p[idx[0]]);
    }

     static INLINE PURE vect_t load(const double * const p) {return _mm256_load_pd(p);}

     static INLINE PURE vect_t loadu(const double * const p) {return _mm256_loadu_pd(p);}

     static INLINE void store(const double * p, vect_t v) {_mm256_store_pd(const_cast<double*>(p), v);}

     static INLINE CONST vect_t add(vect_t a, vect_t b) {return _mm256_add_pd(a, b);}

     static INLINE CONST vect_t sub(const vect_t a, const vect_t b) {return _mm256_sub_pd(a, b);}

     static INLINE CONST vect_t mul(const vect_t a, const vect_t b) {return _mm256_mul_pd(a, b);}

     static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __AVX2__
        return _mm256_madd_pd(a, b, c);
#else
        return _mm256_add_pd(c, _mm256_mul_pd(a, b));
#endif
    }

    static INLINE CONST vect_t nmadd(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __AVX2__
        return _mm256_nmadd_pd(a, b, c);
#else
        return _mm256_sub_pd(c, _mm256_mul_pd(a, b));
#endif
    }

     static INLINE CONST vect_t msub(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __AVX2__
        return _mm256_msub_pd(a, b, c);
#else
        return _mm256_sub_pd(_mm256_mul_pd(a, b), c);
#endif
    }

     static INLINE CONST vect_t eq(const vect_t a, const vect_t b) {return _mm256_cmp_pd(a, b, 0);}

     static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {return _mm256_cmp_pd(a, b, 1);}

     static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) {return _mm256_cmp_pd(a, b, 2);}

     static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {return _mm256_cmp_pd(a, b, 14);}

     static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) {return _mm256_cmp_pd(a, b, 13);}

     static INLINE CONST vect_t vand(const vect_t a, const vect_t b) {return _mm256_and_pd(a, b);}

     static INLINE CONST vect_t vor(const vect_t a, const vect_t b) {return _mm256_or_pd(a, b);}

     static INLINE CONST vect_t floor(const vect_t a) {return _mm256_floor_pd(a);}

     static INLINE CONST vect_t round(const vect_t a) {return _mm256_round_pd(a, _MM_FROUND_TO_NEAREST_INT);}

     static INLINE CONST vect_t hadd(const vect_t a, const vect_t b) {return _mm256_hadd_pd(a, b);}

     static INLINE CONST double hadd_to_scal(const vect_t a) {
         return ((const double*)&a)[0] + ((const double*)&a)[1] + ((const double*)&a)[2] + ((const double*)&a)[3];
     }

#else // __AVX__

     using vect_t = __m128d;

     static const constexpr size_t vect_size = 2;

     static const constexpr size_t alignment = 16;

     static INLINE CONST vect_t zero() {return _mm_setzero_pd();}

     static INLINE CONST vect_t set1(const double x) {return _mm_set1_pd(x);}

     static INLINE CONST vect_t set(const double x1, const double x2) {return _mm_set_pd(x2, x1);}

     template<class T>
     static INLINE PURE vect_t gather(const double * const p, const T * const idx) {
        // TODO AVX2 Gather
        return _mm_set_pd(p[idx[1]], p[idx[0]]);
    }

    static INLINE PURE vect_t load(const double * const p) {return _mm_load_pd(p);}

     static INLINE PURE vect_t loadu(const double * const p) {return _mm_loadu_pd(p);}

     static INLINE void store(const double * p, vect_t v) {_mm_store_pd(const_cast<double*>(p), v);}

     static INLINE CONST vect_t add(vect_t a, vect_t b) {return _mm_add_pd(a, b);}

     static INLINE CONST vect_t sub(const vect_t a, const vect_t b) {return _mm_sub_pd(a, b);}

     static INLINE CONST vect_t mul(const vect_t a, const vect_t b) {return _mm_mul_pd(a, b);}

     static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b) {
        return _mm_fmadd_pd(a, b, c);
    }

    static INLINE CONST vect_t nmadd(const vect_t c, const vect_t a, const vect_t b) {
        return _mm_sub_pd(c, _mm_mul_pd(a, b));

    }

     static INLINE CONST vect_t msub(const vect_t c, const vect_t a, const vect_t b) {
        return _mm_sub_pd(_mm_mul_pd(a, b), c);
    }

     static INLINE CONST vect_t eq(const vect_t a, const vect_t b) {return _mm_cmp_pd(a, b, 0);}

     static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {return _mm_cmp_pd(a, b, 1);}

     static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) {return _mm_cmp_pd(a, b, 2);}

     static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {return _mm_cmp_pd(a, b, 14);}

     static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) {return _mm_cmp_pd(a, b, 13);}

     static INLINE CONST vect_t vand(const vect_t a, const vect_t b) {return _mm_and_pd(a, b);}

     static INLINE CONST vect_t vor(const vect_t a, const vect_t b) {return _mm_or_pd(a, b);}

     static INLINE CONST vect_t floor(const vect_t a) {return _mm_floor_pd(a);}

     static INLINE CONST vect_t round(const vect_t a) {return _mm_round_pd(a, _MM_FROUND_TO_NEAREST_INT);}

     static INLINE CONST vect_t hadd(const vect_t a, const vect_t b) {return _mm_hadd_pd(a, b);}

     static INLINE CONST double hadd_to_scal(const vect_t a) {
         return ((double*)&a)[0] + ((double*)&a)[1];
     }

#endif
};

template<>
struct Simd<float>
{
#ifdef __AVX__
	using vect_t = __m256;

    static const constexpr size_t vect_size = 8;

    static const constexpr size_t alignment = 32;

    static INLINE CONST vect_t zero() {return _mm256_setzero_ps();}

    static INLINE CONST vect_t set1(const float x) {return _mm256_set1_ps(x);}

    static INLINE CONST vect_t set(const float x1, const float x2, const float x3, const float x4, const float x5, const float x6, const float x7, const float x8) {return _mm256_set_ps(x8, x7, x6, x5, x4, x3, x2, x1);}

    template<class T>
    static INLINE PURE vect_t gather(const float * const p, const T * const idx) {
        // TODO AVX2 Gather
        return _mm256_set_ps(p[idx[7]], p[idx[6]], p[idx[5]], p[idx[4]], p[idx[3]], p[idx[2]], p[idx[1]], p[idx[0]]);
    }

    static INLINE PURE vect_t load(const float * const p) {return _mm256_load_ps(p);}

    static INLINE PURE vect_t loadu(const float * const p) {return _mm256_loadu_ps(p);}

    static INLINE void store(const float * p, vect_t v) {_mm256_store_ps(const_cast<float*>(p), v);}

    static INLINE CONST vect_t add(const vect_t a, const vect_t b) {return _mm256_add_ps(a, b);}

    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) {return _mm256_sub_ps(a, b);}

    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) {return _mm256_mul_ps(a, b);}

    static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __AVX2__
        return _mm256_madd_ps(a, b, c);
#else
        return _mm256_add_ps(c, _mm256_mul_ps(a, b));
#endif
    }

    static INLINE CONST vect_t nmadd(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __AVX2__
        return _mm256_nmadd_ps(a, b, c);
#else
        return _mm256_sub_ps(c, _mm256_mul_ps(a, b));
#endif
    }

    static INLINE CONST vect_t msub(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __AVX2__
        return _mm256_msub_ps(a, b, c);
#else
        return _mm256_sub_ps(_mm256_mul_ps(a, b), c);
#endif
    }

    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) {return _mm256_cmp_ps(a, b, 0);}

    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {return _mm256_cmp_ps(a, b, 1);}

    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) {return _mm256_cmp_ps(a, b, 2);}

    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {return _mm256_cmp_ps(a, b, 14);}

    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) {return _mm256_cmp_ps(a, b, 13);}

    static INLINE CONST vect_t vand(const vect_t a, const vect_t b) {return _mm256_and_ps(a, b);}

    static INLINE CONST vect_t vor(const vect_t a, const vect_t b) {return _mm256_or_ps(a, b);}

    static INLINE CONST vect_t floor(const vect_t a) {return _mm256_floor_ps(a);}

    static INLINE CONST vect_t round(const vect_t a) {return _mm256_round_ps(a, _MM_FROUND_TO_NEAREST_INT);}

    static INLINE CONST vect_t hadd(const vect_t a, const vect_t b) {return _mm256_hadd_ps(a, b);}

    static INLINE CONST float hadd_to_scal(const vect_t a) {
        return ((const float*)&a)[0] + ((const float*)&a)[1] + ((const float*)&a)[2] + ((const float*)&a)[3] + ((const float*)&a)[4] + ((const float*)&a)[5] + ((const float*)&a)[6] + ((const float*)&a)[7];
    }

#else // __AVX__

    using vect_t = __m128;

     static const constexpr size_t vect_size = 4;

     static const constexpr size_t alignment = 16;

    static INLINE CONST vect_t zero() {return _mm_setzero_ps();}

    static INLINE CONST vect_t set1(const float x) {return _mm_set1_ps(x);}

    static INLINE CONST vect_t set(const float x1, const float x2, const float x3, const float x4) {return _mm_set_ps(x4, x3, x2, x1);}

    template<class T>
    static INLINE PURE vect_t gather(const float * const p, const T * const idx) {
        // TODO AVX2 Gather
        return _mm256_set_ps(p[idx[3]], p[idx[2]], p[idx[1]], p[idx[0]]);
    }

    static INLINE PURE vect_t load(const float * const p) {return _mm_load_ps(p);}

    static INLINE PURE vect_t loadu(const float * const p) {return _mm_loadu_ps(p);}

    static INLINE void store(const float * p, vect_t v) {_mm_store_ps(const_cast<float*>(p), v);}

    static INLINE CONST vect_t add(const vect_t a, const vect_t b) {return _mm_add_ps(a, b);}

    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) {return _mm_sub_ps(a, b);}

    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) {return _mm_mul_ps(a, b);}

    static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b) {
        return _mm_fadd_ps(c, _mm_mul_ps(a, b));
    }

    static INLINE CONST vect_t nmadd(const vect_t c, const vect_t a, const vect_t b) {
        return _mm_fsub_ps(c, _mm_mul_ps(a, b));
    }

    static INLINE CONST vect_t msub(const vect_t c, const vect_t a, const vect_t b) {
        return _mm_sub_ps(_mm_mul_ps(a, b), c);
    }

    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) {return _mm_cmp_ps(a, b, 0);}

    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {return _mm_cmp_ps(a, b, 1);}

    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) {return _mm_cmp_ps(a, b, 2);}

    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {return _mm_cmp_ps(a, b, 14);}

    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) {return _mm_cmp_ps(a, b, 13);}

    static INLINE CONST vect_t vand(const vect_t a, const vect_t b) {return _mm_and_ps(a, b);}

    static INLINE CONST vect_t vor(const vect_t a, const vect_t b) {return _mm_or_ps(a, b);}

    static INLINE CONST vect_t floor(const vect_t a) {return _mm_floor_ps(a);}

    static INLINE CONST vect_t round(const vect_t a) {return _mm_round_ps(a, _MM_FROUND_TO_NEAREST_INT);}

    static INLINE CONST vect_t hadd(const vect_t a, const vect_t b) {return _mm_hadd_ps(a, b);}

    static INLINE CONST float hadd_to_scal(const vect_t a) {
        return ((float*)&a)[0] + ((float*)&a)[1] + ((float*)&a)[2] + ((float*)&a)[3];
    }

#endif
};

template<>
struct Simd<long>
{
	// TODO
};

template<>
struct Simd<int>
{
	// TODO
};

#undef INLINE
#undef PURE
#undef CONST

#endif /* _FFLASFFPAC_simd_h */
