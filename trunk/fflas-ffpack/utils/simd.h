/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   <bastien.vialla@lirmm.fr>
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
struct Simd;

template<>
struct Simd<double>
{
#ifdef __AVX__	
	using vect_t = __m256d;

	 static const constexpr size_t vect_size = 4;

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

     static INLINE CONST vect_t axpy(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __AVX2__
        return _mm256_madd_pd(a, b, c);
#else
        return _mm256_add_pd(c, _mm256_mul_pd(a, b));
#endif
    }

     static INLINE CONST vect_t saxpy(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __AVX2__
        return _mm256_msub_pd(a, b, c);
#else
        return _mm256_sub_pd(c, _mm256_mul_pd(a, b));
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

     static INLINE CONST vect_t hadd(const vect_t a, const vect_t b) {return _mm256_hadd_pd(a, b);}

     static INLINE CONST double hadd_to_scal(const vect_t a) {
         return ((double*)&a)[0] + ((double*)&a)[1] + ((double*)&a)[2] + ((double*)&a)[3];
     }

#else // __AVX__

     using vect_t = __m128d;

     static const constexpr size_t vect_size = 2;

// TODO double sse

#endif
};

template<>
struct Simd<float>
{
#ifdef __AVX__
	using vect_t = __m256;
    
    static const constexpr size_t vect_size = 8;
    
    static INLINE CONST vect_t zero() {return _mm256_setzero_ps();}
	
    static INLINE CONST vect_t set1(const double x) {return _mm256_set1_ps(x);}
    
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
    
    static INLINE CONST vect_t axpy(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __AVX2__
        return _mm256_madd_ps(a, b, c);
#else
        return _mm256_add_ps(c, _mm256_mul_ps(a, b));
#endif
    }
       
    static INLINE CONST vect_t saxpy(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __AVX2__
        return _mm256_msub_ps(a, b, c);
#else
        return _mm256_sub_ps(c, _mm256_mul_ps(a, b));
#endif
    }
    
    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) {return _mm256_cmp_ps(a, b, 0);}
    
    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {return _mm256_cmp_ps(a, b, 1);}
    
    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) {return _mm256_cmp_ps(a, b, 2);}
    
    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {return _mm256_cmp_ps(a, b, 14);}
    
    static INLINE CONST vect_t greater_eq(vect_t & v, const vect_t & a, const vect_t & b) {return _mm256_cmp_ps(a, b, 13);}
    
    static INLINE CONST vect_t vand(vect_t & v, const vect_t & a, const vect_t & b) {return _mm256_and_ps(a, b);}
    
    static INLINE CONST vect_t vor(vect_t & v, const vect_t & a, const vect_t & b) {return _mm256_or_ps(a, b);}
    
    static INLINE CONST vect_t floor(vect_t & v, const vect_t & a) {return _mm256_floor_ps(a);}
    
    static INLINE CONST vect_t hadd(vect_t & v, const vect_t & a, const vect_t & b) {return _mm256_hadd_ps(a, b);}
    
    static INLINE CONST float hadd_to_scal(const vect_t a) {
        return ((float*)&a)[0] + ((float*)&a)[1] + ((float*)&a)[2] + ((float*)&a)[3] + ((float*)&a)[4] + ((float*)&a)[5] + ((float*)&a)[6] + ((float*)&a)[7];
    }
    
#else // __AVX__
    
    using vect_t = __m128;

     static const constexpr size_t vect_size = 4;

    // TODO double sse
    
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
