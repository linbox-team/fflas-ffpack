/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   Bastien Vialla<bastien.vialla@lirmm.fr>
 * BB <bbboyer@ncsu.edu>
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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd128_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd128_INL

template<bool ArithType, bool Int, bool Signed, int Size>
struct Simd128_impl;

// float
template<>
struct Simd128_impl<true, false, true, 4>{
    using vect_t = __m128;

    using scalar_t = float;

    static const constexpr size_t vect_size = 4;

    static const constexpr size_t alignment = 16;

    static INLINE CONST vect_t zero() {return _mm_setzero_ps();}

    static INLINE CONST vect_t set1(const scalar_t x) {return _mm_set1_ps(x);}

    static INLINE CONST vect_t set(const scalar_t x1, const scalar_t x2, const scalar_t x3, const scalar_t x4) {return _mm_set_ps(x4, x3, x2, x1);}

    template<class T>
    static INLINE PURE vect_t gather(const scalar_t * const p, const T * const idx) {
        // TODO AVX2 Gather
        return _mm256_set_ps(p[idx[3]], p[idx[2]], p[idx[1]], p[idx[0]]);
    }

    static INLINE PURE vect_t load(const scalar_t * const p) {return _mm_load_ps(p);}

    static INLINE PURE vect_t loadu(const scalar_t * const p) {return _mm_loadu_ps(p);}

    static INLINE void store(const scalar_t * p, vect_t v) {_mm_store_ps(const_cast<scalar_t*>(p), v);}

    static INLINE CONST vect_t add(const vect_t a, const vect_t b) {return _mm_add_ps(a, b);}
     static INLINE vect_t addin(vect_t &a, const vect_t b) {return a = add(a,b);}

    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) {return _mm_sub_ps(a, b);}

    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) {return _mm_mul_ps(a, b);}

    static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b) {
        return _mm_add_ps(c, _mm_mul_ps(a, b));
    }

    static INLINE vect_t maddin(vect_t & c, const vect_t a, const vect_t b) { return c = madd(c,a,b); }

    static INLINE CONST vect_t nmadd(const vect_t c, const vect_t a, const vect_t b) {
        return _mm_sub_ps(c, _mm_mul_ps(a, b));
    }

    static INLINE CONST vect_t msub(const vect_t c, const vect_t a, const vect_t b) {
        return _mm_sub_ps(_mm_mul_ps(a, b), c);
    }

    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) {return _mm_cmpeq_ps(a, b);}

    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {return _mm_cmplt_ps(a, b);}

    static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) {return _mm_cmple_ps(a, b);}

    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {return _mm_cmpgt_ps(a, b);}

    static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) {return _mm_cmpge_ps(a, b);}

    static INLINE CONST vect_t vand(const vect_t a, const vect_t b) {return _mm_and_ps(a, b);}

    static INLINE CONST vect_t vor(const vect_t a, const vect_t b) {return _mm_or_ps(a, b);}

    static INLINE CONST vect_t floor(const vect_t a) {
#ifdef __SSE4_1__
        return _mm_floor_ps(a);
#else
        __m128 one = _mm_set1_ps(1.0f);
         __m128 fval = _mm_cvtepi32_ps(_mm_cvttps_epi32(a));
         return _mm_sub_ps(fval, _mm_and_ps(_mm_cmplt_ps(a, fval), one));
#endif
    }

    static INLINE CONST vect_t round(const vect_t a) {
#ifdef __SSE4_1__
        return _mm_round_ps(a, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC);
#else
        return _mm_cvtepi32_ps(_mm_cvtps_epi32(a));
#endif
    }

#ifdef __SSE3__
    static INLINE CONST vect_t hadd(const vect_t a, const vect_t b) {
         return _mm_hadd_ps(a, b);
    }
#endif

    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
        return ((const scalar_t*)&a)[0] + ((const scalar_t*)&a)[1] + ((const scalar_t*)&a)[2] + ((const scalar_t*)&a)[3];
    }
};

// double
template<>
struct Simd128_impl<true, false, true, 8>{
    using vect_t = __m128d;

    using scalar_t = double;

     static const constexpr size_t vect_size = 2;

     static const constexpr size_t alignment = 16;

     static INLINE CONST vect_t zero() {return _mm_setzero_pd();}

     static INLINE CONST vect_t set1(const scalar_t x) {return _mm_set1_pd(x);}

     static INLINE CONST vect_t set(const scalar_t x1, const scalar_t x2) {return _mm_set_pd(x2, x1);}

     template<class T>
     static INLINE PURE vect_t gather(const scalar_t * const p, const T * const idx) {
        // TODO AVX2 Gather
        return _mm_set_pd(p[idx[1]], p[idx[0]]);
    }

    static INLINE PURE vect_t load(const scalar_t * const p) {return _mm_load_pd(p);}

     static INLINE PURE vect_t loadu(const scalar_t * const p) {return _mm_loadu_pd(p);}

     static INLINE void store(const scalar_t * p, vect_t v) {_mm_store_pd(const_cast<scalar_t*>(p), v);}

     static INLINE CONST vect_t add(const vect_t a, const vect_t b) {return _mm_add_pd(a, b);}
     static INLINE vect_t addin(vect_t &a, const vect_t b) {return a = add(a,b);}

     static INLINE CONST vect_t sub(const vect_t a, const vect_t b) {return _mm_sub_pd(a, b);}

     static INLINE CONST vect_t mul(const vect_t a, const vect_t b) {return _mm_mul_pd(a, b);}

     static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b) {
        return _mm_add_pd(c, _mm_mul_pd(a, b));
    }

    static INLINE CONST vect_t nmadd(const vect_t c, const vect_t a, const vect_t b) {
        return _mm_sub_pd(c, _mm_mul_pd(a, b));
    }

    static INLINE CONST vect_t msub(const vect_t c, const vect_t a, const vect_t b) {
        return _mm_sub_pd(_mm_mul_pd(a, b), c);
    }

     static INLINE CONST vect_t eq(const vect_t a, const vect_t b) {return _mm_cmpeq_pd(a, b);}

     static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {return _mm_cmplt_pd(a, b);}

     static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) {return _mm_cmple_pd(a, b);}

     static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {return _mm_cmpgt_pd(a, b);}

     static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) {return _mm_cmpge_pd(a, b);}

     static INLINE CONST vect_t vand(const vect_t a, const vect_t b) {return _mm_and_pd(a, b);}

     static INLINE CONST vect_t vor(const vect_t a, const vect_t b) {return _mm_or_pd(a, b);}

     static INLINE CONST vect_t floor(const vect_t a) {
#ifdef __SSE4_1__
         return _mm_floor_pd(a);
#else
         return _mm_set_pd(std::floor(((const scalar_t*)&a)[1]), std::floor(((const scalar_t*)&a)[0]));
#endif
     }


     static INLINE CONST vect_t round(const vect_t a) {
#ifdef __SSE4_1__
         return _mm_round_pd(a, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC);
#else
         return _mm_set_pd(std::round(((const scalar_t*)&a)[1]), std::round(((const scalar_t*)&a)[0]));
#endif
     }

     static INLINE CONST vect_t hadd(const vect_t a, const vect_t b) {
#ifdef __SSE3__
         return _mm_hadd_pd(a, b);
#else
        return _mm_set_pd(((const scalar_t*)&a)[0] + ((const scalar_t*)&a)[1], ((const scalar_t*)&b)[0] + ((const scalar_t*)&b)[1]);
#endif
     }

     static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
         return ((const scalar_t*)&a)[0] + ((const scalar_t*)&a)[1];
     }

};

#ifdef SIMD_INT
// Trop d'instructions SSE manquantes pour les int8_t

// int8_t
// template<>
// struct Simd128_impl<true, true, true, 1>{
//     using vect_t = __m128i;

//     using scalar_t = int8_t;

//     static const constexpr size_t vect_size = 16;

//     static const constexpr size_t alignment = 16;

//     static INLINE PURE vect_t zero(){return _mm_setzero_si128();}

//     static INLINE PURE vect_t load(const scalar_t * const p) {return _mm_load_si128(p);}

//     static INLINE PURE vect_t loadu(const scalar_t * const p) {return _mm_loadu_si128(p);}

//     static INLINE void store(const scalar_t * p, vect_t v) {_mm_store_si128(const_cast<scalar_t*>(p), v);}

//     static INLINE void storeu(const scalar_t * p, vect_t v) {_mm_storeu_si128(const_cast<scalar_t*>(p), v);}

//     static INLINE CONST vect_t set1(const scalar_t x) {return _mm_set1_epi16(x);} // actually set2

//     static INLINE CONST vect_t set(const scalar_t x1, const scalar_t x2, const scalar_t x3, const scalar_t x4, const scalar_t x5, const scalar_t x6, const scalar_t x7, const scalar_t x8, const scalar_t x9, const scalar_t x10, const scalar_t x11, const scalar_t x12, const scalar_t x13, const scalar_t x14, const scalar_t x15, const scalar_t x16){
//         return _mm_set_epi16(x16,x15,x14,x13,x12,x11,x10,x9,x8, x7, x6, x5, x4, x3, x2, x1);
//     }

//     static INLINE CONST vect_t add(const vect_t a, const vect_t b) {return _mm_add_epi8(a, b);}
//     static INLINE vect_t addin(vect_t &a, const vect_t b) {return a = add(a,b);}

//     static INLINE CONST vect_t sub(const vect_t a, const vect_t b) {return _mm_sub_epi8(a, b);}

//     static INLINE CONST vect_t mul(const vect_t a, const vect_t b) {return _mm_mul_epi8(a,b);}

//     static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b) {
//         return _mm_add_epi8(c,_mm_mul_epi8(a,b));
//     }

//     static INLINE CONST vect_t nmadd(const vect_t c, const vect_t a, const vect_t b) {
//         return _mm_sub_epi8(c,_mm_mul_epi8(a,b));
//     }

//     static INLINE CONST vect_t msub(const vect_t c, const vect_t a, const vect_t b) {
//         return _mm_sub_epi8(_mm_mul_epi8(a,b),c);
//     }

//     static INLINE CONST vect_t eq(const vect_t a, const vect_t b) {return _mm_cmpeq_epi8(a, b);}

//     static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {return _mm_cmplt_epi8(a, b);}

//     static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) {return _mm_cmple_epi8(a, b);}

//     static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {return _mm_cmpgt_epi8(a, b);}

//     static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) {return _mm_cmpge_epi8(a, b);}

//     static INLINE CONST vect_t vand(const vect_t a, const vect_t b) {return _mm_and_epi8(a, b);}

//     static INLINE CONST vect_t vor(const vect_t a, const vect_t b) {return _mm_or_epi8(a, b);}

//     static INLINE CONST vect_t mullo(const vect_t x0, const vect_t x1) {
// #pragma warning "The simd mullo function is emulate, it will induce bad performances."
//         vect_t x2 = x0, x3 = x1;
//         x3 = _mm_unpacklo_epi8(x0, x3);
//         x2 = _mm_unpacklo_epi8(x1, x2);
//         x1 = _mm_unpackhi_epi8(x1, x2);
//         x0 = _mm_unpackhi_epi8(x0, x1);
//     }

//     static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b){
//         return _mm_mulhi_epi(a, b););
//     }

//     static INLINE CONST vect_t mulx(const vect_t a, const vect_t b){
// #pragma warning "The simd mulhx function is emulate, it will induce bad performances."
//         vect_t ah, al;
//         ah = mulhi(a, b);
//         al = mullo(a,b);
//         return set(_mm_extract_epi16(ah,1),_mm_extract_epi16(al,1),_mm_extract_epi16(ah,3),_mm_extract_epi16(al,3),_mm_extract_epi16(ah,5),_mm_extract_epi16(al,5),_mm_extract_epi16(ah,7),_mm_extract_epi16(al,7));
//     }

//     static INLINE CONST vect_t maddx(const vect_t c, const vect_t a, const vect_t b){
// #pragma warning "The simd mulhi function is emulate, it will induce bad performances."
//         return _mm_add_epi32(mulx(a, b),c);
//     }

// };

// int16_t
template<>
struct Simd128_impl<true, true, true, 2>{
    using vect_t = __m128i;

    using scalar_t = int16_t;

    static const constexpr size_t vect_size = 8;

    static const constexpr size_t alignment = 16;

    static INLINE PURE vect_t zero(){return _mm_setzero_si128();}

    static INLINE PURE vect_t load(const scalar_t * const p) {return _mm_load_si128(reinterpret_cast<const vect_t *>(p));}

    static INLINE PURE vect_t loadu(const scalar_t * const p) {return _mm_loadu_si128(reinterpret_cast<const vect_t *>(p));}

    static INLINE void store(const scalar_t * p, vect_t v) {_mm_store_si128(reinterpret_cast<vect_t *>(const_cast<scalar_t*>(p)), v);}

    static INLINE void storeu(const scalar_t * p, vect_t v) {_mm_storeu_si128(reinterpret_cast<vect_t *>(const_cast<scalar_t*>(p)), v);}

    static INLINE CONST vect_t set1(const scalar_t x) {return _mm_set1_epi16(x);} // actually set2

    static INLINE CONST vect_t set(const scalar_t x1, const scalar_t x2, const scalar_t x3, const scalar_t x4, const scalar_t x5, const scalar_t x6, const scalar_t x7, const scalar_t x8){
        return _mm_set_epi16(x8, x7, x6, x5, x4, x3, x2, x1);
    }

    static INLINE CONST vect_t add(const vect_t a, const vect_t b) {return _mm_add_epi16(a, b);}
    static INLINE vect_t addin(vect_t &a, const vect_t b) {return a = add(a,b);}

    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) {return _mm_sub_epi16(a, b);}

	// the mul_epi16 does not exist
	// The following function does not have a proper definition
	// SEE mullo, mulhi, mulx function

	//static INLINE CONST vect_t mul(const vect_t a, const vect_t b) {return _mm_mul_epi16(a,b);}

	//static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b) {
	//	return _mm_add_epi16(c,_mm_mul_epi16(a,b));
	//}
	// static INLINE vect_t maddin(vect_t &c, const vect_t a, const vect_t b) { return c = madd(c,a,b); }

	// static INLINE CONST vect_t nmadd(const vect_t c, const vect_t a, const vect_t b) {
	// 	return _mm_sub_epi16(c,_mm_mul_epi16(a,b));
	// }

	// static INLINE CONST vect_t msub(const vect_t c, const vect_t a, const vect_t b) {
	// 	return _mm_sub_epi16(_mm_mul_epi16(a,b),c);
	// }

    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) {return _mm_cmpeq_epi16(a, b);}

    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {return _mm_cmplt_epi16(a, b);}

    // static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) {return _mm_cmple_epi16(a, b);}

    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {return _mm_cmpgt_epi16(a, b);}

    // static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) {return _mm_cmpge_epi16(a, b);}

    // static INLINE CONST vect_t vand(const vect_t a, const vect_t b) {return _mm_and_epi16(a, b);}

    // static INLINE CONST vect_t vor(const vect_t a, const vect_t b) {return _mm_or_epi16(a, b);}

    static INLINE CONST vect_t mullo(const vect_t a, const vect_t b){return _mm_mullo_epi16(a, b);}

    static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b){
        return _mm_mulhi_epi16(a, b);
    }

    static INLINE CONST vect_t mulx(const vect_t a, const vect_t b){
#pragma warning "The simd mulx function is emulate, it may impact the performances."
        vect_t ah, al;
        ah = mulhi(a, b);
        al = mullo(a,b);
        return set(_mm_extract_epi16(ah,1),_mm_extract_epi16(al,1),_mm_extract_epi16(ah,3),_mm_extract_epi16(al,3),_mm_extract_epi16(ah,5),_mm_extract_epi16(al,5),_mm_extract_epi16(ah,7),_mm_extract_epi16(al,7));
    }

    static INLINE CONST vect_t maddx(const vect_t c, const vect_t a, const vect_t b){
#pragma warning "The simd maddx function is emulate, it may impact the performances."
        return _mm_add_epi32(mulx(a, b),c);
    }

    // static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
    //   return ((const scalar_t*)&a)[0] + ((const scalar_t*)&a)[1] + ((const scalar_t*)&a)[2] + ((const scalar_t*)&a)[3];
    //  }
};

// int32_t
template<>
struct Simd128_impl<true, true, true, 4>{
    using vect_t = __m128i;

    using scalar_t = int32_t;

    static const constexpr size_t vect_size = 4;

    static const constexpr size_t alignment = 16;

    static INLINE PURE vect_t zero(){return _mm_setzero_si128();}

    static INLINE PURE vect_t load(const scalar_t * const p) {return _mm_load_si128(reinterpret_cast<const vect_t *>(p));}

    static INLINE PURE vect_t loadu(const scalar_t * const p) {return _mm_loadu_si128(reinterpret_cast<const vect_t *>(p));}

    static INLINE void store(const scalar_t * p, vect_t v) {_mm_store_si128(reinterpret_cast<vect_t *>(const_cast<scalar_t*>(p)), v);}

    static INLINE void storeu(const scalar_t * p, vect_t v) {_mm_storeu_si128(reinterpret_cast<vect_t *>(const_cast<scalar_t*>(p)), v);}

    static INLINE CONST vect_t set1(const scalar_t x) {return _mm_set1_epi32(x);} // actually set2

    static INLINE CONST vect_t set(const scalar_t x1, const scalar_t x2, const scalar_t x3, const scalar_t x4){
        return _mm_set_epi32(x4, x3, x2, x1);
    }

    static INLINE CONST vect_t add(const vect_t a, const vect_t b) {return _mm_add_epi32(a, b);}
    static INLINE vect_t addin(vect_t &a, const vect_t b) {return a = add(a,b);}

    static INLINE CONST vect_t sub(const vect_t a, const vect_t b) {return _mm_sub_epi32(a, b);}

    static INLINE CONST vect_t mul(const vect_t a, const vect_t b) {return _mm_mul_epi32(a,b);}


	// The following function does not have a proper definition
	// SEE mullo, mulhi, mulx function

	static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b) {
	    return _mm_add_epi32(c,_mm_mul_epi32(a,b));
	}
	static INLINE vect_t maddin(vect_t & c, const vect_t a, const vect_t b) { return c = madd(c,a,b); }
	// static INLINE CONST vect_t nmadd(const vect_t c, const vect_t a, const vect_t b) {
	//     return _mm_sub_epi32(c,_mm_mul_epi32(a,b));
	// }

	// static INLINE CONST vect_t msub(const vect_t c, const vect_t a, const vect_t b) {
	//     return _mm_sub_epi32(_mm_mul_epi32(a,b),c);
	// }

    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) {return _mm_cmpeq_epi32(a, b);}

    static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {return _mm_cmplt_epi32(a, b);}

    // static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) {return _mm_cmple_epi32(a, b);}

    static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {return _mm_cmpgt_epi32(a, b);}

    // static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) {return _mm_cmpge_epi32(a, b);}

    // static INLINE CONST vect_t vand(const vect_t a, const vect_t b) {return _mm_and_epi32(a, b);}

    // static INLINE CONST vect_t vor(const vect_t a, const vect_t b) {return _mm_or_epi32(a, b);}

    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
      return ((const scalar_t*)&a)[0] + ((const scalar_t*)&a)[1] + ((const scalar_t*)&a)[2] + ((const scalar_t*)&a)[3];
     }

    static INLINE CONST vect_t mullo(const vect_t a, const vect_t b){return _mm_mullo_epi32(a, b);}

    static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b){
#pragma warning "The simd mulhi function is emulate, it may impact the performances."
        vect_t a1, a2, b1, b2;
        a1 = set(0, _mm_extract_epi32(a, 0), 0, _mm_extract_epi32(a, 1));
        a2 = set(0, _mm_extract_epi32(a, 1), 0, _mm_extract_epi32(a, 3));
        b1 = set(0, _mm_extract_epi32(b, 0), 0, _mm_extract_epi32(b, 1));
        b2 = set(0, _mm_extract_epi32(b, 1), 0, _mm_extract_epi32(b, 3));
        a1 = _mm_mul_epi32(a1,b1);
        a2 = _mm_mul_epi32(a2,b2);
        return set(_mm_extract_epi32(a1,0),_mm_extract_epi32(a1,2),_mm_extract_epi32(b1,0),_mm_extract_epi32(b2,0));
    }

    static INLINE CONST vect_t mulx(const vect_t a, const vect_t b){return _mm_mul_epi32(a,b);}

    static INLINE CONST vect_t maddx(const vect_t c, const vect_t a, const vect_t b){
        return _mm_add_epi64(_mm_mul_epi32(a, b),c);
    }
};

// int64_t
template<>
struct Simd128_impl<true, true, true, 8>{
    using vect_t = __m128i;
    using half_t = __m128i;

    using scalar_t = int64_t;

    static const constexpr size_t vect_size = 2;

    static const constexpr size_t alignment = 16;

    static INLINE CONST vect_t zero() { return _mm_setzero_si128 (); }

    static INLINE PURE vect_t load(const scalar_t * const p) {return _mm_load_si128(reinterpret_cast<const vect_t *>(p));}

    static INLINE PURE vect_t loadu(const scalar_t * const p) {return _mm_loadu_si128(reinterpret_cast<const vect_t *>(p));}

    static INLINE PURE vect_t loadu_half(const scalar_t * const p) { return loadu(p) ; }

    static INLINE void store(const scalar_t * p, vect_t v) {_mm_store_si128(reinterpret_cast<vect_t *>(const_cast<scalar_t*>(p)), v);}

    static INLINE void storeu(const scalar_t * p, vect_t v) {_mm_storeu_si128(reinterpret_cast<vect_t *>(const_cast<scalar_t*>(p)), v);}

    static INLINE void storeu_half(const scalar_t * p, vect_t v) { storeu(p,v) ; }
    static INLINE void store_half(const scalar_t * p, vect_t v) { store(p,v) ; }

    // static INLINE void set(const scalar_t x1, const scalar_t x2){return _mm_set_epi64(x2, x1);}

    static INLINE CONST vect_t set1(const scalar_t x) {return _mm_set1_epi64x(x);} // actually set2

    static INLINE CONST vect_t add(vect_t a, vect_t b) {return _mm_add_epi64(a, b);}
     static INLINE vect_t addin(vect_t &a, const vect_t b) {return a = add(a,b);}

    static INLINE CONST vect_t sub(vect_t a, vect_t b) {return _mm_sub_epi64(a, b);}

    static INLINE CONST vect_t eq(const vect_t a, const vect_t b) {return _mm_cmpeq_epi64(a, b);}

    /* could do in 4.2 using cmp(eq-gt)_epi64 */

    // static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {return _mm_cmplt_epi64(a, b);}

    // static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) {return _mm_cmple_epi64(a, b);}

    // static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {return _mm_cmpgt_epi64(a, b);}

    // static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) {return _mm_cmpge_epi64(a, b);}

    // static INLINE CONST vect_t vand(const vect_t a, const vect_t b) {return _mm_and_epi64(a, b);}

    // static INLINE CONST vect_t vor(const vect_t a, const vect_t b) {return _mm_or_epi64(a, b);}

    static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
      return ((const scalar_t*)&a)[0] + ((const scalar_t*)&a)[1];
     }

    // XXX bug : x0 read-only
     static INLINE CONST vect_t mullo(const vect_t x0, const vect_t x1){
#pragma warning "The simd mullo function is emulate, it may impact the performances."
        // karatsuba
        vect_t x2, x3, x4;
        x2 = x0;
	vect_t x0_t, x1_t;
        x0_t = _mm_mul_epi32(x1, x0);
        x3 = x1;
        x4 = x2;
        x3 = _mm_srli_epi64(x3, 32);
        x2 = _mm_mul_epi32(x3, x2);
        x4 = _mm_srli_epi64(x4, 32);
        x1_t = _mm_mul_epi32(x4, x1);
        x1_t = _mm_add_epi64(x2, x1_t);
        x1_t = _mm_srli_epi64(x1_t, 32);
        x0_t = _mm_add_epi64(x1_t, x0_t);
        return x0_t;
     }

     static INLINE CONST vect_t mulhi(const vect_t x0, const vect_t x1){
#pragma warning "The simd mulhifunction is emulate, it may impact the performances."
        // karatsuba
        // TODO
     }

     static INLINE CONST vect_t mulx(const vect_t x0, const vect_t x1){
#pragma warning "The simd mulx function does not make sense, rethink your code :)"
     }

     // I NEED THOSE
     static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b) {
	     return _mm_add_epi64(c,_mm_mul_epi32(a,b));
     }
     static INLINE vect_t maddin(vect_t & c, const vect_t a, const vect_t b) { return c = madd(c,a,b); }

};

// uint8_t
template<>
struct Simd128_impl<true, true, false, 1>{
    // static void hello(){std::cout << "uint8_t" << std::endl;}
};

// uint16_t
template<>
struct Simd128_impl<true, true, false, 2>{
    // static void hello(){std::cout << "uint16_t" << std::endl;}
};

// uint32_t
template<>
struct Simd128_impl<true, true, false, 4>{
    // static void hello(){std::cout << "uint32_t" << std::endl;}
};

// uint64_t
template<>
struct Simd128_impl<true, true, false, 8>{
    // static void hello(){std::cout << "uint64_t" << std::endl;}
};


#endif //#ifdef SIMD_INT

template<class T>
using Simd128 = Simd128_impl<std::is_arithmetic<T>::value, std::is_integral<T>::value, std::is_signed<T>::value, sizeof(T)>;

// template<>
// struct Simd128<int64_t>
// {
	// using vect_t = __m128i;

	// static const constexpr size_t vect_size = 2;

	// static const constexpr size_t alignment = 16;

	// static INLINE PURE vect_t load(const int64_t * const p) {return _mm_load_si128(p);}

	// static INLINE PURE vect_t loadu(const int64_t * const p) {return _mm_loadu_si128(p);}

	// static INLINE void store(const int64_t * p, vect_t v) {_mm_store_si128(reinterpret_cast<__m128i*>(p), v);}

	// static INLINE void storeu(const int64_t * p, vect_t v) {_mm_storeu_si128(reinterpret_cast<__m128i*>(p), v);}

	// static INLINE CONST vect_t set1(const int64_t x) {return _mm_set1_epi64x(x);} // actually set2

	// static INLINE CONST vect_t add(vect_t a, vect_t b) {return _mm_add_epi64(a, b);}

	// static INLINE CONST vect_t mul(vect_t a, vect_t b) {return _mm_add_epi32(a, b);}

	// static INLINE CONST vect_t madd(const vect_t c, vect_t a, vect_t b) {
	// 	return _mm_add_epi64(c, __mm_mul_epi32(a,b));
	// }

// };

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd128_INL
