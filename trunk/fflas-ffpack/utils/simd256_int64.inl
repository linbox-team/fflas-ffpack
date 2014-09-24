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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd256_int64_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd256_int64_INL

// int64_t
template<>
struct Simd256_impl<true, true, true, 8>{
    #if defined(__FFLASFFPACK_USE_AVX) or defined(__FFLASFFPACK_USE_AVX2)
	using vect_t = __m256i;

	using half_t = __m128i;

	using scalar_t = int64_t;

	static const constexpr size_t vect_size = 4;

	static const constexpr size_t alignment = 32;

	static INLINE PURE vect_t load(const scalar_t * const p) {return _mm256_load_si256(reinterpret_cast<const vect_t*>(p));}

	static INLINE PURE vect_t loadu(const scalar_t * const p) {return _mm256_loadu_si256(reinterpret_cast<const vect_t*>(p));}

	static INLINE PURE half_t load_half(const scalar_t * const p) {return _mm_load_si128(reinterpret_cast<const half_t*>(p));}

	static INLINE PURE half_t loadu_half(const scalar_t * const p) {return _mm_loadu_si128(reinterpret_cast<const half_t*>(p));}

	static INLINE void store(const scalar_t * p, vect_t v) {_mm256_store_si256(reinterpret_cast<vect_t*>(const_cast<scalar_t*>(p)), v);}

	static INLINE void storeu(const scalar_t * p, vect_t v) {_mm256_storeu_si256(reinterpret_cast<vect_t*>(const_cast<scalar_t*>(p)), v);}

	static INLINE void store_half(const scalar_t * p, half_t v) {_mm_store_si128(reinterpret_cast<half_t*>(const_cast<scalar_t*>(p)), v);}

	static INLINE void storeu_half(const scalar_t * p, half_t v) {_mm_storeu_si128(reinterpret_cast<half_t*>(const_cast<scalar_t*>(p)), v);}

	static INLINE CONST vect_t set1(const scalar_t x) {return _mm256_set1_epi64x(x);} // actually set2

	static INLINE CONST vect_t add(const vect_t a, const vect_t b) {
#ifdef __AVX2__
		return _mm256_add_epi64(a, b);
#else
		half_t a1 = _mm256_extractf128_si256(a,0);
		half_t a2 = _mm256_extractf128_si256(a,1);
		half_t b1 = _mm256_extractf128_si256(b,0);
		half_t b2 = _mm256_extractf128_si256(b,1);
		vect_t res ;
		res=_mm256_insertf128_si256(res,_mm_add_epi64(a1,b1),0);
		res=_mm256_insertf128_si256(res,_mm_add_epi64(a2,b2),1);
		return res;
#endif
	}
	static INLINE vect_t addin(vect_t &a, const vect_t b) {return a = add(a,b);}

	static INLINE CONST vect_t mullo(const vect_t x0, const vect_t x1){
#ifdef __AVX2__
		return _mm256_mullo_epi32(x0,x1);

#else
		// vect_t x2, x3, x4;
		// x2 = _mm256_mul_epi32(x1, x0);
		// x4 = _mm256_srli_epi32(x0, 32);
		// x3 = _mm256_srli_epi32(x1, 32);
		// x1 = _mm256_mul_epi32(x1, x4);
		// x0 = _mm256_mul_epi32(x0, x3);
		// x1 = _mm256_add_epi64(x1, x0);
		// x1 = _mm256_srll_epi32(x1, 32);
		// x0 = _mm256_add_epi64(x2, x1);

		/* use two _mm_mullo_epi32 ? */
		return x0;

#endif
	}

	static INLINE CONST vect_t mul(vect_t a, vect_t b) {
#ifdef __AVX2__
		return _mm256_mul_epi32(a, b);
#else
		half_t aa = _mm256_extractf128_si256(a,0);
		half_t ab = _mm256_extractf128_si256(a,1);
		half_t ba = _mm256_extractf128_si256(b,0);
		half_t bb = _mm256_extractf128_si256(b,1);
		vect_t res ;
		res=_mm256_insertf128_si256(res,_mm_mul_epi32(aa,ab),0);
		res=_mm256_insertf128_si256(res,_mm_mul_epi32(ba,bb),1);
		return res ;
#endif
	}

	static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b) {
#ifdef __AVX2__
		return _mm256_add_epi64(c, __mm256_mul_epi32(a,b));
#else
		return add(mul(a,b),c);
#endif
	}
	static INLINE  vect_t maddin( vect_t & c, const vect_t a, const vect_t b) { return c = madd(c,a,b); }

	static INLINE CONST vect_t zero() {return _mm256_setzero_si256();}
#else
#error "no avx"
#endif
};

// uint64_t
template<>
struct Simd256_impl<true, true, false, 8>{
    // static void hello(){std::cout << "uint64_t" << std::endl;}
};


#endif // __FFLASFFPACK_fflas_ffpack_utils_simd256_int64_INL

