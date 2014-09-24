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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd128_int16_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd128_int16_INL

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

// uint16_t
template<>
struct Simd128_impl<true, true, false, 2>{
	// static void hello(){std::cout << "uint16_t" << std::endl;}
};

#endif
