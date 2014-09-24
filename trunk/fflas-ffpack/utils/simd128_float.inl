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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd128_float_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd128_float_INL

// float
template<>
struct Simd128_impl<true, false, true, 4> {

#if defined(__FFLASFFPACK_USE_SIMD)

	using vect_t = __m128;

	using scalar_t = float;

	static const constexpr size_t vect_size = 4;

	static const constexpr size_t alignment = 16;

	static INLINE CONST vect_t zero()
	{
		return _mm_setzero_ps();
	}


	static INLINE CONST vect_t set1(const scalar_t x)
	{
		return _mm_set1_ps(x);
	}


	static INLINE CONST vect_t set(const scalar_t x1, const scalar_t x2, const scalar_t x3, const scalar_t x4
				       )
	{
		return _mm_set_ps(x4, x3, x2, x1);
	}


	template<class T>
	static INLINE PURE vect_t gather(const scalar_t * const p, const T * const idx)
	{
		// TODO AVX2 Gather
		return _mm256_set_ps(p[idx[3]], p[idx[2]], p[idx[1]], p[idx[0]]);
	}


	static INLINE PURE vect_t load(const scalar_t * const p)
	{
		return _mm_load_ps(p);
	}


	static INLINE PURE vect_t loadu(const scalar_t * const p)
	{
		return _mm_loadu_ps(p);
	}


	static INLINE void store(const scalar_t * p, const vect_t v)
	{
		_mm_store_ps(const_cast<scalar_t*>(p), v);
	}


	static INLINE CONST vect_t add(const vect_t a, const vect_t b)
	{
		return _mm_add_ps(a, b);
	}

	static INLINE vect_t addin(vect_t &a, const vect_t b)
	{
		return a = add(a,b);
	}


	static INLINE CONST vect_t sub(const vect_t a, const vect_t b)
	{
		return _mm_sub_ps(a, b);
	}


	static INLINE CONST vect_t mul(const vect_t a, const vect_t b)
	{
		return _mm_mul_ps(a, b);
	}


	static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b)
	{
		return add(c, mul(a, b));
	}

	static INLINE vect_t maddin(vect_t & c, const vect_t a, const vect_t b)
	{
		return c = madd(c,a,b);
	}

	static INLINE CONST vect_t nmadd(const vect_t c, const vect_t a, const vect_t b)
	{
		return sub(c, mul(a, b));
	}


	static INLINE CONST vect_t msub(const vect_t c, const vect_t a, const vect_t b)
	{
		return sub(mul(a, b), c);
	}


	static INLINE CONST vect_t eq(const vect_t a, const vect_t b)
	{
		return _mm_cmpeq_ps(a, b);
	}


	static INLINE CONST vect_t lesser(const vect_t a, const vect_t b)
	{
		return _mm_cmplt_ps(a, b);
	}


	static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b)
	{
		return _mm_cmple_ps(a, b);
	}


	static INLINE CONST vect_t greater(const vect_t a, const vect_t b)
	{
		return _mm_cmpgt_ps(a, b);
	}


	static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b)
	{
		return _mm_cmpge_ps(a, b);
	}


	static INLINE CONST vect_t vand(const vect_t a, const vect_t b)
	{
		return _mm_and_ps(a, b);
	}


	static INLINE CONST vect_t vor(const vect_t a, const vect_t b)
	{
		return _mm_or_ps(a, b);
	}


	static INLINE CONST vect_t floor(const vect_t a)
	{
#ifdef __SSE4_1__
		return _mm_floor_ps(a);
#else
		__m128 one = _mm_set1_ps(1.0f);
		__m128 fval = _mm_cvtepi32_ps(_mm_cvttps_epi32(a));
		return _mm_sub_ps(fval, _mm_and_ps(_mm_cmplt_ps(a, fval), one));
#endif
	}


	static INLINE CONST vect_t round(const vect_t a)
	{
#ifdef __SSE4_1__
		return _mm_round_ps(a, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC);
#else
		return _mm_cvtepi32_ps(_mm_cvtps_epi32(a));
#endif
	}


	static INLINE CONST vect_t hadd(const vect_t a, const vect_t b)
	{
#ifdef __SSE3__
		return _mm_hadd_ps(a, b);
#endif
	}

	static INLINE CONST scalar_t hadd_to_scal(const vect_t a)
	{
		return ((const scalar_t*)&a)[0] + ((const scalar_t*)&a)[1] + ((const scalar_t*)&a)[2] + ((const scalar_t*)&a)[3];
	}

#else // SIMD
#error "You need SSE instructions to perform 128 bits operations on float"
#endif

} ;

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd128_float_INL
