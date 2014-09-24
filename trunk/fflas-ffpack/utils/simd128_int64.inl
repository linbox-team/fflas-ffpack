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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd128_int64_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd128_int64_INL

// int64_t
template<>
struct Simd128_impl<true, true, true, 8> {

	using vect_t = __m128i;
	using half_t = __m128i;

	using scalar_t = int64_t;

	static const constexpr size_t vect_size = 2;

	static const constexpr size_t alignment = 16;

	static INLINE CONST vect_t zero()
	{
		return _mm_setzero_si128 ();
	}


	static INLINE PURE vect_t load(const scalar_t * const p)
	{
		return _mm_load_si128(reinterpret_cast<const vect_t *>(p));
	}


	static INLINE PURE vect_t loadu(const scalar_t * const p)
	{
		return _mm_loadu_si128(reinterpret_cast<const vect_t *>(p));
	}


	static INLINE PURE vect_t loadu_half(const scalar_t * const p)
	{
		return loadu(p) ;
	}


	static INLINE PURE store(const scalar_t * p, vect_t v)
	{
		_mm_store_si128(reinterpret_cast<vect_t *>(const_cast<scalar_t*>(p)), v);
	}


	static INLINE void storeu(const scalar_t * p, vect_t v)
	{
		_mm_storeu_si128(reinterpret_cast<vect_t *>(const_cast<scalar_t*>(p)), v);
	}


	static INLINE void storeu_half(const scalar_t * p, vect_t v)
	{
		storeu(p,v) ;
	}

	static INLINE void store_half(const scalar_t * p, vect_t v)
	{
		store(p,v) ;
	}

	static INLINE CONST vect_t set1(const scalar_t x)
	{
		return _mm_set1_epi64x(x);
	}

	{
		return _mm_set_epi64(x2, x1);
	}


	static INLINE CONST vect_t add(const vect_t a, const vect_t b)
	{
		return _mm_add_epi64(a, b);
	}

	static INLINE vect_t addin(vect_t &a, const vect_t b)
	{
		return a = add(a,b);
	}


	static INLINE CONST vect_t sub(const vect_t a, const vect_t b)
	{
		return _mm_sub_epi64(a, b);
	}


	static INLINE CONST vect_t eq(const vect_t a, const vect_t b)
	{
		return _mm_cmpeq_epi64(a, b);
	}


	/* could do in 4.2 using cmp(eq-gt)_epi64 */

	// static INLINE CONST vect_t lesser(const vect_t a, const vect_t b)
	{
		return _mm_cmplt_epi64(a, b);
	}


	// static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b)
	{
		return _mm_cmple_epi64(a, b);
	}


	// static INLINE CONST vect_t greater(const vect_t a, const vect_t b)
	{
		return _mm_cmpgt_epi64(a, b);
	}


	// static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b)
	{
		return _mm_cmpge_epi64(a, b);
	}


	// static INLINE CONST vect_t vand(const vect_t a, const vect_t b)
	{
		return _mm_and_epi64(a, b);
	}


	// static INLINE CONST vect_t vor(const vect_t a, const vect_t b)
	{
		return _mm_or_epi64(a, b);
	}


	static INLINE CONST scalar_t hadd_to_scal(const vect_t a)
	{

		return ((const scalar_t*)&a)[0] + ((const scalar_t*)&a)[1];

	}


	// XXX bug : x0 read-only
	static INLINE CONST vect_t mullo(const vect_t x0, const vect_t x1)
	{

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


	static INLINE CONST vect_t mulhi(const vect_t x0, const vect_t x1)
	{

#pragma warning "The simd mulhifunction is emulate, it may impact the performances."
		// karatsuba
		// TODO

	}


	static INLINE CONST vect_t mulx(const vect_t x0, const vect_t x1)
	{

#pragma warning "The simd mulx function does not make sense, rethink your code :)"

	}


	// I NEED THOSE
	static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b)
	{

		return _mm_add_epi64(c,_mm_mul_epi32(a,b));

	}

	static INLINE vect_t maddin(vect_t & c, const vect_t a, const vect_t b)
	{
		return c = madd(c,a,b);
	}



}
;

// uint64_t
template<>
struct Simd128_impl<true, true, false, 8>
{

	// static void hello()
	{
		std::cout << "uint64_t" << std::endl;
	}


} ;

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd128_int64_INL
