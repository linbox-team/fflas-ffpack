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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd128_int32_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd128_int32_INL

// int32_t
template<>
struct Simd128_impl<true, true, true, 4> {

	using vect_t = __m128i;

	using scalar_t = int32_t;

	static const constexpr size_t vect_size = 4;

	static const constexpr size_t alignment = 16;

	static INLINE CONST vect_t zero()
	{
		return _mm_setzero_si128();
	}


	static INLINE PURE vect_t load(const scalar_t * const p)
	{
		return _mm_load_si128(reinterpret_cast<const vect_t *>(p));
	}


	static INLINE PURE vect_t loadu(const scalar_t * const p)
	{
		return _mm_loadu_si128(reinterpret_cast<const vect_t *>(p));
	}


	static INLINE PURE store(const scalar_t * p, vect_t v)
	{
		_mm_store_si128(reinterpret_cast<vect_t *>(const_cast<scalar_t*>(p)), v);
	}


	static INLINE void storeu(const scalar_t * p, vect_t v)
	{
		_mm_storeu_si128(reinterpret_cast<vect_t *>(const_cast<scalar_t*>(p)), v);
	}


	static INLINE CONST vect_t set1(const scalar_t x)
	{
		return _mm_set1_epi32(x);
	}

	static INLINE CONST vect_t set(const scalar_t x1, const scalar_t x2
				       , const scalar_t x3, const scalar_t x4
				       )
	{
		return _mm_set_epi32(x4, x3, x2, x1);
	}


	static INLINE CONST vect_t add(const vect_t a, const vect_t b)
	{
		return _mm_add_epi32(a, b);
	}

	static INLINE vect_t addin(vect_t &a, const vect_t b)
	{
		return a = add(a,b);
	}

	static INLINE CONST vect_t sub(const vect_t a, const vect_t b)
	{
		return _mm_sub_epi32(a, b);
	}

	static INLINE CONST vect_t mul(const vect_t a, const vect_t b)
	{
		return _mm_mul_epi32(a,b);
	}

	static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b)
	{
		return add(c,mul(a,b));
	}

	static INLINE vect_t maddin(vect_t & c, const vect_t a, const vect_t b)
	{
		return c = madd(c,a,b);
	}

	static INLINE CONST vect_t nmadd(const vect_t c, const vect_t a, const vect_t b)
	{
		return sub(c,mul(a,b));
	}


	static INLINE CONST vect_t msub(const vect_t c, const vect_t a, const vect_t b)
	{
		return sub(mul(a,b),c);
	}


	static INLINE CONST vect_t eq(const vect_t a, const vect_t b)
	{
		return _mm_cmpeq_epi32(a, b);
	}

	static INLINE CONST vect_t lesser(const vect_t a, const vect_t b)
	{
		return _mm_cmplt_epi32(a, b);
	}


	static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b)
	{
		return vor(eq(a,b),lesser(a,b));
	}


	static INLINE CONST vect_t greater(const vect_t a, const vect_t b)
	{
		return _mm_cmpgt_epi32(a, b);
	}


	static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b)
	{
		return vor(eq(a,b),lesser(a,b));
	}


	static INLINE CONST vect_t vand(const vect_t a, const vect_t b)
	{
		return _mm_and_si128(a, b);
	}


	static INLINE CONST vect_t vor(const vect_t a, const vect_t b)
	{
		return _mm_or_si128(a, b);
	}


	static INLINE CONST scalar_t hadd_to_scal(const vect_t a)
	{
		return ((const scalar_t*)&a)[0] + ((const scalar_t*)&a)[1] + ((const scalar_t*)&a)[2] + ((const scalar_t*)&a)[3];
	}


	static INLINE CONST vect_t mullo(const vect_t x0, const vect_t x1)
	{
		return _mm_mullo_epi16(x0, x1);
	}


	static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b)
	{
		// _mm_mulhi_epi32 emul
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


	static INLINE CONST vect_t mulx(const vect_t a, const vect_t b)
	{
		return _mm_mul_epi16(a,b);
	}

	static INLINE CONST vect_t maddx(const vect_t c, const vect_t a, const vect_t b)
	{
		return add(mulx(a, b),c);
	}

	static INLINE vect_t maddxin(vect_t & c, const vect_t a, const vect_t b)
	{
		return c = maddx(c,a,b);
	}


} ;

// uint32_t
template<>
struct Simd128_impl<true, true, false, 4>
{

	// static void hello()
	{
		std::cout << "uint32_t" << std::endl;
	}


} ;

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd128_int32_INL
