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
    * Converter from vect_t to a tab.
    * exple:
    *      Converter conv;
    *      conv.v = a;
    *      scalart_t x = conv.t[1]
    */
    union Converter{
        vect_t v;
        scalar_t t[vect_size];
    };

    /*
     *  Return vector of type vect_t with all elements set to zero
     *  Return [0,0] int64_t
     */
    static INLINE CONST vect_t zero()
    {
        return _mm_setzero_si128();
    }

    /*
     *  Broadcast 64-bit integer a to all all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x,x] int64_t
     */
    static INLINE CONST vect_t set1(const scalar_t x)
    {
        return _mm_set1_epi64x(x);
    }

    /*
     *  Broadcast 64-bit integer a to all all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x0,x1] int64_t
     */
    static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1)
    {
        return _mm_set_epi64x(x1,x0);
    }

    /*
     *  Gather 64-bit integer elements with indexes idx[0], ..., idx[1] from the address p in vect_t.
     *  Return [p[idx[0]], p[idx[1]]] int64_t
     */
    template<class T>
    static INLINE PURE vect_t gather(const scalar_t * const p, const T * const idx)
    {
        return set(p[idx[0]], p[idx[1]]);
    }

    /*
     * Load 128-bits of integer data from memory into dst.
     * p must be aligned on a 16-byte boundary or a general-protection exception will be generated.
     * Return [p[0],p[1]] int64_t
     */
    static INLINE PURE vect_t load(const scalar_t * const p)
    {
        return _mm_load_si128(reinterpret_cast<const vect_t*>(p));
    }

    /*
     * Load 128-bits of integer data from memory into dst.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0],p[1]] int64_t
     */
    static INLINE PURE vect_t loadu(const scalar_t * const p)
    {
        return _mm_loadu_si128(reinterpret_cast<const vect_t*>(p));
    }

    /*
     * Store 128-bits of integer data from a into memory.
     * p must be aligned on a 16-byte boundary or a general-protection exception will be generated.
     */
    static INLINE void store(const scalar_t * p, vect_t v)
    {
        _mm_store_si128(reinterpret_cast<vect_t*>(const_cast<scalar_t*>(p)), v);
    }

    /*
     * Store 128-bits of integer data from a into memory.
     * p does not need to be aligned on any particular boundary.
     */
    static INLINE void storeu(const scalar_t * p, vect_t v)
    {
        _mm_storeu_si128(reinterpret_cast<vect_t*>(const_cast<scalar_t*>(p)), v);
    }

    /*
     * Add packed 64-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1] int64_t
                [b0, b1] int64_t
     * Return : [a0+b0, a1+b1]   int64_t
     */
	static INLINE CONST vect_t add(const vect_t a, const vect_t b)
	{
		return _mm_add_epi64(a, b);
	}

	static INLINE vect_t addin(vect_t &a, const vect_t b)
	{
		return a = add(a,b);
	}

	/*
     * Subtract packed 64-bit integers in b from packed 64-bit integers in a, and store the results in vect_t.
     * Args   : [a0, a1] int64_t
                [b0, b1] int64_t
     * Return : [a0-b0, a1-b1]  int64_t
     */
	static INLINE CONST vect_t sub(const vect_t a, const vect_t b)
	{
		return _mm_sub_epi64(a, b);
	}

	/*
     * Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers, and store the low 64 bits of the intermediate integers in vect_t. 
     * Args   : [a0, a1]           int64_t
                [b0, b1]           int64_t
     * Return : [a0*b0 mod 2^16-1, a1*b1 mod 2^16-1] int64_t
     */
	static INLINE CONST vect_t mullo(const vect_t x0, const vect_t x1)
	{
		// _mm_mullo_epi32 emul
#pragma warning "The simd mullo function is emulate, it may impact the performances."

		Converter c0, c1;
		c0.v = x0;
		c1.v = x1;
		return set((__int128(c0.t[0])*c1.t[0]), (__int128(c0.t[1])*c1.t[1]));
	}

	/*
     * Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers, and store the low 64 bits of the intermediate integers in vect_t. 
     * Args   : [a0, a1]           int64_t
                [b0, b1]           int64_t
     * Return : [a0*b0 mod 2^16-1, a1*b1 mod 2^16-1] int64_t
     */
	static INLINE CONST vect_t mul(const vect_t a, const vect_t b)
	{
		return mullo(a, b);
	}

	static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b)
	{
#pragma warning "The simd mulhi function is emulate, it may impact the performances."
		Converter c0, c1;
		c0.v = a;
		c1.v = b;
		return set((__int128(c0.t[0])*c1.t[0]) >> 64, (__int128(c0.t[1])*c1.t[1]) >> 64);
	}

	static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b)
	{
		return add(c,mul(a,b));
	}

	static INLINE vect_t fmaddin(vect_t & c, const vect_t a, const vect_t b)
	{
		return c = fmadd(c,a,b);
	}

	static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b)
	{
		return sub(c,mul(a,b));
	}


	static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b)
	{
		return sub(mul(a,b),c);
	}


	static INLINE CONST vect_t eq(const vect_t a, const vect_t b)
	{
		return _mm_cmpeq_epi64(a, b);
	}

	static INLINE CONST vect_t lesser(const vect_t a, const vect_t b)
	{
#ifdef __SSE4_2__
		return _mm_cmpgt_epi64(b, a);
#else
#pragma warning "The simd lesser function is emulate, it may impact the performances."
		Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set((ca.t[0] < cb.t[0]) ? 0xFFFFFFFFFFFFFFFF : 0, (ca.t[1] < cb.t[1]) ? 0xFFFFFFFFFFFFFFFF : 0);
#endif
	}

	static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b)
	{
		return vor(eq(a,b),lesser(a,b));
	}


	static INLINE CONST vect_t greater(const vect_t a, const vect_t b)
	{
#ifdef __SSE4_2__
		return _mm_cmpgt_epi64(a, b);
#else
#pragma warning "The simd greater function is emulate, it may impact the performances."
        Converter ca, cb;
        ca.v = a;
        cb.v = b;
        return set((ca.t[0] > cb.t[0]) ? 0xFFFFFFFFFFFFFFFF : 0, (ca.t[1] > cb.t[1]) ? 0xFFFFFFFFFFFFFFFF : 0);
#endif
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
        Converter c;
        c.v = a;
		return c.t[0] + c.t[1];
	}

	static INLINE CONST vect_t mulx(const vect_t a, const vect_t b)
	{
		return _mm_mul_epi32(a,b);
	}

	static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b)
	{
		return add(mulx(a, b),c);
	}

	static INLINE vect_t fmaddxin(vect_t & c, const vect_t a, const vect_t b)
	{
		return c = fmaddx(c,a,b);
	}

#else
#error "You need SSE instructions to perform 128 bits operations on int64"
#endif

} ;

// uint64_t
template<>
struct Simd128_impl<true, true, false, 8> {

	// static void hello()
	// {
		// std::cout << "uint64_t" << std::endl;
	// }


} ;

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd128_int64_INL
