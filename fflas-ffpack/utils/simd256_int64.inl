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

/*
 * Simd256 specialized for int64_t
 */
template<>
struct Simd256_impl<true, true, true, 8> {

#if defined(__FFLASFFPACK_USE_AVX2)
	/*
     * alias to 256 bit simd register
     */
	using vect_t = __m256i;

	/*
     * alias to 256 bit simd register
     */
	using half_t = __m128i;

	/*
     * define the scalar type corresponding to the specialization
     */
	using scalar_t = int64_t;

	/*
     * Simd128 for scalar_t, to deal half_t
     */
	using simdHalf = Simd128<scalar_t>;

	/*
     *  number of scalar_t in a simd register
     */
	static const constexpr size_t vect_size = 4;

	/*
     *  alignement required by scalar_t pointer to be loaded in a vect_t
     */
	static const constexpr size_t alignment = 32;

	/*
	 * Converter from vect_t to a tab.
	 * exple: 
	 *		Converter conv;
	 *		conv.v = a;
	 *		scalar_t x = conv.t[i]
	 */
	union Converter{
		vect_t v;
		scalar_t t[vect_size];
	};

	/*
     *  Return vector of type vect_t with all elements set to zero
     *  Return [0,0,0,0] int64_t
     */
	static INLINE CONST vect_t zero() 
	{
		return _mm256_setzero_si256();
	}

	/*
     *  Broadcast 64-bit integer a to all all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x,x,x,x] int64_t
     */
	static INLINE CONST vect_t set1(const scalar_t x) 
	{
		return _mm256_set1_epi64X(x);
	}

	/*
     *  Broadcast 64-bit integer a to all all elements of dst. This intrinsic may generate the vpbroadcastw.
     *  Return [x0,x1,x2,x3] int64_t
     */
	static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1, const scalar_t x2, const scalar_t x3)
	{
		return _mm256_set_epi64x(x3,x2,x1,x0);
	}

	/*
     *  Gather 64-bit integer elements with indexes idx[0], ..., idx[3] from the address p in vect_t.
     *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]]] int64_t
     */
    template<class T>
    static INLINE PURE vect_t gather(const scalar_t * const p, const T * const idx) 
    {
        return set(p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]]);
    }

    /*
     * Load 256-bits of integer data from memory into dst.
     * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
     * Return [p[0],p[1],p[2],p[3]] int32_t
     */
	static INLINE PURE vect_t load(const scalar_t * const p) 
	{
		return _mm256_load_si256(reinterpret_cast<const vect_t*>(p));
	}

	/*
     * Load 256-bits of integer data from memory into dst.
     * p does not need to be aligned on any particular boundary.
     * Return [p[0],p[1],p[2],p[3]] int64_t
     */
	static INLINE PURE vect_t loadu(const scalar_t * const p) 
	{
		return _mm256_loadu_si256(reinterpret_cast<const vect_t*>(p));
	}

	/*
	 * Store 256-bits of integer data from a into memory.
	 * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
	 */
	static INLINE void store(const scalar_t * p, vect_t v) 
	{
		_mm256_store_si256(reinterpret_cast<vect_t*>(const_cast<scalar_t*>(p)), v);
	}

	/*
	 * Store 256-bits of integer data from a into memory.
	 * p does not need to be aligned on any particular boundary.
	 */
	static INLINE void storeu(const scalar_t * p, vect_t v) 
	{
		_mm256_storeu_si256(reinterpret_cast<vect_t*>(const_cast<scalar_t*>(p)), v);
	}

	/*
     * Add packed 64-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3] 						   int64_t
     			[b0, b1, b2, b3] 						   int64_t
     * Return : [a0+b0, a1+b1, a2+b2, a3+b3]   int64_t
     */
	static INLINE CONST vect_t add(const vect_t a, const vect_t b) 
	{
		return _mm256_add_epi64(a, b);
	}

	static INLINE vect_t addin(vect_t &a, const vect_t b) 
	{
		return a = add(a,b);
	}

	/*
     * Subtract packed 64-bits integers in b from packed 64-bits integers in a, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3] 						  int64_t
     			[b0, b1, b2, b3] 						  int64_t
     * Return : [a0-b0, a1-b1, a2-b2, a3-b3]  int64_t
     */
	static INLINE CONST vect_t sub(const vect_t a, const vect_t b) 
	{
		return _mm256_sub_epi64(a, b);
	}

	static INLINE CONST vect_t subin(vect_t &a, const vect_t b) 
	{
		return a = sub(a,b);
	}

	/*
     * Shift packed 64-bit integers in a left by s while shifting in zeros, and store the results in vect_t. 
     * Args   : [a0, a1, a2, a3] int64_t
     * Return : [a0 << s, a1 << s, a2 << s, a3 << s] int64_t
     */
     static INLINE CONST vect_t sll(const vect_t a, const int s)
     {
        return _mm256_slli_epi64(a, s);
     }

     /*
     * Shift packed 64-bit integers in a right by s while shifting in zeros, and store the results in vect_t. 
     * Args   : [a0, a1, a2, a3] int64_t
     * Return : [a0 >> s, a1 >> s, a2 >> s, a3 >> s] int64_t
     */
     static INLINE CONST vect_t srl(const vect_t a, const int s)
     {
        return _mm256_srli_epi64(a, s);
     }

	/*
     * Multiply the packed 64-bits integers in a and b, producing intermediate 128-bit integers, and store the low 64 bits of the intermediate integers in vect_t. 
     * Args   : [a0, a1, a2, a3]           						     int64_t
     			[b0, b1, b2, b3]  		 							 int64_t
     * Return : [a0*b0 mod 2^64-1, a1*b1 mod 2^64-1, a2*b2 mod 2^64-1, a3*b3 mod 2^64-1] int64_t
     */
	static INLINE CONST vect_t mullo(vect_t a, vect_t b) 
	{
		vect_t x2, x3, x4;
		x2 = _mm256_mul_epu32(a, b);
		x4 = srl(a, 32);
		x3 = srl(b, 32);
		b = _mm256_mul_epu32(b, x4);
		a = _mm256_mul_epu32(a, x3);
		a = add(a, x3);
		b = sll(x1, 32);
		a = add(x1, x2);
		return a;
	}

	/*
     * Multiply the packed 64-bits integers in a and b, producing intermediate 128-bit integers, and store the low 64 bits of the intermediate integers in vect_t. 
     * Args   : [a0, a1, a2, a3]           						     int64_t
     			[b0, b1, b2, b3]  		 							 int64_t
     * Return : [a0*b0 mod 2^64-1, a1*b1 mod 2^64-1, a2*b2 mod 2^64-1, a3*b3 mod 2^64-1] int64_t
     */
	static INLINE CONST vect_t mul(const vect_t a, const vect_t b)
	{
		return mullo(a, b);
	}

	/*
     * Multiply packed 64-bit integers in a and b, producing intermediate 128-bit integers, and add the low 64-bits of the intermediate with c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3]           int64_t
                [b0, b1, b2, b3]           int64_t
                [c0, c1, c2, c3]           int64_t
     * Return : [(a0*b0 mod 2^64-1)+c0, (a1*b1 mod 2^64-1)+c1, (a2*b2 mod 2^64-1)+c2, (a3*b3 mod 2^64-1)+c3]
     */
	static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b)
	{
		return add(c,mul(a,b));
	}

	static INLINE vect_t fmaddin(vect_t & c, const vect_t a, const vect_t b)
	{
		return c = fmadd(c,a,b);
	}

	/*
     * Multiply packed 64-bit integers in a and b, producing intermediate 128-bit integers, and substract elements of c to the low 64-bit of the intermiate result, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3]           int64_t
                [b0, b1, b2, b3]           int64_t
                [c0, c1, c2, c3]           int64_t
     * Return : [-(a0*b0 mod 2^64-1)+c0, -(a1*b1 mod 2^64-1)+c1, -(a2*b2 mod 2^64-1)+c2, -(a3*b3 mod 2^64-1)+c3]
     */
	static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b)
	{
		return sub(c,mul(a,b));
	}

	/*
     * Multiply packed 64-bit integers in a and b, producing intermediate 128-bit integers, and substract the low 64-bits of the intermediate with c, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3]           int64_t
                [b0, b1, b2, b3]           int64_t
                [c0, c1, c2, c3]           int64_t
     * Return : [(a0*b0 mod 2^64-1)-c0, (a1*b1 mod 2^64-1)-c1, (a2*b2 mod 2^64-1)-c2, (a3*b3 mod 2^64-1)-c3]
     */
	static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b)
	{
		return sub(mul(a,b),c);
	}

	/*
     * Multiply the packed 64-bits integers in a and b, producing intermediate 128-bit integers, and store the high 64 bits of the intermediate integers in vect_t. 
     * Args   : [a0, a1, a2, a3]           						     int64_t
     			[b0, b1, b2, b3]  		 							 int64_t
     * Return : 
     */
	static INLINE CONST vect_t mulhi(vect_t a, vect_t b) 
	{
		// ugly solution, but it works.
		// tested with gcc, clang, icc
		Converter ca, cb;
		ca.v = a;
		cb.v = b;
		return set((__int128(ca.t[0])*cb.t[0]) >> 64, (__int128(ca.t[1])*cb.t[1]) >> 64, (__int128(ca.t[2])*cb.t[2]) >> 64, (__int128(ca.t[3])*cb.t[3]) >> 64);
	}	

	/*
     * Multiply the low 64-bits integers from each packed 128-bit element in a and b, and store the signed 128-bit results in dst. 
     * Args   : [0, a1, 0, a3]    int64_t
     			[0, b1, 0, b3]    int64_t
     * Return : [a1*b1, a3*b2] int128_t
     */
	static INLINE CONST vect_t mulx(const vect_t a, const vect_t b) 
	{
		// ugly solution, but it works.
		// tested with gcc, clang, icc
		Converter ca, cb;
		ca.v = a;
		cb.v = b;
		return set((__int128(ca.t[1])*cb.t[1]) >> 64, (__int128(ca.t[1])*cb.t[1]), (__int128(ca.t[2])*cb.t[2]) >> 64, (__int128(ca.t[2])*cb.t[2]));
	}

	/*
     * Compare packed 64-bits in a and b for equality, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3]								   int32_t
     			[b0, b1, b2, b3] 								   int32_t
     * Return : [(a0==b0) ? 0xFFFF : 0, (a1==b1) ? 0xFFFF : 0,
                 (a2==b2) ? 0xFFFF : 0, (a3==b3) ? 0xFFFF : 0]                     int32_t
     */
	static INLINE CONST vect_t eq(const vect_t a, const vect_t b) 
	{
		return _mm256_cmpeq_epi64(a, b);
	}

	/*
     * Compare packed 64-bits in a and b for greater-than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3]								   int32_t
     			[b0, b1, b2, b3] 								   int32_t
     * Return : [(a0>b0) ? 0xFFFF : 0, (a1>b1) ? 0xFFFF : 0,
                 (a2>b2) ? 0xFFFF : 0, (a3>b3) ? 0xFFFF : 0]                     	int32_t
     */
	static INLINE CONST vect_t greater(const vect_t a, const vect_t b) 
	{
		return _mm256_cmpgt_epi64(a, b);
	}

	/*
     * Compare packed 64-bits in a and b for lesser-than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3] int32_t
     			[b0, b1, b2, b3] int32_t
     * Return : [(a0<b0) ? 0xFFFF : 0, (a1<b1) ? 0xFFFF : 0,
                 (a2<b2) ? 0xFFFF : 0, (a3<b3) ? 0xFFFF : 0] 					  int32_t
     */
	static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) 
	{
		return _mm256_cmpgt_epi64(b, a);
	}

	/*
     * Compare packed 64-bits in a and b for greater or equal than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3]									 int32_t
     			[b0, b1, b2, b3] 									 int32_t
     * Return : [(a0>=b0) ? 0xFFFF : 0, (a1>=b1) ? 0xFFFF : 0,
                 (a2>=b2) ? 0xFFFF : 0, (a3>=b3) ? 0xFFFF : 0,
                 (a4>=b4) ? 0xFFFF : 0, (a5>=b5) ? 0xFFFF : 0,
                 (a6>=b6) ? 0xFFFF : 0, (a7>=b7) ? 0xFFFF : 0]					  int32_t
     */
	static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) 
	{
		return vor(greater(a, b), eq(a, b));
	}

	/*
     * Compare packed 64-bits in a and b for lesser or equal than, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3, a4, a5, a6, a7] 					int32_t
     			[b0, b1, b2, b3, b4, b5, b6, b7] 					int32_t
     * Return : [(a0<=b0) ? 0xFFFF : 0, (a1<=b1) ? 0xFFFF : 0,
                 (a2<=b2) ? 0xFFFF : 0, (a3<=b3) ? 0xFFFF : 0] 		int32_t
     */
	static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) 
	{
		return vor(lesser(a, b), eq(a, b));
	}

	/*
     * Compute the bitwise AND of packed 64-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3] 
     			[b0, b1, b2, b3]
     * Return : [a0 AND b0, a1 AND b1, a2 AND b2, a3 AND b3]
     */
	static INLINE CONST vect_t vand(const vect_t a, const vect_t b) 
	{
		return _mm256_and_si256(b, a);
	}

	/*
     * Compute the bitwise OR of packed 64-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3]
     			[b0, b1, b2, b3]
     * Return : [a0 OR b0, a1 OR b1, a2 OR b2, a3 OR b3]
     */
	static INLINE CONST vect_t vor(const vect_t a, const vect_t b) 
	{
		return _mm256_or_si256(b, a);
	}

	/*
     * Compute the bitwise XOR of packed 64-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3]
     			[b0, b1, b2, b3]
     * Return : [a0 XOR b0, a1 XOR b1, a2 XOR b2, a3 XOR b3]
     */
	static INLINE CONST vect_t vxor(const vect_t a, const vect_t b) 
	{
		return _mm256_xor_si256(b, a);
	}

	/*
     * Compute the bitwise AND NOT of packed 64-bits integer in a and b, and store the results in vect_t.
     * Args   : [a0, a1, a2, a3]
     			[b0, b1, b2, b3]
     * Return : [a0 ANDNOT b0, a1 ANDNOT b1, a2 ANDNOT b2, a3 ANDNOT b3]
     */
	static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) 
	{
		return _mm256_andnot_si256(b, a);
	}

	/*
     * Horizontally add 64-bits elements of a.
     * Args   : [a0, a1, a2, a3]
     * Return : a0+a1+a2+a3
     */
    static INLINE CONST scalar_t hadd_to_scal(const vect_t a)
    {
        Converter ca;
        ca.v = a;
        return ca.t[0]+ca.t[1]+ca.t[2]+ca.t[3];
    }

	static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b)
	{
		// mulx : 64bits x 64bits -> 128bits
		// add : 64bits x 64bits -> 64bits
		// this operation does not make sense :)
		FFLASFFPACK_abort("not implemented yet");
		return zero();
	}

	static INLINE  vect_t fmaddxin(vect_t & c, const vect_t a, const vect_t b)
	{
		return c = fmaddx(c,a,b);
		return zero();
	}

#else

#error "You need AVX2 instructions to perform 256bits operations on int64_t"

#endif // defined(__FFLASFFPACK_USE_AVX2)

} ;

// uint64_t
template<>
struct Simd256_impl<true, true, false, 8> {

} ;

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd256_int64_INL
