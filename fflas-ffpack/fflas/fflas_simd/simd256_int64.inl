/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   Bastien Vialla<bastien.vialla@lirmm.fr>
 * Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 * Romain Lebreton <romain.lebreton@lirmm.fr>
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

#ifndef __FFLASFFPACK_USE_AVX2
#error "You need AVX2 instructions to perform 256bits operations on int64_t"
#endif

/*
 * Simd256 specialized for int64_t
 */
template <> struct Simd256_impl<true, true, true, 8> {

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
	 * Check if the pointer p is a multiple of alignemnt
	 */
	template <class T> static constexpr bool valid(T *p) { return (int64_t)p % alignment == 0; }

	/*
	 * Check if the number n is a multiple of vect_size
	 */
	template <class T> static constexpr bool compliant(T n) { return n % vect_size == 0; }

	/*
	* Converter from vect_t to a tab.
	* exple:
	*	Converter conv;
	*	conv.v = a;
	*	scalart_t x = conv.t[1]
	*/
	union Converter {
		vect_t v;
		scalar_t t[vect_size];
	};

	/*
	 *  Return vector of type vect_t with all elements set to zero
	 *  Return [0,0,0,0] int64_t
	 */
	static INLINE CONST vect_t zero() { return _mm256_setzero_si256(); }

	/*
	 *  Broadcast 64-bit integer a to all elements of dst. This intrinsic may generate the vpbroadcastw.
	 *  Return [x,x,x,x] int64_t
	 */
	static INLINE CONST vect_t set1(const scalar_t x) { return _mm256_set1_epi64x(x); }

	/*
	 *  Set packed 64-bit integers in dst with the supplied values.
	 *  Return [x0,x1,x2,x3] int64_t
	 */
	static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1, const scalar_t x2, const scalar_t x3) {
		return _mm256_set_epi64x(x3, x2, x1, x0);
	}

	/*
	 *  Gather 64-bit integer elements with indexes idx[0], ..., idx[3] from the address p in vect_t.
	 *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]]] int64_t
	 */
	template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
		return set(p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]]);
	}

	/*
	 * Load 256-bits of integer data from memory into dst.
	 * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
	 * Return [p[0],p[1],p[2],p[3]] int32_t
	 */
	static INLINE PURE vect_t load(const scalar_t *const p) {
		return _mm256_load_si256(reinterpret_cast<const vect_t *>(p));
	}

	/*
	 * Load 256-bits of integer data from memory into dst.
	 * p does not need to be aligned on any particular boundary.
	 * Return [p[0],p[1],p[2],p[3]] int64_t
	 */
	static INLINE PURE vect_t loadu(const scalar_t *const p) {
		return _mm256_loadu_si256(reinterpret_cast<const vect_t *>(p));
	}

	/*
	 * Store 256-bits of integer data from a into memory.
	 * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
	 */
	static INLINE void store(scalar_t *p, vect_t v) {
		_mm256_store_si256(reinterpret_cast<vect_t *>(p), v);
	}

	/*
	 * Store 256-bits of integer data from a into memory.
	 * p does not need to be aligned on any particular boundary.
	 */
	static INLINE void storeu(scalar_t *p, vect_t v) {
		_mm256_storeu_si256(reinterpret_cast<vect_t *>(p), v);
	}

	/*
	 * Store 256-bits of integer data from a into memory using a non-temporal memory hint.
	 * p must be aligned on a 32-byte boundary or a general-protection exception may be generated.
	 */
	static INLINE void stream(scalar_t *p, const vect_t v) {
		_mm256_stream_si256(reinterpret_cast<vect_t *>(p), v);
	}

	/*
	 * Shift packed 64-bit integers in a left by s while shifting in zeros, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 * Return : [a0 << s, a1 << s, a2 << s, a3 << s] int64_t
	 */
	static INLINE CONST vect_t sll(const vect_t a, const int s) { return _mm256_slli_epi64(a, s); }

	/*
	 * Shift packed 64-bit integers in a right by s while shifting in zeros, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 * Return : [a0 >> s, a1 >> s, a2 >> s, a3 >> s] int64_t
	 */
	static INLINE CONST vect_t srl(const vect_t a, const int s) { return _mm256_srli_epi64(a, s); }

	/*
	* Shift packed 64-bit integers in a right by s while shifting in sign bits, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 * Return : [a0 >> s, a1 >> s, a2 >> s, a3 >> s] int64_t
	*/
	static INLINE CONST vect_t sra(const vect_t a, const int s) {
#ifdef __AVX512__
		return _mm256_srai_epi64(a, s);
#else
		const int b = 63 - s;
		vect_t m = sll(set1(1), b);
		vect_t x = srl(a, s);
		vect_t result = sub(vxor(x, m), m); // result = x^m - m
		return result;
#endif
	}

	/*
	 * Add packed 64-bits integer in a and b, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 * Return : [a0+b0, a1+b1, a2+b2, a3+b3]   int64_t
	 */
	static INLINE CONST vect_t add(const vect_t a, const vect_t b) { return _mm256_add_epi64(a, b); }

	static INLINE vect_t addin(vect_t &a, const vect_t b) { return a = add(a, b); }

	/*
	 * Subtract packed 64-bits integers in b from packed 64-bits integers in a, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 * Return : [a0-b0, a1-b1, a2-b2, a3-b3]  int64_t
	 */
	static INLINE CONST vect_t sub(const vect_t a, const vect_t b) { return _mm256_sub_epi64(a, b); }

	static INLINE vect_t subin(vect_t &a, const vect_t b) { return a = sub(a, b); }

	/*
	 * Multiply the packed 64-bits integers in a and b, producing intermediate 128-bit integers, and store the low 64
	 bits of the intermediate integers in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	* Return : [a0*b0 smod 2^32, ..., a3*b3 smod 2^32]	int64_t
	*	   where (a smod p) is the signed representant of a modulo p, that is -p/2 <= (a smod p) < p/2
	 */
	static INLINE CONST vect_t mullo(vect_t a, vect_t b) {
#ifdef __AVX512__
		return _mm256_mullo_epi64(a, b);
#else
#pragma warning "The simd mullo function is emulate, it may impact the performances."
		Converter ca, cb;
		ca.v = a;
		cb.v = b;
		return set(ca.t[0] * cb.t[0], ca.t[1] * cb.t[1], ca.t[2] * cb.t[2], ca.t[3] * cb.t[3]);
#endif
	}

	static INLINE CONST vect_t mul(const vect_t a, const vect_t b) { return mullo(a, b); }

	/*
	 * Multiply the packed 64-bits integers in a and b, producing intermediate 128-bit integers, and store the high 64
	 bits of the intermediate integers in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 * Return : [Floor(a0*b0/2^64), ..., Floor(a3*b3/2^64)] int64_t
	 */
#ifdef __x86_64__
	static INLINE CONST vect_t mulhi(vect_t a, vect_t b) {
#pragma warning "The simd mulhi function is emulate, it may impact the performances."
		// ugly solution, but it works.
		// tested with gcc, clang, icc
		Converter ca, cb;
		ca.v = a;
		cb.v = b;
		return set((scalar_t)((int128_t(ca.t[0]) * cb.t[0]) >> 64), (scalar_t)((int128_t(ca.t[1]) * cb.t[1]) >> 64),
				(scalar_t)((int128_t(ca.t[2]) * cb.t[2]) >> 64), (scalar_t)((int128_t(ca.t[3]) * cb.t[3]) >> 64));
	}
#endif

	/*
	 * Multiply the low 32-bits integers from each packed 64-bit element in a and b, and store the signed 64-bit results
	 in dst.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 * Return : [(a0 smod 2^32)*(b0 smod 2^32), ..., (a3 smod 2^32)*(b3 smod 2^32)]	int64_t
	 *	   where (a smod p) is the signed representant of a modulo p, that is -p/2 <= (a smod p) < p/2
	 */
	static INLINE CONST vect_t mulx(const vect_t a, const vect_t b) { return _mm256_mul_epi32(a, b); }

	/*
	 * Multiply packed 64-bit integers in a and b, producing intermediate 128-bit integers, and add the low 64-bits of
	 the intermediate with c, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 *	    [c0, c1, c2, c3]	int64_t
	 * Return : [(a0*b0+c0) smod 2^64, ..., (a3*b3+c3) smod 2^64]	int64_t
	 */
	static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b) { return add(c, mul(a, b)); }

	static INLINE vect_t fmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fmadd(c, a, b); }

	/*
	 * Multiply the low 32-bit integers from each packed 64-bit element in a and b,
	 * keep the signed 64-bit results and add the low 64-bits of c.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 *	    [c0, c1, c2, c3]	int64_t
	 * Return :	[((a0 smod 2^32)*(b0 smod 2^32)+c0) smod 2^64, ...,
	 *		 ((a3 smod 2^32)*(b3 smod 2^32)+c3) smod 2^64]	int64_t
	 */
	static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b) { return add(c, mulx(a, b)); }

	static INLINE vect_t fmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmaddx(c, a, b); }

	/*
	* Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers,
	* and substract the low 64 bits of the intermediate from elements of c.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 *	    [c0, c1, c2, c3]	int64_t
	 * Return :	[(-a0*b0+c0) smod 2^64, ..., (-a3*b3+c3) smod 2^64]	int64_t
	 */
	static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mul(a, b)); }

	static INLINE vect_t fnmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmadd(c, a, b); }

	/*
	 * Multiply the low 32-bit integers from each packed 64-bit element in a and b,
	 * keep the signed 64-bit results and substract them from elements of c.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 *	    [c0, c1, c2, c3]	int64_t
	 * Return :	[(-(a0 smod 2^32)*(b0 smod 2^32)+c0) smod 2^64, ...,
	 *		 (-(a3 smod 2^32)*(b3 smod 2^32)+c3) smod 2^64]	int64_t
	 */
	static INLINE CONST vect_t fnmaddx(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mulx(a, b)); }

	static INLINE vect_t fnmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmaddx(c, a, b); }

	/*
	 * Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers,
	 * and substract elements of c to the low 64-bits of the intermediate.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 *	    [c0, c1, c2, c3]	int64_t
	 * Return :	[(a0*b0-c0) smod 2^64, ..., (a3*b3-c3) smod 2^64]	int64_t
	 */
	static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b) { return sub(mul(a, b), c); }

	static INLINE vect_t fmsubin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsub(c, a, b); }

	/*
	 * Multiply the low 32-bit integers from each packed 64-bit element in a and b,
	 * keep the signed 64-bit results and substract elements of c from them.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 *	    [c0, c1, c2, c3]	int64_t
	 * Return :	[((a0 smod 2^32)*(b0 smod 2^32)-c0) smod 2^64, ...,
	 *		 ((a3 smod 2^32)*(b3 smod 2^32)-c3) smod 2^64]	int64_t
	 */
	static INLINE CONST vect_t fmsubx(const vect_t c, const vect_t a, const vect_t b) { return sub(mulx(a, b), c); }

	static INLINE vect_t fmsubxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsubx(c, a, b); }

	/*
	 * Compare packed 64-bits in a and b for equality, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 * Return : [(a0==b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1==b1) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a2==b2) ? 0xFFFFFFFFFFFFFFFF : 0, (a3==b3) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
	 */
	static INLINE CONST vect_t eq(const vect_t a, const vect_t b) { return _mm256_cmpeq_epi64(a, b); }

	/*
	 * Compare packed 64-bits in a and b for greater-than, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 * Return : [(a0>b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1>b1) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a2>b2) ? 0xFFFFFFFFFFFFFFFF : 0, (a3>b3) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
	 */
	static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {
		return _mm256_cmpgt_epi64(a, b);
	}

	/*
	 * Compare packed 64-bits in a and b for lesser-than, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 * Return : [(a0<b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1<b1) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a2<b2) ? 0xFFFFFFFFFFFFFFFF : 0, (a3<b3) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
	 */
	static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) { return _mm256_cmpgt_epi64(b, a); }

	/*
	 * Compare packed 64-bits in a and b for greater or equal than, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 * Return : [(a0>=b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1>=b1) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a2>=b2) ? 0xFFFFFFFFFFFFFFFF : 0, (a3>=b3) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a4>=b4) ? 0xFFFFFFFFFFFFFFFF : 0, (a5>=b5) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a6>=b6) ? 0xFFFFFFFFFFFFFFFF : 0, (a7>=b7) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
	 */
	static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

	/*
	 * Compare packed 64-bits in a and b for lesser or equal than, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 * Return : [(a0<=b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1<=b1) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a2<=b2) ? 0xFFFFFFFFFFFFFFFF : 0, (a3<=b3) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
	 */
	static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

	/*
	 * Compute the bitwise AND of packed 64-bits integer in a and b, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 * Return : [a0 AND b0, a1 AND b1, a2 AND b2, a3 AND b3]	int64_t
	 */
	static INLINE CONST vect_t vand(const vect_t a, const vect_t b) { return _mm256_and_si256(b, a); }

	/*
	 * Compute the bitwise OR of packed 64-bits integer in a and b, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 * Return : [a0 OR b0, a1 OR b1, a2 OR b2, a3 OR b3]	int64_t
	 */
	static INLINE CONST vect_t vor(const vect_t a, const vect_t b) { return _mm256_or_si256(b, a); }

	/*
	 * Compute the bitwise XOR of packed 64-bits integer in a and b, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 * Return : [a0 XOR b0, a1 XOR b1, a2 XOR b2, a3 XOR b3]	int64_t
	 */
	static INLINE CONST vect_t vxor(const vect_t a, const vect_t b) { return _mm256_xor_si256(b, a); }

	/*
	 * Compute the bitwise AND NOT of packed 64-bits integer in a and b, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 *	    [b0, b1, b2, b3]	int64_t
	 * Return : [a0 ANDNOT b0, a1 ANDNOT b1, a2 ANDNOT b2, a3 ANDNOT b3]	int64_t
	 */
	static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) { return _mm256_andnot_si256(b, a); }

	/*
	 * Horizontally add 64-bits elements of a.
	 * Args   : [a0, a1, a2, a3]	int64_t
	 * Return : a0+a1+a2+a3	int64_t
	 */
	static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
		Converter ca;
		ca.v = a;
		return scalar_t(ca.t[0] + ca.t[1] + ca.t[2] + ca.t[3]);
	}

	static INLINE CONST vect_t round(const vect_t a) { return a; }

	static INLINE CONST vect_t signbits(const vect_t x) {
		vect_t signBits = sub(zero(), srl(x, 4*sizeof(scalar_t)-1));
		return signBits;
	}

	// mask the high 32 bits of a 64 bits, that is 00000000FFFFFFFF
	static INLINE CONST vect_t mask_high() { return srl(_mm256_set1_epi8(-1), 32); }

	static INLINE CONST vect_t mulhi_fast(vect_t x, vect_t y);

	template <bool overflow, bool poweroftwo>
	static INLINE vect_t mod(vect_t &C, const vect_t &P, const int8_t &shifter, const vect_t &magic, const vect_t &NEGP,
				 const vect_t &MIN, const vect_t &MAX, vect_t &Q, vect_t &T);
}; // Simd256_impl<true, true, true, 8>

/*
 * Simd256 specialized for uint64_t
 */
template <> struct Simd256_impl<true, true, false, 8> : public Simd256_impl<true, true, true, 8> {

	/*
	* define the scalar type corresponding to the specialization
	*/
	using scalar_t = uint64_t;

	/*
	 * Simd128 for scalar_t, to deal half_t
	 */
	using simdHalf = Simd128<scalar_t>;

	/*
	* Converter from vect_t to a tab.
	* exple:
	*	Converter conv;
	*	conv.v = a;
	*	scalart_t x = conv.t[1]
	*/
	union Converter {
		vect_t v;
		scalar_t t[vect_size];
	};

	/*
	 *  Broadcast 64-bit unsigned integer a to all elements of dst. This intrinsic may generate the vpbroadcastw.
	 *  Return [x,x,x,x] uint64_t
	 */
	static INLINE CONST vect_t set1(const scalar_t x) { return _mm256_set1_epi64x(x); }

	/*
	 *  Set packed 64-bit unsigned integers in dst with the supplied values.
	 *  Return [x0,x1,x2,x3] uint64_t
	 */
	static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1, const scalar_t x2, const scalar_t x3) {
		return _mm256_set_epi64x(x3, x2, x1, x0);
	}

	/*
	 *  Gather 64-bit unsigned integer elements with indexes idx[0], ..., idx[3] from the address p in vect_t.
	 *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]]] uint64_t
	 */
	template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
		return set(p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]]);
	}

	/*
	 * Load 256-bits of unsigned integer data from memory into dst.
	 * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
	 * Return [p[0],p[1],p[2],p[3]] uint64_t
	 */
	static INLINE PURE vect_t load(const scalar_t *const p) {
		return _mm256_load_si256(reinterpret_cast<const vect_t *>(p));
	}

	/*
	 * Load 256-bits of unsigned integer data from memory into dst.
	 * p does not need to be aligned on any particular boundary.
	 * Return [p[0],p[1],p[2],p[3]] uint64_t
	 */
	static INLINE PURE vect_t loadu(const scalar_t *const p) {
		return _mm256_loadu_si256(reinterpret_cast<const vect_t *>(p));
	}

	/*
	 * Store 256-bits of unsigned integer data from a into memory.
	 * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
	 */
	static INLINE void store(scalar_t *p, vect_t v) {
		_mm256_store_si256(reinterpret_cast<vect_t *>(p), v);
	}

	/*
	 * Store 256-bits of unsigned integer data from a into memory.
	 * p does not need to be aligned on any particular boundary.
	 */
	static INLINE void storeu(scalar_t *p, vect_t v) {
		_mm256_storeu_si256(reinterpret_cast<vect_t *>(p), v);
	}

	/*
	 * Store 256-bits of unsigned integer data from a into memory using a non-temporal memory hint.
	 * p must be aligned on a 32-byte boundary or a general-protection exception may be generated.
	 */
	static INLINE void stream(scalar_t *p, const vect_t v) {
		_mm256_stream_si256(reinterpret_cast<vect_t *>(p), v);
	}

	/*
	* Shift packed 64-bit unsigned integers in a right by s while shifting in sign bits, and store the results in vect_t.
	 * Args   : [a0, ..., a3]			uint64_t
	 * Return : [Floor(a0/2^s), ..., Floor(a3/2^s)]	uint64_t
	*/
	static INLINE CONST vect_t sra(const vect_t a, const int s) { return _mm256_srli_epi64(a, s); }

	static INLINE CONST vect_t greater(vect_t a, vect_t b) {
		vect_t x;
		x = set1(-(static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1)));
		a = sub(x, a);
		b = sub(x, b);
		return _mm256_cmpgt_epi64(b,a);
	}

	static INLINE CONST vect_t lesser(vect_t a, vect_t b) {
		vect_t x;
		x = set1(-(static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1)));
		a = sub(x, a);
		b = sub(x, b);
		return _mm256_cmpgt_epi64(a, b);
	}

	static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

	static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

	/*
	 * Multiply the packed 64-bits integers in a and b, producing intermediate 128-bit integers, and store the low 64
	 bits of the intermediate integers in vect_t.
	 * Args   : [a0, a1, a2, a3]	   							uint64_t
	 [b0, b1, b2, b3]  		 							uint64_t
	 * Return : [a0*b0 mod 2^64, a1*b1 mod 2^64, a2*b2 mod 2^64, a3*b3 mod 2^64]		uint64_t
	 */
	static INLINE CONST vect_t mullo(vect_t a, vect_t b) {
#pragma warning "The simd mullo function is emulate, it may impact the performances."
		Converter ca, cb;
		ca.v = a;
		cb.v = b;
		return set(ca.t[0] * cb.t[0], ca.t[1] * cb.t[1], ca.t[2] * cb.t[2], ca.t[3] * cb.t[3]);
	}

	/*
	 * Multiply the packed 64-bits integers in a and b, producing intermediate 128-bit integers, and store the high 64
	 bits of the intermediate integers in vect_t.
	 * Args   : [a0, a1, a2, a3]	   							 uint64_t
	 [b0, b1, b2, b3]  		 							 uint64_t
	 * Return :
	 */
#ifdef __x86_64__
	static INLINE CONST vect_t mulhi(vect_t a, vect_t b) {
#pragma warning "The simd mulhi function is emulate, it may impact the performances."
		// ugly solution, but it works.
		// tested with gcc, clang, icc
		Converter c0, c1;
		c0.v = a;
		c1.v = b;
		return set((scalar_t)(((uint128_t)(c0.t[0]) * c1.t[0]) >> 64), (scalar_t)(((uint128_t)(c0.t[1]) * c1.t[1]) >> 64),
			   (scalar_t)(((uint128_t)(c0.t[2]) * c1.t[2]) >> 64), (scalar_t)(((uint128_t)(c0.t[3]) * c1.t[3]) >> 64));
	}
#endif

	/*
	 * Multiply the low 32-bits integers from each packed 64-bit element in a and b, and store the unsigned 64-bit
	 results in dst.
	 * Args   : [a0, a1, a2, a3]	uint64_t
	 [b0, b1, b2, b3]	uint64_t
	 * Return : [a0*b0, a1*b1, a2*b2, a3*b3] uint64_t
	 */
	static INLINE CONST vect_t mulx(const vect_t a, const vect_t b) { return _mm256_mul_epu32(a, b); }

	static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b) { return add(c, mulx(a, b)); }

	static INLINE vect_t fmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmaddx(c, a, b); }

	static INLINE CONST vect_t fnmaddx(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mulx(a, b)); }

	static INLINE vect_t fnmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmaddx(c, a, b); }

	static INLINE CONST vect_t fmsubx(const vect_t c, const vect_t a, const vect_t b) { return sub(mulx(a, b), c); }

	static INLINE vect_t fmsubxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsubx(c, a, b); }

	/*
	 * Horizontally add 64-bits elements of a.
	 * Args   : [a0, a1, a2, a3]
	 * Return : a0+a1+a2+a3
	 */
	static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
		Converter ca;
		ca.v = a;
		return ca.t[0] + ca.t[1] + ca.t[2] + ca.t[3];
	}
}; // Simd256_impl<true, true, false, 8>

#define vect_t Simd256_impl<true, true, true, 8>::vect_t

// warning : may be off by 1 multiple, but we save a mul...
INLINE CONST vect_t Simd256_impl<true, true, true, 8>::mulhi_fast(vect_t x, vect_t y) {
	// unsigned mulhi starts:
	// x1 = xy_high = mulhiu_fast(x,y)
	const vect_t mask = mask_high();

	vect_t x0 = vand(x, mask), x1 = srl(x, 32);
	vect_t y0 = vand(y, mask), y1 = srl(y, 32);

	x0 = Simd256_impl<true, true, false, 8>::mulx(x0, y1); // x0y1
	y0 = Simd256_impl<true, true, false, 8>::mulx(x1, y0); // x1y0
	y1 = Simd256_impl<true, true, false, 8>::mulx(x1, y1); // x1y1

	x1 = vand(y0, mask);
	y0 = srl(y0, 32); // x1y0_lo = x1 // y1yo_hi = y0
	x1 = srl(add(x1, x0), 32);
	y0 = add(y1, y0);

	x1 = add(x1, y0);
	// unsigned mulhi ends

	// fixing signs
	x0 = vand(signbits(x), y);
	x1 = sub(x1, x0);
	x0 = vand(signbits(y), x);
	x1 = sub(x1, x0);
	// end fixing
	return x1;
}

template <bool overflow, bool poweroftwo>
INLINE vect_t Simd256_impl<true, true, true, 8>::mod(vect_t &C, const vect_t &P, const int8_t &shifter, const vect_t &magic, const vect_t &NEGP,
						     const vect_t &MIN, const vect_t &MAX, vect_t &Q, vect_t &T) {
#ifdef __INTEL_COMPILER
	// Works fine with ICC 15.0.1 - A.B.
	C = _mm256_rem_epi64(C, P);
#else
	if (poweroftwo) {
		Q = srl(C, 63);
		vect_t un = set1(1);
		T = sub(sll(un, shifter), un);
		Q = add(C, vand(Q, T));
		Q = sll(srl(Q, shifter), shifter);
		C = sub(C, Q);
		Q = vand(greater(zero(), Q), P);
		C = add(C, Q);
	} else {
		Q = mulhi_fast(C, magic);
		if (overflow) {
			Q = add(Q, C);
		}
		Q = sra(Q, shifter);
		vect_t q1 = Simd256_impl<true, true, false, 8>::mulx(Q, P);
		vect_t q2 = sll(Simd256_impl<true, true, false, 8>::mulx(srl(Q, 32), P), 32);
		C = sub(C, add(q1, q2));
		T = greater_eq(C, P);
		C = sub(C, vand(T, P));
	}
#endif
	NORML_MOD(C, P, NEGP, MIN, MAX, Q, T);
	return C;
}

#undef vect_t

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd256_int64_INL
