/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd128_int64_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd128_int64_INL

#ifndef __FFLASFFPACK_USE_SIMD
#error "You need SSE instructions to perform 128 bits operations on int64"
#endif

/*
 * Simd128 specialized for int64_t
 */
template <> struct Simd128_impl<true, true, true, 8> : public Simd128i_base {

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
	*  Broadcast 64-bit integer a to all elements of dst. This intrinsic may generate the vpbroadcastw.
	*  Return [x,x] int64_t
	*/
	static INLINE CONST vect_t set1(const scalar_t x) { return _mm_set1_epi64x(x); }

	/*
	*  Set packed 64-bit integers in dst with the supplied values.
	*  Return [x0,x1] int64_t
	*/
	static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1) { return _mm_set_epi64x(x1, x0); }

	/*
	*  Gather 64-bit integer elements with indexes idx[0], idx[1] from the address p in vect_t.
	*  Return [p[idx[0]], p[idx[1]]] int64_t
	*/
	template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
		return set(p[idx[0]], p[idx[1]]);
	}

	/*
	* Load 128-bits of integer data from memory into dst.
	* p must be aligned on a 16-byte boundary or a general-protection exception will be generated.
	* Return [p[0],p[1]] int64_t
	*/
	static INLINE PURE vect_t load(const scalar_t *const p) {
		return _mm_load_si128(reinterpret_cast<const vect_t *>(p));
	}

	/*
	* Load 128-bits of integer data from memory into dst.
	* p does not need to be aligned on any particular boundary.
	* Return [p[0],p[1]] int64_t
	*/
	static INLINE PURE vect_t loadu(const scalar_t *const p) {
		return _mm_loadu_si128(reinterpret_cast<const vect_t *>(p));
	}

	/*
	* Store 128-bits of integer data from a into memory.
	* p must be aligned on a 16-byte boundary or a general-protection exception will be generated.
	*/
	static INLINE void store(scalar_t *p, vect_t v) {
		_mm_store_si128(reinterpret_cast<vect_t *>(p), v);
	}

	/*
	* Store 128-bits of integer data from a into memory.
	* p does not need to be aligned on any particular boundary.
	*/
	static INLINE void storeu(scalar_t *p, vect_t v) {
		_mm_storeu_si128(reinterpret_cast<vect_t *>(p), v);
	}

	/*
	* Store 128-bits of integer data from a into memory using a non-temporal memory hint.
	* p must be aligned on a 16-byte boundary or a general-protection exception may be generated.
	*/
	static INLINE void stream(scalar_t *p, const vect_t v) {
		_mm_stream_si128(reinterpret_cast<vect_t *>(p), v);
	}

	/*
	* Shift packed 64-bit integers in a left by s while shifting in zeros, and store the results in vect_t.
	* Args   : [a0, a1] int64_t
	* Return : [a0 << s, a1 << s] int64_t
	*/
	static INLINE CONST vect_t sll(const vect_t a, const int s) { return _mm_slli_epi64(a, s); }

	/*
	* Shift packed 64-bit integers in a right by s while shifting in zeros, and store the results in vect_t.
	* Args   : [a0, a1] int64_t
	* Return : [a0 >> s, a1 >> s] int64_t
	*/
	static INLINE CONST vect_t srl(const vect_t a, const int s) { return _mm_srli_epi64(a, s); }

	/*
	* Shift packed 64-bit integers in a right by s while shifting in sign bits, and store the results in vect_t.
	* Args   : [a0, a1] int64_t
	* Return : [a0 >> s, a1 >> s] int64_t
	*/
	static INLINE CONST vect_t sra(const vect_t a, const int s) {
#ifdef __AVX512__
		return _mm_srai_epi64(a, s);
#else
		const int b = 63 - s;
		vect_t m = sll(set1(1), b);
		vect_t x = srl(a, s);
		vect_t result = sub(vxor(x, m), m); // result = x^m - m
		return result;
#endif // 512
	}

	/*
	* Shuffle 64-bit integers in a using the control in imm8, and store the results in dst.
	* Args   : [a0, a1] int64_t
	* Return : [a[s[0]], a[s[1]]] int64_t
	*/
	template<uint8_t s>
	static INLINE CONST vect_t shuffle(const vect_t a) {
		// Transform s = [d1 d0]_base2 to s1 = [2*d1+1 2*d1 2*d0+1 2*d0]_base4
		constexpr uint8_t s1 = ((s & 1)?(3*4+2):(1*4+0))+16*((s & 2)?(3*4+2):(1*4+0));
		return _mm_shuffle_epi32(a, s1);
	}

	/*
	* Unpack and interleave 64-bit integers from the low half of a and b, and store the results in dst.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
	* Return : [a0, b0] int64_t
	*/
	static INLINE CONST vect_t unpacklo(const vect_t a, const vect_t b) { return _mm_unpacklo_epi64(a, b); }

	/*
	* Unpack and interleave 64-bit integers from the high half of a and b, and store the results in dst.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
	* Return : [a1, b1] int64_t
	*/
	static INLINE CONST vect_t unpackhi(const vect_t a, const vect_t b) { return _mm_unpackhi_epi64(a, b); }

	/*
	* Blend packed 64-bit integers from a and b using control mask imm8, and store the results in dst.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
	* Return : [s[0]?a0:b0, s[1]?a1:b1] int64_t
	*/
	template<uint8_t s>
	static INLINE CONST vect_t blend(const vect_t a, const vect_t b) {
		// _mm_blend_epi16 is faster than _mm_blend_epi32 and require SSE4.1 instead of AVX2
		// We have to transform s = [d1 d0]_base2 to s1 = [d1 d1 d1 d1 d0 d0 d0 d0]_base2
		constexpr uint8_t s1 = (s & 0x1) * 15 + ((s & 0x2) << 3) * 15;
		return _mm_blend_epi16(a, b, s1);
	}

	/*
	* Add packed 64-bits integer in a and b, and store the results in vect_t.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
	* Return : [a0+b0, a1+b1]   int64_t
	*/
	static INLINE CONST vect_t add(const vect_t a, const vect_t b) { return _mm_add_epi64(a, b); }

	static INLINE vect_t addin(vect_t &a, const vect_t b) { return a = add(a, b); }

	/*
	* Subtract packed 64-bit integers in b from packed 64-bit integers in a, and store the results in vect_t.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
	* Return : [a0-b0, a1-b1]  int64_t
	*/
	static INLINE CONST vect_t sub(const vect_t a, const vect_t b) { return _mm_sub_epi64(a, b); }

	static INLINE vect_t subin(vect_t &a, const vect_t b) { return a = sub(a, b); }

	/*
	* Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers, and store the low 64
	bits of the intermediate integers in vect_t.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
	* Return : [a0*b0 smod 2^64, a1*b1 smod 2^64] int64_t
	*	   where (a smod p) is the signed representant of a modulo p, that is -p/2 <= (a smod p) < p/2
	*/
	static INLINE CONST vect_t mullo(const vect_t x0, const vect_t x1) {
#ifdef __AVX512__
		_mm_mullo_epi64(x0, x1);
#else
		// _mm_mullo_epi64 emul
		//#pragma warning "The simd mullo function is emulate, it may impact the performances."
		Converter c0, c1;
		c0.v = x0;
		c1.v = x1;
		return set((scalar_t)(c0.t[0] * c1.t[0]), (scalar_t)(c0.t[1] * c1.t[1]));
#endif // 512
	}

	static INLINE CONST vect_t mul(const vect_t a, const vect_t b) { return mullo(a, b); }

	/*
	* Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers, and store the high 64
	bits of the intermediate integers in vect_t.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
	* Return : [Floor(a0*b0/2^64), Floor(a1*b1/2^64)] int64_t
	*/
#ifdef __FFLASFFPACK_HAVE_INT128
	static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b) {
		//#pragma warning "The simd mulhi function is emulated, it may impact the performances."
		Converter c0, c1;
		c0.v = a;
		c1.v = b;
		return set((scalar_t)((int128_t(c0.t[0]) * c1.t[0]) >> 64), (scalar_t)((int128_t(c0.t[1]) * c1.t[1]) >> 64));
	}
#endif

	/*
	* Multiply the low 32-bit integers from each packed 64-bit element in a and b, and store the signed 64-bit results
	in vect_t.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
	* Return : [(a0 smod 2^32)*(b0 smod 2^32), (a1 smod 2^32)*(b1 smod 2^32)]	int64_t
	*	   where (a smod p) is the signed representant of a modulo p, that is -p/2 <= (a smod p) < p/2
	*/
	static INLINE CONST vect_t mulx(const vect_t a, const vect_t b) { return _mm_mul_epi32(a, b); }

	/*
	* Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers,
	* keep the low 64 bits of the intermediate and add the low 64-bits of c.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
			   [c0, c1] int64_t
	* Return : [(a0*b0+c0) smod 2^64, (a1*b1+c1) smod 2^64]	int64_t
	*/
	static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b) { return add(c, mul(a, b)); }

	static INLINE vect_t fmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fmadd(c, a, b); }

	/*
	* Multiply the low 32-bit integers from each packed 64-bit element in a and b,
	* keep the signed 64-bit results and add the low 64-bits of c.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
			   [c0, c1] int64_t
	* Return : [((a0 smod 2^32)*(b0 smod 2^32)+c0) smod 2^64,
	*		 ((a1 smod 2^32)*(b1 smod 2^32)+c1) smod 2^64]	int64_t
	*/
	static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b) { return add(c, mulx(a, b)); }

	static INLINE vect_t fmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmaddx(c, a, b); }

	/*
	* Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers,
	* and substract the low 64 bits of the intermediate from elements of c.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
			   [c0, c1] int64_t
	* Return : [(-a0*b0+c0) smod 2^64, (-a1*b1+c1) smod 2^64]	int64_t
	*/
	static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mul(a, b)); }

	static INLINE vect_t fnmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmadd(c, a, b); }

	/*
	* Multiply the low 32-bit integers from each packed 64-bit element in a and b,
	* keep the signed 64-bit results and substract them from elements of c.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
			   [c0, c1] int64_t
	* Return : [(-(a0 smod 2^32)*(b0 smod 2^32)+c0) smod 2^64,
	*		 (-(a1 smod 2^32)*(b1 smod 2^32)+c1) smod 2^64]	int64_t
	*/
	static INLINE CONST vect_t fnmaddx(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mulx(a, b)); }

	static INLINE vect_t fnmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmaddx(c, a, b); }

	/*
	* Multiply the packed 64-bit integers in a and b, producing intermediate 128-bit integers,
	* and substract elements of c to the low 64-bits of the intermediate.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
			   [c0, c1] int64_t
	* Return : [(a0*b0-c0) smod 2^64, (a1*b1-c1) smod 2^64]	int64_t
	*/
	static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b) { return sub(mul(a, b), c); }

	static INLINE vect_t fmsubin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsub(c, a, b); }

	/*
	* Multiply the low 32-bit integers from each packed 64-bit element in a and b,
	* keep the signed 64-bit results and substract elements of c from them.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
			   [c0, c1] int64_t
	* Return : [(-(a0 smod 2^32)*(b0 smod 2^32)+c0) smod 2^64,
	*		 (-(a1 smod 2^32)*(b1 smod 2^32)+c1) smod 2^64]	int64_t
	*/
	static INLINE CONST vect_t fmsubx(const vect_t c, const vect_t a, const vect_t b) { return sub(mulx(a, b), c); }

	static INLINE vect_t fmsubxin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsubx(c, a, b); }

	/*
	* Compare packed 64-bits in a and b for equality, and store the results in vect_t.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
	* Return : [(a0==b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1==b1) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
	*/
	static INLINE CONST vect_t eq(const vect_t a, const vect_t b) { return _mm_cmpeq_epi64(a, b); }

	/*
	* Compare packed 64-bits in a and b for greater-than, and store the results in vect_t.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
	* Return : [(a0>b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1>b1) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
	*/
	static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {
#ifdef __SSE4_2__
		return _mm_cmpgt_epi64(a, b);
#else
		//#warning "The simd greater function is emulate, it may impact the performances."
		Converter ca, cb;
		ca.v = a;
		cb.v = b;
		return set((ca.t[0] > cb.t[0]) ? 0xFFFFFFFFFFFFFFFF : 0, (ca.t[1] > cb.t[1]) ? 0xFFFFFFFFFFFFFFFF : 0);
#endif // __SSE4_2__
	}

	/*
	* Compare packed 64-bits in a and b for lesser-than, and store the results in vect_t.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
	* Return : [(a0<b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1<b1) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
	*/
	static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {
#ifdef __SSE4_2__
		return _mm_cmpgt_epi64(b, a);
#else
		//#warning "The simd lesser function is emulate, it may impact the performances."
		Converter ca, cb;
		ca.v = a;
		cb.v = b;
		return set((ca.t[0] < cb.t[0]) ? 0xFFFFFFFFFFFFFFFF : 0, (ca.t[1] < cb.t[1]) ? 0xFFFFFFFFFFFFFFFF : 0);
#endif // __SSE4_2__
	}

	/*
	* Compare packed 64-bits in a and b for greater or equal than, and store the results in vect_t.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
	* Return : [(a0>=b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1>=b1) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
	*/
	static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

	/*
	* Compare packed 64-bits in a and b for lesser or equal than, and store the results in vect_t.
	* Args   : [a0, a1] int64_t
			   [b0, b1] int64_t
	* Return : [(a0<=b0) ? 0xFFFFFFFFFFFFFFFF : 0, (a1<=b1) ? 0xFFFFFFFFFFFFFFFF : 0]	int64_t
	*/
	static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

	/*
	* Horizontally add 64-bits elements of a.
	* Args   : [a0, a1]	int64_t
	* Return : a0+a1	int64_t
	*/
	static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
		Converter conv;
		conv.v = a;
		return scalar_t(conv.t[0] + conv.t[1]);
	}

	static INLINE CONST vect_t round(const vect_t a) { return a; }

	static INLINE CONST vect_t signbits(const vect_t x) {
		vect_t signBits = sub(zero(), srl(x, 4*sizeof(scalar_t)-1));
		return signBits;
	}

	// mask the high 32 bits of a 64 bits, that is 00000000FFFFFFFF
	static INLINE CONST vect_t mask_high() { return srl(_mm_set1_epi8(-1), 32); }

	static INLINE CONST vect_t mulhi_fast(vect_t x, vect_t y);

	template <bool overflow, bool poweroftwo>
	static INLINE vect_t mod(vect_t &C, const vect_t &P, const int8_t &shifter, const vect_t &magic, const vect_t &NEGP,
							 const vect_t &MIN, const vect_t &MAX, vect_t &Q, vect_t &T);
}; // Simd128_impl<true, true, true, 8>

/*
 * Simd128 specialized for uint64_t
 */
template <> struct Simd128_impl<true, true, false, 8> : public Simd128_impl<true, true, true, 8> {

	/*
	* define the scalar type corresponding to the specialization
	*/
	using scalar_t = uint64_t;

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
	 *  Return [x,x] uint64_t
	 */
	static INLINE CONST vect_t set1(const scalar_t x) { return _mm_set1_epi64x(x); }

	/*
	 *  Set packed 64-bit integers in dst with the supplied values.
	 *  Return [x0,x1] uint64_t
	 */
	static INLINE CONST vect_t set(const scalar_t x0, const scalar_t x1) { return _mm_set_epi64x(x1, x0); }

	/*
	 *  Gather 64-bit unsigned integer elements with indexes idx[0], ..., idx[1] from the address p in vect_t.
	 *  Return [p[idx[0]], p[idx[1]]] uint64_t
	 */
	template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
		return set(p[idx[0]], p[idx[1]]);
	}

	/*
	 * Load 128-bits of unsigned integer data from memory into dst.
	 * p must be aligned on a 16-byte boundary or a general-protection exception will be generated.
	 * Return [p[0],p[1]] uint64_t
	 */
	static INLINE PURE vect_t load(const scalar_t *const p) {
		return _mm_load_si128(reinterpret_cast<const vect_t *>(p));
	}

	/*
	 * Load 128-bits of unsigned integer data from memory into dst.
	 * p does not need to be aligned on any particular boundary.
	 * Return [p[0],p[1]] uint64_t
	 */
	static INLINE PURE vect_t loadu(const scalar_t *const p) {
		return _mm_loadu_si128(reinterpret_cast<const vect_t *>(p));
	}

	/*
	 * Store 128-bits of unsigned integer data from a into memory.
	 * p must be aligned on a 32-byte boundary or a general-protection exception will be generated.
	 */
	static INLINE void store(scalar_t *p, vect_t v) {
		_mm_store_si128(reinterpret_cast<vect_t *>(p), v);
	}

	/*
	 * Store 128-bits of unsigned integer data from a into memory.
	 * p does not need to be aligned on any particular boundary.
	 */
	static INLINE void storeu(scalar_t *p, vect_t v) {
		_mm_storeu_si128(reinterpret_cast<vect_t *>(p), v);
	}

	/*
	 * Store 128-bits of unsigned integer data from a into memory using a non-temporal memory hint.
	 * p must be aligned on a 16-byte boundary or a general-protection exception may be generated.
	 */
	static INLINE void stream(scalar_t *p, const vect_t v) {
		_mm_stream_si128(reinterpret_cast<vect_t *>(p), v);
	}

	/*
	* Shift packed 64-bit unsigned integers in a right by s while shifting in sign bits, and store the results in vect_t.
	 * Args   : [a0, a1]				uint64_t
	 * Return : [Floor(a0/2^s), Floor(a1/2^s)]	uint64_t
	*/
	static INLINE CONST vect_t sra(const vect_t a, const int s) { return _mm_srli_epi64(a, s); }

	static INLINE CONST vect_t greater(vect_t a, vect_t b) {
#ifdef __SSE4_2__
		vect_t x;
		x = set1(-(static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1)));
		a = sub(x, a);
		b = sub(x, b);
		return _mm_cmpgt_epi64(b, a);
#else
		//#pragma warning "The simd greater function is emulated, it may impact the performances."
		Converter ca, cb;
		ca.v = a;
		cb.v = b;
		return set((ca.t[0] > cb.t[0]) ? 0xFFFFFFFFFFFFFFFF : 0, (ca.t[1] > cb.t[1]) ? 0xFFFFFFFFFFFFFFFF : 0);
#endif
	}

	static INLINE CONST vect_t lesser(vect_t a, vect_t b) {
#ifdef __SSE4_2__
		vect_t x;
		x = set1(-(static_cast<scalar_t>(1) << (sizeof(scalar_t) * 8 - 1)));
		a = sub(x, a);
		b = sub(x, b);
		return _mm_cmpgt_epi64(a, b);
#else
		//#pragma warning "The simd greater function is emulated, it may impact the performances."
		Converter ca, cb;
		ca.v = a;
		cb.v = b;
		return set((ca.t[0] < cb.t[0]) ? 0xFFFFFFFFFFFFFFFF : 0, (ca.t[1] < cb.t[1]) ? 0xFFFFFFFFFFFFFFFF : 0);
#endif
	}

	static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) { return vor(greater(a, b), eq(a, b)); }

	static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) { return vor(lesser(a, b), eq(a, b)); }

	/*
	* Multiply the packed 64-bit unsigned integers in a and b, producing intermediate 128-bit integers, and store the low 64
	bits of the intermediate integers in vect_t.
	* Args   : [a0, a1] uint64_t
			   [b0, b1] uint64_t
	* Return : [a0*b0 mod 2^64, a1*b1 mod 2^64] uint64_t
	*/
	static INLINE CONST vect_t mullo(const vect_t x0, const vect_t x1) {
		// _mm_mullo_epi32 emul
		//#pragma warning "The simd mullo function is emulated, it may impact the performances."
		Converter c0, c1;
		c0.v = x0;
		c1.v = x1;
		return set((scalar_t)(c0.t[0] * c1.t[0]), (scalar_t)(c0.t[1] * c1.t[1]));
	}

	/*
	* Multiply the packed unsigned 64-bit integers in a and b, producing intermediate 128-bit integers,
	* and store the high 64 bits of the intermediate integers in vect_t.
	* Args   : [a0, a1] uint64_t
			   [b0, b1] uint64_t
	* Return : [Floor(a0*b0/2^16), Floor(a1*b1/2^16)] uint64_t
	*/
#ifdef __FFLASFFPACK_HAVE_INT128
	static INLINE CONST vect_t mulhi(const vect_t a, const vect_t b) {
		//#pragma warning "The simd mulhi function is emulate, it may impact the performances."
		Converter c0, c1;
		c0.v = a;
		c1.v = b;
		return set((scalar_t)((uint128_t(c0.t[0]) * c1.t[0]) >> 64), (scalar_t)((uint128_t(c0.t[1]) * c1.t[1]) >> 64));
	}
#endif

	/*
	* Multiply the low unsigned 32-bit integers from each packed 64-bit element in a and b, and store the signed 64-bit results
	in vect_t.
	* Args   : [a0, a1] uint64_t
			   [b0, b1] uint64_t
	* Return : [(a0 mod 2^32)*(b0 mod 2^32), (a1 mod 2^32)*(b1 mod 2^32)]	uint64_t
	*/
	static INLINE CONST vect_t mulx(const vect_t a, const vect_t b) { return _mm_mul_epu32(a, b); }

	static INLINE CONST vect_t fmaddx(const vect_t c, const vect_t a, const vect_t b) { return add(c, mulx(a, b)); }

	static INLINE CONST vect_t fnmaddx(const vect_t c, const vect_t a, const vect_t b) { return sub(c, mulx(a, b)); }

	static INLINE vect_t fnmaddxin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmaddx(c, a, b); }

	static INLINE CONST vect_t fmsubx(const vect_t c, const vect_t a, const vect_t b) { return sub(mulx(a, b), c); }

	static INLINE CONST vect_t fmsubxin(vect_t c, const vect_t a, const vect_t b) { return c = fmsubx(c, a, b); }

	/*
	* Horizontally add 64-bits elements of a.
	* Args   : [a0, a1, a2, a3]
	* Return : a0+a1+a2+a3
	*/
	static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
		Converter c;
		c.v = a;
		return c.t[0] + c.t[1];
	}
}; //Simd128_impl<true,true,false,8>

#define vect_t Simd128_impl<true,true,true,8>::vect_t

// warning : may be off by 1 multiple, but we save a mul...
INLINE CONST vect_t Simd128_impl<true,true,true,8>::mulhi_fast(vect_t x, vect_t y) {
	// unsigned mulhi starts:
	// x1 = xy_high = mulhiu_fast(x,y)
	const vect_t mask = mask_high();

	vect_t x0 = vand(x, mask), x1 = srl(x, 32);
	vect_t y0 = vand(y, mask), y1 = srl(y, 32);

	x0 = Simd128_impl<true, true, false, 8>::mulx(x0, y1); // x0y1
	y0 = Simd128_impl<true, true, false, 8>::mulx(x1, y0); // x1y0
	y1 = Simd128_impl<true, true, false, 8>::mulx(x1, y1); // x1y1

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

// warning : may be off by 1 multiple, but we save a mul...
template <bool overflow, bool poweroftwo>
INLINE CONST vect_t Simd128_impl<true,true,true,8>::mod(vect_t &C, const vect_t &P, const int8_t &shifter, const vect_t &magic, const vect_t &NEGP,
														const vect_t &MIN, const vect_t &MAX, vect_t &Q, vect_t &T) {
#ifdef __INTEL_COMPILER
	// Works fine with ICC 15.0.1 - A.B.
	// #warning "not tested"
	C = _mm_rem_epi64(C, P);
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
		vect_t q1 = Simd128_impl<true, true, false, 8>::mulx(Q, P);
		vect_t q2 = sll(Simd128_impl<true, true, false, 8>::mulx(srl(Q, 32), P), 32);
		C = sub(C, add(q1, q2));
		T = greater_eq(C, P);
		C = sub(C, vand(T, P));
	}
#endif
	NORML_MOD(C, P, NEGP, MIN, MAX, Q, T);
	return C;
}

#undef vect_t

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd128_int64_INL
