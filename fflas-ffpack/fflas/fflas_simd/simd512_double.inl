/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2018 the FFLAS-FFPACK group
 *
 * Written by   Ozturk Ozden<ozden.ozturk@etu.univ-grenoble-alpes.fr>
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

#ifndef __FFLASFFPACK_simd512_double_INL
#define __FFLASFFPACK_simd512_double_INL

#if not (defined(__FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS))
#error "You need AVX512 instructions to perform 512bits operations on double"
#endif

/*
 * Simd512 specialized for double
 */
template <> struct Simd512_impl<true, false, true, 8> : public Simd512fp_base {
	/*
	 * alias to 512 bit simd register
	 */
	using vect_t = __m512d;

	/*
	 * define the scalar type corresponding to the specialization
	 */
	using scalar_t = double;

	/*
	 *	number of scalar_t in a simd register
	 */
	static const constexpr size_t vect_size = 8;

	/*
	 *	alignement required by scalar_t pointer to be loaded in a vect_t
	 */
	static const constexpr size_t alignment = 64;

	/*
	 * Check if the pointer p is a multiple of alignemnt
	 */
	template <class T> static constexpr bool valid(T *p) { return (int64_t)p % alignment == 0; }

	/*
	 * Check if the number n is a multiple of vect_size
	 */
	template <class T> static constexpr bool compliant(T n) { return n % vect_size == 0; }

	/*
	 *	Return vector of type vect_t with all elements set to zero
	 *  Return [0,0,0,0,0,0,0,0]
	 */
	static INLINE CONST vect_t zero() { return _mm512_setzero_pd(); }

	/*
	 *	Broadcast double-precision (64-bit) floating-point value x to all elements of vect_t.
	 *  Return [x,x,x,x,x,x,x,x]
	 */
	static INLINE CONST vect_t set1(const scalar_t x) { return _mm512_set1_pd(x); }

	/*
	 *	Set packed double-precision (64-bit) floating-point elements in vect_t with the supplied values.
	 *  Return [x1,x2,x3,x4,x5,x6,x7,x8]
	 */
	static INLINE CONST vect_t set(const scalar_t x1, const scalar_t x2, const scalar_t x3, const scalar_t x4, const scalar_t x5, const scalar_t x6, const scalar_t x7, const scalar_t x8) {
		return _mm512_set_pd(x8, x7, x6, x5, x4, x3, x2, x1);
	}

	/*
	 *	Gather double-precision (64-bit) floating-point elements with indexes idx[0], ..., idx[7] from the address p in
	 *vect_t.
	 *  Return [p[idx[0]], p[idx[1]], p[idx[2]], p[idx[3]], p[idx[4]], p[idx[5]], p[idx[6]], p[idx[7]]]
	 */
	template <class T> static INLINE PURE vect_t gather(const scalar_t *const p, const T *const idx) {
		// TODO AVX2 Gather
		return _mm512_set_pd(p[idx[7]], p[idx[6]], p[idx[5]], p[idx[4]] ,p[idx[3]], p[idx[2]], p[idx[1]], p[idx[0]]);
	}

	/*
	 * Load 512-bits (composed of 8 packed double-precision (64-bit) floating-point elements) from memory into vect_t.
	 * p must be aligned on a 64-byte boundary or a general-protection exception will be generated.
	 * Return [p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]]
	 */
	static INLINE PURE vect_t load(const scalar_t *const p) { return _mm512_load_pd(p); }

	/*
	 * Load 512-bits (composed of 8 packed double-precision (64-bit) floating-point elements) from memory into vect_t.
	 * p does not need to be aligned on any particular boundary.
	 * Return [p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]]
	 */
	static INLINE PURE vect_t loadu(const scalar_t *const p) { return _mm512_loadu_pd(p); }

	/*
	 * Store 512-bits (composed of 8 packed double-precision (64-bit) floating-point elements) from p into memory.
	 * p must be aligned on a 64-byte boundary or a general-protection exception will be generated.
	 */
	static INLINE void store(const scalar_t *p, const vect_t v) { _mm512_store_pd(const_cast<scalar_t *>(p), v); }

	/*
	 * Store 512-bits (composed of 8 packed double-precision (64-bit) floating-point elements) from p into memory.
	 * p does not need to be aligned on any particular boundary.
	 */
	static INLINE void storeu(const scalar_t *p, const vect_t v) { _mm512_storeu_pd(const_cast<scalar_t *>(p), v); }

	/*
	 * Store 512-bits (composed of 8 packed double-precision (64-bit) floating-point elements) from a into memory using
	 * a non-temporal memory hint.
	 * p must be aligned on a 64-byte boundary or a general-protection exception may be generated.
	 */
	static INLINE void stream(const scalar_t *p, const vect_t v) { _mm512_stream_pd(const_cast<scalar_t *>(p), v); }

	/*
	* Shuffle double-precision (64-bit) floating-point elements within 128-bit lanes using the control in imm8,
	* and store the results in dst.
	* Args   : [a0, a1, a2, a3] double
			   [b0, b1, b2, b3] double
	* Return : [a[s[0..1]], ..., a[s[6..7]]] double
	*/
/*#if defined(__FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS)
	template<uint8_t s>
	static INLINE CONST vect_t shuffle(const vect_t a) {
		return _mm256_permute4x64_pd(a, s);
	}
#endif
*/
	/*
	* Unpack and interleave double-precision (64-bit) floating-point elements from the low half of each 128-bit lane in a and b,
	* and store the results in dst.
	* Args   : [a0, a1, a2, a3, a4, a5, a6, a7] double
			   [b0, b1, b2, b3, b4, b5, b6, b7] double
	* Return : [a0, b0, a2, b2, a4, b4, a6, b6] double
	*/
	static INLINE CONST vect_t unpacklo_twice(const vect_t a, const vect_t b) { return _mm512_unpacklo_pd(a, b); }

	/*
	* Unpack and interleave double-precision (64-bit) floating-point elements from the high half of each 128-bit lane in a and b,
	* and store the results in dst.
	* Args   : [a0, a1, a2, a3, a4, a5, a6, a7] double
			   [b0, b1, b2, b3, b4, b5, b6, b7] double
	* Return : [a1, b1, a3, b3, a5, b5, a7, b7] double
	*/
	static INLINE CONST vect_t unpackhi_twice(const vect_t a, const vect_t b) { return _mm512_unpackhi_pd(a, b); }

	/*
	* Blend packed double-precision (64-bit) floating-point elements from a and b using control mask s,
	* and store the results in dst.
	* Args   : [a0, a1, a2, a3] double
			   [b0, b1, b2, b3] double
	* Return : [s[0]?a0:b0, ..., s[3]?a3:b3] double
	*/
	/*template<uint8_t s>
	static INLINE CONST vect_t blend(const vect_t a, const vect_t b) {
		return _mm256_blend_pd(a, b, s);
	}
*/
	/*
	* Blend packed double-precision (64-bit) floating-point elements from a and b using mask,
	* and store the results in dst.
	* Args   : [a0, a1, a2, a3] double
			   [b0, b1, b2, b3] double
	* Return : [mask[31]?a0:b0, ..., mask[255]?a3:b3] double
	*/
	/*static INLINE CONST vect_t blendv(const vect_t a, const vect_t b, const vect_t mask) {
		return _mm256_blendv_pd(a, b, mask);
	}
*/
	/*
	 * Add packed double-precision (64-bit) floating-point elements in a and b, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
	 * Return : [a0+b0, a1+b1, a2+b2, a3+b3, a4+b4, a5+b5, a6+b6, a7+b7]
	 */
	static INLINE CONST vect_t add(const vect_t a, const vect_t b) { return _mm512_add_pd(a, b); }

	static INLINE vect_t addin(vect_t &a, const vect_t b) { return a = add(a, b); }

	/*
	 * Subtract packed double-precision (64-bit) floating-point elements in b from packed double-precision (64-bit)
	 * floating-point elements in a, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
	 * Return : [a0-b0, a1-b1, a2-b2, a3-b3, a4-b4, a5-b5, a6-b6, a7-b7]
	 */
	static INLINE CONST vect_t sub(const vect_t a, const vect_t b) { return _mm512_sub_pd(a, b); }

	static INLINE CONST vect_t subin(vect_t &a, const vect_t b) { return a = sub(a, b); }

	/*
	 * Multiply packed double-precision (64-bit) floating-point elements in a and b, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
	 * Return : [a0*b0, a1*b1, a2*b2, a3*b3, a4*b4, a5*b5, a6*b6, a7*b7]
	 */
	static INLINE CONST vect_t mul(const vect_t a, const vect_t b) { return _mm512_mul_pd(a, b); }

	static INLINE CONST vect_t mulin(vect_t &a, const vect_t b) { return a = mul(a, b); }

	/*
	 * Divide packed double-precision (64-bit) floating-point elements in a by packed elements in b,
	 * and store the results in dst.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
	 * Return : [a0/b0, a1/b1, a2/b2, a3/b3, a4/b4, a5/b5, a6/b6, a7/b7]
	 */
	static INLINE CONST vect_t div(const vect_t a, const vect_t b) { return _mm512_div_pd(a, b); }

	/*
	 * Multiply packed double-precision (64-bit) floating-point elements in a and b, add the intermediate result to
	 * packed elements in c, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
	 * Return : [a0*b0+c0, a1*b1+c1, a2*b2+c2, a3*b3+c3, a4*b4+c4, a5*b5+c5, a6*b6+c6, a7*b7+c7]
	 */
	static INLINE CONST vect_t fmadd(const vect_t c, const vect_t a, const vect_t b) {
		return _mm512_fmadd_pd(a, b, c);
	}

	/*
	 * Multiply packed double-precision (64-bit) floating-point elements in a and b, add the intermediate result to
	 * packed elements in c, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
	 * Return : [a0*b0+c0, a1*b1+c1, a2*b2+c2, a3*b3+c3, a4*b4+c4, a5*b5+c5, a6*b6+c6, a7*b7+c7]
	 */
	static INLINE CONST vect_t madd(const vect_t c, const vect_t a, const vect_t b) { return fmadd(c, a, b); }

	/*
	 * Multiply packed double-precision (64-bit) floating-point elements in a and b, add the intermediate result to
	 * packed elements in c, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
	 * Return : [a0*b0+c0, a1*b1+c1, a2*b2+c2, a3*b3+c3, a4*b4+c4, a5*b5+c5, a6*b6+c6, a7*b7+c7]
	 */
	static INLINE CONST vect_t maddx(const vect_t c, const vect_t a, const vect_t b) { return fmadd(c, a, b); }

	static INLINE CONST vect_t fmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fmadd(c, a, b); }

	/*
	 * Multiply packed double-precision (64-bit) floating-point elements in a and b, add the negated intermediate result
	 * to packed elements in c, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
	 * Return : [-(a0*b0)+c0, -(a1*b1)+c1, -(a2*b2)+c2, -(a3*b3)+c3, -(a4*b4)+c4, -(a5*b5)+c5, -(a6*b6)+c6, -(a7*b7)+c7]
	 */
	static INLINE CONST vect_t fnmadd(const vect_t c, const vect_t a, const vect_t b) {
		return _mm512_fnmadd_pd(a, b, c);
	}

	/*
	 * Multiply packed double-precision (64-bit) floating-point elements in a and b, add the negated intermediate result
	 * to packed elements in c, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
	 * Return : [-(a0*b0)+c0, -(a1*b1)+c1, -(a2*b2)+c2, -(a3*b3)+c3, -(a4*b4)+c4, -(a5*b5)+c5, -(a6*b6)+c6, -(a7*b7)+c7]
	 */
	static INLINE CONST vect_t nmadd(const vect_t c, const vect_t a, const vect_t b) { return fnmadd(c, a, b); }

	static INLINE CONST vect_t fnmaddin(vect_t &c, const vect_t a, const vect_t b) { return c = fnmadd(c, a, b); }

	/*
	 * Multiply packed double-precision (64-bit) floating-point elements in a and b, subtract packed elements in c from
	 * the intermediate result, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
	 * Return : [a0*b0-c0, a1*b1-c1, a2*b2-c2, a3*b3-c3, a4*b4-c4, a5*b5-c5, a6*b6-c6, a7*b7-c7]
	 */
	static INLINE CONST vect_t fmsub(const vect_t c, const vect_t a, const vect_t b) {
		return _mm512_fmsub_pd(a, b, c);
	}

	/*
	 * Multiply packed double-precision (64-bit) floating-point elements in a and b, subtract packed elements in c from
	 * the intermediate result, and store the results in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7], [c0, c1, c2, c3, c4, c5, c6, c7]
	 * Return : [a0*b0-c0, a1*b1-c1, a2*b2-c2, a3*b3-c3, a4*b4-c4, a5*b5-c5, a6*b6-c6, a7*b7-c7]
	 */
	static INLINE CONST vect_t msub(const vect_t c, const vect_t a, const vect_t b) { return fmsub(c, a, b); }

	static INLINE CONST vect_t fmsubin(vect_t &c, const vect_t a, const vect_t b) { return c = fmsub(c, a, b); }

	/*
	 * Compare packed double-precision (64-bit) floating-point elements in a and b for equality, and store the results
	 in vect_t.
	 * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
	 * Return : [(a0==b0) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a1==b1) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a2==b2) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a3==b3) ? 0xFFFFFFFFFFFFFFFF : 0]
	 */
	static INLINE CONST vect_t eq(const vect_t a, const vect_t b) {
		int64_t i = 0xFFFFFFFFFFFFFFFF;
		__m512i c = _mm512_set1_epi64(i);
		return _mm512_maskz_expand_ps(_mm512_cmp_pd_mask(a, b, _CMP_EQ_OQ), _mm512_castsi512_ps(c)); 
	}

	/*
	 * Compare packed double-precision (64-bit) floating-point elements in a and b for lesser-than, and store the
	 results in vect_t.
	 * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
	 * Return : [(a0<b0) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a1<b1) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a2<b2) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a3<b3) ? 0xFFFFFFFFFFFFFFFF : 0]
	 */
	static INLINE CONST vect_t lesser(const vect_t a, const vect_t b) {
		int64_t i = 0xFFFFFFFFFFFFFFFF;
		__m512i c = _mm512_set1_epi64(i);
		return _mm512_maskz_expand_ps(_mm512_cmp_pd_mask(a, b, _CMP_LT_OS), _mm512_castsi512_ps(c)); 
	}

	/*
	 * Compare packed double-precision (64-bit) floating-point elements in a and b for lesser or equal than, and store
	 the results in vect_t.
	 * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
	 * Return : [(a0<=b0) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a1<=b1) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a2<=b2) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a3<=b3) ? 0xFFFFFFFFFFFFFFFF : 0]
	 */
	static INLINE CONST vect_t lesser_eq(const vect_t a, const vect_t b) {
		int64_t i = 0xFFFFFFFFFFFFFFFF;
		__m512i c = _mm512_set1_epi64(i);
		return _mm512_maskz_expand_ps(_mm512_cmp_pd_mask(a, b, _CMP_LE_OS), _mm512_castsi512_ps(c));
	}

	/*
	 * Compare packed double-precision (64-bit) floating-point elements in a and b for greater-than, and store the
	 results in vect_t.
	 * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
	 * Return : [(a0>b0) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a1>b1) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a2>b2) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a3>b3) ? 0xFFFFFFFFFFFFFFFF : 0]
	 */
	static INLINE CONST vect_t greater(const vect_t a, const vect_t b) {
		int64_t i = 0xFFFFFFFFFFFFFFFF;
		__m512i c = _mm512_set1_epi64(i);
		return _mm512_maskz_expand_ps(_mm512_cmp_pd_mask(a, b, _CMP_GT_OS), _mm512_castsi512_ps(c)); 
	}

	/*
	 * Compare packed double-precision (64-bit) floating-point elements in a and b for greater or equal than, and store
	 the results in vect_t.
	 * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
	 * Return : [(a0>=b0) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a1>=b1) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a2>=b2) ? 0xFFFFFFFFFFFFFFFF : 0,
	 (a3>=b3) ? 0xFFFFFFFFFFFFFFFF : 0]
	 */
	static INLINE CONST vect_t greater_eq(const vect_t a, const vect_t b) {
		int64_t i = 0xFFFFFFFFFFFFFFFF;
		__m512i c = _mm512_set1_epi64(i);
		return _mm512_maskz_expand_ps(_mm512_cmp_pd_mask(a, b, _CMP_GE_OS), _mm512_castsi512_ps(c)); 
	}

#ifdef __AVX512DQ__
	/*
	 * Compute the bitwise AND of packed double-precision (64-bit) floating-point elements in a and b, and store the
	 * results in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
	 * Return : [a0 AND b0, a1 AND b1, a2 AND b2, a3 AND b3, a4 AND b4, a5 AND b5, a6 AND b6, a7 AND b7]
	 */
	static INLINE CONST vect_t vand(const vect_t a, const vect_t b) { return _mm512_and_pd(a, b); }

	/*
	 * Compute the bitwise OR of packed double-precision (64-bit) floating-point elements in a and b, and store the
	 * results in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
	 * Return : [a0 OR b0, a1 OR b1, a2 OR b2, a3 OR b3, a4 OR b4, a5 OR b5, a6 OR b6, a7 OR b7]
	 */
	static INLINE CONST vect_t vor(const vect_t a, const vect_t b) { return _mm512_or_pd(a, b); }

	/*
	 * Compute the bitwise XOR of packed double-precision (64-bit) floating-point elements in a and b, and store the
	 * results in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
	 * Return : [a0 XOR b0, a1 XOR b1, a2 XOR b2, a3 XOR b3, a4 XOR b4, a5 XOR b5, a6 XOR b6, a7 XOR b7]
	 */
	static INLINE CONST vect_t vxor(const vect_t a, const vect_t b) { return _mm512_xor_pd(a, b); }

	/*
	 * Compute the bitwise AND NOT of packed double-precision (64-bit) floating-point elements in a and b, and store the
	 * results in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7], [b0, b1, b2, b3, b4, b5, b6, b7]
	 * Return : [a0 ANDNOT b0, a1 ANDNOT b1, a2 ANDNOT b2, a3 ANDNOT b3, a4 ANDNOT b4, a5 ANDNOT b5, a6 ANDNOT b6, a7
	 * ANDNOT b7]
	 */
	static INLINE CONST vect_t vandnot(const vect_t a, const vect_t b) { return _mm512_andnot_pd(a, b); }
#else //AVX512DQ
#endif
	/*
	 * Round the packed double-precision (64-bit) floating-point elements in a down to an integer value, and store the
	 * results as packed double-precision floating-point elements in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
	 * Return : [floor(a0), floor(a1), floor(a2), floor(a3), floor(a4), floor(a5), floor(a6), floor(a7)]
	 */
	static INLINE CONST vect_t floor(const vect_t a) { return _mm512_floor_pd(a); }

	/*
	 * Round the packed double-precision (64-bit) floating-point elements in a up to an integer value, and store the
	 * results as packed double-precision floating-point elements in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
	 * Return : [ceil(a0), ceil(a1), ceil(a2), ceil(a3), ceil(a4), ceil(a5), ceil(a6), ceil(a7)]
	 */
	static INLINE CONST vect_t ceil(const vect_t a) { return _mm512_ceil_pd(a); }

	/*
	 * Round the packed double-precision (64-bit) floating-point elements in a, and store the results as packed
	 * double-precision floating-point elements in vect_t.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
	 * Return : [round(a0), round(a1), round(a2), round(a3), round(a4), round(a5), round(a6), round(a7)]
	 */
	static INLINE CONST vect_t round(const vect_t a) {
		int i = 0;
		return _mm512_roundscale_pd(a, i);
	}

	/*
	 * Horizontally add adjacent pairs of double-precision (64-bit) floating-point elements in a and b, and pack the
	 * results in vect_t.
	 * Args   : [a0, a1, a2, a3], [b0, b1, b2, b3]
	 * Return : [a0+a1, b0+b1, a2+a3, b2+b3]
	 */
	//static INLINE CONST vect_t hadd(const vect_t a, const vect_t b) { return _mm256_hadd_pd(a, b); }

	/*
	 * Horizontally add double-precision (64-bit) floating-point elements in a.
	 * Args   : [a0, a1, a2, a3, a4, a5, a6, a7]
	 * Return : a0+a1+a2+a3+a4+a5+a6+a7
	 */
	static INLINE CONST scalar_t hadd_to_scal(const vect_t a) {
		return ((const scalar_t *)&a)[0] + ((const scalar_t *)&a)[1] + ((const scalar_t *)&a)[2] +
				((const scalar_t *)&a)[3] + ((const scalar_t *)&a)[4] + ((const scalar_t *)&a)[5] +
				((const scalar_t *)&a)[6] + ((const scalar_t *)&a)[7];
	}

	static INLINE vect_t mod(vect_t &C, const vect_t &P, const vect_t &INVP, const vect_t &NEGP, const vect_t &MIN,
							 const vect_t &MAX, vect_t &Q, vect_t &T) {
		FLOAT_MOD(C, P, INVP, Q);
		NORML_MOD(C, P, NEGP, MIN, MAX, Q, T);

		return C;
	}

};

#endif // __FFLASFFPACK_simd256_double_INL
