/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* fflas/fflas_finit.inl
 * Copyright (C) 2014 Pascal Giorgi
 *
 * Written by Pascal Giorgi <Pascal.Giorgi@lirmm.fr>
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
#ifndef __FFLASFFPACK_fflas_init_INL
#define __FFLASFFPACK_fflas_init_INL

namespace FFLAS {

#if defined(__AVX2__) || defined(__AVX__)
#ifdef __AVX2__
#define VEC_MODF(C,Q,P,PP,TMP)						\
	{Q= _mm256_mul_pd(C,PP);  Q=_mm256_floor_pd(Q); C = _mm256_fnmadd_pd(Q,P,C); \
		Q= C; Q=_mm256_cmp_pd(Q,P,13); Q=_mm256_and_pd(Q,P); C=_mm256_sub_pd(C,Q);}
#else
#define VEC_MODF(C,Q,P,PP,TMP)						\
	{Q= _mm256_mul_pd(C,PP);  Q=_mm256_floor_pd(Q); TMP = _mm256_mul_pd(Q,P); C= _mm256_sub_pd(C,TMP); \
		Q= C; Q=_mm256_cmp_pd(Q,P,13); Q=_mm256_and_pd(Q,P); C=_mm256_sub_pd(C,Q);}
#endif

	inline void modp( double *T, size_t n, double p, double invp){
		register __m256d C,Q,P,PP,TMP;
		P = _mm256_set1_pd(p);
		PP= _mm256_set1_pd(invp);
		long st=long(T)%32;
		size_t i=0;;
		if (st){ // the array T is not 32 byte aligned (process few elements s.t. (T+i) is 32 bytes aligned)
			for (size_t j=st;j<32;j+=8,i++){
				T[i]=fmod(T[i],p);
				if (T[i]>=p)
					T[i]-=p;
			}
		}
		// perform the loop using 256 bits SIMD
		for (;i<=n-4;i+=4){
			C=_mm256_load_pd(T+i);
			VEC_MODF(C,Q,P,PP,TMP);
			_mm256_store_pd(T+i,C);
		}
		// perform the last elt from T without SIMD
		for (;i<n;i++){
			T[i]=fmod(T[i],p);
			if (T[i]>=p)
				T[i]-=p;
		}
	}
#else

	inline void modp( double *T, size_t n, double p, double invp){
		for(size_t j=0;j<n;j++){
			T[j]= fmod(T[j],p);
			T[j]-= ((T[j]>=p)?p:0);
		}

	}
#endif

	template<>
	void finit (const FFPACK:: Modular<double> & F, const size_t m , const size_t n,
		    double * A, const size_t lda){
		double p, invp;
		p=F.cardinality();
		invp=1./p;
		if(n==lda)
			modp(A,m*n,p,invp);
		else
			for(size_t i=0;i<m;i++)
				modp(A+i*lda,n,p,invp);
#if 1 /* BB: just making the tests pass */
		if(n==lda) {
			for (size_t i = 0 ; i < m*n ; ++i)
				if (A[i] < 0)
					A[i] += p ;
		}
		else {
			for(size_t i=0;i<m;i++)
				for(size_t j=0;j<n;j++)
					if (A[i*lda+j] < 0)
						A[i*lda+j] += p ;
		}
#endif

	}


}

#endif

