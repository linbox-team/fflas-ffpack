/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* fflas/fflas_finit.inl
 * Copyright (C) 2014 Pascal Giorgi
 *
 * Written by Pascal Giorgi <Pascal.Giorgi@lirmm.fr>
 * BB<bboyer@ncsu.edu>
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

// #include <x86intrin.h>
#include <immintrin.h>

namespace FFLAS {

#if defined(__AVX2__) || defined(__AVX__)
#ifdef __AVX2__

#define VEC_MODF_D(C,Q,P,PP,TMP)						\
	{Q= _mm256_mul_pd(C,PP);  Q=_mm256_floor_pd(Q); C = _mm256_fnmadd_pd(Q,P,C); \
		Q= C; Q=_mm256_cmp_pd(Q,P,_CMP_GE_OS); Q=_mm256_and_pd(Q,P); C=_mm256_sub_pd(C,Q);}

#define VEC_MODF_S(C,Q,P,PP,TMP)						\
	{Q= _mm256_mul_ps(C,PP);  Q=_mm256_floor_ps(Q); C = _mm256_fnmadd_ps(Q,P,C); \
		Q= C; Q=_mm256_cmp_ps(Q,P,_CMP_GE_OS); Q=_mm256_and_ps(Q,P); C=_mm256_sub_ps(C,Q);}


#else // __AVX__

#define VEC_MODF_D(C,Q,P,PP,TMP)						\
	{Q= _mm256_mul_pd(C,PP);  Q=_mm256_floor_pd(Q); TMP = _mm256_mul_pd(Q,P); C= _mm256_sub_pd(C,TMP); \
		Q= C; Q=_mm256_cmp_pd(Q,P,_CMP_GE_OS); Q=_mm256_and_pd(Q,P); C=_mm256_sub_pd(C,Q);}

#define VEC_MODF_S(C,Q,P,PP,TMP)						\
	{Q= _mm256_mul_ps(C,PP);  Q=_mm256_floor_ps(Q); TMP = _mm256_mul_ps(Q,P); C= _mm256_sub_ps(C,TMP); \
		Q= C; Q=_mm256_cmp_ps(Q,P,_CMP_GE_OS); Q=_mm256_and_ps(Q,P); C=_mm256_sub_ps(C,Q);}

#endif // AVX or AVX2

	inline void modp( double *T, size_t n, double p, double invp){
		register __m256d C,Q,P,PP,TMP;
		P = _mm256_set1_pd(p);
		PP= _mm256_set1_pd(invp);
		long st=long(T)%32;
		size_t i=0;;
		if (st){ // the array T is not 32 byte aligned (process few elements s.t. (T+i) is 32 bytes aligned)
			for (size_t j=st;j<32;j+=8,i++){
				T[i]=fmod(T[i],p);
				// if (T[i]>=p)
					// T[i]-=p;
			}
		}
		// perform the loop using 256 bits SIMD
		for (;i<=n-4;i+=4){
			C=_mm256_load_pd(T+i);
			VEC_MODF_D(C,Q,P,PP,TMP);
			_mm256_store_pd(T+i,C);
		}
		// perform the last elt from T without SIMD
		for (;i<n;i++){
			T[i]=fmod(T[i],p);
			// if (T[i]>=p)
				// T[i]-=p;
		}
	}

	inline void modp( float *T, size_t n, float p, float invp){
		register __m256  C,Q,P,PP,TMP;
		P = _mm256_set1_ps(p);
		PP= _mm256_set1_ps(invp);
		long st=long(T)%32;
		size_t i=0;;
		if (st){ // the array T is not 32 byte aligned (process few elements s.t. (T+i) is 32 bytes aligned)
			for (size_t j=st;j<32;j+=4,i++){
				// std::cout << "avx2" << std::endl;
				T[i]=fmodf(T[i],p);
				// if (T[i]>=p)
					// T[i]-=p;
			}
		}
		// perform the loop using 256 bits SIMD
		for (;i<=n-8;i+=8){
			C=_mm256_load_ps(T+i);
			VEC_MODF_S(C,Q,P,PP,TMP);
			_mm256_store_ps(T+i,C);
		}
		// perform the last elt from T without SIMD
		for (;i<n;i++){
			T[i]=fmodf(T[i],p);
			// if (T[i]>=p)
				// T[i]-=p;
		}
	}

#else // no AVX

	inline void modp( double *T, size_t n, double p, double invp){
		for(size_t j=0;j<n;j++){
			T[j]= fmod(T[j],p);
			// T[j]-= ((T[j]>=p)?p:0);
		}

	}

	inline void modp( float *T, size_t n, float p, float invp){
		for(size_t j=0;j<n;j++){
			T[j]= fmodf(T[j],p);
			// T[j]-= ((T[j]>=p)?p:0);
		}

	}

#endif // AVX(2) oupa

	template<>
	void finit (const FFPACK:: Modular<double> & F, const size_t m,
		    double * A, const size_t incX)
	{
		if(incX == 1) {
			double p, invp;
			p=(double)F.cardinality();
			invp=1./p;
			modp(A,m,p,invp);
			for(size_t i=0;i<m;i++)
				if (A[i] < 0)
					A[i] += p ;

		}
		else { /*  faster with copy, use incX=1, copy back ? */
			double * Xi = A ;
			for (; Xi < A+m*incX; Xi+=incX )
				F.init( *Xi , *Xi);

		}
	}

	template<>
	void finit (const FFPACK:: ModularBalanced<double> & F, const size_t m,
		    double * A, const size_t incX)
	{
		if(incX == 1) {
			double p, invp;
			p=(double)F.cardinality();
			invp=1./p;
			double pmax = (p-1)/2 ;
			double pmin = pmax-p+1;
			modp(A,m,p,invp);
			for(size_t i=0;i<m;i++)
				if (A[i] < pmin)
					A[i] += p ;
				else if (A[i] > pmax)
					A[i] -= p ;

		}
		else { /*  faster with copy, use incX=1, copy back ? */
			double * Xi = A ;
			for (; Xi < A+m*incX; Xi+=incX )
				F.init( *Xi , *Xi);
		}

	}

	template<>
	void finit (const FFPACK:: Modular<float> & F, const size_t m,
		    float * A, const size_t incX)
	{
		if(incX == 1) {
			float p, invp;
			p=(float)F.cardinality();
			invp=1.f/p;
			modp(A,m,p,invp);
			for(size_t i=0;i<m;i++)
				if (A[i] < 0)
					A[i] += p ;

		}
		else { /*  faster with copy, use incX=1, copy back ? */
			float * Xi = A ;
			for (; Xi < A+m*incX; Xi+=incX )
				F.init( *Xi , *Xi);

		}

	}

	template<>
	void finit (const FFPACK:: ModularBalanced<float> & F, const size_t m,
		    float * A, const size_t incX)
	{
		if(incX == 1) {
			float p, invp;
			p=(float)F.cardinality();
			invp=1.f/p;
			float pmax = (p-1)/2 ;
			float pmin = pmax-p+1;
			modp(A,m,p,invp);
			for(size_t i=0;i<m;i++) {
				if (A[i] < pmin)
					A[i] += p ;
				else if (A[i] > pmax)
					A[i] -= p ;
			}

		}
		else { /*  faster with copy, use incX=1, copy back ? */
			float * Xi = A ;
			for (; Xi < A+m*incX; Xi+=incX )
				F.init( *Xi , *Xi);
		}

	}


#if 0
	template<>
	void finit (const FFPACK:: Modular<double> & F, const size_t m , const size_t n,
		    double * A, const size_t lda){
		double p, invp;
		p=(double)F.cardinality();
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

	template<>
	void finit (const FFPACK:: ModularBalanced<double> & F, const size_t m , const size_t n,
		    double * A, const size_t lda){
		double p, invp;
		p=(double)F.cardinality();
		invp=1./p;
		if(n==lda)
			modp(A,m*n,p,invp);
		else
			for(size_t i=0;i<m;i++)
				modp(A+i*lda,n,p,invp);
#if 1 /* BB: just making the tests pass */
		double pmax = (p-1)/2 ;
		double pmin = pmax-p+1;

		if(n==lda) {
			for (size_t i = 0 ; i < m*n ; ++i)
				if (A[i] < pmin)
					A[i] += p ;
				else if (A[i] > pmax)
					A[i] -= p ;
		}
		else {
			for(size_t i=0;i<m;i++)
				for(size_t j=0;j<n;j++)
					if (A[i*lda+j] < pmin)
						A[i*lda+j] += p ;
					else if (A[i*lda+j] > pmax)
						A[i*lda+j] -= p ;

		}
#endif

	}
#endif

}

#endif // __FFLASFFPACK_fflas_init_INL

