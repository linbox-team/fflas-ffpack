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
#include "fflas-ffpack/field/unparametric.h"

#define FFLASFFPACK_COPY_INIT 100

namespace FFLAS {

#define DEBUG_M256(M) \
	{                    \
		float ssssssss[8] ; \
		_mm256_store_ps(ssssssss,M); \
		for (size_t ii = 0 ; ii < 8 ; ++ii) \
		std::cout << ssssssss[ii] << ", " ; \
		std::cout << std::endl; \
	}


#if defined(__AVX2__) || defined(__AVX__)
#ifdef __AVX2__
	// compute C modulo P in the range [MIN, MAX]
#define VEC_MODF_D(C,Q,P,NEGP,INVP,T,MIN,MAX)				\
	{								\
		Q = _mm256_mul_pd(C,INVP);  Q=_mm256_floor_pd(Q);	\
		C = _mm256_fnmadd_pd(Q,P,C);				\
		Q = _mm256_cmp_pd(C,MAX,_CMP_GT_OS);			\
		T = _mm256_cmp_pd(C,MIN,_CMP_LT_OS);			\
		Q = _mm256_and_pd(Q,NEGP);				\
		T = _mm256_and_pd(T,P);					\
		Q = _mm256_or_pd(Q,T);				\
		C = _mm256_add_pd(C,Q);					\
	}
	// compute C modulo P in the range [MIN, MAX]
#define VEC_MODF_S(C,Q,P,NEGP,INVP,T,MIN,MAX)				\
	{								\
		Q = _mm256_mul_ps(C,INVP);  Q=_mm256_floor_ps(Q);	\
		C = _mm256_fnmadd_ps(Q,P,C);				\
		Q = _mm256_cmp_ps(C,MAX,_CMP_GT_OS);			\
		T = _mm256_cmp_ps(C,MIN,_CMP_LT_OS);			\
		Q = _mm256_and_ps(Q,NEGP);				\
		T = _mm256_and_ps(T,P);					\
		Q = _mm256_or_ps(Q,T);				\
		C = _mm256_add_ps(C,Q);					\
	}

#else // __AVX__
	// compute C modulo P in the range [MIN, MAX]
#define VEC_MODF_D(C,Q,P,NEGP,INVP,T,MIN,MAX)				\
	{								\
		Q = _mm256_mul_pd(C,INVP);  Q=_mm256_floor_pd(Q);	\
		T = _mm256_mul_pd(Q,P); C= _mm256_sub_pd(C,T);		\
		Q = _mm256_cmp_pd(C,MAX,_CMP_GT_OS);			\
		T = _mm256_cmp_pd(C,MIN,_CMP_LT_OS);			\
		Q = _mm256_and_pd(Q,NEGP);				\
		T = _mm256_and_pd(T,P);					\
		Q = _mm256_or_pd(Q,T);				\
		C = _mm256_add_pd(C,Q);					\
	}
	// compute C modulo P in the range [MIN, MAX]
#define VEC_MODF_S(C,Q,P,NEGP,INVP,T,MIN,MAX)				\
	{  							\
		Q = _mm256_mul_ps(C,INVP); Q=_mm256_floor_ps(Q);	\
		T = _mm256_mul_ps(Q,P); C= _mm256_sub_ps(C,T);		\
		Q = _mm256_cmp_ps(C,MAX,_CMP_GT_OS);			\
		T = _mm256_cmp_ps(C,MIN,_CMP_LT_OS);			\
		Q = _mm256_and_ps(Q,NEGP);				\
		T = _mm256_and_ps(T,P);					\
		Q = _mm256_or_ps(Q,T);				\
		C = _mm256_add_ps(C,Q);					\
	}

#endif // AVX or AVX2

	inline void modp( double *T, size_t n, double p, double invp, double min, double max){
		register __m256d C,Q,P,NEGP,INVP,TMP,MIN,MAX;
		P   = _mm256_set1_pd(p);
		NEGP= _mm256_set1_pd(-p);
		INVP= _mm256_set1_pd(invp);
		MIN = _mm256_set1_pd(min);
		MAX = _mm256_set1_pd(max);
		long st=long(T)%32;
		size_t i=0;
		if (n < 4) {
			for (;i<n;i++){
				T[i]=fmod(T[i],p);
				T[i]-=(T[i]>max)?p:0;
				T[i]+=(T[i]<min)?p:0;
			}
				return;

		}
		if (st){ // the array T is not 32 byte aligned (process few elements s.t. (T+i) is 32 bytes aligned)
			for (size_t j=(size_t)st;j<32;j+=8,i++){
				T[i]=fmod(T[i],p);
				T[i]-=(T[i]>max)?p:0;
				T[i]+=(T[i]<min)?p:0;
			}
		}
		FFLASFFPACK_check((long(T+i)%32==0));
		// perform the loop using 256 bits SIMD
		for (;i<=n-4;i+=4){
			C=_mm256_load_pd(T+i);
			VEC_MODF_D(C,Q,P,NEGP,INVP,TMP,MIN,MAX);
			_mm256_store_pd(T+i,C);
		}
		// perform the last elt from T without SIMD
		for (;i<n;i++){
			T[i]=fmod(T[i],p);
			T[i]-=(T[i]>max)?p:0;
			T[i]+=(T[i]<min)?p:0;
		}
	}

	inline void modp( float *T,  size_t n, float p, float invp, float min, float max){
		register __m256 C,Q,P,NEGP,INVP,TMP,MIN,MAX;
		P   = _mm256_set1_ps(p);
		NEGP= _mm256_set1_ps(-p);
		INVP= _mm256_set1_ps(invp);
		MIN= _mm256_set1_ps(min);
		MAX= _mm256_set1_ps(max);
		long st=long(T)%32;
		size_t i=0;;
		if (n < 8) {
			for (;i<n;i++){
				T[i]=fmodf(T[i],p);
				T[i]-=(T[i]>max)?p:0;
				T[i]+=(T[i]<min)?p:0;
			}
				return;
		}
		if (st){ // the array T is not 32 byte aligned (process few elements s.t. (T+i) is 32 bytes aligned)
			for (size_t j=(size_t)st;j<32;j+=4,i++){
				T[i]=fmodf(T[i],p);
				T[i]-=(T[i]>max)?p:0;
				T[i]+=(T[i]<min)?p:0;
			}
		}
		// perform the loop using 256 bits SIMD
		for (;i<=n-8;i+=8){
			C=_mm256_load_ps(T+i);
			VEC_MODF_S(C,Q,P,NEGP,INVP,TMP,MIN,MAX);
			_mm256_store_ps(T+i,C);
		}
		// perform the last elt from T without SIMD
		for (;i<n;i++){
			T[i]=fmodf(T[i],p);
			T[i]-=(T[i]>max)?p:0;
			T[i]+=(T[i]<min)?p:0;
		}
	}

#else // no AVX

	inline void modp( double *T, size_t n, double p, double invp, double min, double max){
		for(size_t j=0;j<n;j++){
			T[j]= fmod(T[j],p);
			T[j]-=(T[j]>max)?p:0;
			T[j]+=(T[j]<min)?p:0;
		}

	}

	inline void modp( float *T, size_t n, float p, float invp, float min, float max){
		for(size_t j=0;j<n;j++){
			T[j]= fmodf(T[j],p);
			T[j]-=(T[j]>max)?p:0;
			T[j]+=(T[j]<min)?p:0;
		}

	}

#endif // AVX(2)

	template<>
	void finit (const FFPACK:: Modular<double> & F, const size_t m,
		    double * A, const size_t incX)
	{
		if(incX == 1) {
			double p, invp;
			p=(double)F.cardinality();
			invp=1./p;
			modp(A,m,p,invp,0,p-1);
		}
		else { /*  faster with copy, use incX=1, copy back ? */
			if (m < FFLASFFPACK_COPY_INIT) {
				double * Xi = A ;
				for (; Xi < A+m*incX; Xi+=incX )
					F.init( *Xi , *Xi);
			}
			else {
				double * Ac = new double[m] ;
				fcopy(F,m,A,incX,Ac,1);
				finit(F,m,Ac,1);
				fcopy(F,m,Ac,1,A,incX);
				delete[] Ac;
			}
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
			modp(A,m,p,invp,0,p-1);
		}
		else { /*  faster with copy, use incX=1, copy back ? */
			if (m < FFLASFFPACK_COPY_INIT) {
				float * Xi = A ;
				for (; Xi < A+m*incX; Xi+=incX )
					F.init( *Xi , *Xi);
			}
			else {
				float * Ac = new float[m] ;
				fcopy(F,m,A,incX,Ac,1);
				finit(F,m,Ac,1);
				fcopy(F,m,Ac,1,A,incX);
				delete[] Ac;
			}

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
			modp(A,m,p,invp,pmin,pmax);
		}
		else { /*  faster with copy, use incX=1, copy back ? */
			if (m < FFLASFFPACK_COPY_INIT) {
				double * Xi = A ;
				for (; Xi < A+m*incX; Xi+=incX )
					F.init( *Xi , *Xi);
			}
			else {
				double * Ac = new double[m] ;
				fcopy(F,m,A,incX,Ac,1);
				finit(F,m,Ac,1);
				fcopy(F,m,Ac,1,A,incX);
				delete[] Ac;
			}
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
			modp(A,m,p,invp,pmin,pmax);
		}
		else { /*  faster with copy, use incX=1, copy back ? */
			if (m < FFLASFFPACK_COPY_INIT) {
				float * Xi = A ;
				for (; Xi < A+m*incX; Xi+=incX )
					F.init( *Xi , *Xi);
			}
			else {
				float * Ac = new float[m] ;
				fcopy(F,m,A,incX,Ac,1);
				finit(F,m,Ac,1);
				fcopy(F,m,Ac,1,A,incX);
				delete[] Ac;
			}
		}

	}

	template<class Field, class OtherElement>
	void
	finit (const Field& F, const size_t n,
	       const OtherElement * Y, const size_t incY,
	       typename Field::Element * X, const size_t incX)
	{
		typename Field::Element * Xi = X ;
		const OtherElement * Yi = Y ;

		if (incX == 1 && incY == 1)
			for (; Xi < X + n ; ++Xi, ++Yi)
				F.init( *Xi , *Yi);
		else
			for (; Xi < X+n*incX; Xi+=incX, Yi += incX )
				F.init( *Xi , *Yi);
	}


	template<class Field>
	void
	finit (const Field& F, const size_t m , const size_t n,
	       typename Field::Element * A, const size_t lda)
	{
		if (n == lda)
			finit(F,n*m,A,1);
		else
			for (size_t i = 0 ; i < m ; ++i)
				finit(F,n,A+i*lda,1);
		return;
	}

	template<class Field, class OtherElement>
	void
	finit (const Field& F, const size_t m , const size_t n,
	       const OtherElement * B, const size_t ldb,
	       typename Field::Element * A, const size_t lda)
	{
		if (n == lda && lda == ldb)
			finit(F,n*m,B,1,A,1);
		else
			for (size_t i = 0 ; i < m ; ++i)
				finit(F,n,B+i*ldb,1,A+i*lda,1);
		return;
	}

} // end of namespace FFLAS
#endif // __FFLASFFPACK_fflas_init_INL

