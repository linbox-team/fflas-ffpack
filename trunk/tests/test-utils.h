/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Brice Boyer <bboyer@imag.fr>
 * This file is Free Software and part of FFLAS-FFPACK.
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

/*! @file tests/test-utils.h
 * @ingroup tests
 * @brief Utilities to create matrices with prescribed shapes, properties,...
 * To be used in the tests
 */

#ifndef __FFLASFFPACK_tests_test_utils_H
#define __FFLASFFPACK_tests_test_utils_H

#include "fflas-ffpack/field/nonzero-randiter.h"
#include "fflas-ffpack/utils/debug.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "givaro/givinteger.h"
#include <givaro/givintprime.h>
#include <chrono>
#include <random>

namespace FFPACK {

	/*! @brief  Random Matrix.
	 * Creates a \c m x \c n matrix with random entries.
	 * @param F field
	 * @param A pointer to the matrix (preallocated to at least \c m x \c lda field elements)
	 * @param m number of rows in \p A
	 * @param n number of cols in \p A
	 * @param lda leading dimension of \p A
	 * @return pointer to \c A.
	 */
	template<class Field>
	typename Field::Element * RandomMatrix(const Field & F,
					       typename Field::Element * A,
					       size_t m, size_t n, size_t lda)
	{
		typedef typename Field::RandIter Randiter ;
		Randiter R(F);
		for (size_t i=0 ; i<m ; ++i)
			for (size_t j= 0; j<n ;++j)
				R.random( A[i*lda+j] );
		return A;

	}

	template<class T >
	T * RandomMatrix(const Givaro::UnparametricRing< T > & F,
			 T * A,
			 size_t m, size_t n, size_t lda)
	{
		Givaro::Modular<T> G(101);
		RandomMatrix(G,A,m,n,lda);
		return A;

	}
	/*! Random integer in range.
	 * @param a min bound
	 * @param b max bound
	 * @return a random integer in [a,b[  */
	int RandInt(int a, int b)
	{
		int x = a ;
		x += rand()%(b-a+1);
		FFLASFFPACK_check(x<b && x>=a);
		return x ;
	}

	/*! Random integer in range.
	 * @param a min bound
	 * @param b max bound
	 * @return a random integer in [a,b[  */
	size_t RandInt(size_t a, size_t b)
	{
		size_t x = a ;
		x += (size_t)rand()%(b-a);
		FFLASFFPACK_check(x<b && x>=a);
		return x ;
	}

	/*! @brief  Random Matrix with prescribed rank.
	 * Creates an \c m x \c n matrix with random entries and rank \c r.
	 * @param F field
	 * @param A pointer to the matrix (preallocated to at least \c m x \c lda field elements)
	 * @param r rank of the matrix to build
	 * @param m number of rows in \p A
	 * @param n number of cols in \p A
	 * @param lda leading dimension of \p A
	 * @return pointer to \c A.
	 */
	template<class Field>
	typename Field::Element_ptr RandomMatrixWithRank (const Field & F,
							  typename Field::Element_ptr A, size_t lda,
							  size_t r, size_t m, size_t n)
	{
		FFLASFFPACK_check(r <= std::min(m,n));
		FFLASFFPACK_check(n <= lda);
		typedef typename Field::RandIter Randiter ;
		//typedef typename Field::Element  Element ;
		typedef typename Field::Element_ptr  Element_ptr;
		Randiter R(F);
		NonzeroRandIter<Field,Randiter> nzR(F,R);

		size_t * P = FFLAS::fflas_new<size_t>(n);
		size_t * Q = FFLAS::fflas_new<size_t>(m);
		for (size_t i = 0 ; i < m ; ++i ) Q[i] = 0;
		for (size_t i = 0 ; i < n ; ++i ) P[i] = 0;

		Element_ptr U = FFLAS::fflas_new(F,m,n);
		Element_ptr L = FFLAS::fflas_new(F,m,m);


		/*  Create L, lower invertible */
		for (size_t i=0 ; i<m ; ++i){
			for (size_t j= 0; j<i ;++j) R.random( L[i*m+j] );
			nzR.random( L[i*m+i] );
			for (size_t j= i+1; j<m ;++j) F.init(L[i*m+j],F.zero);
		}

		/*  Create U, upper or rank r */
		for (size_t i=0 ; i<r ; ++i){
			for (size_t j= 0 ; j<i ;++j) F.init(U[i*n+j],0UL);
			nzR.random( U[i*n+i] );
			for (size_t j= i+1; j<n ;++j) R.random( U[i*n+j] );
		}
		for (size_t i=r ; i<m ; ++i)
			for (size_t j= 0 ; j<n ;++j)
				F.init(U[i*n+j],F.zero);

		/*  Create a random P,Q */
		for (size_t i = 0 ; i < n ; ++i)
			P[i] = i + RandInt(0UL,n-i);
		for (size_t i = 0 ; i < m ; ++i)
			Q[i] = i + RandInt(0UL,m-i);

		/*  compute product */

		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				m,0,(int)n, U, n, P);
		FFPACK::applyP (F, FFLAS::FflasLeft,  FFLAS::FflasNoTrans,
				m,0,(int)m, L, m, Q);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
			      m, n, m, F.one, L, m, U, n, F.zero, A, lda);
		//! @todo compute LU with ftrtr

		FFLAS::fflas_delete(P);
		FFLAS::fflas_delete(L);
		FFLAS::fflas_delete(U);
		FFLAS::fflas_delete(Q);

		return A;

	}

	void RandomRankProfile (size_t N, size_t R, size_t* rkp){
		size_t curr = 0;
		std::vector<bool> rows(N,false);
		while (curr<R){
			size_t i;
			while (rows [i = rand() % N]);
			rows[i] = true;
			rkp [curr] = i;
			curr++;
		}
	}

	
	template<class Field>
	void RandomMatrixWithRankandRPM (const Field& F, typename Field::Element_ptr A, size_t lda,
					 size_t R, size_t M, size_t N,
					 const size_t * RRP, const size_t * CRP){
	
	typedef typename Field::RandIter Randiter ;
	Randiter RI(F);
	FFPACK::NonzeroRandIter<Field,Randiter> nzR(F,RI);
	typename Field::Element_ptr L= FFLAS::fflas_new(F,M,N);
	FFLAS::fzero(F, M, N, L, N);
	for (size_t k = 0; k < R; ++k){
		size_t i = RRP[k];
		size_t j = CRP[k];
		nzR.random (L [i*N+j]);
		for (size_t l=i+1; l < M; ++l)
			RI.random (L [l*N+j]);
	}
	typename Field::Element_ptr U= FFLAS::fflas_new(F,N,N);
	FFLAS::fzero(F, N, N, U, N);
	for (size_t i = 0; i < N; ++i){
		nzR.random (U [i*N+i]);
		for (size_t j=i+1; j < N; ++j)
			RI.random (U [i*N+j]);
	}


	typename Field::Element alpha, beta;
	F.init(alpha,1.0);
	F.init(beta,0.0);
	PAR_BLOCK{
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M,N,N, alpha, L, N, U, N, beta, A, lda, FFLAS::ParSeqHelper::Parallel());
	}
	FFLAS::fflas_delete(L);
	FFLAS::fflas_delete(U);

	}

        /*! @brief  Random Matrix with prescribed rank, with random  rank profile matrix
	 * Creates an \c m x \c n matrix with random entries, rank \c r and with a rank profile matrix
	 * chosen uniformly at random.
	 * @param F field
	 * @param A pointer to the matrix (preallocated to at least \c m x \c lda field elements)
	 * @param r rank of the matrix to build
	 * @param m number of rows in \p A
	 * @param n number of cols in \p A
	 * @param lda leading dimension of \p A
	 * @return pointer to \c A.
	 */
	template<class Field>
	void RandomMatrixWithRankandRandomRPM (const Field& F, typename Field::Element_ptr A, size_t lda,
					       size_t R, size_t M, size_t N){
		    // generate the r pivots in the rank profile matrix E
		size_t pivot_r[R];
		size_t pivot_c[R];
			
		RandomRankProfile (M, R, pivot_r);
		RandomRankProfile (N, R, pivot_c);
		RandomMatrixWithRankandRPM (F, A, lda, R, M, N, pivot_r, pivot_c);
	}
	
	/*! @brief  Random Matrix with prescribed det.
	 * @bug duplicate with linbox
	 * Creates a \c m x \c n matrix with random entries and rank \c r.
	 * @param F field
	 * @param A pointer to the matrix (preallocated to at least \c m x \c lda field elements)
	 * @param r rank of the matrix to build
	 * @param m number of rows in \p A
	 * @param n number of cols in \p A
	 * @param lda leading dimension of \p A
	 * @return pointer to \c A.
	 */
	template<class Field>
	typename Field::Element * RandomMatrixWithDet(const Field & F,
						      typename Field::Element * A,
						      typename Field::Element d,
						      size_t n, size_t lda)
	{
		FFLASFFPACK_check(n <= lda);
		typedef typename Field::RandIter Randiter ;
		typedef typename Field::Element  Element ;
		Randiter R(F);
		NonzeroRandIter<Field,Randiter> nzR(F,R);

		size_t * P = FFLAS::fflas_new<size_t>(n);
		size_t * Q = FFLAS::fflas_new<size_t>(n);
		for (size_t i = 0 ; i < n ; ++i ) Q[i] = 0;
		for (size_t i = 0 ; i < n ; ++i ) P[i] = 0;

		Element * U = FFLAS::fflas_new<Element>(n*lda);
		Element * L = FFLAS::fflas_new<Element>(n*n);

		/*  Create a random P,Q */

		for (size_t i = 0 ; i < n ; ++i)
			P[i] = i + RandInt(0UL,n-i);
		for (size_t i = 0 ; i < n ; ++i)
			Q[i] = i + RandInt(0UL,n-i);

		/*  det of P,Q */
		int d1 =1 ;
		for (size_t i = 0 ; i < n ; ++i)
			if (P[i] != i)
				d1 = -d1;
		for (size_t i = 0 ; i < n ; ++i)
			if (Q[i] != i)
				d1 = -d1;



		/*  Create L, lower det d */
		for (size_t i=0 ; i<n ; ++i)
			for (size_t j= 0; j<i ;++j)
				R.random( L[i*n+j] );

		Element dd = F.one;
		for (size_t i=0 ; i<n-1 ; ++i) {
			nzR.random( L[i*n+i] );
			F.mulin(dd,L[i*n+i]);
		}

		F.div(dd,d,dd);
		if (d1<0) F.negin(dd);
		L[n*n-1] = dd ;

		for (size_t i=0 ; i<n ; ++i)
			for (size_t j= i+1; j<n ;++j)
				F.init(L[i*n+j],0UL);


		/*  Create U, upper or rank r */
		for (size_t i=0 ; i<n ; ++i) {
			for (size_t j= 0; j<i ;++j)
				U[i*lda+j] = F.zero;
			U[i*lda+i] = F.one;
			for (size_t j= i+1; j<n ;++j)
				R.random( U[i*lda+j] );
		}

		/*  compute product */

		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				n,0,(int)n, U, lda, P);
		FFPACK::applyP (F, FFLAS::FflasLeft,  FFLAS::FflasNoTrans,
				n,0,(int)n, L, n, Q);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
			      n,n,n, 1.0, L,n, U,lda, 0.0, A,lda);
		//! @todo compute LU with ftrtr

		FFLAS::fflas_delete( P);
		FFLAS::fflas_delete( L);
		FFLAS::fflas_delete( U);
		FFLAS::fflas_delete( Q);

		return A;

	}


	template<typename Field>
	Givaro::Integer maxFieldElt() {return (Givaro::Integer)Field::getMaxModulus();}
	template<>
	Givaro::Integer maxFieldElt<Givaro::ZRing<Givaro::Integer>>() {return (Givaro::Integer)-1;}

	/*** Field chooser for test according to characteristic q and bitsize b ***/
	/* if q=-1 -> field is chosen randomly with a charateristic of b bits
	   if b=0 -> bitsize is chosen randomly according to maxFieldElt
	 */
	template<typename Field>
	Field* chooseField(Givaro::Integer q, unsigned long b){
		Givaro::Integer maxV= maxFieldElt<Field>();
		auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
		std::mt19937 mt_rand(seed);
		if (maxV>0 && (q> maxV || b> maxV.bitsize()))
			return nullptr;
		if (b<=1){
			//srand((double)std::chrono::high_resolution_clock::now());
			auto bitrand = std::bind(std::uniform_int_distribution<unsigned long>(2,maxV.bitsize()-1),
						 mt_rand);
			b = bitrand();
		}
		Givaro::IntPrimeDom IPD;
		Givaro::Integer tmp,p;
		if (q==-1){
			// Choose characteristic as a random prime of b bits
			do{
				Givaro::Integer _p;
				Givaro::Integer::seeding(Givaro::Integer(mt_rand()));
				Givaro::Integer::random_exact_2exp(_p,b);
				IPD.prevprime( tmp, _p+1 );
				p =  tmp;
			}while( (p < 2) );
		}
		else p=q;

		return new Field(p);
	}



} // FFPACK
#endif
