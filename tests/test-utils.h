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
	typename Field::Element * RandomMatrixWithRank(const Field & F,
						       typename Field::Element * A,
						       size_t r,
						       size_t m, size_t n, size_t lda)
	{
		FFLASFFPACK_check(r <= std::min(m,n));
		FFLASFFPACK_check(n <= lda);
		typedef typename Field::RandIter Randiter ;
		typedef typename Field::Element  Element ;
		Randiter R(F);
		NonzeroRandIter<Field,Randiter> nzR(F,R);

		size_t * P = new size_t[n];
		size_t * Q = new size_t[m];
		for (size_t i = 0 ; i < m ; ++i ) Q[i] = 0;
		for (size_t i = 0 ; i < n ; ++i ) P[i] = 0;

		Element * U = new Element[m*lda];
		Element * L = new Element[m*m];


		/*  Create L, lower invertible */
		for (size_t i=0 ; i<m ; ++i)
			for (size_t j= 0; j<i ;++j)
				R.random( L[i*m+j] );

		for (size_t i=0 ; i<m ; ++i)
			nzR.random( L[i*m+i] );

		for (size_t i=0 ; i<m ; ++i)
			for (size_t j= i+1; j<m ;++j)
				F.init(L[i*m+j],0UL);


		/*  Create U, upper or rank r */
		for (size_t i=0 ; i<r ; ++i)
			for (size_t j= i+1; j<r ;++j)
				R.random( U[i*lda+j] );
		for (size_t i=0 ; i<r ; ++i)
			nzR.random( U[i*lda+i] );
		for (size_t i=0 ; i<r ; ++i)
			for (size_t j= 0 ; j<i ;++j)
				F.init(U[i*lda+j],0UL);

		for (size_t i=r ; i<m ; ++i)
			for (size_t j= 0 ; j<n ;++j)
				F.init(U[i*lda+j],0UL);

		for (size_t i=0 ; i<r ; ++i)
			for (size_t j= r ; j<n ;++j)
				R.random( U[i*lda+j] );

		/*  Create a random P,Q */

		for (size_t i = 0 ; i < n ; ++i)
			P[i] = i + RandInt(0UL,n-i);
		for (size_t i = 0 ; i < m ; ++i)
			Q[i] = i + RandInt(0UL,m-i);

		/*  compute product */

		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				m,0,(int)n, U, lda, P);
		FFPACK::applyP (F, FFLAS::FflasLeft,  FFLAS::FflasNoTrans,
				m,0,(int)m, L, m, Q);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
			      m,n,m, 1.0, L,m, U,lda, 0.0, A,lda);
		//! @todo compute LU with ftrtr

		delete[] P;
		delete[] L;
		delete[] U;
		delete[] Q;

		return A;

	}

} // FFPACK
#endif
