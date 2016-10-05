/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

/*! @file utils/fflas_randommatrix.h
 * @ingroup tests
 * @brief Utilities to create matrices with prescribed shapes, properties,...
 * To be used in benchmarks/tests
 */

#ifndef __FFLASFFPACK_randommatrix_H
#define __FFLASFFPACK_randommatrix_H

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/utils/debug.h"
#include <givaro/givinteger.h>
#include <givaro/givintprime.h>
#include <givaro/givranditer.h>
#include <chrono>
#include <random>

namespace FFPACK {

    /**! @brief  Random non-zero Matrix.
     * Creates a \c m x \c n matrix with random entries, and at least one of them is non zero.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param [out] A the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
	 * @param G a random iterator
     * @return \c A.
     */
    template<class Field, class RandIter>
    inline typename Field::Element_ptr
	NonZeroRandomMatrix(const Field & F, size_t m, size_t n, typename Field::Element_ptr A, size_t lda, RandIter& G) {
		bool ok=false;
		while (!ok)
			for (size_t i=0 ; i<m ; ++i)
				for (size_t j= 0; j<n ;++j)
					if (!F.isZero(G.random (A[i*lda+j])))
						ok = true;
		return A;
	}

	/**! @brief  Random non-zero Matrix.
     * Creates a \c m x \c n matrix with random entries, and at least one of them is non zero.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param [out] A the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @return \c A.
     */
    template<class Field, class RandIter>
    inline typename Field::Element_ptr
	NonZeroRandomMatrix(const Field & F, size_t m, size_t n,
						typename Field::Element_ptr A, size_t lda) {
		typename Field::RandIter G(F);
		return NonZeroRandomMatrix(F, m, n, A, lda, G);
	}

	/**! @brief  Random Matrix.
     * Creates a \c m x \c n matrix with random entries.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param [out] A the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
	 * @param G a random iterator
     * @return \c A.
     */
    template<class Field, class RandIter>
    inline typename Field::Element_ptr
	RandomMatrix(const Field & F, size_t m, size_t n, typename Field::Element_ptr A, size_t lda, RandIter& G) {
		for (size_t i=0 ; i<m ; ++i)
			for (size_t j= 0; j<n ;++j)
				G.random (A[i*lda+j]);
		return A;
	}

	/**! @brief  Random Matrix.
	 * Creates a \c m x \c n matrix with random entries.
	 * @param F field
	 * @param m number of rows in \p A
	 * @param n number of cols in \p A
	 * @param [out] A the matrix (preallocated to at least \c m x \c lda field elements)
	 * @param lda leading dimension of \p A
	 * @return \c A.
	 */
    template<class Field>
    inline typename Field::Element_ptr
	RandomMatrix(const Field & F, size_t m, size_t n, typename Field::Element_ptr A, size_t lda) {
		typename Field::RandIter G(F);
		return RandomMatrix (F, m, n, A, lda, G);
	}

	/*! Random integer in range.
     * @param a min bound
     * @param b max bound
     * @return a random integer in [a,b[  */
    inline size_t RandInt(size_t a, size_t b)
    {
        size_t x = a ;
        x += (size_t)rand()%(b-a);
        FFLASFFPACK_check(x<b && x>=a);
        return x ;
    }
} // FFPACK

#include "fflas-ffpack/ffpack/ffpack.h"

namespace FFPACK{
    /*! @brief  Random Matrix with prescribed rank.
     * Creates an \c m x \c n matrix with random entries and rank \c r.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param r rank of the matrix to build
     * @param A pointer to the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @return pointer to \c A.
     */
    template<class Field>
	inline typename Field::Element_ptr
	RandomMatrixWithRank (const Field & F, size_t m, size_t n, size_t r,
						  typename Field::Element_ptr A, size_t lda){
		typename Field::RandIter G(F);
		return RandomMatrixWithRank(F, m, n, r, A, lda);
	}

	/*! @brief  Random Matrix with prescribed rank.
     * Creates an \c m x \c n matrix with random entries and rank \c r.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param r rank of the matrix to build
     * @param A pointer to the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
	 * @param G a random iterator
     * @return pointer to \c A.
     */
    template<class Field, class RandIter>
	inline typename Field::Element_ptr
	RandomMatrixWithRank (const Field & F, size_t m, size_t n, size_t r,
						  typename Field::Element_ptr A, size_t lda, RandIter& G){
		FFLASFFPACK_check(r <= std::min(m,n));
		FFLASFFPACK_check(n <= lda);
		typedef typename Field::Element_ptr  Element_ptr;

		Givaro::GeneralRingNonZeroRandIter<Field,RandIter> nzG (G);

		size_t * P = FFLAS::fflas_new<size_t>(n);
		size_t * Q = FFLAS::fflas_new<size_t>(m);
		for (size_t i = 0 ; i < m ; ++i ) Q[i] = 0;
		for (size_t i = 0 ; i < n ; ++i ) P[i] = 0;

		Element_ptr U = FFLAS::fflas_new(F,m,n);
		Element_ptr L = FFLAS::fflas_new(F,m,m);

            /*  Create L, lower invertible */
		for (size_t i=0 ; i<m ; ++i){
			for (size_t j= 0; j<i ;++j) G.random( L[i*m+j] );
			nzG.random( L[i*m+i] );
			for (size_t j= i+1; j<m ;++j) F.init(L[i*m+j],F.zero);
		}

            /*  Create U, upper or rank r */
		for (size_t i=0 ; i<r ; ++i){
			for (size_t j= 0 ; j<i ;++j) F.init(U[i*n+j],0U);
			nzG.random( U[i*n+i] );
			for (size_t j= i+1; j<n ;++j) G.random( U[i*n+j] );
		}
		for (size_t i=r ; i<m ; ++i)
			for (size_t j= 0 ; j<n ;++j)
				F.init(U[i*n+j],F.zero);

            /*  Create a random P,Q */
		for (size_t i = 0 ; i < n ; ++i)
			P[i] = i + RandInt(0U,n-i);
		for (size_t i = 0 ; i < m ; ++i)
			Q[i] = i + RandInt(0U,m-i);

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

    inline void RandomRankProfile (size_t N, size_t R, size_t* rkp){
        size_t curr = 0;
        std::vector<bool> rows(N,false);
        while (curr<R){
            size_t i;
            while (rows [i = RandInt(0U, N)]);
            rows[i] = true;
            rkp [curr] = i;
            curr++;
        }
    }

/*! @brief  Random Matrix with prescribed rank and rank profile matrix
     * Creates an \c m x \c n matrix with random entries and rank \c r.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param r rank of the matrix to build
     * @param A pointer to the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
	 * @param RRP the R dimensional array with row positions of the rank profile matrix' pivots
	 * @param CRP the R dimensional array with column positions of the rank profile matrix' pivots
	 * @param G a random iterator
     * @return pointer to \c A.
     */
	template<class Field,class RandIter>
	inline typename Field::Element_ptr
	RandomMatrixWithRankandRPM (const Field& F, size_t M, size_t N, size_t R,
								typename Field::Element_ptr A, size_t lda,
								const size_t * RRP, const size_t * CRP, RandIter& G){

		Givaro::GeneralRingNonZeroRandIter<Field,RandIter> nzG(G);

		typename Field::Element_ptr L= FFLAS::fflas_new(F,M,N);

		FFLAS::pfzero(F, M, N, L, N);
			// Disabling the  parallel loop, as there is no way to declare G as SHARED in paladin
			//FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads> H;
			//SYNCH_GROUP (FOR1D(k, R, H,
		for (size_t k=0; k<R; ++k){
			size_t i = RRP[k];
			size_t j = CRP[k];
			nzG.random (L [i*N+j]);
			for (size_t l=i+1; l < M; ++l)
				G.random (L [l*N+j]);
		}
			//));

		typename Field::Element_ptr U= FFLAS::fflas_new(F,N,N);
		FFLAS::pfzero(F, N, N, U, N);
			//SYNCH_GROUP ( FOR1D(i, N, H,
		for (size_t i=0; i<N; ++i){
			nzG.random (U [i*N+i]);
			for (size_t j=i+1; j < N; ++j)
				G.random (U [i*N+j]);
		}
			//));

		typename Field::Element alpha, beta;
		F.init(alpha,1.0);
		F.init(beta,0.0);
                // auto sp=SPLITTER(); //CP: broken with Modular<Integer>. Need to reorganize  the helper behaviour with ParSeq and ModeTraits
            auto sp=NOSPLIT();
            FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M,N,N, alpha, L, N, U, N, beta, A, lda, sp);
            FFLAS::fflas_delete(L);
            FFLAS::fflas_delete(U);
			return A;
        }

/*! @brief  Random Matrix with prescribed rank and rank profile matrix
     * Creates an \c m x \c n matrix with random entries and rank \c r.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param r rank of the matrix to build
     * @param A pointer to the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
	 * @param RRP the R dimensional array with row positions of the rank profile matrix' pivots
	 * @param CRP the R dimensional array with column positions of the rank profile matrix' pivots
     * @return pointer to \c A.
     */
    template<class Field>
	inline typename Field::Element_ptr
	RandomMatrixWithRankandRPM (const Field& F, size_t M, size_t N, size_t R,
								typename Field::Element_ptr A, size_t lda,
								const size_t * RRP, const size_t * CRP){
		typename Field::RandIter G(F);
		return RandomMatrixWithRankandRPM (F, M, N, R, A, lda, RRP, CRP, G);
	}

	/*! @brief  Random Matrix with prescribed rank, with random  rank profile matrix
     * Creates an \c m x \c n matrix with random entries, rank \c r and with a 
	 * rank profile matrix chosen uniformly at random.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param r rank of the matrix to build
     * @param A pointer to the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @return pointer to \c A.
     */
    template<class Field, class RandIter>
	inline typename Field::Element_ptr
	RandomMatrixWithRankandRandomRPM (const Field& F, size_t M, size_t N, size_t R,
									  typename Field::Element_ptr A, size_t lda, RandIter& G){
            // generate the r pivots in the rank profile matrix E
            size_t pivot_r[R];
            size_t pivot_c[R];

            RandomRankProfile (M, R, pivot_r);
            RandomRankProfile (N, R, pivot_c);
            return RandomMatrixWithRankandRPM (F, M, N, R, A, lda, pivot_r, pivot_c, G);
        }

    /*! @brief  Random Matrix with prescribed rank, with random  rank profile matrix
     * Creates an \c m x \c n matrix with random entries, rank \c r and with a 
	 * rank profile matrix chosen uniformly at random.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param r rank of the matrix to build
     * @param A pointer to the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @return pointer to \c A.
     */
    template<class Field>
	inline typename Field::Element_ptr
	RandomMatrixWithRankandRandomRPM (const Field& F, size_t M, size_t N, size_t R,
									  typename Field::Element_ptr A, size_t lda){
		typename Field::RandIter G(F);
		return RandomMatrixWithRankandRandomRPM (F, M, N, R, A, lda, G);
	}

	/*! @brief  Random Matrix with prescribed det.
     * Creates a \c m x \c n matrix with random entries and rank \c r.
     * @param F field
	 * @param d the prescribed value for the determinant of A
     * @param n number of cols in \p A
     * @param A pointer to the matrix to be generated (preallocated to at least \c n x \c lda field elements)
     * @param lda leading dimension of \p A
     * @return pointer to \c A.
     */
    template<class Field>
	inline typename Field::Element_ptr
	RandomMatrixWithDet(const Field & F, size_t n, const typename Field::Element d,
						typename Field::Element_ptr A, size_t lda) {
		typename Field::RandIter G(F);
		return RandomMatrixWithDet (F, n, d, A, lda, G);
	}
	/*! @brief  Random Matrix with prescribed det.
     * Creates a \c m x \c n matrix with random entries and rank \c r.
     * @param F field
	 * @param d the prescribed value for the determinant of A
     * @param n number of cols in \p A
     * @param A pointer to the matrix to be generated (preallocated to at least \c n x \c lda field elements)
     * @param lda leading dimension of \p A
     * @return pointer to \c A.
     */
	template<class Field, class RandIter>
	inline typename Field::Element_ptr
	RandomMatrixWithDet(const Field & F, size_t n, const typename Field::Element d,
						typename Field::Element_ptr A, size_t lda, RandIter& G){
		FFLASFFPACK_check(n <= lda);
		typedef typename Field::Element  Element ;
		Givaro::GeneralRingNonZeroRandIter<Field,RandIter> nzG (G);

		size_t * P = FFLAS::fflas_new<size_t>(n);
		size_t * Q = FFLAS::fflas_new<size_t>(n);
		for (size_t i = 0 ; i < n ; ++i ) Q[i] = 0;
		for (size_t i = 0 ; i < n ; ++i ) P[i] = 0;

		Element * U = FFLAS::fflas_new<Element>(n*lda);
		Element * L = FFLAS::fflas_new<Element>(n*n);

            /*  Create a random P,Q */

		for (size_t i = 0 ; i < n ; ++i)
			P[i] = i + RandInt(0U,n-i);
		for (size_t i = 0 ; i < n ; ++i)
			Q[i] = i + RandInt(0U,n-i);

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
				G.random( L[i*n+j] );

		Element dd = F.one;
		for (size_t i=0 ; i<n-1 ; ++i) {
			nzG.random( L[i*n+i] );
			F.mulin(dd,L[i*n+i]);
		}

		F.div(dd,d,dd);
		if (d1<0) F.negin(dd);
		L[n*n-1] = dd ;

		for (size_t i=0 ; i<n ; ++i)
			for (size_t j= i+1; j<n ;++j)
				F.init(L[i*n+j],0U);


            /*  Create U, upper or rank r */
		for (size_t i=0 ; i<n ; ++i) {
			for (size_t j= 0; j<i ;++j)
				U[i*lda+j] = F.zero;
			U[i*lda+i] = F.one;
			for (size_t j= i+1; j<n ;++j)
				G.random( U[i*lda+j] );
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
} // FFPACK
#endif
