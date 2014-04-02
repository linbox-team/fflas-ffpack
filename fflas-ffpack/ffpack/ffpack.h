/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* ffpack.h
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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

/** @file ffpack.h
 * \brief Set of elimination based routines for dense linear algebra.
 * Matrices are supposed over finite prime field of characteristic less than 2^26.

*/

#ifndef __FFLASFFPACK_ffpack_H
#define __FFLASFFPACK_ffpack_H

#ifdef __FFLASFFPACK_USE_OPENMP
#include <omp.h>
#endif

#include "fflas-ffpack/fflas/fflas.h"

#include <list>
#include <vector>
#include <iostream> // std::cout
#include <algorithm>

// The use of the small size LQUP is currently disabled:
// need for a better handling of element base (double, float, generic) combined
// with different thresholds.
// TransPosed version has to be implemented too.
#ifndef __FFPACK_LUDIVINE_CUTOFF
#define __FFPACK_LUDIVINE_CUTOFF 0
#endif

#ifndef __FFPACK_CHARPOLY_THRESHOLD
#define __FFPACK_CHARPOLY_THRESHOLD 30
#endif
/** @brief <b>F</b>inite <b>F</b>ield <b>PACK</b>
 * Set of elimination based routines for dense linear algebra.
 *
 * This namespace enlarges the set of BLAS routines of the class FFLAS, with higher
 * level routines based on elimination.
 \ingroup ffpack
 */
namespace FFPACK  {

	// public:
	enum FFPACK_LUDIVINE_TAG
	{
		FfpackLQUP=1,
		FfpackSingular=2
	};

	enum FFPACK_CHARPOLY_TAG
	{
		FfpackLUK=1,
		FfpackKG=2,
		FfpackHybrid=3,
		FfpackKGFast=4,
		FfpackDanilevski=5,
		FfpackArithProg=6,
		FfpackKGFastG=7
	};

	class CharpolyFailed{};

	enum FFPACK_MINPOLY_TAG
	{
		FfpackDense=1,
		FfpackKGF=2
	};
	void RankProfilesFromPLUQ (size_t* RowRankProfile, size_t* ColumnRankProfile,
				   const size_t * P, const size_t * Q,
				   const size_t M, const size_t N, const size_t R);

	void LAPACKPerm2MathPerm (size_t * MathP, const size_t * LapackP,
				  const size_t N);

	void MathPerm2LAPACKPerm (size_t * LapackP, const size_t * MathP,
				  const size_t N);
	template <class Element>
	void applyS (Element* A, const size_t lda, const size_t width,
		     const size_t M2,
		     const size_t R1, const size_t R2,
		     const size_t R3, const size_t R4);

	template <class Element>
	void applyT (Element* A, const size_t lda, const size_t width,
		     const size_t N2,
		     const size_t R1, const size_t R2,
		     const size_t R3, const size_t R4);

	void composePermutationsP (size_t * MathP,
				   const size_t * P1,
				   const size_t * P2,
				   const size_t R, const size_t N);
	void composePermutationsQ (size_t * MathP,
				   const size_t * Q1,
				   const size_t * Q2,
				   const size_t R, const size_t N);

	void cyclic_shift_mathPerm (size_t * P,  const size_t s);
	template<typename Base_t>
	void cyclic_shift_row_col(Base_t * A, size_t m, size_t n, size_t lda);
	template<typename Base_t>
	void cyclic_shift_row(Base_t * A, size_t m, size_t n, size_t lda);
	template<typename Base_t>
	void cyclic_shift_col(Base_t * A, size_t m, size_t n, size_t lda);


	/** Apply a permutation submatrix of P (between ibeg and iend) to a matrix
	 * to (iend-ibeg) vectors of size M stored in A (as column for NoTrans
	 * and rows for Trans).
	 * Side==FFLAS::FflasLeft for row permutation Side==FFLAS::FflasRight for a column
	 * permutation
	 * Trans==FFLAS::FflasTrans for the inverse permutation of P
	 * @param F
	 * @param Side
	 * @param Trans
	 * @param M
	 * @param ibeg
	 * @param iend
	 * @param A
	 * @param lda
	 * @param P
	 * @warning not sure the submatrix is still a permutation and the one we expect in all cases... examples for iend=2, ibeg=1 and P=[2,2,2]
	 */
	template<class Field>
	void
	applyP( const Field& F,
		const FFLAS::FFLAS_SIDE Side,
		const FFLAS::FFLAS_TRANSPOSE Trans,
		const size_t M, const int ibeg, const int iend,
		typename Field::Element * A, const size_t lda, const size_t * P )
	{

		if ( Side == FFLAS::FflasRight ) {
			typename Field::Element tmp;
			if ( Trans == FFLAS::FflasTrans )
				for (size_t j = 0 ; j < M ; ++j){
					for ( size_t i=(size_t)ibeg; i<(size_t) iend; ++i)
						if ( P[i]!= i ) {
							F.assign(tmp,A[j*lda+P[i]]);
							F.assign(A[j*lda+P[i]],A[j*lda+i]);
							F.assign(A[j*lda+i],tmp);
							// std::swap(A[j*lda+P[i]],A[j*lda+i]);
						}
					//FFLAS::fswap( F, M, A + P[i]*1, lda, A + i*1, lda );
				}
			else // Trans == FFLAS::FflasNoTrans
				for (size_t j = 0 ; j < M ; ++j){
					for (int i=iend; i-->ibeg; )
						if ( P[i]!=(size_t)i ) {
							F.assign(tmp,A[j*lda+P[i]]);
							F.assign(A[j*lda+P[i]],A[j*lda+(size_t)i]);
							F.assign(A[j*lda+(size_t)i],tmp);
							// std::swap(A[j*lda+P[i]],A[j*lda+(size_t)i]);
						}
					//FFLAS::fswap( F, M, A + P[i]*1, lda, A + i*1, lda );
				}
		}
		else { // Side == FFLAS::FflasLeft
			if ( Trans == FFLAS::FflasNoTrans )
				for (size_t i=(size_t)ibeg; i<(size_t)iend; ++i){
					if ( P[i]!= (size_t) i )
						FFLAS::fswap( F, M,
							      A + P[i]*lda, 1,
							      A + i*lda, 1 );
				}
			else // Trans == FFLAS::FflasTrans
				for (int i=iend; i-->ibeg; ){
					if ( P[i]!= (size_t) i ){
						FFLAS::fswap( F, M,
							      A + P[i]*lda, 1,
							      A + (size_t)i*lda, 1 );
					}
				}
		}

	}



#ifdef __FFLASFFPACK_USE_OPENMP


	// Parallel applyP with OPENMP tasks
	template<class Field>
	void
	papplyP( const Field& F,
		 const FFLAS::FFLAS_SIDE Side,
		 const FFLAS::FFLAS_TRANSPOSE Trans,
		 const size_t m, const int ibeg, const int iend,
		 typename Field::Element * A, const size_t lda, const size_t * P )
	{
		int numthreads = omp_get_max_threads();
		size_t BLOCKSIZE=std::max(2*m/numthreads,(size_t)1); // Assume that there is at least 2 ApplyP taking place in parallel
		size_t NBlocks = m/BLOCKSIZE;
		size_t LastBlockSize = m % BLOCKSIZE;
		if (LastBlockSize)
			NBlocks++;
		else
			LastBlockSize=BLOCKSIZE;
		for (size_t t = 0; t < NBlocks; ++t)
		{
			size_t BlockDim = BLOCKSIZE;
			if (t == NBlocks-1)
				BlockDim = LastBlockSize;
#pragma omp task shared (A, P, F) firstprivate(BlockDim)
			applyP(F, Side, Trans, BlockDim, ibeg, iend, A+BLOCKSIZE*t*((Side == FFLAS::FflasRight)?lda:1), lda, P);
		}
#pragma omp taskwait
	}

	// Parallel applyT with OPENMP tasks
	template <class Element>
	void papplyT (Element* A, const size_t lda, const size_t width,
		      const size_t N2,
		      const size_t R1, const size_t R2,
		      const size_t R3, const size_t R4)
	{
		int numthreads = omp_get_max_threads();
		size_t BLOCKSIZE=std::max(width/numthreads,(size_t)1);
		size_t NBlocks = width/BLOCKSIZE;
		size_t LastBlockSize = width % BLOCKSIZE;
		if (LastBlockSize)
			NBlocks++;
		else
			LastBlockSize=BLOCKSIZE;
		for (size_t t = 0; t < NBlocks; ++t)
		{
			size_t BlockDim = BLOCKSIZE;
			if (t == NBlocks-1)
				BlockDim = LastBlockSize;
#pragma omp task shared (A) firstprivate(BlockDim,R1,R2,R3,R4,N2)
			applyT (A+BLOCKSIZE*t*lda, lda, BlockDim, N2, R1, R2, R3, R4);
		}
#pragma omp taskwait
	}


	// Parallel applyS tasks with OPENMP tasks
	template <class Element>
	void papplyS (Element* A, const size_t lda, const size_t width,
		      const size_t M2,
		      const size_t R1, const size_t R2,
		      const size_t R3, const size_t R4)
	{
		int numthreads = omp_get_max_threads();
		size_t BLOCKSIZE=std::max(width/numthreads,(size_t)1);
		size_t NBlocks = width/BLOCKSIZE;
		size_t LastBlockSize = width % BLOCKSIZE;
		if (LastBlockSize)
			NBlocks++;
		else
			LastBlockSize=BLOCKSIZE;
		for (size_t t = 0; t < NBlocks; ++t)
		{
			size_t BlockDim = BLOCKSIZE;
			if (t == NBlocks-1)
				BlockDim = LastBlockSize;
#pragma omp task shared (A) firstprivate(BlockDim,R1,R2,R3,R4,M2)
			applyS(A+BLOCKSIZE*t, lda, BlockDim, M2, R1, R2, R3, R4);
		}
#pragma omp taskwait
	}


#endif


	/** Computes the rank of the given matrix using a LQUP factorization.
	 * The input matrix is modified.
	 * @param F field
	 * @param M row dimension of the matrix
	 * @param N column dimension of the matrix
	 * @param A input matrix
	 * @param lda leading dimension of A
	 */
	template <class Field>
	size_t
	Rank( const Field& F, const size_t M, const size_t N,
	      typename Field::Element * A, const size_t lda)
	{
		if (M == 0 and  N  == 0)
			return 0 ;

		size_t *P = new size_t[N];
		size_t *Q = new size_t[M];
		size_t R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, N,
				     A, lda, P, Q, FfpackLQUP);
		delete[] Q;
		delete[] P;
		return R;
	}

	/**  Returns true if the given matrix is singular.
	 * The method is a block elimination with early termination
	 *
	 * using LQUP factorization  with early termination.
	 * If <code>M != N</code>,
	 * then the matrix is virtually padded with zeros to make it square and
	 * it's determinant is zero.
	 * @warning The input matrix is modified.
	 * @param F field
	 * @param M row dimension of the matrix
	 * @param N column dimension of the matrix.
	 * @param [in,out] A input matrix
	 * @param lda leading dimension of A
	 */
	template <class Field>
	bool
	IsSingular( const Field& F, const size_t M, const size_t N,
		    typename Field::Element * A, const size_t lda)
	{
		if ( (M==0) and (N==0) )
			return  false ;
		if ( (M==0) or (N==0) )
			return  true ;
		if ( M != N )
			return  true ;


		size_t *P = new size_t[N];
		size_t *Q = new size_t[M];
		bool singular  = !LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, N,
					    A, lda, P, Q, FfpackSingular);

		delete[] P;
		delete[] Q;
		return singular;
	}

	/** @brief Returns the determinant of the given matrix.
	 * @details The method is a block elimination with early termination
	 * using LQUP factorization  with early termination.
	 * If <code>M != N</code>,
	 * then the matrix is virtually padded with zeros to make it square and
	 * it's determinant is zero.
	 * @warning The input matrix is modified.
	 * @param F field
	 * @param M row dimension of the matrix
	 * @param N column dimension of the matrix.
	 * @param [in,out] A input matrix
	 * @param lda leading dimension of A
	 */
	template <class Field>
	typename Field::Element
	Det( const Field& F, const size_t M, const size_t N,
	     typename Field::Element * A, const size_t lda)
	{
		if ( (M==0) and (N==0) )
			return  F.one ;
		if ( (M==0) or (N==0) )
			return  F.zero ;
		if ( M != N )
			return  F.zero ;

		typename Field::Element det; F.init(det);
		bool singular;
		size_t *P = new size_t[N];
		size_t *Q = new size_t[M];
		singular  = !LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans,  M, N,
				       A, lda, P, Q, FfpackSingular);
		if (singular){
			F.assign(det,F.zero);
			delete[] P;
			delete[] Q;
			return det;
		}
		else{
			F.assign(det,F.one);
			typename Field::Element *Ai=A;
			for (; Ai < A+ M*lda+N; Ai+=lda+1 )
				F.mulin( det, *Ai );
			int count=0;
			for (size_t i=0;i<N;++i)
				if (P[i] != i) ++count;

			if ((count&1) == 1)
				F.negin(det);
		}
		delete[] P;
		delete[] Q;
		return det;
	}

	// forward declaration
	template<class Field>
	void
	solveLB2( const Field& F, const FFLAS::FFLAS_SIDE Side,
		  const size_t M, const size_t N, const size_t R,
		  typename Field::Element * L, const size_t ldl,
		  const size_t * Q,
		  typename Field::Element * B, const size_t ldb ) ;



	/** Solve the system \f$A X = B\f$ or \f$X A = B\f$.
	 * Solving using the \c LQUP decomposition of \p A
	 * already computed inplace with \c LUdivine(FFLAS::FflasNoTrans, FFLAS::FflasNonUnit).
	 * Version for A square.
	 * If A is rank deficient, a solution is returned if the system is consistent,
	 * Otherwise an info is 1
	 *
	 * @param F field
	 * @param Side Determine wheter the resolution is left or right looking.
	 * @param M row dimension of \p B
	 * @param N col dimension of \p B
	 * @param R rank of \p A
	 * @param A input matrix
	 * @param lda leading dimension of \p A
	 * @param P column permutation of the \c LQUP decomposition of \p A
	 * @param Q column permutation of the \c LQUP decomposition of \p A
	 * @param B Right/Left hand side matrix. Initially stores \p B, finally stores the solution \p X.
	 * @param ldb leading dimension of \p B
	 * @param info Success of the computation: 0 if successfull, >0 if system is inconsistent
	 */
	template <class Field>
	void
	fgetrs (const Field& F,
		const FFLAS::FFLAS_SIDE Side,
		const size_t M, const size_t N, const size_t R,
		typename Field::Element *A, const size_t lda,
		const size_t *P, const size_t *Q,
		typename Field::Element *B, const size_t ldb,
		int * info)
	{

		*info =0;
		if (Side == FFLAS::FflasLeft) { // Left looking solve A X = B

			solveLB2 (F, FFLAS::FflasLeft, M, N, R, A, lda, Q, B, ldb);

			applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
				N, 0,(int) R, B, ldb, Q);

			bool consistent = true;
			for (size_t i = R; i < M; ++i)
				for (size_t j = 0; j < N; ++j)
					if (!F.isZero (*(B + i*ldb + j)))
						consistent = false;
			if (!consistent) {
				std::cerr<<"System is inconsistent"<<std::endl;
				*info = 1;
			}
			// The last rows of B are now supposed to be 0
#if 0
			for (size_t i = R; i < M; ++i)
				for (size_t j = 0; j < N; ++j)
					*(B + i*ldb + j) = F.zero;
#endif

			ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
			       R, N, F.one, A, lda , B, ldb);

			applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
				N, 0,(int) R, B, ldb, P);

		}
		else { // Right Looking X A = B

			applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
				M, 0,(int) R, B, ldb, P);

			ftrsm (F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
			       M, R, F.one, A, lda , B, ldb);

			fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M, N-R, R, F.one,
			       B, ldb, A+R, lda, F.mOne, B+R, ldb);

			bool consistent = true;
			for (size_t i = 0; i < M; ++i)
				for (size_t j = R; j < N; ++j)
					if (!F.isZero (*(B + i*ldb + j)))
						consistent = false;
			if (!consistent) {
				std::cerr<<"System is inconsistent"<<std::endl;
				*info = 1;
			}
			// The last cols of B are now supposed to be 0

			applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				M, 0,(int) R, B, ldb, Q);

			solveLB2 (F, FFLAS::FflasRight, M, N, R, A, lda, Q, B, ldb);
		}
	}

	/** Solve the system A X = B or X A = B.
	 * Solving using the LQUP decomposition of A
	 * already computed inplace with LUdivine(FFLAS::FflasNoTrans, FFLAS::FflasNonUnit).
	 * Version for A rectangular.
	 * If A is rank deficient, a solution is returned if the system is consistent,
	 * Otherwise an info is 1
	 *
	 * @param F field
	 * @param Side Determine wheter the resolution is left or right looking.
	 * @param M row dimension of A
	 * @param N col dimension of A
	 * @param NRHS number of columns (if Side = FFLAS::FflasLeft) or row (if Side = FFLAS::FflasRight) of the matrices X and B
	 * @param R rank of A
	 * @param A input matrix
	 * @param lda leading dimension of A
	 * @param P column permutation of the LQUP decomposition of A
	 * @param Q column permutation of the LQUP decomposition of A
	 * @param X solution matrix
	 * @param ldx leading dimension of X
	 * @param B Right/Left hand side matrix.
	 * @param ldb leading dimension of B
	 * @param info Succes of the computation: 0 if successfull, >0 if system is inconsistent
	 */
	template <class Field>
	typename Field::Element *
	fgetrs (const Field& F,
		const FFLAS::FFLAS_SIDE Side,
		const size_t M, const size_t N, const size_t NRHS, const size_t R,
		typename Field::Element *A, const size_t lda,
		const size_t *P, const size_t *Q,
		typename Field::Element *X, const size_t ldx,
		const typename Field::Element *B, const size_t ldb,
		int * info)
	{

		*info =0;

		typename Field::Element* W;
		size_t ldw;

		if (Side == FFLAS::FflasLeft) { // Left looking solve A X = B

			// Initializing X to 0 (to be optimized)
			FFLAS::fzero(F,N,NRHS,X,ldx);
			// for (size_t i = 0; i <N; ++i)
			// for (size_t j=0; j< NRHS; ++j)
			// F.assign (*(X+i*ldx+j), F.zero);

			if (M > N){ // Cannot copy B into X
				W = new typename Field::Element [M*NRHS];
				ldw = NRHS;
				FFLAS::fcopy(F,M,NRHS,W,ldw,B,ldb);
				// for (size_t i=0; i < M; ++i)
				// FFLAS::fcopy (F, NRHS, W + i*ldw, 1, B + i*ldb, 1);

				solveLB2 (F, FFLAS::FflasLeft, M, NRHS, R, A, lda, Q, W, ldw);

				applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
					NRHS, 0,(int) R, W, ldw, Q);

				bool consistent = true;
				for (size_t i = R; i < M; ++i)
					for (size_t j = 0; j < NRHS; ++j)
						if (!F.isZero (*(W + i*ldw + j)))
							consistent = false;
				if (!consistent) {
					std::cerr<<"System is inconsistent"<<std::endl;
					*info = 1;
					delete[] W;
					return X;
				}
				// Here the last rows of W are supposed to be 0

				ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
				       R, NRHS, F.one, A, lda , W, ldw);

				FFLAS::fcopy(F,R,NRHS,X,ldx,W,ldw);
				// for (size_t i=0; i < R; ++i)
				// FFLAS::fcopy (F, NRHS, X + i*ldx, 1, W + i*ldw, 1);

				delete[] W;
				applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
					NRHS, 0,(int) R, X, ldx, P);

			}
			else { // Copy B to X directly

				FFLAS::fcopy(F,M,NRHS,X,ldx,B,ldb);
				// for (size_t i=0; i < M; ++i)
				// FFLAS::fcopy (F, NRHS, X + i*ldx, 1, B + i*ldb, 1);

				solveLB2 (F, FFLAS::FflasLeft, M, NRHS, R, A, lda, Q, X, ldx);

				applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
					NRHS, 0,(int) R, X, ldx, Q);

				bool consistent = true;
				for (size_t i = R; i < M; ++i)
					for (size_t j = 0; j < NRHS; ++j)
						if (!F.isZero (*(X + i*ldx + j)))
							consistent = false;
				if (!consistent) {
					std::cerr<<"System is inconsistent"<<std::endl;
					*info = 1;
					return X;
				}
				// Here the last rows of W are supposed to be 0

				ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
				       R, NRHS, F.one, A, lda , X, ldx);

				applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
					NRHS, 0,(int) R, X, ldx, P);
			}
			return X;

		}
		else { // Right Looking X A = B

			FFLAS::fzero(F,NRHS,M,X,ldx);
			// for (size_t i = 0; i <NRHS; ++i)
			// for (size_t j=0; j< M; ++j)
			// F.assign (*(X+i*ldx+j), F.zero);

			if (M < N) {
				W = new typename Field::Element [NRHS*N];
				ldw = N;
				FFLAS::fcopy (F,NRHS, N, W, ldw, B, ldb);
				// for (size_t i=0; i < NRHS; ++i)
				// FFLAS::fcopy (F, N, W + i*ldw, 1, B + i*ldb, 1);

				applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
					NRHS, 0,(int) R, W, ldw, P);

				ftrsm (F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
				       NRHS, R, F.one, A, lda , W, ldw);

				fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, NRHS, N-R, R, F.one,
				       W, ldw, A+R, lda, F.mOne, W+R, ldw);

				bool consistent = true;
				for (size_t i = 0; i < NRHS; ++i)
					for (size_t j = R; j < N; ++j)
						if (!F.isZero (*(W + i*ldw + j)))
							consistent = false;
				if (!consistent) {
					std::cerr<<"System is inconsistent"<<std::endl;
					*info = 1;
					delete[] W;
					return X;
				}
				// The last N-R cols of W are now supposed to be 0
				FFLAS::fcopy (F, NRHS,R, X ,ldx,  W , ldb);
				// for (size_t i=0; i < NRHS; ++i)
				// FFLAS::fcopy (F, R, X + i*ldx, 1, W + i*ldb, 1);
				delete[] W;
				applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
					NRHS, 0,(int) R, X, ldx, Q);

				solveLB2 (F, FFLAS::FflasRight, NRHS, M, R, A, lda, Q, X, ldx);

			}
			else { // M >=N
				FFLAS::fcopy(F,NRHS,N,X,ldx,B,ldb);
				// for (size_t i=0; i < NRHS; ++i)
				// FFLAS::fcopy (F, N, X + i*ldx, 1, B + i*ldb, 1);

				applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
					NRHS, 0,(int) R, X, ldx, P);

				ftrsm (F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
				       NRHS, R, F.one, A, lda , X, ldx);

				fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, NRHS, N-R, R, F.one,
				       X, ldx, A+R, lda, F.mOne, X+R, ldx);

				bool consistent = true;
				for (size_t i = 0; i < NRHS; ++i)
					for (size_t j = R; j < N; ++j)
						if (!F.isZero (*(X + i*ldx + j)))
							consistent = false;
				if (!consistent) {
					std::cerr<<"System is inconsistent"<<std::endl;
					*info = 1;
					return X;
				}
				// The last N-R cols of W are now supposed to be 0

				applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
					NRHS, 0,(int) R, X, ldx, Q);

				solveLB2 (F, FFLAS::FflasRight, NRHS, M, R, A, lda, Q, X, ldx);

			}
			return X;
		}
	}

	/** @brief Square system solver
	 * @param F The computation domain
	 * @param Side Determine wheter the resolution is left or right looking
	 * @param M row dimension of B
	 * @param N col dimension of B
	 * @param A input matrix
	 * @param lda leading dimension of A
	 * @param B Right/Left hand side matrix. Initially contains B, finally contains the solution X.
	 * @param ldb leading dimension of B
	 * @param info Success of the computation: 0 if successfull, >0 if system is inconsistent
	 * @return the rank of the system
	 *
	 * Solve the system A X = B or X A = B.
	 * Version for A square.
	 * If A is rank deficient, a solution is returned if the system is consistent,
	 * Otherwise an info is 1
	 */
	template <class Field>
	size_t
	fgesv (const Field& F,
	       const FFLAS::FFLAS_SIDE Side,
	       const size_t M, const size_t N,
	       typename Field::Element *A, const size_t lda,
	       typename Field::Element *B, const size_t ldb,
	       int * info)
	{

		size_t Na;
		if (Side == FFLAS::FflasLeft)
			Na = M;
		else
			Na = N;

		size_t* P = new size_t[Na];
		size_t* Q = new size_t[Na];

		size_t R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, Na, Na, A, lda, P, Q, FfpackLQUP);

		fgetrs (F, Side, M, N, R, A, lda, P, Q, B, ldb, info);

		delete[] P;
		delete[] Q;

		return R;
	}

	/**  @brief Rectangular system solver
	 * @param F The computation domain
	 * @param Side Determine wheter the resolution is left or right looking
	 * @param M row dimension of A
	 * @param N col dimension of A
	 * @param NRHS number of columns (if Side = FFLAS::FflasLeft) or row (if Side = FFLAS::FflasRight) of the matrices X and B
	 * @param A input matrix
	 * @param lda leading dimension of A
	 * @param B Right/Left hand side matrix. Initially contains B, finally contains the solution X.
	 * @param ldb leading dimension of B
	 * @param X
	 * @param ldx
	 * @param info Success of the computation: 0 if successfull, >0 if system is inconsistent
	 * @return the rank of the system
	 *
	 * Solve the system A X = B or X A = B.
	 * Version for A square.
	 * If A is rank deficient, a solution is returned if the system is consistent,
	 * Otherwise an info is 1
	 */
	template <class Field>
	size_t
	fgesv (const Field& F,
	       const FFLAS::FFLAS_SIDE Side,
	       const size_t M, const size_t N, const size_t NRHS,
	       typename Field::Element *A, const size_t lda,
	       typename Field::Element *X, const size_t ldx,
	       const typename Field::Element *B, const size_t ldb,
	       int * info)
	{

		size_t* P = new size_t[N];
		size_t* Q = new size_t[M];

		size_t R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, N, A, lda, P, Q, FfpackLQUP);

		fgetrs (F, Side, M, N, NRHS, R, A, lda, P, Q, X, ldx, B, ldb, info);

		delete[] P;
		delete[] Q;

		return R;
	}

	/**  Solve the system Ax=b.
	 * Solving using LQUP factorization and
	 * two triangular system resolutions.
	 * The input matrix is modified.
	 * @param F The computation domain
	 * @param M row dimension of the matrix
	 * @param A input matrix
	 * @param lda leading dimension of A
	 * @param x solution vector
	 * @param incx increment of x
	 * @param b right hand side vector
	 * @param incb increment of b
	 */
	/// Solve linear system using LQUP factorization.
	template <class Field>
	typename Field::Element*
	Solve( const Field& F, const size_t M,
	       typename Field::Element * A, const size_t lda,
	       typename Field::Element * x, const int incx,
	       const typename Field::Element * b, const int incb )
	{

		size_t *P = new size_t[M];
		size_t *rowP = new size_t[M];

		if (LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, M, A, lda, P, rowP, FfpackLQUP) < M){
			std::cerr<<"SINGULAR MATRIX"<<std::endl;
			delete[] P;
			delete[] rowP;
			return x;
		}
		else{
			FFLAS::fcopy( F, M, x, incx, b, incb );

			ftrsv(F,  FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, M,
			      A, lda , x, incx);
			ftrsv(F,  FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, M,
			      A, lda , x, incx);
			applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
				M, 0,(int) M, x, incx, P );
			delete[] rowP;
			delete[] P;

			return x;

		}
	}



	/**  Computes a basis of the Left/Right nullspace of the matrix A.
	 * return the dimension of the nullspace.
	 *
	 * @param F The computation domain
	 * @param Side
	 * @param M
	 * @param N
	 * @param[in,out] A input matrix of dimension M x N, A is modified
	 * @param lda
	 * @param[out] NS output matrix of dimension N x NSdim (allocated here)
	 * @param[out] ldn
	 * @param[out] NSdim the dimension of the Nullspace (N-rank(A))
	 *
	 */
	template <class Field>
	size_t NullSpaceBasis (const Field& F, const FFLAS::FFLAS_SIDE Side,
			       const size_t M, const size_t N,
			       typename Field::Element* A, const size_t lda,
			       typename Field::Element*& NS, size_t& ldn,
			       size_t& NSdim)
	{
		if (Side == FFLAS::FflasRight) { // Right NullSpace
			size_t* P = new size_t[N];
			size_t* Qt = new size_t[M];

			size_t R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, N, A, lda, P, Qt);
			delete [] Qt;

			ldn = N-R;
			NSdim = ldn;

			if (NSdim == 0) {
				delete[] P;
				NS = NULL ;
				return NSdim ;
			}

			NS = new typename Field::Element [N*ldn];

			if (R == 0) {
				delete[] P;
				FFLAS::fidentity(F,N,ldn,NS,ldn);
				return NSdim;
			}

			FFLAS::fcopy (F, R, ldn, NS , ldn,  A + R,  lda );

			ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, R, ldn,
			       F.mOne, A, lda, NS, ldn);

			FFLAS::fidentity(F,NSdim,NSdim,NS+R*ldn,ldn);

			applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
				NSdim, 0,(int) R, NS, ldn, P);

			delete [] P;

			return NSdim;
		}
		else { // Left NullSpace
			size_t* P = new size_t[M];
			size_t* Qt = new size_t[N];

			size_t R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasTrans, M, N, A, lda, P, Qt);
			delete [] Qt;

			ldn = M;
			NSdim = M-R;

			if (NSdim == 0) {
				delete[] P;
				NS = NULL;
				return NSdim;
			}

			NS = new typename Field::Element [NSdim*ldn];


			if (R == 0) {
				delete[] P;
				FFLAS::fidentity(F,NSdim,ldn,NS,ldn);
				return NSdim;
			}


			FFLAS::fcopy (F, NSdim, R, NS, ldn, A + R *lda, lda);

			ftrsm (F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, NSdim, R,
			       F.mOne, A, lda, NS, ldn);

			FFLAS::fidentity(F,NSdim,NSdim,NS+R,ldn);

			applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				NSdim, 0,(int) R, NS, ldn, P);

			delete [] P;

			return NSdim;
		}
	}

	/**  Computes the row rank profile of A.
	 *
	 * @param F
	 * @param M
	 * @param N
	 * @param A input matrix of dimension M x N
	 * @param lda
	 * @param rkprofile return the rank profile as an array of row indexes, of dimension r=rank(A)
	 *
	 * rkprofile is allocated during the computation.
	 * @returns R
	 */
	template <class Field>
	size_t RowRankProfile (const Field& F, const size_t M, const size_t N,
			       typename Field::Element* A, const size_t lda,
			       size_t* &rkprofile)
	{
		size_t *P = new size_t[N];
		size_t *Q = new size_t[M];
		size_t R;

		R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, N, A, lda, P, Q);
		rkprofile = new size_t[R];

		for (size_t i=0; i<R; ++i)
			rkprofile[i] = Q[i];
		delete[] P;
		delete[] Q;
		return R;
	}

	/**  Computes the column rank profile of A.
	 *
	 * @param F
	 * @param M
	 * @param N
	 * @param A input matrix of dimension
	 * @param lda
	 * @param rkprofile return the rank profile as an array of row indexes, of dimension r=rank(A)
	 *
	 * A is modified
	 * rkprofile is allocated during the computation.
	 * @returns R
	 */
	template <class Field>
	size_t ColumnRankProfile (const Field& F, const size_t M, const size_t N,
				  typename Field::Element* A, const size_t lda,
				  size_t* &rkprofile)
	{
		size_t *P = new size_t[N];
		size_t *Q = new size_t[M];
		size_t R;

		R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasTrans, M, N, A, lda, P, Q);
		rkprofile = new size_t[R];

		for (size_t i=0; i<R; ++i)
			rkprofile[i] = Q[i];
		delete[] P;
		delete[] Q;
		return R;
	}

	/** RowRankProfileSubmatrixIndices.
	 * Computes the indices of the submatrix r*r X of A whose rows correspond to
	 * the row rank profile of A.
	 *
	 * @param F
	 * @param M
	 * @param N
	 * @param A input matrix of dimension
	 * @param rowindices array of the row indices of X in A
	 * @param colindices array of the col indices of X in A
	 * @param lda
	 * @param[out] R
	 *
	 * rowindices and colindices are allocated during the computation.
	 * A is modified
	 * @returns R
	 */
	template <class Field>
	size_t RowRankProfileSubmatrixIndices (const Field& F,
					       const size_t M, const size_t N,
					       typename Field::Element* A,
					       const size_t lda,
					       size_t*& rowindices,
					       size_t*& colindices,
					       size_t& R)
	{
		size_t *P = new size_t[N];
		size_t *Q = new size_t[M];

		R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, N, A, lda, P, Q);
		rowindices = new size_t[M];
		colindices = new size_t[N];
		for (size_t i=0; i<R; ++i){
			rowindices [i] = Q [i];
		}
		for (size_t i=0; i<N; ++i)
			colindices [i] = i;
		size_t tmp;
		for (size_t i=0; i<R; ++i){
			if (i != P[i]){
				tmp = colindices[i];
				colindices[i] = colindices[P[i]];
				colindices[P[i]] = tmp;
			}
		}

		delete[] P;
		delete[] Q;

		return R;
	}

	/** Computes the indices of the submatrix r*r X of A whose columns correspond to
	 * the column rank profile of A.
	 *
	 * @param F
	 * @param M
	 * @param N
	 * @param A input matrix of dimension
	 * @param rowindices array of the row indices of X in A
	 * @param colindices array of the col indices of X in A
	 * @param lda
	 * @param[out] R
	 *
	 * rowindices and colindices are allocated during the computation.
	 * @warning A is modified
	 * \return R
	 */
	template <class Field>
	size_t ColRankProfileSubmatrixIndices (const Field& F,
					       const size_t M, const size_t N,
					       typename Field::Element* A,
					       const size_t lda,
					       size_t*& rowindices,
					       size_t*& colindices,
					       size_t& R)
	{
		size_t *P = new size_t[M];
		size_t *Q = new size_t[N];

		R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasTrans, M, N, A, lda, P, Q);
		rowindices = new size_t[M];
		colindices = new size_t[N];
		for (size_t i=0; i<R; ++i)
			colindices [i] = Q [i];

		for (size_t i=0; i<N; ++i)
			rowindices [i] = i;

		size_t tmp;
		for (size_t i=0; i<R; ++i){
			if (i != P[i]){
				tmp = rowindices[i];
				rowindices[i] = rowindices[P[i]];
				rowindices[P[i]] = tmp;
			}
		}
		delete[] P;
		delete[] Q;

		return R;
	}

	/** Compute the r*r submatrix X of A, by picking the row rank profile rows of A.
	 *
	 * @param F
	 * @param M
	 * @param N
	 * @param A input matrix of dimension M x N
	 * @param lda
	 * @param X the output matrix
	 * @param[out] R
	 *
	 * A is not modified
	 * X is allocated during the computation.
	 * @return R
	 */
	template <class Field>
	size_t RowRankProfileSubmatrix (const Field& F,
					const size_t M, const size_t N,
					typename Field::Element* A,
					const size_t lda,
					typename Field::Element*& X, size_t& R)
	{

		size_t * rowindices, * colindices;
		typedef typename Field::Element Element ;

		Element * A2 = new Element[M*N] ;
		FFLAS::fcopy(F,M,N,A2,N,A,lda);

		RowRankProfileSubmatrixIndices (F, M, N, A2, N, rowindices, colindices, R);

		X = new Element[R*R];
		for (size_t i=0; i<R; ++i)
			for (size_t j=0; j<R; ++j)
				F.assign (*(X + i*R + j), *(A + rowindices[i]*lda + colindices[j]));
		delete[] A2;
		delete[] rowindices;
		delete[] colindices;
		return R;
	}

	/** Compute the \f$ r\times r\f$ submatrix X of A, by picking the row rank profile rows of A.
	 *
	 *
	 * @param F
	 * @param M
	 * @param N
	 * @param A input matrix of dimension M x N
	 * @param lda
	 * @param X the output matrix
	 * @param[out] R
	 *
	 * A is not modified
	 * X is allocated during the computation.
	 * \returns R
	 */
	template <class Field>
	size_t ColRankProfileSubmatrix (const Field& F, const size_t M, const size_t N,
					typename Field::Element* A, const size_t lda,
					typename Field::Element*& X, size_t& R)
	{

		size_t * rowindices, * colindices;
		typedef typename Field::Element Element ;

		Element * A2 = new Element[M*N];
		FFLAS::fcopy(F,M,N,A2,N,A,lda);

		ColRankProfileSubmatrixIndices (F, M, N, A2, N, rowindices, colindices, R);

		X = new Element[R*R];
		for (size_t i=0; i<R; ++i)
			for (size_t j=0; j<R; ++j)
				F.assign (*(X + i*R + j), *(A + rowindices[i]*lda + colindices[j]));
		delete[] A2;
		delete[] colindices;
		delete[] rowindices;
		return R;
	}

	/** LQUPtoInverseOfFullRankMinor.
	 * Suppose A has been factorized as L.Q.U.P, with rank r.
	 * Then Qt.A.Pt has an invertible leading principal r x r submatrix
	 * This procedure efficiently computes the inverse of this minor and puts it into X.
	 * @note It changes the lower entries of A_factors in the process (NB: unless A was nonsingular and square)
	 *
	 * @param F
	 * @param rank       rank of the matrix.
	 * @param A_factors  matrix containing the L and U entries of the factorization
	 * @param lda
	 * @param QtPointer  theLQUP->getQ()->getPointer() (note: getQ returns Qt!)
	 * @param X          desired location for output
	 * @param ldx
	 */
	template <class Field>
	typename Field::Element*
	LQUPtoInverseOfFullRankMinor( const Field& F, const size_t rank,
				      typename Field::Element * A_factors, const size_t lda,
				      const size_t* QtPointer,
				      typename Field::Element * X, const size_t ldx)
	{

		// upper entries are okay, just need to move up bottom ones
		const size_t* srcRow = QtPointer;
		for (size_t row=0; row<rank; row++, srcRow++)
			if (*srcRow != row) {
				typename Field::Element* oldRow = A_factors + (*srcRow) * lda;
				typename Field::Element* newRow = A_factors + row * lda;
				for (size_t col=0; col<row; col++, oldRow++, newRow++)
					F.assign(*newRow, *oldRow);
			}

		// X <- (Qt.L.Q)^(-1)
		//invL( F, rank, A_factors, lda, X, ldx);
		ftrtri (F, FFLAS::FflasLower, FFLAS::FflasUnit, rank, A_factors, lda);
		FFLAS::fcopy(F,rank,rank,A_factors,lda,X,ldx);
		// for (size_t i=0; i<rank; ++i)
		// FFLAS::fcopy (F, rank, A_factors+i*lda, 1, X+i*ldx,1);

		// X = U^-1.X
		ftrsm( F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans,
		       FFLAS::FflasNonUnit, rank, rank, F.one, A_factors, lda, X, ldx);

		return X;

	}

	//---------------------------------------------------------------------
	// TURBO: rank computation algorithm
	//---------------------------------------------------------------------
	template <class Field>
	size_t
	TURBO (const Field& F, const size_t M, const size_t N,
	       typename Field::Element* A, const size_t lda, size_t * P, size_t * Q, const size_t cutoff);

	template<class Field>
	size_t
	PLUQ (const Field& Fi, const FFLAS::FFLAS_DIAG Diag,
	      const size_t M, const size_t N,
	      typename Field::Element * A, const size_t lda,
	      size_t*P, size_t *Q);
	template<class Field>
	size_t
	PLUQ_basecase (const Field& Fi, const FFLAS::FFLAS_DIAG Diag,
		       const size_t M, const size_t N,
		       typename Field::Element * A, const size_t lda,
		       size_t*P, size_t *Q);


	/** @brief Compute the LQUP factorization of the given matrix.
	 * Using
	 * a block algorithm and return its rank.
	 * The permutations P and Q are represented
	 * using LAPACK's convention.
	 * @param F field
	 * @param Diag  precise whether U should have a unit diagonal or not
	 * @param trans UNKOWN TAG, probably the \c LU of \f$A^t\f$
	 * @param M matrix row dimension
	 * @param N matrix column dimension
	 * @param A input matrix
	 * @param lda leading dimension of \p A
	 * @param P the column permutation
	 * @param Qt the transpose of the row permutation \p Q
	 * @param LuTag flag for setting the earling termination if the matrix
	 * is singular
	 * @param cutoff UNKOWN TAG, probably a switch to a faster algo below \c cutoff
	 *
	 * @return the rank of \p A
	 * @bib
	 * - Jeannerod CP, <i>\c LSP Matrix Decomposition Revisited</i>, 2006
	 * - Pernet C, Brassel M <i>\c LUdivine, une divine factorisation \c LU</i>, 2002
	 * .
	 */
	template <class Field>
	size_t
	LUdivine (const Field& F, const FFLAS::FFLAS_DIAG Diag,  const FFLAS::FFLAS_TRANSPOSE trans,
		  const size_t M, const size_t N,
		  typename Field::Element * A, const size_t lda,
		  size_t* P, size_t* Qt
		  , const FFPACK_LUDIVINE_TAG LuTag=FfpackLQUP
		  , const size_t cutoff=__FFPACK_LUDIVINE_CUTOFF
		 );

	//! LUpdate
	template <class Field>
	size_t LUpdate (const Field& F,
			const FFLAS::FFLAS_DIAG Diag, const FFLAS::FFLAS_TRANSPOSE trans,
			const size_t M, const size_t N,
			typename Field::Element * A, const size_t lda,
			const size_t R,
			const size_t K,
			typename Field::Element * B, const size_t ldb,
			size_t*P, size_t *Q
			, const FFPACK::FFPACK_LUDIVINE_TAG LuTag  =FFPACK::FfpackLQUP
			, const size_t cutoff  =__FFPACK_LUDIVINE_CUTOFF
		       );

	template<class Element>
	class callLUdivine_small;

	//! LUdivine small case
	template <class Field>
	size_t
	LUdivine_small (const Field& F, const FFLAS::FFLAS_DIAG Diag,  const FFLAS::FFLAS_TRANSPOSE trans,
			const size_t M, const size_t N,
			typename Field::Element * A, const size_t lda,
			size_t* P, size_t* Q,
			const FFPACK_LUDIVINE_TAG LuTag=FfpackLQUP);

	//! LUdivine gauss
	template <class Field>
	size_t
	LUdivine_gauss (const Field& F, const FFLAS::FFLAS_DIAG Diag,
			const size_t M, const size_t N,
			typename Field::Element * A, const size_t lda,
			size_t* P, size_t* Q,
			const FFPACK_LUDIVINE_TAG LuTag=FfpackLQUP);




	/** Compute the inverse of a triangular matrix.
	 * @param F
	 * @param Uplo whether the matrix is upper of lower triangular
	 * @param Diag whether the matrix if unit diagonal
	 * @param N
	 * @param A
	 * @param lda
	 *
	 */
	template<class Field>
	void
	ftrtri (const Field& F, const FFLAS::FFLAS_UPLO Uplo, const FFLAS::FFLAS_DIAG Diag,
		const size_t N, typename Field::Element * A, const size_t lda)
	{
		if (N == 1){
			if (Diag == FFLAS::FflasNonUnit)
				F.invin (*A);
		}
		else {
			size_t N1 = N/2;
			size_t N2 = N - N1;
			ftrtri (F, Uplo, Diag, N1, A, lda);
			ftrtri (F, Uplo, Diag, N2, A + N1*(lda+1), lda);
			if (Uplo == FFLAS::FflasUpper){
				ftrmm (F, FFLAS::FflasLeft, Uplo, FFLAS::FflasNoTrans, Diag, N1, N2,
				       F.one, A, lda, A + N1, lda);
				ftrmm (F, FFLAS::FflasRight, Uplo, FFLAS::FflasNoTrans, Diag, N1, N2,
				       F.mOne, A + N1*(lda+1), lda, A + N1, lda);
			}
			else {
				ftrmm (F, FFLAS::FflasLeft, Uplo, FFLAS::FflasNoTrans, Diag, N2, N1,
				       F.one, A + N1*(lda+1), lda, A + N1*lda, lda);
				ftrmm (F, FFLAS::FflasRight, Uplo, FFLAS::FflasNoTrans, Diag, N2, N1,
				       F.mOne, A, lda, A + N1*lda, lda);
			}
		}
	}


	/**  Compute the product UL.
	 * Product UL of the upper, resp lower triangular matrices U and L
	 * stored one above the other in the square matrix A.
	 * Diag == Unit if the matrix U is unit diagonal
	 * @param F
	 * @param diag
	 * @param N
	 * @param A
	 * @param lda
	 *
	 */
	template<class Field>
	void
	ftrtrm (const Field& F, const FFLAS::FFLAS_DIAG diag, const size_t N,
		typename Field::Element * A, const size_t lda)
	{

		if (N == 1)
			return;
		size_t N1 = N/2;
		size_t N2 = N-N1;

		ftrtrm (F, diag, N1, A, lda);

		fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, N1, N1, N2, F.one,
		       A+N1, lda, A+N1*lda, lda, F.one, A, lda);

		ftrmm (F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans,
		       (diag == FFLAS::FflasUnit) ? FFLAS::FflasNonUnit : FFLAS::FflasUnit,
		       N1, N2, F.one, A + N1*(lda+1), lda, A + N1, lda);

		ftrmm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, diag, N2, N1,
		       F.one, A + N1*(lda+1), lda, A + N1*lda, lda);

		ftrtrm (F, diag, N2, A + N1*(lda+1), lda);

	}

	/** Extracts a triangular matrix from a compact storage A=L\U of rank R.
	 * row and column dimension of T are greater or equal to that of A and
	 * if UpLo = FflasLower, rowdim(T) = rowdim(A), else coldim(T) = coldim (A)
	 * @param F: base field
	 * @param UpLo: selects if the upper or lower triangular matrix is returned
	 * @param diag: selects if the triangular matrix unit-diagonal
	 * @param M: row dimension of T
	 * @param N: column dimension of T
	 * @param R: rank of the triangular matrix (how many rows/columns need to be copied)
	 * @param T: output matrix
	 * @param ldt: leading dimension of T
	 * @param A: input matrix
	 * @param lda: leading dimension of A
	 */
	template <class Field>
	void
	TriangularFromLU (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
			  const FFLAS::FFLAS_DIAG diag,
			  const size_t M, const size_t N, const size_t R,
			  typename Field::Element * T, const size_t ldt,
			  const typename Field::Element * A, const size_t lda){
		const typename Field::Element * Ai = A;
		typename Field::Element * Ti = T;
		if (Uplo == FFLAS::FflasUpper){
			for (size_t i=0; i<R; i++, Ai += lda, Ti += ldt){
				//!@todo just one triangular fzero+fcopy ?
				if (diag == FFLAS::FflasNonUnit){
					FFLAS::fzero(F,i,Ti,1);
					FFLAS::fcopy (F, N-i, Ti+i, 1, Ai+i, 1);
				}
				else {
					FFLAS::fzero(F,i,Ti,1);
					F.assign (*(Ti+i), F.one);
					FFLAS::fcopy (F, N-i-1, Ti+i+1, 1, Ai+i+1, 1);
				}
			}
			Ti = T+R*ldt;
			// FFLAS::fzero(F,M-R,N,Ti+(M-R)*ldt,ldt);
			for (size_t i=R; i<M; i++, Ti+=ldt) {
				FFLAS::fzero(F,N,Ti,1);
			}
		}
		else {
			for (size_t i=0; i<R; i++, Ai += lda, Ti += ldt){
				if (diag == FFLAS::FflasNonUnit){
					FFLAS::fcopy (F, i+1, Ti, 1, Ai, 1);
					// FFLAS::fzero(F,N-i-1,Ti+i+1,1);
					for (size_t j=i+1; j<N; j++)
						F.assign (Ti[j], F.zero);
				}
				else {
					FFLAS::fcopy (F, i, Ti, 1, Ai, 1);
					F.assign (Ti[i], F.one);
					// FFLAS::fzero(F,N-i-1,Ti+i+1,1);
					for (size_t j=i+1; j<N; j++)
						F.assign (Ti[j], F.zero);
				}
			}
			Ti = T+R*ldt;
			for (size_t i=R; i<M; i++, Ti+=ldt)
				FFLAS::fcopy(F, i, Ti, 1, Ai, 1);
			// FFLAS::fzero(F,N-R,Ti+R,1);
			for (size_t j=R; j<N; j++)
				F.assign (Ti[j], F.zero);
		}
	}

	/** Extracts an echelon matrix from a compact storage A=L\U of rank R.
	 * Either L or U is in Echelon form (depending on Uplo
	 * The echelon structure is defined by the first R values of the array P.
	 * row and column dimension of T are equal to that of A
	 * @param F: base field
	 * @param UpLo: selects if the upper or lower triangular matrix is returned
	 * @param diag: selects if the echelon matrix has unit pivots
	 * @param M: row dimension of T
	 * @param N: column dimension of T
	 * @param R: rank of the triangular matrix (how many rows/columns need to be copied)
	 * @param P: positions of the R pivots
	 * @param T: output matrix
	 * @param ldt: leading dimension of T
	 * @param A: input matrix
	 * @param lda: leading dimension of A
	 */
	template <class Field>
	void
	EchelonFromLU (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
		       const FFLAS::FFLAS_DIAG diag,
		       const size_t M, const size_t N, const size_t R, const size_t* P,
		       typename Field::Element * T, const size_t ldt,
		       const typename Field::Element * A, const size_t lda){

		const typename Field::Element * Ai = A;
		typename Field::Element * Ti = T;
		if (Uplo == FFLAS::FflasUpper){
			for (size_t i=0; i<R; i++, Ai += lda, Ti += ldt){
				size_t piv = P[i];
				if (diag == FFLAS::FflasNonUnit){
					FFLAS::fzero(F,piv,Ti,1);
					// for (size_t j=0; j<piv; j++)
					// F.assign (Ti[j], F.zero);
					FFLAS::fcopy (F, N-piv, Ti+piv, 1, Ai+piv, 1);
				}
				else {
					FFLAS::fzero(F,piv,Ti,1);
					// for (size_t j=0; j<piv; j++)
					// F.assign (Ti[j], F.zero);
					F.assign (Ti[piv], F.one);
					FFLAS::fcopy (F, N-piv-1, Ti+piv+1, 1, Ai+piv+1, 1);
				}
			}
			Ti = T+R*ldt;
			// FFLAS::fzero(F,M-R,N,Ti+R*ldt,ldt);
			for (size_t i=R; i<M; i++, Ti+=ldt)
				for (size_t j=0; j<N; j++)
					F.assign (Ti[j], F.zero);
		}
		else {
			for (size_t i=0; i<R; i++, Ai++, Ti++){
				size_t piv = P[i];
				if (diag == FFLAS::FflasNonUnit){
					FFLAS::fzero(F,piv,Ti,ldt);
					// for (size_t j=0; j<piv; j++)
					// F.assign (*(Ti+j*ldt), F.zero);
					FFLAS::fcopy (F, M-piv, Ti+piv*ldt, ldt, Ai+piv*lda, lda);
				}
				else {
					FFLAS::fzero(F,piv,Ti,ldt);
					// for (size_t j=0; j<piv; j++)
					// F.assign (*(Ti+j*ldt), F.zero);
					F.assign (*(Ti+piv*ldt), F.one);
					FFLAS::fcopy (F, M-piv-1, Ti+(piv+1)*ldt, ldt, Ai+(piv+1)*lda, lda);
				}
			}
			Ti = T+R;
			FFLAS::fzero(F,M,N-R,Ti,ldt);
			// for (size_t i=0; i<M; i++, Ti+=ldt)
			// for (size_t j=0; j<N-R; j++)
			// F.assign (Ti[j], F.zero);
		}
	}
	/*****************/
	/* ECHELON FORMS */
	/*****************/

	/** Compute the Column Echelon form of the input matrix in-place.
	 *
	 * After the computation A = [ M \ V ] such that AU = C is a column echelon
	 * decomposition of A, with U = P^T [   V    ] and C = M + Q [ Ir ]
	 *                                  [ 0 In-r ]               [ 0  ]
	 * Qt = Q^T
	 * If transform=false, the matrix U is not computed.
	 * See also test-colechelon for an example of use
	 * @param F
	 * @param M
	 * @param N
	 * @param A
	 * @param lda
	 * @param P
	 * @param Qt
	 * @param transform
	 */
	template <class Field>
	size_t
	ColumnEchelonForm (const Field& F, const size_t M, const size_t N,
			   typename Field::Element * A, const size_t lda,
			   size_t* P, size_t* Qt, bool transform = true);

	/**  Compute the Row Echelon form of the input matrix in-place.
	 *
	 * After the computation A = [ L \ M ] such that L A = R is a row echelon
	 * decomposition of A, with L =  [ L  0   ] P  and R = M + [Ir 0] Q^T
	 *                               [    In-r]
	 * Qt = Q^T
	 * If transform=false, the matrix L is not computed.
	 * See also test-rowechelon for an example of use
	 * @param F
	 * @param M
	 * @param N
	 * @param A
	 * @param lda
	 * @param P
	 * @param Qt
	 * @param transform
	 */
	template <class Field>
	size_t
	RowEchelonForm (const Field& F, const size_t M, const size_t N,
			typename Field::Element * A, const size_t lda,
			size_t* P, size_t* Qt, const bool transform = false);

	/** Compute the Reduced Column Echelon form of the input matrix in-place.
	 *
	 * After the computation A = [ V   ] such that AU = R is a reduced col echelon
	 *                           [ M 0 ]
	 * decomposition of A, where U = P^T [ V      ] and R = Q [ Ir   ]
	 *                                   [ 0 In-r ]           [ M  0 ]
	 * Qt = Q^T
	 * If transform=false, the matrix U is not computed and the matrix A = R
	 *
	 * @param F
	 * @param M
	 * @param N
	 * @param A
	 * @param lda
	 * @param P
	 * @param Qt
	 * @param transform
	 */
	template <class Field>
	size_t
	ReducedColumnEchelonForm (const Field& F, const size_t M, const size_t N,
				  typename Field::Element * A, const size_t lda,
				  size_t* P, size_t* Qt, const bool transform = true);

	/** Compute the Reduced Row Echelon form of the input matrix in-place.
	 *
	 * After the computation A = [ V1 M ] such that L A = R is a reduced row echelon
	 *                           [ V2 0 ]
	 * decomposition of A, where L =  [ V1  0   ] P and R =  [ Ir M  ] Q^T
	 *                                [ V2 In-r ]            [ 0     ]
	 * Qt = Q^T
	 * If transform=false, the matrix U is not computed and the matrix A = R
	 * @param F
	 * @param M
	 * @param N
	 * @param A
	 * @param lda
	 * @param P
	 * @param Qt
	 * @param transform
	 */
	template <class Field>
	size_t
	ReducedRowEchelonForm (const Field& F, const size_t M, const size_t N,
			       typename Field::Element * A, const size_t lda,
			       size_t* P, size_t* Qt, const bool transform = true);

	/**  Variant by the block recursive algorithm.
	 * (See A. Storjohann Thesis 2000)
	 * !!!!!! Warning !!!!!!
	 * This code is NOT WORKING properly for some echelon structures.
	 * This is due to a limitation of the way we represent permutation matrices
	 * (LAPACK's standard):
	 *  - a composition of transpositions Tij of the form
	 *    P = T_{1,j1} o T_{2,j2] o...oT_{r,jr}, with jk>k for all 0 < k <= r <= n
	 *  - The permutation on the columns, performed by this block recursive algorithm
	 *  cannot be represented with such a composition.
	 * Consequently this function should only be used for benchmarks
	 */
	template <class Field>
	size_t
	ReducedRowEchelonForm2 (const Field& F, const size_t M, const size_t N,
				typename Field::Element * A, const size_t lda,
				size_t* P, size_t* Qt, const bool transform = true){
		for (size_t i=0; i<N; ++i)
			Qt[i] = i;
		return REF (F, M, N, A, lda, 0, 0, N, P, Qt);

	}

	//! REF
	template <class Field>
	size_t
	REF (const Field& F, const size_t M, const size_t N,
	     typename Field::Element * A, const size_t lda,
	     const size_t colbeg, const size_t rowbeg, const size_t colsize,
	     size_t* Qt, size_t* P);

	/*****************/
	/*   INVERSION   */
	/*****************/
	/**  @brief Invert the given matrix in place
	 * or computes its nullity if it is singular.
	 *
	 * An inplace \f$2n^3\f$ algorithm is used.
	 * @param F The computation domain
	 * @param M order of the matrix
	 * @param [in,out] A input matrix (\f$M \times M\f$)
	 * @param lda leading dimension of A
	 * @param nullity dimension of the kernel of A
	 * @return pointer to \f$A\f$ and \f$A \gets A^{-1}\f$
	 */
	template <class Field>
	typename Field::Element*
	Invert (const Field& F, const size_t M,
		typename Field::Element * A, const size_t lda,
		int& nullity)
	{
		FFLASFFPACK_check(lda >= M);

		if (M == 0) {
			nullity = 0 ;
			return NULL ;
		}

		size_t * P = new size_t[M];
		size_t * Q = new size_t[M];
		size_t R =  ReducedColumnEchelonForm (F, M, M, A, lda, P, Q);
		nullity = (int)(M - R);
		applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
			M, 0, (int)R, A, lda, P);
		delete [] P;
		delete [] Q;
		return A;
	}

	/** @brief Invert the given matrix in place
	 * or computes its nullity if it is singular.
	 *
	 * @pre \p X is preallocated and should be large enough to store the
	 * \f$ m \times m\f$ matrix \p A.
	 *
	 * @param F The computation domain
	 * @param M order of the matrix
	 * @param [in] A input matrix (\f$M \times M\f$)
	 * @param lda leading dimension of \p A
	 * @param [out] X this is the inverse of \p A if \p A is invertible
	 * (non \c NULL and \f$ \mathtt{nullity} = 0\f$). It is untouched
	 * otherwise.
	 * @param ldx leading dimension of \p X
	 * @param nullity dimension of the kernel of \p A
	 * @return pointer to \f$X = A^{-1}\f$
	 */
	template <class Field>
	typename Field::Element*
	Invert (const Field& F, const size_t M,
		const typename Field::Element * A, const size_t lda,
		typename Field::Element * X, const size_t ldx,
		int& nullity)
	{
		FFLASFFPACK_check(lda >= M);
		FFLASFFPACK_check(ldx >= M);
		if (M == 0) {
			nullity = 0 ;
			return NULL ;
		}


		FFLAS::fcopy(F,M,M,X,ldx,A,lda);
		Invert (F,  M, X, lda, nullity);
		return X;
	}

	/** @brief Invert the given matrix or computes its nullity if it is singular.
	 *
	 * An \f$2n^3f\f$ algorithm is used.
	 * This routine can be \% faster than FFPACK::Invert but is not totally inplace.
	 *
	 * @pre \p X is preallocated and should be large enough to store the
	 * \f$ m \times m\f$ matrix \p A.
	 *
	 * @warning A is overwritten here !
	 * @bug not tested.
	 * @param F
	 * @param M order of the matrix
	 * @param [in,out] A input matrix (\f$M \times M\f$). On output, \p A
	 * is modified and represents a "psycological" factorisation \c LU.
	 * @param lda leading dimension of A
	 * @param [out] X this is the inverse of \p A if \p A is invertible
	 * (non \c NULL and \f$ \mathtt{nullity} = 0\f$). It is untouched
	 * otherwise.
	 * @param ldx leading dimension of \p X
	 * @param nullity dimension of the kernel of \p A
	 * @return pointer to \f$X = A^{-1}\f$
	 */
	template <class Field>
	typename Field::Element*
	Invert2( const Field& F, const size_t M,
		 typename Field::Element * A, const size_t lda,
		 typename Field::Element * X, const size_t ldx,
		 int& nullity)
	{
		FFLASFFPACK_check(lda >= M);
		FFLASFFPACK_check(ldx >= M);

		if (M == 0) {
			nullity = 0 ;
			return NULL ;
		}

		size_t *P = new size_t[M];
		size_t *rowP = new size_t[M];


		nullity = int(M - LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, M, A, lda, P, rowP, FfpackLQUP));

		if (nullity > 0){
			delete[] P;
			delete[] rowP;
			return NULL;
		}
		else {
			// Initializing X to 0
#if 0/*  timer remnants */
			t1.clear();
			t1.start();
#endif
			//! @todo this init is not all necessary (done after ftrtri)
			FFLAS::fzero(F,M,M,X,ldx);

			// X = L^-1 in n^3/3
			ftrtri (F, FFLAS::FflasLower, FFLAS::FflasUnit, M, A, lda);
			for (size_t i=0; i<M; ++i){
				for (size_t j=i; j<M; ++j)
					F.assign(*(X +i*ldx+j), F.zero);
				F.assign (*(X+i*(ldx+1)), F.one);
			}
			for (size_t i=1; i<M; ++i)
				FFLAS::fcopy (F, i, (X+i*ldx), 1, (A+i*lda), 1);

			ftrsm( F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
			       M, M, F.one, A, lda , X, ldx);

			// X = P^-1.X
			applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans,
				M, 0,(int) M, X, ldx, P );

			delete[] P;
			delete[] rowP;
			return X;
		}
	}


	/*****************************/
	/* CHARACTERISTIC POLYNOMIAL */
	/*****************************/


	/**
	 * Compute the characteristic polynomial of A using Krylov
	 * Method, and LUP factorization of the Krylov matrix
	 */
	template <class Field, class Polynomial>
	std::list<Polynomial>&
	CharPoly( const Field& F, std::list<Polynomial>& charp, const size_t N,
		  typename Field::Element * A, const size_t lda,
		  const FFPACK_CHARPOLY_TAG CharpTag= FfpackArithProg);

	template<class Polynomial, class Field>
	Polynomial & mulpoly(const Field& F, Polynomial &res, const Polynomial & P1, const Polynomial & P2)
	{
		size_t i,j;
		// Warning: assumes that res is allocated to the size of the product
		res.resize(P1.size()+P2.size()-1);
		FFLAS::fzero(F,res.size(),&res[0],1);
		for ( i=0;i<P1.size();i++)
			for ( j=0;j<P2.size();j++)
				F.axpyin(res[i+j],P1[i],P2[j]);
		return res;
	}

	template <class Field, class Polynomial>
	Polynomial&
	CharPoly( const Field& F, Polynomial& charp, const size_t N,
		  typename Field::Element * A, const size_t lda,
		  const FFPACK_CHARPOLY_TAG CharpTag= FfpackArithProg)
	{

		std::list<Polynomial> factor_list;
		CharPoly (F, factor_list, N, A, lda, CharpTag);
		typename std::list<Polynomial >::const_iterator it;
		it = factor_list.begin();

		charp.resize(N+1);

		Polynomial P = charp = *(it++);

		while( it!=factor_list.end() ){
			mulpoly (F,charp, P, *it);
			P = charp;
			++it;
		}

		return charp;
	}

	/**********************/
	/* MINIMAL POLYNOMIAL */
	/**********************/

	/** Compute the minimal polynomial.
	 * Minpoly of (A,v) using an LUP
	 * factorization of the Krylov Base (v, Av, .., A^kv)
	 * U,X must be (n+1)*n
	 * U contains the Krylov matrix and X, its LSP factorization
	 */
	template <class Field, class Polynomial>
	Polynomial&
	MinPoly( const Field& F, Polynomial& minP, const size_t N,
		 const typename Field::Element *A, const size_t lda,
		 typename Field::Element* X, const size_t ldx, size_t* P,
		 const FFPACK::FFPACK_MINPOLY_TAG MinTag= FFPACK::FfpackDense,
		 const size_t kg_mc=0, const size_t kg_mb=0, const size_t kg_j=0 );


	//! Solve L X = B or X L = B in place.
	//! L is M*M if Side == FFLAS::FflasLeft and N*N if Side == FFLAS::FflasRight, B is M*N.
	//! Only the R non trivial column of L are stored in the M*R matrix L
	//! Requirement :  so that L could  be expanded in-place
	template<class Field>
	void
	solveLB( const Field& F, const FFLAS::FFLAS_SIDE Side,
		 const size_t M, const size_t N, const size_t R,
		 typename Field::Element * L, const size_t ldl,
		 const size_t * Q,
		 typename Field::Element * B, const size_t ldb )
	{

		size_t LM = (Side == FFLAS::FflasRight)?N:M;
		int i = (int)R ;
		for (; i--; ){ // much faster for
			if (  Q[i] > (size_t) i){
				//for (size_t j=0; j<=Q[i]; ++j)
				//F.init( *(L+Q[i]+j*ldl), 0 );
				//std::cerr<<"1 deplacement "<<i<<"<-->"<<Q[i]<<endl;
				FFLAS::fcopy( F, LM-Q[i]-1, L+Q[i]*(ldl+1)+ldl,ldl, L+(Q[i]+1)*ldl+i, ldl );
				for ( size_t j=Q[i]*ldl; j<LM*ldl; j+=ldl)
					F.assign( *(L+i+j), F.zero );
			}
		}
		ftrsm( F, Side, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, M, N, F.one, L, ldl , B, ldb);
		//write_field(F,std::cerr<<"dans solveLB "<<endl,L,N,N,ldl);
		// Undo the permutation of L
		for (size_t ii=0; ii<R; ++ii){
			if ( Q[ii] > (size_t) ii){
				//for (size_t j=0; j<=Q[ii]; ++j)
				//F.init( *(L+Q[ii]+j*ldl), 0 );
				FFLAS::fcopy( F, LM-Q[ii]-1, L+(Q[ii]+1)*ldl+ii, ldl, L+Q[ii]*(ldl+1)+ldl,ldl );
				for ( size_t j=Q[ii]*ldl; j<LM*ldl; j+=ldl)
					F.assign( *(L+Q[ii]+j), F.zero );
			}
		}
	}

	//! Solve L X = B in place.
	//! L is M*M or N*N, B is M*N.
	//! Only the R non trivial column of L are stored in the M*R matrix L
	template<class Field>
	void
	solveLB2( const Field& F, const FFLAS::FFLAS_SIDE Side,
		  const size_t M, const size_t N, const size_t R,
		  typename Field::Element * L, const size_t ldl,
		  const size_t * Q,
		  typename Field::Element * B, const size_t ldb )
	{
		typename Field::Element * Lcurr,* Rcurr,* Bcurr;
		size_t ib,  Ldim;
		int k;
		if ( Side == FFLAS::FflasLeft ){
			size_t j = 0;
			while ( j<R ) {
				ib = Q[j];
				k = (int)ib ;
				while ((j<R) && ( (int) Q[j] == k)  ) {k++;j++;}
				Ldim = (size_t)k-ib;
				Lcurr = L + j-Ldim + ib*ldl;
				Bcurr = B + ib*ldb;
				Rcurr = Lcurr + Ldim*ldl;

				ftrsm( F, Side, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, Ldim, N, F.one,
				       Lcurr, ldl , Bcurr, ldb );

				fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-(size_t)k, N, Ldim, F.mOne,
				       Rcurr , ldl, Bcurr, ldb, F.one, Bcurr+Ldim*ldb, ldb);
			}
		}
		else{ // Side == FFLAS::FflasRight
			int j=(int)R-1;
			while ( j >= 0 ) {
				ib = Q[j];
				k = (int) ib;
				while ( (j >= 0) &&  ( (int)Q[j] == k)  ) {--k;--j;}
				Ldim = ib-(size_t)k;
				Lcurr = L + j+1 + (k+1)*(int)ldl;
				Bcurr = B + ib+1;
				Rcurr = Lcurr + Ldim*ldl;

				fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M,  Ldim, N-ib-1, F.mOne,
				       Bcurr, ldb, Rcurr, ldl,  F.one, Bcurr-Ldim, ldb);

				ftrsm (F, Side, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, M, Ldim, F.one,
				       Lcurr, ldl , Bcurr-Ldim, ldb );
			}
		}
	}


	template<class Field>
	void trinv_left( const Field& F, const size_t N, const typename Field::Element * L, const size_t ldl,
			 typename Field::Element * X, const size_t ldx )
	{
		FFLAS::fcopy(F,N,N,X,ldx,L,ldl);
		// for (size_t i=0; i<N; ++i)
		// FFLAS::fcopy (F, N, X+i*ldx, 1, L+i*ldl, 1);
		ftrtri (F, FFLAS::FflasLower, FFLAS::FflasUnit, N, X, ldx);
		//invL(F,N,L,ldl,X,ldx);
	}

	template <class Field>
	size_t KrylovElim( const Field& F, const size_t M, const size_t N,
			   typename Field::Element * A, const size_t lda, size_t*P,
			   size_t *Q, const size_t deg, size_t *iterates, size_t * inviterates, const size_t maxit,size_t virt);

	template <class Field>
	size_t  SpecRankProfile (const Field& F, const size_t M, const size_t N,
				 typename Field::Element * A, const size_t lda, const size_t deg, size_t *rankProfile);
	template <class Field, class Polynomial>
	std::list<Polynomial>&
	CharpolyArithProg (const Field& F, std::list<Polynomial>& frobeniusForm,
			   const size_t N, typename Field::Element * A, const size_t lda, const size_t c);

	template <class Field>
	void CompressRows (Field& F, const size_t M,
			   typename Field::Element * A, const size_t lda,
			   typename Field::Element * tmp, const size_t ldtmp,
			   const size_t * d, const size_t nb_blocs);

	template <class Field>
	void CompressRowsQK (Field& F, const size_t M,
			     typename Field::Element * A, const size_t lda,
			     typename Field::Element * tmp, const size_t ldtmp,
			     const size_t * d,const size_t deg, const size_t nb_blocs);

	template <class Field>
	void DeCompressRows (Field& F, const size_t M, const size_t N,
			     typename Field::Element * A, const size_t lda,
			     typename Field::Element * tmp, const size_t ldtmp,
			     const size_t * d, const size_t nb_blocs);
	template <class Field>
	void DeCompressRowsQK (Field& F, const size_t M, const size_t N,
			       typename Field::Element * A, const size_t lda,
			       typename Field::Element * tmp, const size_t ldtmp,
			       const size_t * d, const size_t deg, const size_t nb_blocs);

	template <class Field>
	void CompressRowsQA (Field& F, const size_t M,
			     typename Field::Element * A, const size_t lda,
			     typename Field::Element * tmp, const size_t ldtmp,
			     const size_t * d, const size_t nb_blocs);
	template <class Field>
	void DeCompressRowsQA (Field& F, const size_t M, const size_t N,
			       typename Field::Element * A, const size_t lda,
			       typename Field::Element * tmp, const size_t ldtmp,
			       const size_t * d, const size_t nb_blocs);


	namespace Protected {


		// Subroutine for Keller-Gehrig charpoly algorithm
		// Compute the new d after a LSP ( d[i] can be zero )
		template<class Field>
		size_t
		newD( const Field& F, size_t * d, bool& KeepOn,
		      const size_t l, const size_t N,
		      typename Field::Element * X,
		      const size_t* Q,
		      std::vector<std::vector<typename Field::Element> >& minpt);

		template<class Field>
		size_t
		updateD(const Field& F, size_t * d, size_t k,
			std::vector<std::vector<typename Field::Element> >& minpt );

		//---------------------------------------------------------------------
		// RectangleCopyTURBO: Copy A to T, with respect to the row permutation
		//                     defined by the lsp factorization of located in
		//                     A-dist2pivot
		//---------------------------------------------------------------------
		template <class Field>
		void
		RectangleCopyTURBO( const Field& F, const size_t M, const size_t N,
				    const size_t dist2pivot, const size_t rank,
				    typename Field::Element * T, const size_t ldt,
				    const typename Field::Element * A, const size_t lda )
		{

			const typename Field::Element * Ai = A;
			typename Field::Element * T1i = T, T2i = T + rank*ldt;
			size_t x = dist2pivot;
			for (; Ai<A+M*lda; Ai+=lda){
				while ( F.isZero(*(Ai-x)) ) { // test if the pivot is 0
					FFLAS::fcopy( F, N, T2i, 1, Ai, 1);
					Ai += lda;
					T2i += ldt;
				}
				FFLAS::fcopy( F, N, T1i, 1, Ai, 1);
				T1i += ldt;
				x--;
			}
		}



		//---------------------------------------------------------------------
		// LUdivine_construct: (Specialisation of LUdivine)
		// LUP factorisation of X, the Krylov base matrix of A^t and v, in A.
		// X contains the nRowX first vectors v, vA, .., vA^{nRowX-1}
		// A contains the LUP factorisation of the nUsedRowX first row of X.
		// When all rows of X have been factorized in A, and rank is full,
		// then X is updated by the following scheme: X <= ( X; X.B ), where
		// B = A^2^i.
		// This enables to make use of Matrix multiplication, and stop computing
		// Krylov vector, when the rank is not longer full.
		// P is the permutation matrix stored in an array of indexes
		//---------------------------------------------------------------------

		template <class Field>
		size_t
		LUdivine_construct( const Field& F, const FFLAS::FFLAS_DIAG Diag,
				    const size_t M, const size_t N,
				    const typename Field::Element * A, const size_t lda,
				    typename Field::Element * X, const size_t ldx,
				    typename Field::Element * u, size_t* P,
				    bool computeX, const FFPACK_MINPOLY_TAG MinTag= FFPACK::FfpackDense
				    , const size_t kg_mc =0
				    , const size_t kg_mb =0
				    , const size_t kg_j  =0
				  );

		template <class Field, class Polynomial>
		std::list<Polynomial>&
		KellerGehrig( const Field& F, std::list<Polynomial>& charp, const size_t N,
			      const typename Field::Element * A, const size_t lda );

		template <class Field, class Polynomial>
		int
		KGFast ( const Field& F, std::list<Polynomial>& charp, const size_t N,
			 typename Field::Element * A, const size_t lda,
			 size_t * kg_mc, size_t* kg_mb, size_t* kg_j );

		template <class Field, class Polynomial>
		std::list<Polynomial>&
		KGFast_generalized (const Field& F, std::list<Polynomial>& charp,
				    const size_t N,
				    typename Field::Element * A, const size_t lda);


		template<class Field>
		void
		fgemv_kgf( const Field& F,  const size_t N,
			   const typename Field::Element * A, const size_t lda,
			   const typename Field::Element * X, const size_t incX,
			   typename Field::Element * Y, const size_t incY,
			   const size_t kg_mc, const size_t kg_mb, const size_t kg_j );

		template <class Field, class Polynomial>
		std::list<Polynomial>&
		LUKrylov( const Field& F, std::list<Polynomial>& charp, const size_t N,
			  typename Field::Element * A, const size_t lda,
			  typename Field::Element * U, const size_t ldu);

		template <class Field, class Polynomial>
		std::list<Polynomial>&
		Danilevski (const Field& F, std::list<Polynomial>& charp,
			    const size_t N, typename Field::Element * A, const size_t lda);

		template <class Field, class Polynomial>
		std::list<Polynomial>&
		LUKrylov_KGFast( const Field& F, std::list<Polynomial>& charp, const size_t N,
				 typename Field::Element * A, const size_t lda,
				 typename Field::Element * X, const size_t ldx);
	} // Protected
} // FFPACK

#include "ffpack_ludivine.inl"
#include "ffpack_minpoly.inl"
#include "ffpack_charpoly_kglu.inl"
#include "ffpack_charpoly_kgfast.inl"
#include "ffpack_charpoly_kgfastgeneralized.inl"
#include "ffpack_charpoly_danilevski.inl"
#include "ffpack_charpoly.inl"
#include "ffpack_krylovelim.inl"
#include "ffpack_frobenius.inl"
#include "ffpack_echelonforms.inl"
//#include "ffpack_pluq-copy.inl"
#include "ffpack_pluq.inl"
#endif // __FFLASFFPACK_ffpack_H

