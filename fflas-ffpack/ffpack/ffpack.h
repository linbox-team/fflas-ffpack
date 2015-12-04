/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* ffpack.h
 * Copyright (C) 2005 Clement Pernet
 *               2014 FFLAS-FFPACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 * Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

#include <fflas-ffpack/fflas-ffpack-config.h>

#ifdef __FFLASFFPACK_USE_OPENMP
#include <omp.h>
#endif

#include "fflas-ffpack/fflas/fflas.h"
//#include "parallel.h"
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
namespace FFPACK  { /* tags */

	enum FFPACK_LU_TAG
	{
		FfpackSlabRecursive = 1,
		FfpackTileRecursive = 2,
		FfpackSingular = 3
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

}
namespace FFPACK { /* Permutations */

	/*****************/
	/* PERMUTATIONS  */
	/*****************/


	void LAPACKPerm2MathPerm (size_t * MathP, const size_t * LapackP,
				  const size_t N);

	void MathPerm2LAPACKPerm (size_t * LapackP, const size_t * MathP,
				  const size_t N);

	template <class Field>
	void MatrixApplyS (const Field& F, typename Field::Element_ptr A, const size_t lda, const size_t width,
		     const size_t M2,
		     const size_t R1, const size_t R2,
		     const size_t R3, const size_t R4);

	template <class Element>
	void PermApplyS (Element* A, const size_t lda, const size_t width,
			 const size_t M2,
			 const size_t R1, const size_t R2,
			 const size_t R3, const size_t R4);

	template <class Field>
	void MatrixApplyT (const Field& F, typename Field::Element_ptr A, const size_t lda, const size_t width,
			   const size_t N2,
			   const size_t R1, const size_t R2,
			   const size_t R3, const size_t R4);

	template <class Element>
	void PermApplyT (Element* A, const size_t lda, const size_t width,
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
	template<class Field>
	void cyclic_shift_row(const Field& F, typename Field::Element_ptr A, size_t m, size_t n, size_t lda);
	template<class Field>
	void cyclic_shift_col(const Field& F, typename Field::Element_ptr A, size_t m, size_t n, size_t lda);


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
		const size_t M, const size_t ibeg, const size_t iend,
		typename Field::Element_ptr A, const size_t lda, const size_t * P );



//#ifdef __FFLASFFPACK_USE_OPENMP


	//! Parallel applyP with OPENMP tasks
	template<class Field>
	void
	papplyP( const Field& F,
		 const FFLAS::FFLAS_SIDE Side,
		 const FFLAS::FFLAS_TRANSPOSE Trans,
		 const size_t m, const size_t ibeg, const size_t iend,
		 typename Field::Element_ptr A, const size_t lda, const size_t * P );

	//! Parallel applyT with OPENMP tasks
	template <class Field>
	void pMatrixApplyT (const Field& F, typename Field::Element_ptr A, const size_t lda,
		      const size_t width, const size_t N2,
		      const size_t R1, const size_t R2,
		      const size_t R3, const size_t R4) ;


	//! Parallel applyS tasks with OPENMP tasks
	template <class Field>
	void pMatrixApplyS (const Field& F, typename Field::Element_ptr A, const size_t lda,
			    const size_t width, const size_t M2,
			    const size_t R1, const size_t R2,
			    const size_t R3, const size_t R4) ;

	template<class Field>
	size_t
	pPLUQ(const Field& Fi, const FFLAS::FFLAS_DIAG Diag,
	      const size_t M, const size_t N,
	      typename Field::Element_ptr A, const size_t lda,
	      size_t* P, size_t* Q, int nt);


//#endif

} // FFPACK permutations
// #include "ffpack_permutation.inl"

namespace FFPACK { /* fgetrs, fgesv */

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
		typename Field::Element_ptr A, const size_t lda,
		const size_t *P, const size_t *Q,
		typename Field::Element_ptr B, const size_t ldb,
		int * info);

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
	typename Field::Element_ptr
	fgetrs (const Field& F,
		const FFLAS::FFLAS_SIDE Side,
		const size_t M, const size_t N, const size_t NRHS, const size_t R,
		typename Field::Element_ptr A, const size_t lda,
		const size_t *P, const size_t *Q,
		typename Field::Element_ptr X, const size_t ldx,
		typename Field::ConstElement_ptr B, const size_t ldb,
		int * info);

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
	       typename Field::Element_ptr A, const size_t lda,
	       typename Field::Element_ptr B, const size_t ldb,
	       int * info);

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
	       typename Field::Element_ptr A, const size_t lda,
	       typename Field::Element_ptr X, const size_t ldx,
	       typename Field::ConstElement_ptr B, const size_t ldb,
	       int * info);

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

} // FFPACK fgesv, fgetrs
// #include "ffpack_fgesv.inl"
// #include "ffpack_fgetrs.inl"

namespace FFPACK { /* ftrtr */



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
		const size_t N, typename Field::Element_ptr A, const size_t lda);


	template<class Field>
	void trinv_left( const Field& F, const size_t N, typename Field::ConstElement_ptr L, const size_t ldl,
			 typename Field::Element_ptr X, const size_t ldx );

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
			typename Field::Element_ptr A, const size_t lda);

} // FFPACK ftrtr
// #include "ffpack_ftrtr.inl"

namespace FFPACK { /* PLUQ */

	/** @brief Compute the PLUQ factorization of the given matrix.
	 * Using a block algorithm and return its rank.
	 * The permutations P and Q are represented
	 * using LAPACK's convention.
	 * @param F field
	 * @param Diag   whether U should have a unit diagonal or not
	 * @param trans, \c LU of \f$A^t\f$
	 * @param M matrix row dimension
	 * @param N matrix column dimension
	 * @param A input matrix
	 * @param lda leading dimension of \p A
	 * @param P the row permutation
	 * @param Q the column permutation

	 * @return the rank of \p A
	 * @bib
	 * - Dumas J-G.,  Pernet C., and Sultan Z. <i>\c Simultaneous computation of the row and column rank profiles </i>, ISSAC'13, 2013
	 * .
	 */
	template<class Field>
	size_t
	PLUQ (const Field& F, const FFLAS::FFLAS_DIAG Diag,
	      const size_t M, const size_t N,
	      typename Field::Element_ptr A, const size_t lda,
	      size_t*P, size_t *Q);

} // FFPACK PLUQ
// #include "ffpack_pluq.inl"

namespace FFPACK { /* ludivine */

	/** @brief Compute the CUP factorization of the given matrix.
	 * Using
	 * a block algorithm and return its rank.
	 * The permutations P and Q are represented
	 * using LAPACK's convention.
	 * @param F field
	 * @param Diag  whether U should have a unit diagonal or not
	 * @param trans  \c LU of \f$A^t\f$
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
	 * - Jeannerod C-P, Pernet, C. and Storjohann, A. <i>\c Rank-profile revealing Gaussian elimination and the CUP matrix decomposition  </i>, J. of Symbolic Comp., 2013
	 * - Pernet C, Brassel M <i>\c LUdivine, une divine factorisation \c LU</i>, 2002
	 * .
	 */
	template <class Field>
	size_t
	LUdivine (const Field& F, const FFLAS::FFLAS_DIAG Diag,  const FFLAS::FFLAS_TRANSPOSE trans,
			  const size_t M, const size_t N,
			  typename Field::Element_ptr A, const size_t lda,
			  size_t* P, size_t* Qt,
			  const FFPACK_LU_TAG LuTag = FfpackSlabRecursive,
			  const size_t cutoff=__FFPACK_LUDIVINE_CUTOFF);

	template<class Element>
	class callLUdivine_small;

	//! LUdivine small case
	template <class Field>
	size_t
	LUdivine_small (const Field& F, const FFLAS::FFLAS_DIAG Diag,  const FFLAS::FFLAS_TRANSPOSE trans,
					const size_t M, const size_t N,
					typename Field::Element_ptr A, const size_t lda,
					size_t* P, size_t* Q,
					const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);

	//! LUdivine gauss
	template <class Field>
	size_t
	LUdivine_gauss (const Field& F, const FFLAS::FFLAS_DIAG Diag,
					const size_t M, const size_t N,
					typename Field::Element_ptr A, const size_t lda,
					size_t* P, size_t* Q,
					const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);

	namespace Protected {



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
				    typename Field::ConstElement_ptr A, const size_t lda,
				    typename Field::Element_ptr X, const size_t ldx,
				    typename Field::Element_ptr u, size_t* P,
				    bool computeX, const FFPACK_MINPOLY_TAG MinTag= FfpackDense
				    , const size_t kg_mc =0
				    , const size_t kg_mb =0
				    , const size_t kg_j  =0
				  );

	} // Protected

} //FFPACK ludivine, turbo
// #include "ffpack_ludivine.inl"

namespace FFPACK { /* echelon */
	/*****************/
	/* ECHELON FORMS */
	/*****************/

	/** Compute the Column Echelon form of the input matrix in-place.
	 *
	 * If LuTag == FfpackTileRecursive, then after the computation A = [ M \ V ]
	 * such that AU = C is a column echelon decomposition of A,
	 * with U = P^T [   V    ] and C = M + Q [ Ir ]
	 *              [ 0 In-r ]               [ 0  ]
	 * If LuTag == FfpackTileRecursive then A = [ N \ V ] such that the same holds with M = Q N
	 *
	 * Qt = Q^T
	 * If transform=false, the matrix V is not computed.
	 * See also test-colechelon for an example of use
	 * @param F
	 * @param M
	 * @param N
	 * @param A
	 * @param lda
	 * @param P the column permutation
	 * @param Qt the row position of the pivots in the echelon form
	 * @param transform
	 */
	template <class Field>
	size_t
	ColumnEchelonForm (const Field& F, const size_t M, const size_t N,
					   typename Field::Element_ptr A, const size_t lda,
					   size_t* P, size_t* Qt, bool transform = false,
					   const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);

	/**  Compute the Row Echelon form of the input matrix in-place.
	 *
	 * If LuTag == FfpackTileRecursive, then after the computation A = [ L \ M ]
	 * such that X A = R is a row echelon decomposition of A,
	 * with X =  [ L  0   ] P  and R = M + [Ir 0] Q^T
	 *           [    In-r]
	 * If LuTag == FfpackTileRecursive then A = [ L \ N ] such that the same holds with M =  N Q^T
	 * Qt = Q^T
	 * If transform=false, the matrix L is not computed.
	 * See also test-rowechelon for an example of use
	 * @param F
	 * @param M
	 * @param N
	 * @param A
	 * @param lda
	 * @param P the row permutation
	 * @param Qt the column position of the pivots in the echelon form
	 * @param transform
	 */
	template <class Field>
	size_t
	RowEchelonForm (const Field& F, const size_t M, const size_t N,
					typename Field::Element_ptr A, const size_t lda,
					size_t* P, size_t* Qt, const bool transform = false,
					const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);

	/** Compute the Reduced Column Echelon form of the input matrix in-place.
	 *
	 * After the computation A = [ V   ] such that AX = R is a reduced col echelon
	 *                           [ M 0 ]
	 * decomposition of A, where X = P^T [ V      ] and R = Q [ Ir   ]
	 *                                   [ 0 In-r ]           [ M  0 ]
	 * Qt = Q^T
	 * If transform=false, the matrix X is not computed and the matrix A = R
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
							  typename Field::Element_ptr A, const size_t lda,
							  size_t* P, size_t* Qt, const bool transform = false,
							  const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);

	/** Compute the Reduced Row Echelon form of the input matrix in-place.
	 *
	 * After the computation A = [ V1 M ] such that X A = R is a reduced row echelon
	 *                           [ V2 0 ]
	 * decomposition of A, where X =  [ V1  0   ] P and R =  [ Ir M  ] Q^T
	 *                                [ V2 In-r ]            [ 0     ]
	 * Qt = Q^T
	 * If transform=false, the matrix X is not computed and the matrix A = R
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
			       typename Field::Element_ptr A, const size_t lda,
						   size_t* P, size_t* Qt, const bool transform = false,
						   const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);

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
				typename Field::Element_ptr A, const size_t lda,
				size_t* P, size_t* Qt, const bool transform = true);

	//! REF
	template <class Field>
	size_t
	REF (const Field& F, const size_t M, const size_t N,
	     typename Field::Element_ptr A, const size_t lda,
	     const size_t colbeg, const size_t rowbeg, const size_t colsize,
	     size_t* Qt, size_t* P);

} // FFPACK
// #include "ffpack_echelonforms.inl"

namespace FFPACK { /* invert */
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
	typename Field::Element_ptr
	Invert (const Field& F, const size_t M,
		typename Field::Element_ptr A, const size_t lda,
		int& nullity);

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
	typename Field::Element_ptr
	Invert (const Field& F, const size_t M,
			typename Field::ConstElement_ptr A, const size_t lda,
			typename Field::Element_ptr X, const size_t ldx,
			int& nullity);

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
	typename Field::Element_ptr
	Invert2( const Field& F, const size_t M,
			 typename Field::Element_ptr A, const size_t lda,
			 typename Field::Element_ptr X, const size_t ldx,
			 int& nullity);

} // FFPACK invert
// #include "ffpack_invert.inl"

namespace FFPACK { /* charpoly */
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
			  typename Field::Element_ptr A, const size_t lda,
			  const FFPACK_CHARPOLY_TAG CharpTag= FfpackArithProg);

	template<class Polynomial, class Field>
	Polynomial & mulpoly(const Field& F, Polynomial &res, const Polynomial & P1, const Polynomial & P2);

	template <class Field, class Polynomial>
	Polynomial&
	CharPoly( const Field& F, Polynomial& charp, const size_t N,
		  typename Field::Element_ptr A, const size_t lda,
		  const FFPACK_CHARPOLY_TAG CharpTag= FfpackArithProg);


	namespace Protected {
		template <class Field, class Polynomial>
		std::list<Polynomial>&
		KellerGehrig( const Field& F, std::list<Polynomial>& charp, const size_t N,
			      typename Field::ConstElement_ptr A, const size_t lda );

		template <class Field, class Polynomial>
		int
		KGFast ( const Field& F, std::list<Polynomial>& charp, const size_t N,
			 typename Field::Element_ptr A, const size_t lda,
			 size_t * kg_mc, size_t* kg_mb, size_t* kg_j );

		template <class Field, class Polynomial>
		std::list<Polynomial>&
		KGFast_generalized (const Field& F, std::list<Polynomial>& charp,
				    const size_t N,
				    typename Field::Element_ptr A, const size_t lda);


		template<class Field>
		void
		fgemv_kgf( const Field& F,  const size_t N,
			   typename Field::ConstElement_ptr A, const size_t lda,
			   typename Field::ConstElement_ptr X, const size_t incX,
			   typename Field::Element_ptr Y, const size_t incY,
			   const size_t kg_mc, const size_t kg_mb, const size_t kg_j );

		template <class Field, class Polynomial>
		std::list<Polynomial>&
		LUKrylov( const Field& F, std::list<Polynomial>& charp, const size_t N,
			  typename Field::Element_ptr A, const size_t lda,
			  typename Field::Element_ptr U, const size_t ldu);

		template <class Field, class Polynomial>
		std::list<Polynomial>&
		Danilevski (const Field& F, std::list<Polynomial>& charp,
			    const size_t N, typename Field::Element_ptr A, const size_t lda);

		template <class Field, class Polynomial>
		std::list<Polynomial>&
		LUKrylov_KGFast( const Field& F, std::list<Polynomial>& charp, const size_t N,
				 typename Field::Element_ptr A, const size_t lda,
				 typename Field::Element_ptr X, const size_t ldx);
	} // Protected
} // FFPACK charpoly
// #include "ffpack_charpoly_kglu.inl"
// #include "ffpack_charpoly_kgfast.inl"
// #include "ffpack_charpoly_kgfastgeneralized.inl"
// #include "ffpack_charpoly_danilevski.inl"
// #include "ffpack_charpoly.inl"

namespace FFPACK { /* frobenius, charpoly */

	template <class Field, class Polynomial>
	std::list<Polynomial>&
	CharpolyArithProg (const Field& F, std::list<Polynomial>& frobeniusForm,
			   const size_t N, typename Field::Element_ptr A, const size_t lda, const size_t c);


} // FFPACK frobenius
// #include "ffpack_frobenius.inl"

namespace FFPACK { /* minpoly */


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
		 typename Field::ConstElement_ptr A, const size_t lda,
		 typename Field::Element_ptr X, const size_t ldx, size_t* P,
		 const FFPACK_MINPOLY_TAG MinTag= FFPACK::FfpackDense,
		 const size_t kg_mc=0, const size_t kg_mb=0, const size_t kg_j=0 );

} // FFPACK minpoly
// #include "ffpack_minpoly.inl"

namespace FFPACK { /* Krylov Elim */

	template <class Field>
	size_t KrylovElim( const Field& F, const size_t M, const size_t N,
			   typename Field::Element_ptr A, const size_t lda, size_t*P,
			   size_t *Q, const size_t deg, size_t *iterates, size_t * inviterates, const size_t maxit,size_t virt);

	template <class Field>
	size_t  SpecRankProfile (const Field& F, const size_t M, const size_t N,
				 typename Field::Element_ptr A, const size_t lda, const size_t deg, size_t *rankProfile);

} // FFPACK KrylovElim
// #include "ffpack_krylovelim.inl"

namespace FFPACK { /* Solutions */
	/********/
	/* RANK */
	/********/



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
	      typename Field::Element_ptr A, const size_t lda) ;

	/********/
	/* DET  */
	/********/


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
		    typename Field::Element_ptr A, const size_t lda);

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
	     typename Field::Element_ptr A, const size_t lda);

	/*********/
	/* SOLVE */
	/*********/


	/// Solve linear system using LQUP factorization.
	template <class Field>
	typename Field::Element_ptr
	Solve( const Field& F, const size_t M,
	       typename Field::Element_ptr A, const size_t lda,
	       typename Field::Element_ptr x, const int incx,
	       typename Field::ConstElement_ptr b, const int incb );


	//! Solve L X = B or X L = B in place.
	//! L is M*M if Side == FFLAS::FflasLeft and N*N if Side == FFLAS::FflasRight, B is M*N.
	//! Only the R non trivial column of L are stored in the M*R matrix L
	//! Requirement :  so that L could  be expanded in-place
	template<class Field>
	void
	solveLB( const Field& F, const FFLAS::FFLAS_SIDE Side,
		 const size_t M, const size_t N, const size_t R,
		 typename Field::Element_ptr L, const size_t ldl,
		 const size_t * Q,
		 typename Field::Element_ptr B, const size_t ldb );

	//! Solve L X = B in place.
	//! L is M*M or N*N, B is M*N.
	//! Only the R non trivial column of L are stored in the M*R matrix L
	template<class Field>
	void
	solveLB2( const Field& F, const FFLAS::FFLAS_SIDE Side,
		  const size_t M, const size_t N, const size_t R,
		  typename Field::Element_ptr L, const size_t ldl,
		  const size_t * Q,
		  typename Field::Element_ptr B, const size_t ldb );


	/*************/
	/* NULLSPACE */
	/*************/

	/**  Computes a vector of the Left/Right nullspace of the matrix A.
	 *
	 * @param F The computation domain
	 * @param Side
	 * @param M
	 * @param N
	 * @param[in,out] A input matrix of dimension M x N, A is modified to its LU version
	 * @param lda
	 * @param[out] X output vector
	 * @param incX
	 *
	 */
	template <class Field>
	void RandomNullSpaceVector (const Field& F, const FFLAS::FFLAS_SIDE Side,
                                const size_t M, const size_t N,
                                typename Field::Element_ptr A, const size_t lda,
                                typename Field::Element_ptr X, const size_t incX);

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
			       typename Field::Element_ptr A, const size_t lda,
			       typename Field::Element_ptr& NS, size_t& ldn,
			       size_t& NSdim);

	/*****************/
	/* RANK PROFILES */
	/*****************/

	/** @brief Computes the row rank profile of A.
	 *
	 * @param F
	 * @param M
	 * @param N
	 * @param A input matrix of dimension M x N
	 * @param lda
	 * @param rkprofile return the rank profile as an array of row indexes, of dimension r=rank(A)
	 * @param LuTag: chooses the elimination algorithm. SlabRecursive for LUdivine, TileRecursive for PLUQ
	 *
	 * A is modified
	 * rkprofile is allocated during the computation.
	 * @returns R
	 */
	template <class Field>
	size_t RowRankProfile (const Field& F, const size_t M, const size_t N,
						   typename Field::Element_ptr A, const size_t lda,
						   size_t* &rkprofile,
						   const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);


	/**  @brief Computes the column rank profile of A.
	 *
	 * @param F
	 * @param M
	 * @param N
	 * @param A input matrix of dimension
	 * @param lda
	 * @param rkprofile return the rank profile as an array of row indexes, of dimension r=rank(A)
	 * @param LuTag: chooses the elimination algorithm. SlabRecursive for LUdivine, TileRecursive for PLUQ
	 *
	 * A is modified
	 * rkprofile is allocated during the computation.
	 * @returns R
	 */
	template <class Field>
	size_t ColumnRankProfile (const Field& F, const size_t M, const size_t N,
							  typename Field::Element_ptr A, const size_t lda,
							  size_t* &rkprofile,
							  const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);

	/**  @brief Recovers the column/row rank profile from the permutation of an LU decomposition.
	 *
	 * Works with both the CUP/PLE decompositions (obtained by LUdivine) or the PLUQ decomposition
	 * Assumes that the output vector containing the rank profile is already allocated.
	 * @param P the permutation carrying the rank profile information
	 * @param N the row/col dimension for a row/column rank profile
	 * @param R the rank of the matrix (
	 * @param rkprofile return the rank profile as an array of indices
	 * @param LuTag: chooses the elimination algorithm. SlabRecursive for LUdivine, TileRecursive for PLUQ
	 *
	 * A is modified
	 *
	 */
	void RankProfileFromLU (const size_t* P, const size_t N, const size_t R,
							size_t* rkprofile, const FFPACK_LU_TAG LuTag);

	/**  @brief Recovers the row and column rank profiles of any leading submatrix from the PLUQ decomposition.
	 *
	 * Only works with the PLUQ decomposition
	 * Assumes that the output vectors containing the rank profiles are already allocated.
	 *
	 * @param P the permutation carrying the rank profile information
	 * @param M the row dimension of the initial matrix
	 * @param N the column dimension of the initial matrix
	 * @param R the rank of the initial matrix
	 * @param LSm the row dimension of the leading submatrix considered
	 * @param LSn the column dimension of the leading submatrix considered
	 * @param P the row permutation of the PLUQ decomposition
	 * @param Q the column permutation of the PLUQ decomposition
	 * @param RRP return the row rank profile of the leading
	 * @param LuTag: chooses the elimination algorithm. SlabRecursive for LUdivine, TileRecursive for PLUQ
	 * @return the rank of the LSm x LSn leading submatrix
	 *
	 * A is modified
	 * @bib
	 * - Dumas J-G., Pernet C., and Sultan Z. <i>\c Simultaneous computation of the row and column rank profiles </i>, ISSAC'13.
	 */
	size_t LeadingSubmatrixRankProfiles (const size_t M, const size_t N, const size_t R,
										 const size_t LSm, const size_t LSn,
										 const size_t* P, const size_t* Q,
										 size_t* RRP, size_t* CRP);
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
										   typename Field::Element_ptr A,
										   const size_t lda,
										   size_t*& rowindices,
										   size_t*& colindices,
										   size_t& R);

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
					       typename Field::Element_ptr A,
					       const size_t lda,
					       size_t*& rowindices,
					       size_t*& colindices,
					       size_t& R);

	/** Computes the r*r submatrix X of A, by picking the row rank profile rows of A.
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
									typename Field::Element_ptr A,
									const size_t lda,
									typename Field::Element_ptr& X, size_t& R);

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
					typename Field::Element_ptr A, const size_t lda,
					typename Field::Element_ptr& X, size_t& R);

	/*********************************************/
	/* Accessors to Triangular and Echelon forms */
	/*********************************************/

	/** Extracts a triangular matrix from a compact storage A=L\U of rank R.
	 * if OnlyNonZeroVectors is false, then T and A have the same dimensions
	 * Otherwise, T is R x N if UpLo = FflasUpper, else T is  M x R
	 * @param F: base field
	 * @param UpLo: selects if the upper or lower triangular matrix is returned
	 * @param diag: selects if the triangular matrix unit-diagonal
	 * @param M: row dimension of T
	 * @param N: column dimension of T
	 * @param R: rank of the triangular matrix (how many rows/columns need to be copied)
	 * @param A: input matrix
	 * @param lda: leading dimension of A
	 * @param T: output matrix
	 * @param ldt: leading dimension of T
	 * @param OnlyNonZeroVectors: decides whether the last zero rows/columns should be ignored
	 */
	template <class Field>
	void
	getTriangular (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
				   const FFLAS::FFLAS_DIAG diag,
				   const size_t M, const size_t N, const size_t R,
				   typename Field::ConstElement_ptr A, const size_t lda,
				   typename Field::Element_ptr T, const size_t ldt,
				   const bool OnlyNonZeroVectors = false);

	/** Cleans up a compact storage A=L\U to reveal a triangular matrix of rank R.
	 * @param F: base field
	 * @param UpLo: selects if the upper or lower triangular matrix is revealed
	 * @param diag: selects if the triangular matrix unit-diagonal
	 * @param M: row dimension of A
	 * @param N: column dimension of A
	 * @param R: rank of the triangular matrix
	 * @param A: input/output matrix
	 * @param lda: leading dimension of A
	 */
	template <class Field>
	void
	getTriangular (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
				   const FFLAS::FFLAS_DIAG diag,
				   const size_t M, const size_t N, const size_t R,
				   typename Field::Element_ptr A, const size_t lda);

	/** Extracts a matrix in echelon form from a compact storage A=L\U of rank R obtained by
	 * RowEchelonForm or ColumnEchelonForm.
	 * Either L or U is in Echelon form (depending on Uplo)
	 * The echelon structure is defined by the first R values of the array P.
	 * row and column dimension of T are greater or equal to that of A
	 * @param F: base field
	 * @param UpLo: selects if the upper or lower triangular matrix is returned
	 * @param diag: selects if the echelon matrix has unit pivots
	 * @param M: row dimension of T
	 * @param N: column dimension of T
	 * @param R: rank of the triangular matrix (how many rows/columns need to be copied)
	 * @param P: positions of the R pivots
	 * @param A: input matrix
	 * @param lda: leading dimension of A
	 * @param T: output matrix
	 * @param ldt: leading dimension of T
	 * @param OnlyNonZeroVectors: decides whether the last zero rows/columns should be ignored
	 * @param LuTag: which factorized form (CUP/PLE if FfpackSlabRecursive, PLUQ if FfpackTileRecursive)
	 */
	template <class Field>
	void
	getEchelonForm (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
					const FFLAS::FFLAS_DIAG diag,
					const size_t M, const size_t N, const size_t R, const size_t* P,
					typename Field::ConstElement_ptr A, const size_t lda,
					typename Field::Element_ptr T, const size_t ldt,
					const bool OnlyNonZeroVectors = false,
					const FFPACK_LU_TAG LuTag = FfpackSlabRecursive);

    /** Cleans up a compact storage A=L\U obtained by RowEchelonForm or ColumnEchelonForm
	 * to reveal an echelon form of rank R.
	 * Either L or U is in Echelon form (depending on Uplo)
	 * The echelon structure is defined by the first R values of the array P.
	 * @param F: base field
	 * @param UpLo: selects if the upper or lower triangular matrix is returned
	 * @param diag: selects if the echelon matrix has unit pivots
	 * @param M: row dimension of A
	 * @param N: column dimension of A
	 * @param R: rank of the triangular matrix (how many rows/columns need to be copied)
	 * @param P: positions of the R pivots
	 * @param A: input/output matrix
	 * @param lda: leading dimension of A
	 * @param LuTag: which factorized form (CUP/PLE if FfpackSlabRecursive, PLUQ if FfpackTileRecursive)
	 */
	template <class Field>
	void
	getEchelonForm (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
					const FFLAS::FFLAS_DIAG diag,
					const size_t M, const size_t N, const size_t R, const size_t* P,
					typename Field::Element_ptr A, const size_t lda,
					const FFPACK_LU_TAG LuTag = FfpackSlabRecursive);

	/** Extracts a transformation matrix to echelon form from a compact storage A=L\U
	 * of rank R obtained by RowEchelonForm or ColumnEchelonForm.
	 * If Uplo == FflasLower:
	 *   T is N x N (already allocated) such that A T = C is a transformation of A in
	 *   Column echelon form
	 * Else
	 *   T is M x M (already allocated) such that T A = E is a transformation of A in
	 *   Row Echelon form
	 * @param F: base field
	 * @param UpLo: Lower means Transformation to Column Echelon Form, Upper, to Row Echelon Form
	 * @param diag: selects if the echelon matrix has unit pivots
	 * @param M: row dimension of A
	 * @param N: column dimension of A
	 * @param R: rank of the triangular matrix
	 * @param P: permutation matrix
	 * @param A: input matrix
	 * @param lda: leading dimension of A
	 * @param T: output matrix
	 * @param ldt: leading dimension of T
	 * @param LuTag: which factorized form (CUP/PLE if FfpackSlabRecursive, PLUQ if FfpackTileRecursive)
	 */
	template <class Field>
	void
	getEchelonTransform (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
						 const FFLAS::FFLAS_DIAG diag,
						 const size_t M, const size_t N, const size_t R, const size_t* P, const size_t* Q,
						 typename Field::ConstElement_ptr A, const size_t lda,
						 typename Field::Element_ptr T, const size_t ldt,
						 const FFPACK_LU_TAG LuTag = FfpackSlabRecursive);
	/** Extracts a matrix in echelon form from a compact storage A=L\U of rank R obtained by
	 * ReducedRowEchelonForm or ReducedColumnEchelonForm with transform = true.
	 * Either L or U is in Echelon form (depending on Uplo)
	 * The echelon structure is defined by the first R values of the array P.
	 * row and column dimension of T are greater or equal to that of A
	 * @param F: base field
	 * @param UpLo: selects if the upper or lower triangular matrix is returned
	 * @param diag: selects if the echelon matrix has unit pivots
	 * @param M: row dimension of T
	 * @param N: column dimension of T
	 * @param R: rank of the triangular matrix (how many rows/columns need to be copied)
	 * @param P: positions of the R pivots
	 * @param A: input matrix
	 * @param lda: leading dimension of A
	 * @param ldt: leading dimension of T
	 * @param LuTag: which factorized form (CUP/PLE if FfpackSlabRecursive, PLUQ if FfpackTileRecursive)
	 * @param OnlyNonZeroVectors: decides whether the last zero rows/columns should be ignored
	 */
	template <class Field>
	void
	getReducedEchelonForm (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
						   const size_t M, const size_t N, const size_t R, const size_t* P,
						   typename Field::ConstElement_ptr A, const size_t lda,
						   typename Field::Element_ptr T, const size_t ldt,
						   const bool OnlyNonZeroVectors = false,
						   const FFPACK_LU_TAG LuTag = FfpackSlabRecursive);

    /** Cleans up a compact storage A=L\U of rank R obtained by ReducedRowEchelonForm or
	 * ReducedColumnEchelonForm with transform = true.
	 * Either L or U is in Echelon form (depending on Uplo)
	 * The echelon structure is defined by the first R values of the array P.
	 * @param F: base field
	 * @param UpLo: selects if the upper or lower triangular matrix is returned
	 * @param diag: selects if the echelon matrix has unit pivots
	 * @param M: row dimension of A
	 * @param N: column dimension of A
	 * @param R: rank of the triangular matrix (how many rows/columns need to be copied)
	 * @param P: positions of the R pivots
	 * @param A: input/output matrix
	 * @param lda: leading dimension of A
	 * @param LuTag: which factorized form (CUP/PLE if FfpackSlabRecursive, PLUQ if FfpackTileRecursive)
	 */
	template <class Field>
	void
	getReducedEchelonForm (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
						   const size_t M, const size_t N, const size_t R, const size_t* P,
						   typename Field::Element_ptr A, const size_t lda,
						   const FFPACK_LU_TAG LuTag = FfpackSlabRecursive);

	/** Extracts a transformation matrix to echelon form from a compact storage A=L\U
	 * of rank R obtained by RowEchelonForm or ColumnEchelonForm.
	 * If Uplo == FflasLower:
	 *   T is N x N (already allocated) such that A T = C is a transformation of A in
	 *   Column echelon form
	 * Else
	 *   T is M x M (already allocated) such that T A = E is a transformation of A in
	 *   Row Echelon form
	 * @param F: base field
	 * @param UpLo: selects Col or Row Echelon Form
	 * @param diag: selects if the echelon matrix has unit pivots
	 * @param M: row dimension of A
	 * @param N: column dimension of A
	 * @param R: rank of the triangular matrix
	 * @param P: permutation matrix
	 * @param A: input matrix
	 * @param lda: leading dimension of A
	 * @param T: output matrix
	 * @param ldt: leading dimension of T
	 * @param LuTag: which factorized form (CUP/PLE if FfpackSlabRecursive, PLUQ if FfpackTileRecursive)
	 */
	template <class Field>
	void
	getReducedEchelonTransform (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
								const size_t M, const size_t N, const size_t R, const size_t* P, const size_t* Q,
								typename Field::ConstElement_ptr A, const size_t lda,
								typename Field::Element_ptr T, const size_t ldt,
								const FFPACK_LU_TAG LuTag = FfpackSlabRecursive);
		/** Auxiliary routine: determines the permutation that changes a PLUQ decomposition
		 * into a echelon form revealing PLUQ decomposition
		 */
	void
	PLUQtoEchelonPermutation (const size_t N, const size_t R, const size_t * P, size_t * outPerm);

} // FFPACK
// #include "ffpack.inl"

namespace FFPACK { /* not used */

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
	typename Field::Element_ptr
	LQUPtoInverseOfFullRankMinor( const Field& F, const size_t rank,
				      typename Field::Element_ptr A_factors, const size_t lda,
				      const size_t* QtPointer,
				      typename Field::Element_ptr X, const size_t ldx);

} // FFPACK
// include precompiled instantiation headers (avoiding to recompile them)
#ifdef FFPACK_COMPILED 
#include "fflas-ffpack/interfaces/libs/ffpack_inst.h"
#endif

#include "ffpack_fgesv.inl"
#include "ffpack_fgetrs.inl"
#include "ffpack_ftrtr.inl"
#include "ffpack_pluq.inl"
#include "ffpack_pluq_mp.inl"
#include "ffpack_ppluq.inl"
#include "ffpack_ludivine.inl"
#include "ffpack_ludivine_mp.inl"
#include "ffpack_echelonforms.inl"
#include "ffpack_invert.inl"
#include "ffpack_charpoly_kglu.inl"
#include "ffpack_charpoly_kgfast.inl"
#include "ffpack_charpoly_kgfastgeneralized.inl"
#include "ffpack_charpoly_danilevski.inl"
#include "ffpack_charpoly.inl"
#include "ffpack_frobenius.inl"
#include "ffpack_minpoly.inl"
#include "ffpack_krylovelim.inl"
#include "ffpack_permutation.inl"
#include "ffpack_rankprofiles.inl"
#include "ffpack.inl"
#endif // __FFLASFFPACK_ffpack_H

