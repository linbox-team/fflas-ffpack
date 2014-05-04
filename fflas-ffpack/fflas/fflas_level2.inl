/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by BB <bbboyer@ncsu.edu>
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

/** @file fflas/fflas_level2.h
 * @brief  Matrix-Vector operations
 * or anything of \f$n^2\f$ complexity
 */

#ifndef __FFLASFFPACK_fflas_fflas_level2_INL
#define __FFLASFFPACK_fflas_fflas_level2_INL

namespace FFLAS {

	//---------------------------------------------------------------------
	// Level 2 routines
	//---------------------------------------------------------------------

	/** \brief fcopy : \f$A \gets B \f$.
	 * @param F field
	 * @param m number of rows to copy
	 * @param n number of cols to copy
	 * \param A matrix in \p F
	 * \param lda stride of \p A
	 * \param B vector in \p F
	 * \param ldb stride of \p B
	 */
	template<class Field>
	void
	fcopy (const Field& F, const size_t m, const size_t n,
	       const typename Field::Element * B, const size_t ldb ,
	       typename Field::Element * A, const size_t lda );

	/** \brief fzero : \f$A \gets 0 \f$.
	 * @param F field
	 * @param m number of rows to zero
	 * @param n number of cols to zero
	 * \param A matrix in \p F
	 * \param lda stride of \p A
	 * @warning may be buggy if Element is larger than int
	 */

	template<class Field>
	void
	fzero (const Field& F, const size_t m, const size_t n,
	       typename Field::Element * A, const size_t lda)
	{
		/*  use memset only with Elements that are ok */
		if (n == lda) { // contigous data
			// memset(A,(int) F.zero,m*n); // might be bogus ?
			fzero(F,m*n,A,1);
		}
		else { // not contiguous (strided)
			for (size_t i = 0 ; i < m ; ++i)
				// memset(A+i*lda,(int) F.zero,n) ; // might be bogus ?
				fzero(F,n,A+i*lda,1);
		}
	}

	//! creates a diagonal matrix
	template<class Field>
	void
	fidentity (const Field& F, const size_t m, const size_t n,
		   typename Field::Element * A, const size_t lda, const typename Field::Element & d) // =F.one...
	{
		fzero(F,m,n,A,lda);
		for (size_t i = 0 ; i < std::min(m,n) ; ++i)
			F.assign(A[i*lda+i],d);
	}

	//! creates a diagonal matrix
	template<class Field>
	void
	fidentity (const Field& F, const size_t m, const size_t n,
		   typename Field::Element * A, const size_t lda)
	{
		fzero(F,m,n,A,lda);
		for (size_t i = 0 ; i < std::min(m,n) ; ++i)
			F.assign(A[i*lda+i],F.one);
	}

	/** finit
	 * \f$A \gets  A mod F\f$.
	 * @param F field
	 * @param m number of rows
	 * @param n number of cols
	 * \param A matrix in \p F
	 * \param lda stride of \p A
	 * @internal
	 */
	template<class Field>
	void
	finit (const Field& F, const size_t m , const size_t n,
	       typename Field::Element * A, const size_t lda);

	/** finit
	 * \f$A \gets  B mod F\f$.
	 * @param F field
	 * @param m number of rows
	 * @param n number of cols
	 * \param A matrix in \p F
	 * \param lda stride of \p A
	 * \param B matrix in \p OtherElement
	 * \param ldb stride of \p B
	 * @internal
	 */
	template<class Field, class OtherElement>
	void
	finit (const Field& F, const size_t m , const size_t n,
	       const OtherElement * B, const size_t ldb,
	       typename Field::Element * A, const size_t lda);

	/** fconvert
	 * \f$A \gets  B mod F\f$.
	 * @param F field
	 * @param m number of rows
	 * @param n number of cols
	 * \param A matrix in \p OtherElement
	 * \param lda stride of \p A
	 * \param B matrix in \p F
	 * \param ldb stride of \p B
	 * @internal
	 */
	template<class Field, class OtherElement>
	void
	fconvert (const Field& F, const size_t m , const size_t n,
	        OtherElement * A, const size_t lda,
	       const typename Field::Element* B, const size_t ldb)
	{
		//!@todo check if n == lda
		for (size_t i = 0 ; i < m ; ++i)
			fconvert(F,n,A+i*lda,1,B+i*ldb,1);
		return;
	}

	/** fnegin
	 * \f$A \gets - A\f$.
	 * @param F field
	 * @param m number of rows
	 * @param n number of cols
	 * \param A matrix in \p F
	 * \param lda stride of \p A
	 * @internal
	 */
	template<class Field>
	void
	fnegin (const Field& F, const size_t m , const size_t n,
	       typename Field::Element * A, const size_t lda)
	{
		//!@todo check if n == lda
		for (size_t i = 0 ; i < m ; ++i)
			fnegin(F,n,A+i*lda,1);
		return;
	}

	/** fneg
	 * \f$A \gets  - B\f$.
	 * @param F field
	 * @param m number of rows
	 * @param n number of cols
	 * \param A matrix in \p F
	 * \param lda stride of \p A
	 * @internal
	 */
	template<class Field>
	void
	fneg (const Field& F, const size_t m , const size_t n,
	       const typename Field::Element * B, const size_t ldb,
	       typename Field::Element * A, const size_t lda)
	{
		//!@todo check if n == lda
		for (size_t i = 0 ; i < m ; ++i)
			fneg(F,n,B+i*ldb,1,A+i*lda,1);
		return;
	}

	/** fscalin
	 * \f$A \gets a \cdot A\f$.
	 * @param F field
	 * @param m number of rows
	 * @param n number of cols
	 * @param alpha homotecie scalar
	 * \param A matrix in \p F
	 * \param lda stride of \p A
	 * @internal
	 */
	template<class Field>
	void
	fscalin (const Field& F, const size_t m , const size_t n,
	       const typename Field::Element alpha,
	       typename Field::Element * A, const size_t lda);

	/** fscal
	 * \f$B \gets a \cdot A\f$.
	 * @param F field
	 * @param m number of rows
	 * @param n number of cols
	 * @param alpha homotecie scalar
	 * \param[in] A matrix in \p F
	 * \param lda stride of \p A
	 * \param[out] B matrix in \p F
	 * \param ldb stride of \p B
	 * @internal
	 */
	template<class Field>
	void
	fscal (const Field& F, const size_t m , const size_t n,
	       const typename Field::Element alpha,
	       const typename Field::Element * A, const size_t lda,
	       typename Field::Element * B, const size_t ldb);

	/** \brief faxpy : \f$y \gets \alpha \cdot x + y\f$.
	 * @param F field
	 * @param m row dimension
	 * @param n column dimension
	 * @param alpha scalar
	 * \param[in] X vector in \p F
	 * \param ldx leading dimension of \p X
	 * \param[in,out] Y vector in \p F
	 * \param ldy leading dimension of \p Y
	 */
	template<class Field>
	void
	faxpy (const Field& F, const size_t m, const size_t n
	       , const typename Field::Element alpha,
	       const typename Field::Element * X, const size_t ldx,
	       typename Field::Element * Y, const size_t ldy );

	/** \brief faxpby : \f$y \gets \alpha \cdot x + \beta \cdot y\f$.
	 * @param F field
	 * @param m row dimension
	 * @param n column dimension
	 * @param alpha scalar
	 * \param[in] X vector in \p F
	 * \param ldx leading dimension of \p X
	 * \param beta scalar
	 * \param[in,out] Y vector in \p F
	 * \param ldy leading dimension of \p Y
	 * \note this is a catlas function
	 */
	template<class Field>
	void
	faxpby (const Field& F, const size_t m, const size_t n,
	       const typename Field::Element alpha,
	       const typename Field::Element * X, const size_t ldx,
	       const typename Field::Element beta,
	       typename Field::Element * Y, const size_t ldy );

	/** \brief fmove : \f$A \gets B \f$ and \f$ B \gets 0\f$.
	 * @param F field
	 * @param m number of rows to copy
	 * @param n number of cols to copy
	 * \param A matrix in \p F
	 * \param lda stride of \p A
	 * \param B vector in \p F
	 * \param ldb stride of \p B
	 */
	template<class Field>
	void
	fmove (const Field& F, const size_t m, const size_t n,
	       typename Field::Element * A, const size_t lda,
	       typename Field::Element * B, const size_t ldb )
	{
		fcopy(F,m,n,A,lda,B,ldb);
		fzero(F,m,n,B,ldb);
	}

	/** fadd : matrix addition.
	 * Computes \p C = \p A + \p B.
	 * @param F field
	 * @param M rows
	 * @param N cols
	 * @param A dense matrix of size \c MxN
	 * @param lda leading dimension of \p A
	 * @param B dense matrix of size \c MxN
	 * @param ldb leading dimension of \p B
	 * @param C dense matrix of size \c MxN
	 * @param ldc leading dimension of \p C
	 */
	template <class Field>
	void
	fadd (const Field& F, const size_t M, const size_t N,
	      const typename Field::Element* A, const size_t lda,
	      const typename Field::Element* B, const size_t ldb,
	      typename Field::Element* C, const size_t ldc);



	/** fsub : matrix subtraction.
	 * Computes \p C = \p A - \p B.
	 * @param F field
	 * @param M rows
	 * @param N cols
	 * @param A dense matrix of size \c MxN
	 * @param lda leading dimension of \p A
	 * @param B dense matrix of size \c MxN
	 * @param ldb leading dimension of \p B
	 * @param C dense matrix of size \c MxN
	 * @param ldc leading dimension of \p C
	 */
	template <class Field>
	void
	fsub (const Field& F, const size_t M, const size_t N,
	      const typename Field::Element* A, const size_t lda,
	      const typename Field::Element* B, const size_t ldb,
	      typename Field::Element* C, const size_t ldc);

	//! fsubin
	//! C = C - B
	template <class Field>
	void
	fsubin (const Field& F, const size_t M, const size_t N,
		const typename Field::Element* B, const size_t ldb,
		typename Field::Element* C, const size_t ldc);

	/** fadd : matrix addition with scaling.
	 * Computes \p C = \p A + alpha \p B.
	 * @param F field
	 * @param M rows
	 * @param N cols
	 * @param A dense matrix of size \c MxN
	 * @param lda leading dimension of \p A
	 * @param alpha some scalar
	 * @param B dense matrix of size \c MxN
	 * @param ldb leading dimension of \p B
	 * @param C dense matrix of size \c MxN
	 * @param ldc leading dimension of \p C
	 */
	template <class Field>
	void
	fadd (const Field& F, const size_t M, const size_t N,
	      const typename Field::Element* A, const size_t lda,
	      const typename Field::Element alpha,
	      const typename Field::Element* B, const size_t ldb,
	      typename Field::Element* C, const size_t ldc);

	//! faddin
	template <class Field>
	void
	faddin (const Field& F, const size_t M, const size_t N,
		const typename Field::Element* B, const size_t ldb,
		typename Field::Element* C, const size_t ldc);


	/**  @brief finite prime Field GEneral Matrix Vector multiplication.
	 *
	 *  Computes  \f$Y \gets \alpha \mathrm{op}(A) X + \beta Y \f$.
	 * @param F field
	 * \param TransA if \c TransA==FflasTrans then \f$\mathrm{op}(A)=A^t\f$.
	 * @param M rows
	 * @param N cols
	 * @param alpha scalar
	 * @param A dense matrix of size \c MxN
	 * @param lda leading dimension of \p A
	 * @param X dense vector of size \c N
	 * @param incX stride of \p X
	 * @param beta scalar
	 * @param[out] Y dense vector of size \c M
	 * @param incY stride of \p Y
	 */
	template<class Field>
	void
	fgemv (const Field& F, const FFLAS_TRANSPOSE TransA,
	       const size_t M, const size_t N,
	       const typename Field::Element alpha,
	       const typename Field::Element * A, const size_t lda,
	       const typename Field::Element * X, const size_t incX,
	       const  typename Field::Element beta,
	       typename Field::Element * Y, const size_t incY);

	/**  @brief fger: rank one update of a general matrix
	 *
	 *  Computes  \f$A \gets \alpha x . y^T + A\f$
	 * @param F field
	 * @param M rows
	 * @param N cols
	 * @param alpha scalar
	 * @param[in,out] A dense matrix of size \c MxN and leading dimension \p lda
	 * @param lda leading dimension of \p A
	 * @param x dense vector of size \c M
	 * @param incx stride of \p X
	 * @param y dense vector of size \c N
	 * @param incy stride of \p Y
	 */
	template<class Field>
	void
	fger (const Field& F, const size_t M, const size_t N,
	      const typename Field::Element alpha,
	      const typename Field::Element * x, const size_t incx,
	      const typename Field::Element * y, const size_t incy,
	      typename Field::Element * A, const size_t lda);

	/** @brief ftrsv: TRiangular System solve with Vector
	 *  Computes  \f$ X \gets \mathrm{op}(A^{-1}) X\f$
	 *  @param F field
	 * @param X vector of size \p N on a field \p F
	 * @param incX stride of \p  X
	 * @param A a matrix of leading dimension \p lda and size \p N
	 * @param lda leading dimension of \p A
	 * @param N number of rows or columns of \p A according to \p TransA
	 * \param TransA if \c TransA==FflasTrans then \f$\mathrm{op}(A)=A^t\f$.
	 * \param Diag if \c Diag==FflasUnit then \p A is unit.
	 * \param Uplo if \c Uplo==FflasUpper then \p A is upper triangular
	 */
	template<class Field>
	void
	ftrsv (const Field& F, const FFLAS_UPLO Uplo,
	       const FFLAS_TRANSPOSE TransA, const FFLAS_DIAG Diag,
	       const size_t N,const typename Field::Element * A, const size_t lda,
	       typename Field::Element * X, int incX);

} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_level2_INL

