/*
 * Copyright (C) 2014,2015 the FFLAS-FFPACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *            Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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


namespace FFLAS {
    //---------------------------------------------------------------------
    // Level 2 routines
    //---------------------------------------------------------------------

    /** \brief fassign : \f$A \gets B \f$.
     * @param F field
     * @param m number of rows to copy
     * @param n number of cols to copy
     * \param A matrix in \p F
     * \param lda stride of \p A
     * \param B vector in \p F
     * \param ldb stride of \p B
     */
    template INST_OR_DECL
    void
    fassign (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m, const size_t n,
             const FFLAS_ELT* B, const size_t ldb ,
             FFLAS_ELT* A, const size_t lda );

    /** \brief fzero : \f$A \gets 0 \f$.
     * @param F field
     * @param m number of rows to zero
     * @param n number of cols to zero
     * \param A matrix in \p F
     * \param lda stride of \p A
     * @warning may be buggy if Element is larger than int
     */

    template INST_OR_DECL
    void
    fzero (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m, const size_t n,
           FFLAS_ELT* A, const size_t lda);
    // {
    // 	/*  use memset only with Elements that are ok */
    // 	if (n == lda) { // contigous data
    // 		// memset(A,(int) F.zero,m*n); // might be bogus ?
    // 		fzero(F,m*n,A,1);
    // 	}
    // 	else { // not contiguous (strided)
    // 		for (size_t i = 0 ; i < m ; ++i)
    // 			// memset(A+i*lda,(int) F.zero,n) ; // might be bogus ?
    // 			fzero(F,n,A+i*lda,1);
    // 	}
    // }
    /** \brief fequal : test \f$A = B \f$.
     * @param F field
     * @param m row dimension
     * @param n column dimension
     * \param A m x n matrix in \p F
     * \param lda leading dimension of A
     * \param B m x n matrix in \p F
     * \param ldb leading dimension of B
     */
    template INST_OR_DECL
    bool
    fequal (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m, const size_t n,
            const FFLAS_ELT* A, const size_t lda,
            const FFLAS_ELT* B, const size_t ldb);
    // {
    // 	bool res=true;
    // 	for (size_t i = 0 ; i < m ; ++i)
    // 		res &= fequal (F, n, A + i*lda, 1, B + i*ldb, 1);
    // 	return res;
    // }
    /** \brief fiszero : test \f$A = 0 \f$.
     * @param F field
     * @param m row dimension
     * @param n column dimension
     * \param A m x n matrix in \p F
     * \param lda leading dimension of A
     */
    template INST_OR_DECL
    bool
    fiszero (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m, const size_t n,
             const FFLAS_ELT* A, const size_t lda);
    // {
    // 	bool res=true;
    // 	for (size_t i = 0 ; i < m ; ++i)
    // 		res &= fiszero (F, n, A + i*lda, 1);
    // 	return res;
    // }

    //! creates a diagonal matrix
    template INST_OR_DECL
    void
    fidentity (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m, const size_t n,
               FFLAS_ELT* A, const size_t lda, const FFLAS_ELT & d);
    // {
    // 	fzero(F,m,n,A,lda);
    // 	for (size_t i = 0 ; i < std::min(m,n) ; ++i)
    // 		F.assign(A[i*lda+i],d);
    // }

    //! creates a diagonal matrix
    template INST_OR_DECL
    void
    fidentity (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m, const size_t n,
               FFLAS_ELT* A, const size_t lda);
    // {
    // 	fzero(F,m,n,A,lda);
    // 	for (size_t i = 0 ; i < std::min(m,n) ; ++i)
    // 		F.assign(A[i*lda+i],F.one);
    // }

    /** freduce
     * \f$A \gets  A mod F\f$.
     * @param F field
     * @param m number of rows
     * @param n number of cols
     * \param A matrix in \p F
     * \param lda stride of \p A
     * @internal
     */
    template INST_OR_DECL
    void
    freduce (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m , const size_t n,
             FFLAS_ELT* A, const size_t lda);

    /** freduce
     * \f$A \gets  B mod F\f$.
     * @param F field
     * @param m number of rows
     * @param n number of cols
     * \param A matrix in \p F
     * \param lda stride of \p A
     * \param B matrix in \p Element
     * \param ldb stride of \p B
     * @internal
     */
    template INST_OR_DECL
    void
    freduce (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m , const size_t n,
             const FFLAS_ELT* B, const size_t ldb,
             FFLAS_ELT* A, const size_t lda);

    /** finit
     * \f$A \gets  B mod F\f$.
     * @param F field
     * @param m number of rows
     * @param n number of cols
     * \param A matrix in \p F
     * \param lda stride of \p A
     * \param B matrix in \p F
     * \param ldb stride of \p B
     * @internal
     */
    template INST_OR_DECL
    void
    finit (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m , const size_t n,
           const FFLAS_ELT* B, const size_t ldb,
           FFLAS_ELT* A, const size_t lda);


    /** fnegin
     * \f$A \gets - A\f$.
     * @param F field
     * @param m number of rows
     * @param n number of cols
     * \param A matrix in \p F
     * \param lda stride of \p A
     * @internal
     */
    template INST_OR_DECL
    void
    fnegin (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m , const size_t n,
            FFLAS_ELT* A, const size_t lda);
    // {
    // 	//!@todo check if n == lda
    // 	for (size_t i = 0 ; i < m ; ++i)
    // 		fnegin(F,n,A+i*lda,1);
    // 	return;
    // }

    /** fneg
     * \f$A \gets  - B\f$.
     * @param F field
     * @param m number of rows
     * @param n number of cols
     * \param A matrix in \p F
     * \param lda stride of \p A
     * @internal
     */
    template INST_OR_DECL
    void
    fneg (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m , const size_t n,
          const FFLAS_ELT* B, const size_t ldb,
          FFLAS_ELT* A, const size_t lda);
    // {
    // 	//!@todo check if n == lda
    // 	for (size_t i = 0 ; i < m ; ++i)
    // 		fneg(F,n,B+i*ldb,1,A+i*lda,1);
    // 	return;
    // }

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
    template INST_OR_DECL
    void
    fscalin (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m , const size_t n,
             const FFLAS_ELT alpha,
             FFLAS_ELT* A, const size_t lda);

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
    template INST_OR_DECL
    void
    fscal (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m , const size_t n,
           const FFLAS_ELT alpha,
           const FFLAS_ELT* A, const size_t lda,
           FFLAS_ELT* B, const size_t ldb);

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
    template INST_OR_DECL
    void
    faxpy (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m, const size_t n
           , const FFLAS_ELT alpha,
           const FFLAS_ELT* X, const size_t ldx,
           FFLAS_ELT* Y, const size_t ldy );

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
    // template INST_OR_DECL
    // void
    // faxpby (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m, const size_t n,
    // 	const FFLAS_ELT alpha,
    // 	const FFLAS_ELT* X, const size_t ldx,
    // 	const FFLAS_ELT beta,
    // 	FFLAS_ELT* Y, const size_t ldy );

    /** \brief fmove : \f$A \gets B \f$ and \f$ B \gets 0\f$.
     * @param F field
     * @param m number of rows to copy
     * @param n number of cols to copy
     * \param A matrix in \p F
     * \param lda stride of \p A
     * \param B vector in \p F
     * \param ldb stride of \p B
     */
    template INST_OR_DECL
    void
    fmove (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t m, const size_t n,
           FFLAS_ELT* A, const size_t lda,
           FFLAS_ELT* B, const size_t ldb );
    // {
    // 	fassign(F,m,n,A,lda,B,ldb);
    // 	fzero(F,m,n,B,ldb);
    // }

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
    template  INST_OR_DECL
    void
    fadd (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
          const FFLAS_ELT* A, const size_t lda,
          const FFLAS_ELT* B, const size_t ldb,
          FFLAS_ELT* C, const size_t ldc);



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
    template  INST_OR_DECL
    void
    fsub (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
          const FFLAS_ELT* A, const size_t lda,
          const FFLAS_ELT* B, const size_t ldb,
          FFLAS_ELT* C, const size_t ldc);

    //! fsubin
    //! C = C - B
    template  INST_OR_DECL
    void
    fsubin (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
            const FFLAS_ELT* B, const size_t ldb,
            FFLAS_ELT* C, const size_t ldc);

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
    template  INST_OR_DECL
    void
    fadd (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
          const FFLAS_ELT* A, const size_t lda,
          const FFLAS_ELT alpha,
          const FFLAS_ELT* B, const size_t ldb,
          FFLAS_ELT* C, const size_t ldc);

    //! faddin
    template  INST_OR_DECL
    void
    faddin (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
            const FFLAS_ELT* B, const size_t ldb,
            FFLAS_ELT* C, const size_t ldc);


    /**  @brief finite prime FFLAS_FIELD<FFLAS_ELT> GEneral Matrix Vector multiplication.
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
    template INST_OR_DECL
    FFLAS_ELT*
    fgemv (const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS_TRANSPOSE TransA,
           const size_t M, const size_t N,
           const FFLAS_ELT alpha,
           const FFLAS_ELT* A, const size_t lda,
           const FFLAS_ELT* X, const size_t incX,
           const  FFLAS_ELT beta,
           FFLAS_ELT* Y, const size_t incY);

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
    template INST_OR_DECL
    void
    fger (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
          const FFLAS_ELT alpha,
          const FFLAS_ELT* x, const size_t incx,
          const FFLAS_ELT* y, const size_t incy,
          FFLAS_ELT* A, const size_t lda);

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
    template INST_OR_DECL
    void
    ftrsv (const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS_UPLO Uplo,
           const FFLAS_TRANSPOSE TransA, const FFLAS_DIAG Diag,
           const size_t N,const FFLAS_ELT* A, const size_t lda,
           FFLAS_ELT* X, int incX);



} // FFLAS

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
