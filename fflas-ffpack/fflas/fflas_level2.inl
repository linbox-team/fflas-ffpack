/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

#include "givaro/zring.h"
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
    template<class Field>
    void
    fassign (const Field& F, const size_t m, const size_t n,
             typename Field::ConstElement_ptr B, const size_t ldb ,
             typename Field::Element_ptr A, const size_t lda );

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
           typename Field::Element_ptr A, const size_t lda)
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
    /** \brief fzero : \f$A \gets 0 \f$ for a triangular matrix.
     * @param F field
     * @param shape shape of the triangular matrix
     * @param m number of rows to zero
     * @param n number of cols to zero
     * \param A matrix in \p F
     * \param lda stride of \p A
     * @warning may be buggy if Element is larger than int
     */

    template<class Field>
    void
    fzero (const Field& F, const FFLAS_UPLO shape, const FFLAS_DIAG diag,
           const size_t n, typename Field::Element_ptr A, const size_t lda)
    {
        ptrdiff_t inc_row;
        typename Field::Element_ptr Ai;
        size_t size = n;
        switch (shape){
            case FflasUpper:{inc_row = lda+1; Ai = A; break;}
            case FflasLower:{inc_row = -lda; Ai = A+(n-1)*lda; break;}
            case FflasLeftTri:{inc_row = lda; Ai = A; break;}
            case FflasRightTri:{inc_row = -lda+1; Ai = A+(n-1)*lda; break;}
        }
        if (diag == FflasUnit){
            size--;
            if (shape == FflasUpper || shape == FflasRightTri) Ai++;
        }
        for (; size > 0; Ai += inc_row, size--)
            fzero (F, size, Ai, 1);
    }

    /** \brief frand : \f$A \gets random \f$.
     * @param F field
     * @param G randomiterator
     * @param m number of rows to randomize
     * @param n number of cols to randomize
     * \param A matrix in \p F
     * \param lda stride of \p A
     */
    template<class Field, class RandIter>
    void
    frand (const Field& F, RandIter& G, const size_t m, const size_t n,
           typename Field::Element_ptr A, const size_t lda)
    {
        /*  use memset only with Elements that are ok */
        if (n == lda) { // contigous data
            // memset(A,(int) F.zero,m*n); // might be bogus ?
            frand(F,G,m*n,A,1);
        }
        else { // not contiguous (strided)
            for (size_t i = 0 ; i < m ; ++i)
                // memset(A+i*lda,(int) F.zero,n) ; // might be bogus ?
                frand(F,G,n,A+i*lda,1);
        }
    }
    /** \brief fequal : test \f$A = B \f$.
     * @param F field
     * @param m row dimension
     * @param n column dimension
     * \param A m x n matrix in \p F
     * \param lda leading dimension of A
     * \param B m x n matrix in \p F
     * \param ldb leading dimension of B
     */
    template<class Field>
    bool
    fequal (const Field& F, const size_t m, const size_t n,
            typename Field::ConstElement_ptr A, const size_t lda,
            typename Field::ConstElement_ptr B, const size_t ldb)
    {
        bool res=true;
        for (size_t i = 0 ; i < m ; ++i)
            res &= fequal (F, n, A + i*lda, 1, B + i*ldb, 1);
        return res;
    }
    /** \brief fiszero : test \f$A = 0 \f$.
     * @param F field
     * @param m row dimension
     * @param n column dimension
     * \param A m x n matrix in \p F
     * \param lda leading dimension of A
     */
    template<class Field>
    bool
    fiszero (const Field& F, const size_t m, const size_t n,
             typename Field::ConstElement_ptr A, const size_t lda)
    {
        bool res=true;
        for (size_t i = 0 ; i < m ; ++i)
            res &= fiszero (F, n, A + i*lda, 1);
        return res;
    }

    //! creates a diagonal matrix
    template<class Field>
    void
    fidentity (const Field& F, const size_t m, const size_t n,
               typename Field::Element_ptr A, const size_t lda, const typename Field::Element & d) // =F.one...
    {
        fzero(F,m,n,A,lda);
        for (size_t i = 0 ; i < std::min(m,n) ; ++i)
            F.assign(A[i*lda+i],d);
    }

    //! creates a diagonal matrix
    template<class Field>
    void
    fidentity (const Field& F, const size_t m, const size_t n,
               typename Field::Element_ptr A, const size_t lda)
    {
        fzero(F,m,n,A,lda);
        for (size_t i = 0 ; i < std::min(m,n) ; ++i)
            F.assign(A[i*lda+i],F.one);
    }

    /** freduce
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
    freduce (const Field& F, const size_t m , const size_t n,
             typename Field::Element_ptr A, const size_t lda);

        //! freduce for square symmetric matrices
    template<class Field>
    void
    freduce (const Field& F,  const FFLAS_UPLO uplo, const size_t N,
             typename Field::Element_ptr A, const size_t lda);

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
    template<class Field>
    void
    freduce (const Field& F, const size_t m , const size_t n,
             typename Field::ConstElement_ptr B, const size_t ldb,
             typename Field::Element_ptr A, const size_t lda);

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
    template<class Field, class OtherElement_ptr>
    void
    finit (const Field& F, const size_t m , const size_t n,
           const OtherElement_ptr B, const size_t ldb,
           typename Field::Element_ptr A, const size_t lda);

    /** finit
     * Initializes \p A in \p F$.
     * @param F field
     * @param m number of rows
     * @param n number of cols
     * \param A matrix in \p F
     * \param lda stride of \p A
     * @internal
     */
    template<class Field, class OtherElement_ptr>
    void
    finit (const Field& F, const size_t m , const size_t n,
           typename Field::Element_ptr A, const size_t lda);

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
    template<class Field, class OtherElement_ptr>
    void
    fconvert (const Field& F, const size_t m , const size_t n,
              OtherElement_ptr A, const size_t lda,
              typename Field::ConstElement_ptr B, const size_t ldb)
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
            typename Field::Element_ptr A, const size_t lda)
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
          typename Field::ConstElement_ptr B, const size_t ldb,
          typename Field::Element_ptr A, const size_t lda)
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
             typename Field::Element_ptr A, const size_t lda);

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
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::Element_ptr B, const size_t ldb);

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
           typename Field::ConstElement_ptr X, const size_t ldx,
           typename Field::Element_ptr Y, const size_t ldy );

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
            typename Field::ConstElement_ptr X, const size_t ldx,
            const typename Field::Element beta,
            typename Field::Element_ptr Y, const size_t ldy );

    /** \brief fmove : \f$A \gets B \f$ and \f$ B \gets 0\f$.
     * @param F field
     * @param m number of rows to copy
     * @param n number of cols to copy
     * \param A matrix in \p F
     * \param lda stride of \p A
     * \param B matrix in \p F
     * \param ldb stride of \p B
     */
    template<class Field>
    void
    fmove (const Field& F, const size_t m, const size_t n,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::Element_ptr B, const size_t ldb )
    {
        fassign(F,m,n,A,lda,B,ldb);
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
          typename Field::ConstElement_ptr A, const size_t lda,
          typename Field::ConstElement_ptr B, const size_t ldb,
          typename Field::Element_ptr C, const size_t ldc);



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
          typename Field::ConstElement_ptr A, const size_t lda,
          typename Field::ConstElement_ptr B, const size_t ldb,
          typename Field::Element_ptr C, const size_t ldc);

    //! fsubin
    //! C = C - B
    template <class Field>
    void
    fsubin (const Field& F, const size_t M, const size_t N,
            typename Field::ConstElement_ptr B, const size_t ldb,
            typename Field::Element_ptr C, const size_t ldc);

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
          typename Field::ConstElement_ptr A, const size_t lda,
          const typename Field::Element alpha,
          typename Field::ConstElement_ptr B, const size_t ldb,
          typename Field::Element_ptr C, const size_t ldc);

    //! faddin
    template <class Field>
    void
    faddin (const Field& F, const size_t M, const size_t N,
            typename Field::ConstElement_ptr B, const size_t ldb,
            typename Field::Element_ptr C, const size_t ldc);

    //! fadding for symmetric matrices
    template <class Field>
    void
    faddin (const Field& F,
            const FFLAS_UPLO uplo,
            const size_t N,
            typename Field::ConstElement_ptr B, const size_t ldb,
            typename Field::Element_ptr C, const size_t ldc);


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
    typename Field::Element_ptr
    fgemv (const Field& F, const FFLAS_TRANSPOSE TransA,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr X, const size_t incX,
           const  typename Field::Element beta,
           typename Field::Element_ptr Y, const size_t incY);

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
          typename Field::ConstElement_ptr x, const size_t incx,
          typename Field::ConstElement_ptr y, const size_t incy,
          typename Field::Element_ptr A, const size_t lda);

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
           const size_t N,typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::Element_ptr X, int incX);

    /** @brief bitsize:
     *  Computes  the largest bitsize of the matrix' coefficients.
     *  If the matrix is over a modular prime field, it returns the bitsize of the largest
     *  element (in a bsolute value)
     * @param F field
     * @param M rows
     * @param N cols
     * @param incX stride of \p  X
     * @param A a matrix of leading dimension \p lda and size \p MxN
     * @param lda leading dimension of \p A
     */
    template<class Field>
    inline size_t bitsize(const Field& F, size_t M, size_t N, const typename Field::ConstElement_ptr A, size_t lda){
        Givaro::Integer min = F.minElement() ,max = F.maxElement();
        return std::max(max,-min).bitsize();
    }

    template<>
    inline size_t bitsize<Givaro::ZRing<Givaro::Integer> >(const Givaro::ZRing<Givaro::Integer>& F, size_t M, size_t N, const Givaro::Integer* A, size_t lda){
        size_t bs = 1;
        for (size_t i=0; i<M; ++i)
            for (size_t j=0; j<N; ++j)
                bs = std::max(bs, A[i*lda+j].bitsize());
        return bs;
    }
    /** @brief ftrsm: TRiangular Matrix Vector prodcut
     * Computes  \f$ X \gets \mathrm{op}(A) X\f$
     * @param F field
     * @param X vector of size \p N on a field \p F
     * @param incX stride of \p  X
     * @param A a matrix of leading dimension \p lda and size \p N
     * @param lda leading dimension of \p A
     * @param N number of rows and columns of \p A
     * \param TransA if \c TransA==FflasTrans then \f$\mathrm{op}(A)=A^T\f$.
     * \param Diag if \c Diag==FflasUnit then \p A is unit diagonal.
     * \param Uplo if \c Uplo==FflasUpper then \p A is upper triangular
     */
    template<class Field>
    void
    ftrmv (const Field& F, const FFLAS_UPLO Uplo,
           const FFLAS_TRANSPOSE TransA, const FFLAS_DIAG Diag,
           const size_t N,typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::Element_ptr X, int incX){
        // Defaulting to ftrmm for the moment, waiting for specialized implem.
        ftrmm (F, FFLAS::FflasLeft, Uplo, TransA, Diag, N, 1, F.one, A, lda, X, incX);
    }

} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_level2_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
