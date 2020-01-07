/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
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

/** @file fflas/fflas_level3.h
 * @brief  Matrix-Matrix operations
 * or anything of \f$>n^2\f$ complexity.
 */

#ifndef __FFLASFFPACK_fflas_fflas_level3_INL
#define __FFLASFFPACK_fflas_fflas_level3_INL

//#include <givaro/zring.h>

#include "fflas_bounds.inl"
#include "fflas_helpers.inl"

namespace FFLAS { namespace Protected {
    //-----------------------------------------------------------------------------
    // Some conversion functions
    //-----------------------------------------------------------------------------


    //---------------------------------------------------------------------
    // Finite Field matrix => double matrix
    // Special design for upper-triangular matrices
    //---------------------------------------------------------------------
    template<class Field>
    void MatF2MatD_Triangular (const Field& F,
                               Givaro::DoubleDomain::Element_ptr S, const size_t lds,
                               typename Field::ConstElement_ptr const E,
                               const size_t lde,
                               const size_t m, const size_t n)
    {

        typename Field::ConstElement_ptr Ei = E;
        Givaro::DoubleDomain::Element_ptr Si = S;
        size_t i=0, j;
        for ( ; i<m;++i, Ei+=lde, Si+=lds)
            for ( j=i; j<n;++j)
                F.convert(*(Si+j),*(Ei+j));
    }

    //---------------------------------------------------------------------
    // Finite Field matrix => float matrix
    // Special design for upper-triangular matrices
    //---------------------------------------------------------------------
    //! @todo do finit(...,FFLAS_TRANS,FFLAS_DIAG)
    //! @todo do fconvert(...,FFLAS_TRANS,FFLAS_DIAG)
    template<class Field>
    void MatF2MatFl_Triangular (const Field& F,
                                Givaro::FloatDomain::Element_ptr S, const size_t lds,
                                typename Field::ConstElement_ptr const E,
                                const size_t lde,
                                const size_t m, const size_t n)
    {

        typename Field::ConstElement_ptr Ei = E;
        Givaro::FloatDomain::Element_ptr Si = S;
        size_t i=0, j;
        for ( ; i<m;++i, Ei+=lde, Si+=lds)
            for ( j=i; j<n;++j)
                F.convert(*(Si+j),*(Ei+j));
    }

    /**
     * Computes the maximal size for delaying the modular reduction
     *         in a triangular system resolution.
     *
     *  Compute the maximal dimension k, such that a unit diagonal triangular
     *  system of dimension k can be solved over Z without overflow of the
     *  underlying floating point representation.
     *
     *  @bib
     *  - Dumas, Giorgi, Pernet 06, arXiv:cs/0601133.
     *
     * \param F Finite Field/Ring of the computation
     *
     */
    // Specialized routines for ftrsm
    template <class Element> class ftrsmLeftUpperNoTransNonUnit;
    template <class Element> class ftrsmLeftUpperNoTransUnit;
    template <class Element> class ftrsmLeftUpperTransNonUnit;
    template <class Element> class ftrsmLeftUpperTransUnit;
    template <class Element> class ftrsmLeftLowerNoTransNonUnit;
    template <class Element> class ftrsmLeftLowerNoTransUnit;
    template <class Element> class ftrsmLeftLowerTransNonUnit;
    template <class Element> class ftrsmLeftLowerTransUnit;
    template <class Element> class ftrsmRightUpperNoTransNonUnit;
    template <class Element> class ftrsmRightUpperNoTransUnit;
    template <class Element> class ftrsmRightUpperTransNonUnit;
    template <class Element> class ftrsmRightUpperTransUnit;
    template <class Element> class ftrsmRightLowerNoTransNonUnit;
    template <class Element> class ftrsmRightLowerNoTransUnit;
    template <class Element> class ftrsmRightLowerTransNonUnit;
    template <class Element> class ftrsmRightLowerTransUnit;

    // Specialized routines for ftrmm
    template <class Element> class ftrmmLeftUpperNoTransNonUnit;
    template <class Element> class ftrmmLeftUpperNoTransUnit;
    template <class Element> class ftrmmLeftUpperTransNonUnit;
    template <class Element> class ftrmmLeftUpperTransUnit;
    template <class Element> class ftrmmLeftLowerNoTransNonUnit;
    template <class Element> class ftrmmLeftLowerNoTransUnit;
    template <class Element> class ftrmmLeftLowerTransNonUnit;
    template <class Element> class ftrmmLeftLowerTransUnit;
    template <class Element> class ftrmmRightUpperNoTransNonUnit;
    template <class Element> class ftrmmRightUpperNoTransUnit;
    template <class Element> class ftrmmRightUpperTransNonUnit;
    template <class Element> class ftrmmRightUpperTransUnit;
    template <class Element> class ftrmmRightLowerNoTransNonUnit;
    template <class Element> class ftrmmRightLowerNoTransUnit;
    template <class Element> class ftrmmRightLowerTransNonUnit;
    template <class Element> class ftrmmRightLowerTransUnit;

} // protected
} // FFLAS

namespace FFLAS {

    //---------------------------------------------------------------------
    // Level 3 routines
    //---------------------------------------------------------------------
    // set by default for ftrsm to be thread safe
    // undef it at your own risk, and only if you run it in sequential
#define __FFLAS__TRSM_READONLY

    /** @brief ftrsm: <b>TR</b>iangular <b>S</b>ystem solve with <b>M</b>atrix.
     * Computes  \f$ B \gets \alpha \mathrm{op}(A^{-1}) B\f$ or  \f$B \gets \alpha B \mathrm{op}(A^{-1})\f$.
     * \param F field
     * \param Side if \c Side==FflasLeft then  \f$ B \gets \alpha \mathrm{op}(A^{-1}) B\f$ is computed.
     * \param Uplo if \c Uplo==FflasUpper then \p A is upper triangular
     * \param TransA if \c TransA==FflasTrans then \f$\mathrm{op}(A)=A^t\f$.
     * \param Diag if \c Diag==FflasUnit then \p A is unit.
     * \param M rows of \p B
     * \param N cols of \p B
     * @param alpha scalar
     * \param A triangular invertible matrix. If \c Side==FflasLeft then \p A is \f$N\times N\f$, otherwise \p A is \f$M\times M\f$
     * @param lda leading dim of \p A
     * @param B matrix of size \p MxN
     * @param ldb leading dim of \p B
     * @bug \f$\alpha\f$ must be non zero.
     */
    template<class Field>
    void
    ftrsm (const Field& F, const FFLAS_SIDE Side,
           const FFLAS_UPLO Uplo,
           const FFLAS_TRANSPOSE TransA,
           const FFLAS_DIAG Diag,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
#ifdef __FFLAS__TRSM_READONLY
           typename Field::ConstElement_ptr A,
#else
           typename Field::Element_ptr A,
#endif
           const size_t lda,
           typename Field::Element_ptr B, const size_t ldb);

    /** @brief ftrmm: <b>TR</b>iangular <b>M</b>atrix <b>M</b>ultiply.
     * Computes  \f$ B \gets \alpha \mathrm{op}(A) B\f$ or  \f$B \gets \alpha B \mathrm{op}(A)\f$.
     * @param F field
     * \param Side if \c Side==FflasLeft then  \f$ B \gets \alpha \mathrm{op}(A) B\f$ is computed.
     * \param Uplo if \c Uplo==FflasUpper then \p A is upper triangular
     * \param TransA if \c TransA==FflasTrans then \f$\mathrm{op}(A)=A^t\f$.
     * \param Diag if \c Diag==FflasUnit then \p A is implicitly unit.
     * \param M rows of \p B
     * \param N cols of \p B
     * @param alpha scalar
     * \param A triangular matrix. If \c Side==FflasLeft then \p A is \f$N\times N\f$, otherwise \p A is \f$M\times M\f$
     * @param lda leading dim of \p A
     * @param B matrix of size \p MxN
     * @param ldb leading dim of \p B
     */
    template<class Field>
    void
    ftrmm (const Field& F, const FFLAS_SIDE Side,
           const FFLAS_UPLO Uplo,
           const FFLAS_TRANSPOSE TransA,
           const FFLAS_DIAG Diag,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::Element_ptr B, const size_t ldb);

    /** @brief ftrmm: <b>TR</b>iangular <b>M</b>atrix <b>M</b>ultiply with 3 operands
     * Computes  \f$ C \gets \alpha \mathrm{op}(A) B + beta C\f$ or  \f$C \gets \alpha B \mathrm{op}(A) + beta C\f$.
     * @param F field
     * \param Side if \c Side==FflasLeft then  \f$ B \gets \alpha \mathrm{op}(A) B\f$ is computed.
     * \param Uplo if \c Uplo==FflasUpper then \p A is upper triangular
     * \param TransA if \c TransA==FflasTrans then \f$\mathrm{op}(A)=A^t\f$.
     * \param Diag if \c Diag==FflasUnit then \p A is implicitly unit.
     * \param M rows of \p B
     * \param N cols of \p B
     * @param alpha scalar
     * \param A triangular matrix. If \c Side==FflasLeft then \p A is \f$N\times N\f$, otherwise \p A is \f$M\times M\f$
     * @param lda leading dim of \p A
     * @param B matrix of size \p MxN
     * @param ldb leading dim of \p B
     * @param beta scalar
     * @param C matrix of size \p MxN
     * @param ldc leading dim of \p C
     */
    template<class Field>
    void
    ftrmm (const Field& F, const FFLAS_SIDE Side,
           const FFLAS_UPLO Uplo,
           const FFLAS_TRANSPOSE TransA,
           const FFLAS_DIAG Diag,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr B, const size_t ldb,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc);

    /** @brief  fsyrk: Symmetric Rank K update
     *
     * Computes the Lower or Upper triangular part of \f$C = \alpha A \times A^T + \beta C\f$ or \f$C = \alpha A^T \times A + \beta C\f$
     * \param F field.
     * \param UpLo whether to compute the upper or the lower triangular part of the symmetric matrix \p C
     * \param trans if \c ta==FflasNoTrans then comput \f$C = \alpha A \times A^T + \beta C\f$, else  \f$C = \alpha A^T \times A + \beta C\f$
     * \param n order of matrix \p C
     * \param k see \p A
     * \param alpha scalar
     * \param A \f$A\f$ is \f$n \times k\f$ or \f$A\f$ is \f$k \times n\f$
     * \param lda leading dimension of \p A
     * \param beta scalar
     * \param C \f$C\f$ is \f$n \times n\f$
     * \param ldc leading dimension of \p C
     * @warning \f$\alpha\f$ \e must be invertible
     */
    template<class Field>
    typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc);

    template<class Field, typename  FieldTrait>
    inline typename Field::Element_ptr
    fsyrk_strassen (const Field& F,
                    const FFLAS_UPLO UpLo,
                    const FFLAS_TRANSPOSE trans,
                    const size_t N,
                    const size_t K,
                    const typename Field::Element y1,
                    const typename Field::Element y2,
                    const typename Field::Element alpha,
                    typename Field::ConstElement_ptr A, const size_t lda,
                    const typename Field::Element beta,
                    typename Field::Element_ptr C, const size_t ldc,
                    MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait> & H);

    template<class Field>
    typename Field::Element_ptr
    fsyr2k (const Field& F,
            const FFLAS_UPLO UpLo,
            const FFLAS_TRANSPOSE trans,
            const size_t n,
            const size_t k,
            const typename Field::Element alpha,
            typename Field::ConstElement_ptr A, const size_t lda,
            typename Field::ConstElement_ptr B, const size_t ldb,
            const typename Field::Element beta,
            typename Field::Element_ptr C, const size_t ldc);

    /** @brief  fsyrk: Symmetric Rank K update with diagonal scaling
     *
     * Computes the Lower or Upper triangular part of
     * \f$C = \alpha A \times D \times A^T + \beta C\f$ or
     * \f$C = \alpha A^T \times D \times A + \beta C\f$ where \p D is a diagonal matrix.
     * Matrix \p A is updated into \f$ D\times A\f$ (if trans = FflasTrans) or
     * \f$ A\times D\f$ (if trans = FflasNoTrans).
     * \param F field.
     * \param UpLo whether to compute the upper or the lower triangular part of the symmetric
     *        matrix \p C
     * \param trans if \c ta==FflasNoTrans then compute \f$C = \alpha A \times A^T + \beta C\f$,
     *              else  \f$C = \alpha A^T \times A + \beta C\f$
     * \param n order of matrix \p C
     * \param k see \p A
     * \param alpha scalar
     * \param A \f$A\f$ is \f$n \times k\f$ or \f$A\f$ is \f$k \times n\f$
     * \param lda leading dimension of \p A
     * \param D \f$D\f$ is \f$k \times k\f$ diagonal matrix, stored as a vector of k coefficients
     * \param lda leading dimension of \p A
     * \param beta scalar
     * \param C \f$C\f$ is \f$n \times n\f$
     * \param ldc leading dimension of \p C
     * @warning \f$\alpha\f$ \e must be invertible
     */
    template<class Field>
    typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::ConstElement_ptr D, const size_t incD,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc, const size_t threshold=__FFLASFFPACK_FSYRK_THRESHOLD);
    template<class Field>
    typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::ConstElement_ptr D, const size_t incD,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           const ParSeqHelper::Sequential seq,
           const size_t threshold=__FFLASFFPACK_FSYRK_THRESHOLD);
    template<class Field, class Cut, class Param>
    typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::ConstElement_ptr D, const size_t incD,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           const ParSeqHelper::Parallel<Cut,Param> par,
           const size_t threshold=__FFLASFFPACK_FSYRK_THRESHOLD);
    /** @brief  fsyrk: Symmetric Rank K update with diagonal scaling
     *
     * Computes the Lower or Upper triangular part of
     * \f$C = \alpha A \times Delta D \times A^T + \beta C\f$ or
     * \f$C = \alpha A^T \times Delta D \times A + \beta C\f$ where \p D is a diagonal matrix
     * and \p Delta is a block diagonal with either 1 on the diagonal or 2x2 swap blocks
     * Matrix \p A is updated into \f$ D\times A\f$ (if trans = FflasTrans) or
     * \f$ A\times D\f$ (if trans = FflasNoTrans).
     * \param F field.
     * \param UpLo whether to compute the upper or the lower triangular part of the symmetric
     *        matrix \p C
     * \param trans if \c ta==FflasNoTrans then compute \f$C = \alpha A Delta D \times A^T + \beta C\f$,
     *              else  \f$C = \alpha A^T Delta D \times A + \beta C\f$
     * \param n see \p B
     * \param k see \p A
     * \param alpha scalar
     * \param A \f$A\f$ is \f$n \times k\f$ or \f$A\f$ is \f$k \times n\f$
     * \param lda leading dimension of \p A
     * \param D \f$D\f$ is \f$k \times k\f$ diagonal matrix, stored as a vector of k coefficients
     * \param twoBlocks a vector boolean indicating the beginning of each 2x2 blocs in Delta
     * \param lda leading dimension of \p A
     * \param beta scalar
     * \param C \f$C\f$ is \f$n \times n\f$
     * \param ldc leading dimension of \p C
     * @warning \f$\alpha\f$ \e must be invertible
     */
    template<class Field>
    typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::ConstElement_ptr D, const size_t incD,
           const std::vector<bool>& twoBlock,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc, const size_t threshold=__FFLASFFPACK_FSYRK_THRESHOLD);

    /** @brief  fsyr2k: Symmetric Rank 2K update
     *
     * Computes the Lower or Upper triangular part of \f$C = \alpha ( A \times B^T + B \times A^T) + \beta C\f$ or \f$C = \alpha ( A^T \times B + B^T \times A ) + \beta C\f$
     * \param F field.
     * \param UpLo whether to compute the upper or the lower triangular part of the symmetric matrix \p C
     * \param trans if \c ta==FflasNoTrans then compute \f$C = \alpha ( A \times B^T + B \times A^T ) + \beta C\f$, else  \f$C = \alpha ( A^T \times B + B^T \times A) + \beta C\f$
     * \param n order of matrix \p C
     * \param k see \p A
     * \param alpha scalar
     * \param A \f$A\f$ is \f$n \times k\f$ (FflasNoTrans) or \f$A\f$ is \f$k \times n\f$ (FflasTrans)
     * \param lda leading dimension of \p A
     * \param beta scalar
     * \param C \f$C\f$ is \f$n \times n\f$
     * \param ldc leading dimension of \p C
     * @warning \f$\alpha\f$ \e must be invertible
     */
    template<class Field>
    typename Field::Element_ptr
    fsyr2k (const Field& F,
            const FFLAS_UPLO UpLo,
            const FFLAS_TRANSPOSE trans,
            const size_t n,
            const size_t k,
            const typename Field::Element alpha,
            typename Field::ConstElement_ptr A, const size_t lda,
            typename Field::ConstElement_ptr B, const size_t ldb,
            const typename Field::Element beta,
            typename Field::Element_ptr C, const size_t ldc);

    /** @brief  fgemm: <b>F</b>ield <b>GE</b>neral <b>M</b>atrix <b>M</b>ultiply.
     *
     * Computes \f$C = \alpha \mathrm{op}(A) \times \mathrm{op}(B) + \beta C\f$
     * Automatically set Winograd recursion level
     * \param F field.
     * \param ta if \c ta==FflasTrans then \f$\mathrm{op}(A)=A^t\f$, else \f$\mathrm{op}(A)=A\f$,
     * \param tb same for matrix \p B
     * \param m see \p A
     * \param n see \p B
     * \param k see \p A
     * \param alpha scalar
     * \param beta scalar
     * \param A \f$\mathrm{op}(A)\f$ is \f$m \times k\f$
     * \param B \f$\mathrm{op}(B)\f$ is \f$k \times n\f$
     * \param C \f$C\f$ is \f$m \times n\f$
     * \param lda leading dimension of \p A
     * \param ldb leading dimension of \p B
     * \param ldc leading dimension of \p C
     * \param w recursive levels of Winograd's algorithm are used. No argument (or -1) does auto computation of \p w.
     * @warning \f$\alpha\f$ \e must be invertible
     */
    template<class Field>
    typename Field::Element_ptr
    fgemm( const Field& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr B, const size_t ldb,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc);

    template<typename Field>
    typename Field::Element_ptr
    fgemm( const Field& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr B, const size_t ldb,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           const ParSeqHelper::Sequential seq);

    template<typename Field, class Cut, class Param>
    typename Field::Element_ptr
    fgemm( const Field& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr B, const size_t ldb,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           const ParSeqHelper::Parallel<Cut,Param> par);

    template<typename Field>
    typename Field::Element_ptr
    pfgemm (const Field& F,
            const FFLAS_TRANSPOSE ta,
            const FFLAS_TRANSPOSE tb,
            const size_t m,
            const size_t n,
            const size_t k,
            const typename Field::Element alpha,
            typename Field::ConstElement_ptr A, const size_t lda,
            typename Field::ConstElement_ptr B, const size_t ldb,
            const typename Field::Element beta,
            typename Field::Element_ptr C, const size_t ldc,
            size_t numthreads = 0){
        PAR_BLOCK{
            size_t nt = numthreads ? numthreads : NUM_THREADS;
            ParSeqHelper::Parallel<CuttingStrategy::Block,StrategyParameter::Threads> par(nt);
            fgemm(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,par);
        }
        return C;
    }

    template<class Field>
    typename Field::Element*
    pfgemm_1D_rec( const Field& F,
                   const FFLAS_TRANSPOSE ta,
                   const FFLAS_TRANSPOSE tb,
                   const size_t m,
                   const size_t n,
                   const size_t k,
                   const typename Field::Element alpha,
                   const typename Field::Element_ptr A, const size_t lda,
                   const typename Field::Element_ptr B, const size_t ldb,
                   const typename Field::Element beta,
                   typename Field::Element * C, const size_t ldc, size_t seuil);

    template<class Field>
    typename Field::Element*
    pfgemm_2D_rec( const Field& F,
                   const FFLAS_TRANSPOSE ta,
                   const FFLAS_TRANSPOSE tb,
                   const size_t m,
                   const size_t n,
                   const size_t k,
                   const typename Field::Element alpha,
                   const typename Field::Element_ptr A, const size_t lda,
                   const typename Field::Element_ptr B, const size_t ldb,
                   const typename Field::Element beta,
                   typename Field::Element * C, const size_t ldc, size_t seuil);

    template<class Field>
    typename Field::Element*
    pfgemm_3D_rec( const Field& F,
                   const FFLAS_TRANSPOSE ta,
                   const FFLAS_TRANSPOSE tb,
                   const size_t m,
                   const size_t n,
                   const size_t k,
                   const typename Field::Element alpha,
                   const typename Field::Element_ptr A, const size_t lda,
                   const typename Field::Element_ptr B, const size_t ldb,
                   const typename Field::Element beta,
                   typename Field::Element_ptr C, const size_t ldc, size_t seuil, size_t * x);

    template<class Field>
    typename Field::Element_ptr
    pfgemm_3D_rec2( const Field& F,
                    const FFLAS_TRANSPOSE ta,
                    const FFLAS_TRANSPOSE tb,
                    const size_t m,
                    const size_t n,
                    const size_t k,
                    const typename Field::Element alpha,
                    const typename Field::Element_ptr A, const size_t lda,
                    const typename Field::Element_ptr B, const size_t ldb,
                    const typename Field::Element beta,
                    typename Field::Element_ptr C, const size_t ldc, size_t seuil, size_t *x);

    /** @brief  fgemm: <b>F</b>ield <b>GE</b>neral <b>M</b>atrix <b>M</b>ultiply.
     *
     * Computes \f$C = \alpha \mathrm{op}(A) \times \mathrm{op}(B) + \beta C\f$
     * Version with Helper. Input and Output are not supposed to be reduced.
     * \param F field.
     * \param ta if \c ta==FflasTrans then \f$\mathrm{op}(A)=A^t\f$, else \f$\mathrm{op}(A)=A\f$,
     * \param tb same for matrix \p B
     * \param m see \p A
     * \param n see \p B
     * \param k see \p A
     * \param alpha scalar
     * \param beta scalar
     * \param A \f$\mathrm{op}(A)\f$ is \f$m \times k\f$
     * \param B \f$\mathrm{op}(B)\f$ is \f$k \times n\f$
     * \param C \f$C\f$ is \f$m \times n\f$
     * \param lda leading dimension of \p A
     * \param ldb leading dimension of \p B
     * \param ldc leading dimension of \p C
     * \param H helper, driving the computation (algorithm, delayed modular reduction, switch of base type, etc)
     * @warning \f$\alpha\f$ \e must be invertible
     */
    // template<class Field, class AlgoT, class FieldTrait, class ParSeqTrait>
    // inline  typename Field::Element_ptr
    // fgemm (const Field& F,
    //        const FFLAS_TRANSPOSE ta,
    //        const FFLAS_TRANSPOSE tb,
    //        const size_t m, const size_t n, const size_t k,
    //        const typename Field::Element alpha,
    //        typename Field::Element_ptr A, const size_t lda,
    //        typename Field::Element_ptr B, const size_t ldb,
    //        const typename Field::Element beta,
    //        typename Field::Element_ptr C, const size_t ldc,
    //        MMHelper<Field, AlgoT, FieldTrait, ParSeqTrait> & H);

} // FFLAS

#include "fflas-ffpack/paladin/parallel.h"

namespace FFLAS {

    /** @brief fsquare: Squares a matrix.
     * compute \f$ C \gets \alpha \mathrm{op}(A) \mathrm{op}(A) + \beta C\f$ over a Field \p F
     * Avoid the conversion of B
     * @param ta  if \c ta==FflasTrans, \f$\mathrm{op}(A)=A^T\f$.
     * @param F field
     * @param n size of \p A
     * @param alpha scalar
     * @param beta scalar
     * @param A dense matrix of size \c nxn
     * @param lda leading dimension of \p A
     * @param C dense matrix of size \c nxn
     * @param ldc leading dimension of \p C
     */
    template<class Field>
    typename Field::Element_ptr fsquare (const Field& F,
                                         const FFLAS_TRANSPOSE ta,
                                         const size_t n,
                                         const typename Field::Element alpha,
                                         typename Field::ConstElement_ptr A,
                                         const size_t lda,
                                         const typename Field::Element beta,
                                         typename Field::Element_ptr C,
                                         const size_t ldc);


} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_level3_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
