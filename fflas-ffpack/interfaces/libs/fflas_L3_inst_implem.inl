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
    template INST_OR_DECL
    void
    ftrsm (const FFLAS_FIELD <FFLAS_ELT>& F, const FFLAS_SIDE Side,
           const FFLAS_UPLO Uplo,
           const FFLAS_TRANSPOSE TransA,
           const FFLAS_DIAG Diag,
           const size_t M, const size_t N,
           const FFLAS_ELT alpha,
#ifdef __FFLAS__TRSM_READONLY
           const FFLAS_ELT* A,
#else
           FFLAS_ELT* A,
#endif
           const size_t lda,
           FFLAS_ELT* B, const size_t ldb);

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
    template INST_OR_DECL
    void
    ftrmm (const FFLAS_FIELD <FFLAS_ELT>& F, const FFLAS_SIDE Side,
           const FFLAS_UPLO Uplo,
           const FFLAS_TRANSPOSE TransA,
           const FFLAS_DIAG Diag,
           const size_t M, const size_t N,
           const FFLAS_ELT alpha,
           const FFLAS_ELT* A, const size_t lda,
           FFLAS_ELT* B, const size_t ldb);

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
    template INST_OR_DECL
    FFLAS_ELT* fgemm( const FFLAS_FIELD <FFLAS_ELT>& F,
                      const FFLAS_TRANSPOSE ta,
                      const FFLAS_TRANSPOSE tb,
                      const size_t m, const size_t n, const size_t k,
                      const FFLAS_ELT alpha,
                      const FFLAS_ELT* A, const size_t lda,
                      const FFLAS_ELT* B, const size_t ldb,
                      const FFLAS_ELT beta,
                      FFLAS_ELT* C, const size_t ldc);

    template INST_OR_DECL
    FFLAS_ELT*
    fgemm( const FFLAS_FIELD <FFLAS_ELT>& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m,
           const size_t n,
           const size_t k,
           const FFLAS_ELT alpha,
           const FFLAS_ELT* A, const size_t lda,
           const FFLAS_ELT* B, const size_t ldb,
           const FFLAS_ELT beta,
           FFLAS_ELT* C, const size_t ldc,
           const ParSeqHelper::Sequential seq);

    template INST_OR_DECL
    FFLAS_ELT*
    fgemm( const FFLAS_FIELD <FFLAS_ELT>& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m,
           const size_t n,
           const size_t k,
           const FFLAS_ELT alpha,
           const FFLAS_ELT* A, const size_t lda,
           const FFLAS_ELT* B, const size_t ldb,
           const FFLAS_ELT beta,
           FFLAS_ELT* C, const size_t ldc,
           const ParSeqHelper::Parallel<CuttingStrategy::Recursive,StrategyParameter::TwoDAdaptive> par);

    template INST_OR_DECL
    FFLAS_ELT*
    fgemm( const FFLAS_FIELD <FFLAS_ELT>& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m,
           const size_t n,
           const size_t k,
           const FFLAS_ELT alpha,
           const FFLAS_ELT* A, const size_t lda,
           const FFLAS_ELT* B, const size_t ldb,
           const FFLAS_ELT beta,
           FFLAS_ELT* C, const size_t ldc,
           const ParSeqHelper::Parallel<CuttingStrategy::Block,StrategyParameter::Threads> par);


    /** @brief fsquare: Squares a matrix.
     * compute \f$ C \gets \alpha \mathrm{op}(A) \mathrm{op}(A) + \beta C\f$ over a FFLAS_FIELD <FFLAS_ELT> \p F
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
    template INST_OR_DECL
    FFLAS_ELT* fsquare (const FFLAS_FIELD <FFLAS_ELT>& F,
                        const FFLAS_TRANSPOSE ta,
                        const size_t n,
                        const FFLAS_ELT alpha,
                        const FFLAS_ELT* A, const size_t lda,
                        const FFLAS_ELT beta,
                        FFLAS_ELT* C, const size_t ldc);


} // FFLAS

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
