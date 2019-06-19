/* ffpack_inst_implem.inl
 * Copyright (C) 2005 Clement Pernet
 *               2014 FFLAS-FFPACK group
 *               2015 FFLAS-FFPACK group
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

namespace FFPACK {

    void composePermutationsLLM (size_t * MathP,
                                 const size_t * P1,
                                 const size_t * P2,
                                 const size_t R, const size_t N);

    void composePermutationsLLL (size_t * P1,
                                 const size_t * P2,
                                 const size_t R, const size_t N);

    void composePermutationsMLM (size_t * MathP1,
                                 const size_t * P2,
                                 const size_t R, const size_t N);

    void cyclic_shift_mathPerm (size_t * P,  const size_t s);
    template<typename Base_t>
    void cyclic_shift_row_col(Base_t * A, size_t m, size_t n, size_t lda);
    template INST_OR_DECL
    void cyclic_shift_row(const FFLAS_FIELD<FFLAS_ELT>& F, FFLAS_ELT* A, size_t m, size_t n, size_t lda);
    template INST_OR_DECL
    void cyclic_shift_col(const FFLAS_FIELD<FFLAS_ELT>& F, FFLAS_ELT* A, size_t m, size_t n, size_t lda);


    template INST_OR_DECL
    void applyP( const FFLAS_FIELD<FFLAS_ELT>& F,
                 const FFLAS::FFLAS_SIDE Side,
                 const FFLAS::FFLAS_TRANSPOSE Trans,
                 const size_t M, const size_t ibeg, const size_t iend,
                 FFLAS_ELT* A, const size_t lda, const size_t * P );

    template INST_OR_DECL
    void fgetrs (const FFLAS_FIELD<FFLAS_ELT>& F,
                 const FFLAS::FFLAS_SIDE Side,
                 const size_t M, const size_t N, const size_t R,
                 FFLAS_ELT* A, const size_t lda,
                 const size_t *P, const size_t *Q,
                 FFLAS_ELT* B, const size_t ldb,
                 int * info);

    template INST_OR_DECL
    FFLAS_ELT* fgetrs (const FFLAS_FIELD<FFLAS_ELT>& F,
                       const FFLAS::FFLAS_SIDE Side,
                       const size_t M, const size_t N, const size_t NRHS, const size_t R,
                       FFLAS_ELT* A, const size_t lda,
                       const size_t *P, const size_t *Q,
                       FFLAS_ELT* X, const size_t ldx,
                       const FFLAS_ELT* B, const size_t ldb,
                       int * info);
    template INST_OR_DECL
    size_t fgesv (const FFLAS_FIELD<FFLAS_ELT>& F,
                  const FFLAS::FFLAS_SIDE Side,
                  const size_t M, const size_t N,
                  FFLAS_ELT* A, const size_t lda,
                  FFLAS_ELT* B, const size_t ldb,
                  int * info);

    template INST_OR_DECL
    size_t fgesv (const FFLAS_FIELD<FFLAS_ELT>& F,
                  const FFLAS::FFLAS_SIDE Side,
                  const size_t M, const size_t N, const size_t NRHS,
                  FFLAS_ELT* A, const size_t lda,
                  FFLAS_ELT* X, const size_t ldx,
                  const FFLAS_ELT* B, const size_t ldb,
                  int * info);

    template INST_OR_DECL
    void ftrtri (const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_UPLO Uplo, const FFLAS::FFLAS_DIAG Diag,
                 const size_t N, FFLAS_ELT* A, const size_t lda, const size_t threshold);


    template INST_OR_DECL
    void trinv_left( const FFLAS_FIELD<FFLAS_ELT>& F, const size_t N, const FFLAS_ELT* L, const size_t ldl,
                     FFLAS_ELT* X, const size_t ldx );

    template INST_OR_DECL
    void ftrtrm (const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_SIDE side, const FFLAS::FFLAS_DIAG diag, const size_t N,
                 FFLAS_ELT* A, const size_t lda);

    template INST_OR_DECL
    size_t PLUQ (const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_DIAG Diag,
                 const size_t M, const size_t N,
                 FFLAS_ELT* A, const size_t lda,
                 size_t*P, size_t *Q);

    template INST_OR_DECL
    size_t LUdivine (const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_DIAG Diag,  const FFLAS::FFLAS_TRANSPOSE trans,
                     const size_t M, const size_t N,
                     FFLAS_ELT* A, const size_t lda,
                     size_t* P, size_t* Qt,
                     const FFPACK_LU_TAG LuTag,
                     const size_t cutoff);

    template INST_OR_DECL
    size_t LUdivine_small (const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_DIAG Diag,  const FFLAS::FFLAS_TRANSPOSE trans,
                           const size_t M, const size_t N,
                           FFLAS_ELT* A, const size_t lda,
                           size_t* P, size_t* Q,
                           const FFPACK_LU_TAG LuTag);

    template INST_OR_DECL
    size_t LUdivine_gauss (const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_DIAG Diag,
                           const size_t M, const size_t N,
                           FFLAS_ELT* A, const size_t lda,
                           size_t* P, size_t* Q,
                           const FFPACK_LU_TAG LuTag);

    template INST_OR_DECL
    size_t RowEchelonForm (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
                           FFLAS_ELT* A, const size_t lda,
                           size_t* P, size_t* Qt, const bool transform,
                           const FFPACK_LU_TAG LuTag);

    template INST_OR_DECL
    size_t ReducedRowEchelonForm (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
                                  FFLAS_ELT* A, const size_t lda,
                                  size_t* P, size_t* Qt, const bool transform,
                                  const FFPACK_LU_TAG LuTag);

    template INST_OR_DECL
    size_t ColumnEchelonForm (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
                              FFLAS_ELT* A, const size_t lda,
                              size_t* P, size_t* Qt, const bool transform,
                              const FFPACK_LU_TAG LuTag);
    template INST_OR_DECL
    size_t ReducedColumnEchelonForm (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
                                     FFLAS_ELT* A, const size_t lda,
                                     size_t* P, size_t* Qt, const bool transform,
                                     const FFPACK_LU_TAG LuTag);

    template INST_OR_DECL
    FFLAS_ELT* Invert (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M,
                       FFLAS_ELT* A, const size_t lda,
                       int& nullity);

    template INST_OR_DECL
    FFLAS_ELT* Invert (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M,
                       const FFLAS_ELT* A, const size_t lda,
                       FFLAS_ELT* X, const size_t ldx,
                       int& nullity);

    template INST_OR_DECL
    FFLAS_ELT* Invert2( const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M,
                        FFLAS_ELT* A, const size_t lda,
                        FFLAS_ELT* X, const size_t ldx,
                        int& nullity);

    template INST_OR_DECL
    std::list<Givaro::Poly1Dom<FFLAS_FIELD<FFLAS_ELT> >::Element>&
    CharPoly (const Givaro::Poly1Dom<FFLAS_FIELD<FFLAS_ELT> >& R,
              std::list<Givaro::Poly1Dom<FFLAS_FIELD<FFLAS_ELT> >::Element> & charp,
              const size_t N, FFLAS_ELT* A, const size_t lda,
              FFLAS_FIELD<FFLAS_ELT>::RandIter& G,
              const FFPACK_CHARPOLY_TAG CharpTag, const size_t degree);

    // template INST_OR_DECL
    // std::list<Givaro::Poly1Dom<FFLAS_FIELD<FFLAS_ELT> >::Element> &
    // CharPoly<FFLAS_FIELD,Givaro::Poly1Dom<FFLAS_FIELD<FFLAS_ELT> > >
    //         (const FFLAS_FIELD<FFLAS_ELT>& F,
    // 		 std::list<Givaro::Poly1Dom<FFLAS_FIELD<FFLAS_ELT> >::Element> & charp,
    // 		 const size_t N, FFLAS_ELT* A, const size_t lda,
    // 		 const FFPACK_CHARPOLY_TAG CharpTag);

    template INST_OR_DECL
    Givaro::Poly1Dom<FFLAS_FIELD<FFLAS_ELT> >::Element&
    CharPoly(const Givaro::Poly1Dom<FFLAS_FIELD<FFLAS_ELT> >& R,
             Givaro::Poly1Dom<FFLAS_FIELD<FFLAS_ELT> >::Element& charp,
             const size_t N, FFLAS_ELT* A, const size_t lda,
             FFLAS_FIELD<FFLAS_ELT>::RandIter& G, const FFPACK_CHARPOLY_TAG CharpTag, const size_t degree);

    template INST_OR_DECL
    Givaro::Poly1Dom<FFLAS_FIELD<FFLAS_ELT> >::Element&
    CharPoly(const Givaro::Poly1Dom<FFLAS_FIELD<FFLAS_ELT> >& R,
             Givaro::Poly1Dom<FFLAS_FIELD<FFLAS_ELT> >::Element& charp,
             const size_t N, FFLAS_ELT* A, const size_t lda,
             const FFPACK_CHARPOLY_TAG CharpTag, const size_t degree);

    template INST_OR_DECL
    std::vector<FFLAS_ELT>& MinPoly( const FFLAS_FIELD<FFLAS_ELT>& F, std::vector<FFLAS_ELT>& minP, const size_t N,
                                     const FFLAS_ELT* A, const size_t lda,
                                     FFLAS_FIELD<FFLAS_ELT>::RandIter& G);
    template INST_OR_DECL
    std::vector<FFLAS_ELT>& MinPoly( const FFLAS_FIELD<FFLAS_ELT>& F, std::vector<FFLAS_ELT>& minP, const size_t N,
                                     const FFLAS_ELT* A, const size_t lda);

    template INST_OR_DECL
    std::vector<FFLAS_ELT>& MatVecMinPoly (const FFLAS_FIELD<FFLAS_ELT>& F, std::vector<FFLAS_ELT>& minP,
                                           const size_t N, const FFLAS_ELT* A, const size_t lda,
                                           const FFLAS_ELT* V, const size_t incv);

    template INST_OR_DECL
    size_t KrylovElim( const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
                       FFLAS_ELT* A, const size_t lda, size_t*P,
                       size_t *Q, const size_t deg, size_t *iterates, size_t * inviterates, const size_t maxit,size_t virt);

    template INST_OR_DECL
    size_t SpecRankProfile (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
                            FFLAS_ELT* A, const size_t lda, const size_t deg, size_t *rankProfile);

    template INST_OR_DECL
    size_t Rank (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
                 FFLAS_ELT* A, const size_t lda);

    template INST_OR_DECL
    bool IsSingular (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
                     FFLAS_ELT* A, const size_t lda);
    template INST_OR_DECL
    FFLAS_ELT& Det (const FFLAS_FIELD<FFLAS_ELT>& F, FFLAS_ELT& det, const size_t N,
                   FFLAS_ELT* A, const size_t lda, size_t *P, size_t *Q);
    template INST_OR_DECL
    FFLAS_ELT& Det (const FFLAS_FIELD<FFLAS_ELT>& F, FFLAS_ELT& det, const size_t N,
                   FFLAS_ELT* A, const size_t lda,
                   const FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive, FFLAS::StrategyParameter::Threads>& parH, size_t *P, size_t *Q);

    template INST_OR_DECL
    FFLAS_ELT* Solve( const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M,
                      FFLAS_ELT* A, const size_t lda,
                      FFLAS_ELT* x, const int incx,
                      const FFLAS_ELT* b, const int incb );
    template INST_OR_DECL
    void solveLB( const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_SIDE Side,
                  const size_t M, const size_t N, const size_t R,
                  FFLAS_ELT* L, const size_t ldl,
                  const size_t * Q,
                  FFLAS_ELT* B, const size_t ldb );

    template INST_OR_DECL
    void solveLB2( const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_SIDE Side,
                   const size_t M, const size_t N, const size_t R,
                   FFLAS_ELT* L, const size_t ldl,
                   const size_t * Q,
                   FFLAS_ELT* B, const size_t ldb );

    template INST_OR_DECL
    void RandomNullSpaceVector (const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_SIDE Side,
                                const size_t M, const size_t N,
                                FFLAS_ELT* A, const size_t lda,
                                FFLAS_ELT* X, const size_t incX);
    template INST_OR_DECL
    size_t NullSpaceBasis (const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_SIDE Side,
                           const size_t M, const size_t N,
                           FFLAS_ELT* A, const size_t lda,
                           FFLAS_ELT*& NS, size_t& ldn,
                           size_t& NSdim);
    template INST_OR_DECL
    size_t RowRankProfile (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
                           FFLAS_ELT* A, const size_t lda,
                           size_t* &rkprofile,
                           const FFPACK_LU_TAG LuTag);

    template INST_OR_DECL
    size_t ColumnRankProfile (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
                              FFLAS_ELT* A, const size_t lda,
                              size_t* &rkprofile,
                              const FFPACK_LU_TAG LuTag);

    void RankProfileFromLU (const size_t* P, const size_t N, const size_t R,
                            size_t* rkprofile, const FFPACK_LU_TAG LuTag);

    size_t LeadingSubmatrixRankProfiles (const size_t M, const size_t N, const size_t R,
                                         const size_t LSm, const size_t LSn,
                                         const size_t* P, const size_t* Q,
                                         size_t* RRP, size_t* CRP);
    template INST_OR_DECL
    size_t RowRankProfileSubmatrixIndices (const FFLAS_FIELD<FFLAS_ELT>& F,
                                           const size_t M, const size_t N,
                                           FFLAS_ELT* A, const size_t lda,
                                           size_t*& rowindices, size_t*& colindices,
                                           size_t& R);

    template INST_OR_DECL
    size_t ColRankProfileSubmatrixIndices (const FFLAS_FIELD<FFLAS_ELT>& F,
                                           const size_t M, const size_t N,
                                           FFLAS_ELT* A, const size_t lda,
                                           size_t*& rowindices, size_t*& colindices,
                                           size_t& R);
    template INST_OR_DECL
    size_t RowRankProfileSubmatrix (const FFLAS_FIELD<FFLAS_ELT>& F,
                                    const size_t M, const size_t N,
                                    FFLAS_ELT* A, const size_t lda,
                                    FFLAS_ELT*& X, size_t& R);
    template INST_OR_DECL
    size_t ColRankProfileSubmatrix (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t M, const size_t N,
                                    FFLAS_ELT* A, const size_t lda,
                                    FFLAS_ELT*& X, size_t& R);

    template INST_OR_DECL
    void getTriangular <FFLAS_FIELD<FFLAS_ELT> > (const FFLAS_FIELD<FFLAS_ELT> & F, const FFLAS::FFLAS_UPLO Uplo,
                                                  const FFLAS::FFLAS_DIAG diag,
                                                  const size_t M, const size_t N, const size_t R,
                                                  const FFLAS_ELT* A, const size_t lda,
                                                  FFLAS_ELT* T, const size_t ldt,
                                                  const bool OnlyNonZeroVectors);

    template INST_OR_DECL
    void getTriangular <FFLAS_FIELD<FFLAS_ELT> >(const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_UPLO Uplo,
                                                 const FFLAS::FFLAS_DIAG diag,
                                                 const size_t M, const size_t N, const size_t R,
                                                 FFLAS_ELT* A, const size_t lda);

    template INST_OR_DECL
    void getEchelonForm <FFLAS_FIELD<FFLAS_ELT> > (const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_UPLO Uplo,
                                                   const FFLAS::FFLAS_DIAG diag,
                                                   const size_t M, const size_t N, const size_t R, const size_t* P,
                                                   const FFLAS_ELT* A, const size_t lda,
                                                   FFLAS_ELT* T, const size_t ldt,
                                                   const bool OnlyNonZeroVectors,
                                                   const FFPACK_LU_TAG LuTag);
    template INST_OR_DECL
    void getEchelonForm <FFLAS_FIELD<FFLAS_ELT> > (const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_UPLO Uplo,
                                                   const FFLAS::FFLAS_DIAG diag,
                                                   const size_t M, const size_t N, const size_t R, const size_t* P,
                                                   FFLAS_ELT* A, const size_t lda,
                                                   const FFPACK_LU_TAG LuTag);

    template INST_OR_DECL
    void getEchelonTransform <FFLAS_FIELD<FFLAS_ELT> > (const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_UPLO Uplo,
                                                        const FFLAS::FFLAS_DIAG diag,
                                                        const size_t M, const size_t N, const size_t R, const size_t* P, const size_t* Q,
                                                        const FFLAS_ELT* A, const size_t lda,
                                                        FFLAS_ELT* T, const size_t ldt,
                                                        const FFPACK_LU_TAG LuTag);
    template INST_OR_DECL
    void getReducedEchelonForm<FFLAS_FIELD<FFLAS_ELT> > (const FFLAS_FIELD<FFLAS_ELT> & F, const FFLAS::FFLAS_UPLO Uplo,
                                                         const size_t M, const size_t N, const size_t R, const size_t* P,
                                                         const FFLAS_ELT* A, const size_t lda,
                                                         FFLAS_ELT* T, const size_t ldt,
                                                         const bool OnlyNonZeroVectors,
                                                         const FFPACK_LU_TAG LuTag);

    template INST_OR_DECL
    void getReducedEchelonForm<FFLAS_FIELD<FFLAS_ELT> > (const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_UPLO Uplo,
                                                         const size_t M, const size_t N, const size_t R, const size_t* P,
                                                         FFLAS_ELT* A, const size_t lda,
                                                         const FFPACK_LU_TAG LuTag);
    template INST_OR_DECL
    void getReducedEchelonTransform<FFLAS_FIELD<FFLAS_ELT> > (const FFLAS_FIELD<FFLAS_ELT>& F, const FFLAS::FFLAS_UPLO Uplo,
                                                              const size_t M, const size_t N, const size_t R, const size_t* P, const size_t* Q,
                                                              const FFLAS_ELT* A, const size_t lda,
                                                              FFLAS_ELT* T, const size_t ldt,
                                                              const FFPACK_LU_TAG LuTag);
    void PLUQtoEchelonPermutation (const size_t N, const size_t R, const size_t * P, size_t * outPerm);

    template INST_OR_DECL
    FFLAS_ELT* LQUPtoInverseOfFullRankMinor( const FFLAS_FIELD<FFLAS_ELT>& F, const size_t rank,
                                             FFLAS_ELT* A_factors, const size_t lda,
                                             const size_t* QtPointer,
                                             FFLAS_ELT* X, const size_t ldx);
} // FFPACK
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
