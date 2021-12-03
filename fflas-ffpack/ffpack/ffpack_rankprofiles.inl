/* ffpack_rankprofiles.inl
 * Copyright (C) 2015 FFLAS-FFACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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

#ifndef __FFLASFFPACK_ffpack_rank_profiles_INL
#define __FFLASFFPACK_ffpack_rank_profiles_INL

namespace FFPACK{
    template <class Field>
    inline size_t RowRankProfile (const Field& F, const size_t M, const size_t N,
                                  typename Field::Element_ptr A, const size_t lda,
                                  size_t* &rkprofile, const FFPACK_LU_TAG LuTag){
        FFLAS::ParSeqHelper::Sequential seqH;
        return FFPACK::RowRankProfile (F, M, N, A, lda, rkprofile, LuTag, seqH);
    }

    template <class Field>
    inline size_t pRowRankProfile (const Field& F, const size_t M, const size_t N,
                                  typename Field::Element_ptr A, const size_t lda,
                                   size_t* &rkprofile, size_t numthreads, const FFPACK_LU_TAG LuTag){
        size_t r;
        PAR_BLOCK{
            size_t nt = numthreads ? numthreads : NUM_THREADS;
            FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::Threads> parH(nt);
            r = FFPACK::RowRankProfile (F, M, N, A, lda, rkprofile, LuTag, parH);
        }
        return r;
    }

    template <class Field, class PSHelper>
    inline size_t RowRankProfile (const Field& F, const size_t M, const size_t N,
                                  typename Field::Element_ptr A, const size_t lda,
                                  size_t* &rkprofile, const FFPACK_LU_TAG LuTag, PSHelper& psH){


        size_t *P = FFLAS::fflas_new<size_t>((LuTag==FfpackSlabRecursive)?N:M);
        size_t *Q = FFLAS::fflas_new<size_t>((LuTag==FfpackSlabRecursive)?M:N);
        size_t R;

        if (LuTag == FfpackSlabRecursive){
            R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, N, A, lda, P, Q);
            std::swap(P,Q);
        } else
            R = PLUQ (F, FFLAS::FflasNonUnit, M, N, A, lda, P, Q, psH);

        rkprofile = FFLAS::fflas_new<size_t> (R);

        RankProfileFromLU (P, M, R, rkprofile, LuTag);

        FFLAS::fflas_delete (Q);
        FFLAS::fflas_delete (P);
        return R;
    }

    template <class Field>
    inline size_t ColumnRankProfile (const Field& F, const size_t M, const size_t N,
                                     typename Field::Element_ptr A, const size_t lda,
                                     size_t* &rkprofile, const FFPACK_LU_TAG LuTag){
        FFLAS::ParSeqHelper::Sequential seqH;
        return FFPACK::ColumnRankProfile (F, M, N, A, lda, rkprofile, LuTag, seqH);
    }

    template <class Field>
    inline size_t pColumnRankProfile (const Field& F, const size_t M, const size_t N,
                                     typename Field::Element_ptr A, const size_t lda,
                                      size_t* &rkprofile, size_t numthreads, const FFPACK_LU_TAG LuTag){
        size_t r;
        PAR_BLOCK{
            size_t nt = numthreads ? numthreads : NUM_THREADS;
            FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::Threads> parH(nt);
            r = FFPACK::ColumnRankProfile (F, M, N, A, lda, rkprofile, LuTag, parH);
        }
        return r;
    }

    template <class Field, class PSHelper>
    inline size_t ColumnRankProfile (const Field& F, const size_t M, const size_t N,
                              typename Field::Element_ptr A, const size_t lda,
                              size_t* &rkprofile, const FFPACK_LU_TAG LuTag, PSHelper& psH){

        size_t *P = FFLAS::fflas_new<size_t>(M);
        size_t *Q = FFLAS::fflas_new<size_t>(N);
        size_t R;

        if (LuTag == FfpackSlabRecursive){
            R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasTrans, M, N, A, lda, P, Q);
        } else
            R = PLUQ (F, FFLAS::FflasNonUnit, M, N, A, lda, P, Q, psH);

        rkprofile = FFLAS::fflas_new<size_t> (R);

        RankProfileFromLU (Q, N, R, rkprofile, LuTag);

        FFLAS::fflas_delete (P);
        FFLAS::fflas_delete (Q);
        return R;
    }


    inline void RankProfileFromLU (const size_t* Q, const size_t N, const size_t R,
                                   size_t* rkprofile, const FFPACK_LU_TAG LuTag){

        if (LuTag == FfpackSlabRecursive)
            std::copy(Q, Q+R, rkprofile);
        else {
            size_t * RP = FFLAS::fflas_new<size_t>(N);
            for (size_t i=0;i < N; ++i)
                RP [i] = i;
            for (size_t i=0; i<N; ++i)
                if (Q[i] != i)
                    std::swap (RP [i], RP [Q [i]]);

            std::copy(RP, RP+R, rkprofile);
            std::sort (rkprofile, rkprofile + R);
            FFLAS::fflas_delete(RP);
        }
    }

    inline size_t LeadingSubmatrixRankProfiles (const size_t M, const size_t N, const size_t R,
                                                const size_t LSm, const size_t LSn,
                                                const size_t* P, const size_t* Q,
                                                size_t* RRP, size_t* CRP){
        size_t LSr=0; // rank of the LSm x LSn leading submatrix

        size_t* MathP = FFLAS::fflas_new<size_t>(M);
        size_t* MathQ = FFLAS::fflas_new<size_t>(N);

        LAPACKPerm2MathPerm (MathP, P, M);
        LAPACKPerm2MathPerm (MathQ, Q, N);
        for (size_t i = 0; i < R; i++)
            if (MathP[i] < LSm && MathQ[i] < LSn){
                RRP [LSr] = MathP[i];
                CRP [LSr] = MathQ[i];
                LSr++;
            }
        std::sort (RRP, RRP+LSr);
        std::sort (CRP, CRP+LSr);
        FFLAS::fflas_delete(MathP);
        FFLAS::fflas_delete(MathQ);
        return LSr;

    }


    template <class Field>
    size_t RowRankProfileSubmatrixIndices (const Field& F,
                                           const size_t M, const size_t N,
                                           typename Field::Element_ptr A,
                                           const size_t lda,
                                           size_t*& rowindices,
                                           size_t*& colindices,
                                           size_t& R)
    {
        size_t *P = FFLAS::fflas_new<size_t>(N);
        size_t *Q = FFLAS::fflas_new<size_t>(M);

        R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, N, A, lda, P, Q);
        rowindices = FFLAS::fflas_new<size_t>(M);
        colindices = FFLAS::fflas_new<size_t>(N);
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

        FFLAS::fflas_delete( P);
        FFLAS::fflas_delete( Q);

        return R;
    }

    template <class Field>
    size_t ColRankProfileSubmatrixIndices (const Field& F,
                                           const size_t M, const size_t N,
                                           typename Field::Element_ptr A,
                                           const size_t lda,
                                           size_t*& rowindices,
                                           size_t*& colindices,
                                           size_t& R)
    {
        size_t *P = FFLAS::fflas_new<size_t>(M);
        size_t *Q = FFLAS::fflas_new<size_t>(N);

        R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasTrans, M, N, A, lda, P, Q);
        rowindices = FFLAS::fflas_new<size_t>(M);
        colindices = FFLAS::fflas_new<size_t>(N);
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
        FFLAS::fflas_delete( P);
        FFLAS::fflas_delete( Q);

        return R;
    }

    template <class Field>
    size_t RowRankProfileSubmatrix (const Field& F,
                                    const size_t M, const size_t N,
                                    typename Field::Element_ptr A,
                                    const size_t lda,
                                    typename Field::Element_ptr& X, size_t& R)
    {

        size_t * rowindices, * colindices;

        typename Field::Element_ptr A2 = FFLAS::fflas_new (F, M, N) ;
        FFLAS::fassign(F,M,N,A,lda,A2,N);

        RowRankProfileSubmatrixIndices (F, M, N, A2, N, rowindices, colindices, R);

        X = FFLAS::fflas_new (F, R, R);
        for (size_t i=0; i<R; ++i)
            for (size_t j=0; j<R; ++j)
                F.assign (*(X + i*R + j), *(A + rowindices[i]*lda + colindices[j]));
        FFLAS::fflas_delete (A2);
        FFLAS::fflas_delete( rowindices);
        FFLAS::fflas_delete( colindices);
        return R;
    }

    template <class Field>
    size_t ColRankProfileSubmatrix (const Field& F, const size_t M, const size_t N,
                                    typename Field::Element_ptr A, const size_t lda,
                                    typename Field::Element_ptr& X, size_t& R)
    {

        size_t * rowindices, * colindices;

        typename Field::Element_ptr A2 = FFLAS::fflas_new (F, M, N);
        FFLAS::fassign(F,M,N,A,lda,A2,N);

        ColRankProfileSubmatrixIndices (F, M, N, A2, N, rowindices, colindices, R);

        X = FFLAS::fflas_new (F, R, R);
        for (size_t i=0; i<R; ++i)
            for (size_t j=0; j<R; ++j)
                F.assign (*(X + i*R + j), *(A + rowindices[i]*lda + colindices[j]));
        FFLAS::fflas_delete (A2);
        FFLAS::fflas_delete( colindices);
        FFLAS::fflas_delete( rowindices);
        return R;
    }

    template <class Field>
    typename Field::Element_ptr
    LQUPtoInverseOfFullRankMinor( const Field& F, const size_t rank,
                                  typename Field::Element_ptr A_factors, const size_t lda,
                                  const size_t* QtPointer,
                                  typename Field::Element_ptr X, const size_t ldx)
    {

        // upper entries are okay, just need to move up bottom ones
        const size_t* srcRow = QtPointer;
        for (size_t row=0; row<rank; row++, srcRow++)
            if (*srcRow != row) {
                typename Field::Element_ptr oldRow = A_factors + (*srcRow) * lda;
                typename Field::Element_ptr newRow = A_factors + row * lda;
                for (size_t col=0; col<row; col++, oldRow++, newRow++)
                    F.assign(*newRow, *oldRow);
            }

        // X <- (Qt.L.Q)^(-1)
        //invL( F, rank, A_factors, lda, X, ldx);
        ftrtri (F, FFLAS::FflasLower, FFLAS::FflasUnit, rank, A_factors, lda);
        FFLAS::fassign(F,rank,rank,A_factors,lda,X,ldx);

        // X = U^-1.X
        ftrsm( F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans,
               FFLAS::FflasNonUnit, rank, rank, F.one, A_factors, lda, X, ldx);

        return X;

    }

} // namespace FFPACK

#endif // __FFLASFFPACK_ffpack_rank_profiles_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
