/*
 * Copyright (C) 2016 the FFLAS-FFPACK group
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

#ifndef __FFLASFFPACK_fflas_fsyrk_INL
#define __FFLASFFPACK_fflas_fsyrk_INL

namespace FFLAS {
    template<class Field>
    inline typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t N,
           const size_t K,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc){

        //@TODO: write an optimized iterative basecase
        if (N==1){ // Base case
            F.mulin (*C, beta);
            size_t incA = (trans==FFLAS::FflasNoTrans)?1:lda;
            F.axpyin (*C, alpha, fdot (F, K, A, incA, A, incA));
            return C;
        } else {
            size_t N1 = N>>1;
            size_t N2 = N - N1;
            // Comments written for the case UpLo==FflasUpper, trans==FflasNoTrans
            FFLAS_TRANSPOSE oppTrans;
            if (trans==FflasNoTrans) {oppTrans=FflasTrans;}
            else {oppTrans=FflasNoTrans;}

            typename Field::ConstElement_ptr A2 = A + N1*(trans==FflasNoTrans?lda:1);
            typename Field::Element_ptr C12 = C + N1;
            typename Field::Element_ptr C21 = C + N1*ldc;
            typename Field::Element_ptr C22 = C12 + N1*ldc;
            // C11 <- alpha A1 x A1^T + beta C11
            fsyrk (F, UpLo, trans, N1, K, alpha, A, lda, beta, C, ldc);
            // C22 <- alpha A2 x A2^T + beta C22
            fsyrk (F, UpLo, trans, N2, K, alpha, A2, lda, beta, C22, ldc);

            if (UpLo == FflasUpper) {
                // C12 <- alpha A1 * A2^T + beta C12
                fgemm (F, trans, oppTrans, N1, N2, K, alpha, A, lda, A2, lda, beta, C12, ldc);
            } else {
                // C21 <- alpha A2 * A1^T + beta C21
                fgemm (F, trans, oppTrans, N2, N1, K, alpha, A2, lda, A, lda, beta, C21, ldc);
            }
            return C;
        }
    }
    template<class Field>
    inline typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t N,
           const size_t K,
           const typename Field::Element alpha,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::ConstElement_ptr D, const size_t incD,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc, const size_t threshold){
        return fsyrk (F, UpLo, trans, N, K, alpha, A, lda, D, incD, beta, C, ldc, ParSeqHelper::Sequential(), threshold);
    }
    template<class Field>
    inline typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t N,
           const size_t K,
           const typename Field::Element alpha,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::ConstElement_ptr D, const size_t incD,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           const ParSeqHelper::Sequential seq,
           const size_t threshold){

        size_t incRow,incCol;
        FFLAS_TRANSPOSE oppTrans;
        if (trans==FflasNoTrans) {incRow=lda;incCol=1;oppTrans=FflasTrans;}
        else {incRow = 1; incCol = lda;oppTrans=FflasNoTrans;}

        if (N <= threshold){
            typename Field::Element_ptr temp = fflas_new (F, N, K);
            size_t ldt;
            if (trans==FFLAS::FflasNoTrans){
                ldt =K;
                fassign(F, N, K, A, lda, temp, ldt);
            } else{
                ldt = N;
                fassign(F, K, N, A, lda, temp, ldt);
            }
            typename Field::Element_ptr Ai = A;
            typename Field::ConstElement_ptr Di = D;
            for (; Ai != A + K*incCol; Ai += incCol, Di+=incD){
                fscalin (F, N, *Di, Ai, incRow);
            }
            FFLAS::fgemm (F, trans, oppTrans, N, N, K, alpha, A, lda, temp, ldt, beta, C, ldc);
            FFLAS::fflas_delete(temp);
            return C;

        }else {
            size_t N1 = N>>1;
            size_t N2 = N - N1;
            // Comments written for the case UpLo==FflasUpper, trans==FflasNoTrans

            typename Field::Element_ptr A2 = A + N1*incRow;
            typename Field::Element_ptr C12 = C + N1;
            typename Field::Element_ptr C21 = C + N1*ldc;
            typename Field::Element_ptr C22 = C12 + N1*ldc;

            typename Field::Element_ptr temp = fflas_new (F, std::max(N2,N1),K);
            size_t ldt, incRowT,incColT;
            if (trans==FflasNoTrans) {ldt=K; incRowT=ldt; incColT=1;}
            else {ldt = N2; incRowT=1; incColT=ldt;}

            // temp <- A2 x D1
            typename Field::Element_ptr Ai = A2, Ti = temp;
            typename Field::ConstElement_ptr Di = D;
            for (; Ai != A2 + K*incCol; Ai += incCol, Ti += incColT, Di+=incD){
                fscal (F, N2, *Di, Ai, incRow, Ti, incRowT);
            }
            if (UpLo == FflasUpper) {
                // C12 <- alpha A1 x temp^T + beta C12
                fgemm (F, trans, oppTrans, N1, N2, K, alpha, A, lda, temp, ldt, beta, C12, ldc);
            } else {
                // C21 <- alpha temp x A11^T + beta C21
                fgemm (F, trans, oppTrans, N2, N1, K, alpha, temp, ldt, A, lda, beta, C21, ldc);
            }
            fflas_delete (temp);

            // C11 <- alpha A1 x D1 x A1^T + beta C11 and A1 <- A1 x D1
            fsyrk (F, UpLo, trans, N1, K, alpha, A, lda, D, incD, beta, C, ldc, threshold);
            // C22 <- alpha A2 x D1 x A2^T + beta C22 and A2 <- A2 x D1
            fsyrk (F, UpLo, trans, N2, K, alpha, A2, lda, D, incD, beta, C22, ldc, threshold);

            return C;
        }
    }
    template<class Field, class Cut, class Param>
    inline typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t N,
           const size_t K,
           const typename Field::Element alpha,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::ConstElement_ptr D, const size_t incD,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           const ParSeqHelper::Parallel<Cut,Param> par,
           const size_t threshold){

        size_t incRow,incCol;
        FFLAS_TRANSPOSE oppTrans;
        if (trans==FflasNoTrans) {incRow=lda;incCol=1;oppTrans=FflasTrans;}
        else {incRow = 1; incCol = lda;oppTrans=FflasNoTrans;}

        if (N <= threshold){
            typename Field::Element_ptr temp = fflas_new (F, N, K);
            size_t ldt;
            if (trans==FFLAS::FflasNoTrans){
                ldt =K;
                fassign(F, N, K, A, lda, temp, ldt);
            } else{
                ldt = N;
                fassign(F, K, N, A, lda, temp, ldt);
            }
            typename Field::Element_ptr Ai = A;
            typename Field::ConstElement_ptr Di = D;
            for (; Ai != A + K*incCol; Ai += incCol, Di+=incD){
                fscalin (F, N, *Di, Ai, incRow);
            }
            FFLAS::fgemm (F, trans, oppTrans, N, N, K, alpha, A, lda, temp, ldt, beta, C, ldc, par);
            FFLAS::fflas_delete(temp);
            return C;

        }else {
            size_t N1 = N>>1;
            size_t N2 = N - N1;
            // Comments written for the case UpLo==FflasUpper, trans==FflasNoTrans

            typename Field::Element_ptr A2 = A + N1*incRow;
            typename Field::Element_ptr C12 = C + N1;
            typename Field::Element_ptr C21 = C + N1*ldc;
            typename Field::Element_ptr C22 = C12 + N1*ldc;

            size_t nt = par.numthreads();
            if (nt == 1)
                return fsyrk(F, UpLo, trans, N, K, alpha, A, lda, D, incD, beta, C, ldc, ParSeqHelper::Sequential(), threshold);
            size_t nt2 = nt >> 1;
            size_t ntr = nt - nt2;
            ParSeqHelper::Parallel<Cut, Param> ps_rec1(nt2);
            ParSeqHelper::Parallel<Cut, Param> ps_rec2(ntr);
            ParSeqHelper::Parallel<CuttingStrategy::Block,StrategyParameter::Threads> ps_fgemm (nt);

            typename Field::Element_ptr temp = fflas_new (F, std::max(N2,N1),K);
            size_t ldt, incRowT,incColT;
            if (trans==FflasNoTrans) {ldt=K; incRowT=ldt; incColT=1;}
            else {ldt = N2; incRowT=1; incColT=ldt;}

            // temp <- A2 x D1
            typename Field::Element_ptr Ai = A2, Ti = temp;
            typename Field::ConstElement_ptr Di = D;
            //TODO parallelize fscal
            //TODO more cache efficient fscal
            for (; Ai != A2 + K*incCol; Ai += incCol, Ti += incColT, Di+=incD){
                fscal (F, N2, *Di, Ai, incRow, Ti, incRowT);
            }
            //TODO parallelize fgemm
            if (UpLo == FflasUpper) {
                // C12 <- alpha A1 x temp^T + beta C12
                fgemm (F, trans, oppTrans, N1, N2, K, alpha, A, lda, temp, ldt, beta, C12, ldc, ps_fgemm);
            } else {
                // C21 <- alpha temp x A11^T + beta C21
                fgemm (F, trans, oppTrans, N2, N1, K, alpha, temp, ldt, A, lda, beta, C21, ldc, ps_fgemm);
            }
            fflas_delete (temp);

            SYNCH_GROUP(
                        // C11 <- alpha A1 x D1 x A1^T + beta C11 and A1 <- A1 x D1
                        TASK(MODE(READ(D[0]) READWRITE(A[0]) WRITE(C[0]) CONSTREFERENCE(A, D, C, F, ps_rec1)),
                             fsyrk (F, UpLo, trans, N1, K, alpha, A, lda, D, incD, beta, C, ldc, ps_rec1, threshold));
                        // C22 <- alpha A2 x D1 x A2^T + beta C22 and A2 <- A2 x D1
                        TASK(MODE(READ(D[0]) READWRITE(A2[0]) WRITE(C22[0]) CONSTREFERENCE(A2, D, C22, F, ps_rec2)),
                             fsyrk (F, UpLo, trans, N2, K, alpha, A2, lda, D, incD, beta, C22, ldc, ps_rec2, threshold));
                       );
            return C;
        }
    }

    template<class Field>
    inline typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t N,
           const size_t K,
           const typename Field::Element alpha,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::ConstElement_ptr D, const size_t incD,
           const std::vector<bool>& twoBlocks,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc, const size_t threshold){

        size_t incRow,incCol;
        FFLAS_TRANSPOSE oppTrans;
        if (trans==FflasNoTrans) {incRow=lda;incCol=1;oppTrans=FflasTrans;}
        else {incRow = 1; incCol = lda;oppTrans=FflasNoTrans;}

        if (N <= threshold){
            typename Field::Element_ptr temp = fflas_new (F, N, K);
            size_t ldt,incRowT,incColT;
            if (trans==FFLAS::FflasNoTrans){
                ldt =K;
                incRowT=ldt; incColT=1;
                fassign(F, N, K, A, lda, temp, ldt);
            } else{
                ldt = N;
                incRowT=1; incColT=ldt;
                fassign(F, K, N, A, lda, temp, ldt);
            }
            typename Field::Element_ptr Ai = A;
            typename Field::Element_ptr tempi = temp;
            typename Field::ConstElement_ptr Di = D;
            for (size_t i=0; i<K;i++, Ai += incCol, Di+=incD, tempi+=incColT){
                if (!twoBlocks[i])
                    fscalin (F, N, *Di, Ai, incRow);
                else{
                    fscal (F, N, *Di, Ai+incCol,incRow, Ai,incRow);
                    fscal (F, N, *Di, tempi,incRowT, Ai+incCol,incRow);
                    Ai+=incCol; Di+=incD; tempi+=incColT;i++;
                }
            }
            FFLAS::fgemm (F, trans, oppTrans, N, N, K, alpha, A, lda, temp, ldt, beta, C, ldc);
            FFLAS::fflas_delete(temp);
            return C;

        }else {
            size_t N1 = N>>1;
            if (twoBlocks[N1-1]) N1++; // don't split a 2x2 block
            size_t N2 = N - N1;
            // Comments written for the case UpLo==FflasUpper, trans==FflasNoTrans

            typename Field::Element_ptr A2 = A + N1*incRow;
            typename Field::Element_ptr C12 = C + N1;
            typename Field::Element_ptr C21 = C + N1*ldc;
            typename Field::Element_ptr C22 = C12 + N1*ldc;

            typename Field::Element_ptr temp = fflas_new (F, std::max(N2,N1),K);
            size_t ldt, incRowT,incColT;
            if (trans==FflasNoTrans) {ldt=K; incRowT=ldt; incColT=1;}
            else {ldt = N2; incRowT=1; incColT=ldt;}

            // temp <- A2 x D1
            typename Field::Element_ptr Ai = A2, Ti = temp;
            typename Field::ConstElement_ptr Di = D;
            for (size_t i=0; i<K; Ai += incCol, Ti += incColT, Di+=incD,i++){
                if (!twoBlocks[i])
                    fscal (F, N2, *Di, Ai, incRow, Ti, incRowT);
                else {
                    fscal (F, N2, *Di, Ai, incRow, Ti+incColT, incRowT);
                    fscal (F, N2, *Di, Ai+incCol, incRow, Ti, incRowT);
                    Ti+=incColT; Ai+=incCol; Di+=incD; i++;
                }

            }
            if (UpLo == FflasUpper) {
                // C12 <- alpha A1 x temp^T + beta C12
                fgemm (F, trans, oppTrans, N1, N2, K, alpha, A, lda, temp, ldt, beta, C12, ldc);
            } else {
                // C21 <- alpha temp x A11^T + beta C21
                fgemm (F, trans, oppTrans, N2, N1, K, alpha, temp, ldt, A, lda, beta, C21, ldc);
            }
            fflas_delete (temp);

            // C11 <- alpha A1 x D1 x A1^T + beta C11 and A1 <- A1 x D1
            fsyrk (F, UpLo, trans, N1, K, alpha, A, lda, D, incD, twoBlocks, beta, C, ldc, threshold);
            // C22 <- alpha A2 x D1 x A2^T + beta C22 and A2 <- A2 x D1
            fsyrk (F, UpLo, trans, N2, K, alpha, A2, lda, D, incD, twoBlocks, beta, C22, ldc, threshold);

            return C;
        }
    }
}

#endif //__FFLASFFPACK_fflas_fsyrk_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
