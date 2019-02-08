/*
 * Copyright (C) 2017 FFLAS-FFACK group
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
#ifndef __FFLASFFPACK_ffpack_ftrssyr2k_INL
#define __FFLASFFPACK_ffpack_ftrssyr2k_INL

#include "fflas-ffpack/utils/fflas_io.h"

namespace FFPACK {
    template<class Field>
    inline void ftrssyr2k (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
                           const FFLAS::FFLAS_DIAG diagA, const size_t N,
                           typename Field::ConstElement_ptr A, const size_t lda,
                           typename Field::Element_ptr B, const size_t ldb, const size_t threshold){

        if (!N) return;
        if (N <= 1){ // base case TODO: write a basecase with dimension = threshold
            typename Field::Element two;
            F.init(two, 2);
            if(diagA == FFLAS::FflasNonUnit){
                F.mulin(two,*A);
            }
            F.divin (*B,two);
            return;
        } else { // recursive case
            size_t N1 = N>>1;
            size_t N2 = N - N1;
            size_t A2rowdim = (Uplo == FFLAS::FflasUpper)? N1 : N2;
            size_t A2coldim = (Uplo == FFLAS::FflasUpper)? N2 : N1;
            typename Field::ConstElement_ptr A2;
            typename Field::Element_ptr B2;
            if (Uplo == FFLAS::FflasUpper){ A2 = A + N1; B2 = B + N1;}
            else { A2 = A + N1*lda; B2 = B + N1*ldb;}
            typename Field::ConstElement_ptr A3 = A + N1*(lda+1);
            typename Field::Element_ptr B3 = B + N1*(ldb+1);

            /* Solving [ A1^T      ] [ X1 X2 ] + [ X1^T      ] [ A1 A2 ] = [ B1   B2 ]
             *         [ A2^T A3^T ] [    X3 ]   [ X2^T X3^T ] [    A3 ]   [ B2^T B3 ]
             */

            // Solving  A1^T X1 + X1^T A1 = B1 in B1
            ftrssyr2k (F, Uplo, diagA, N1, A, lda, B, ldb, threshold);
            if (Uplo == FFLAS::FflasUpper){
                // B2 <- B2 - X1^T . A2
                ftrmm (F, FFLAS::FflasLeft, Uplo, FFLAS::FflasTrans, FFLAS::FflasNonUnit, A2rowdim, A2coldim, F.mOne, B, ldb, A2, lda, F.one, B2, ldb);
                // B2 <- A1^-T B2
                ftrsm (F, FFLAS::FflasLeft, Uplo, FFLAS::FflasTrans, diagA, A2rowdim, A2coldim, F.one, A, lda, B2, ldb);
            } else {
                // B2 <- B2 - A2 . X1^T
                ftrmm (F, FFLAS::FflasRight, Uplo, FFLAS::FflasTrans, FFLAS::FflasNonUnit, A2rowdim, A2coldim, F.mOne, B, ldb, A2, lda, F.one, B2, ldb); // To be checked
                // B2 <- B2 A1^-T
                ftrsm (F, FFLAS::FflasRight, Uplo, FFLAS::FflasTrans, diagA, A2rowdim, A2coldim, F.one, A, lda, B2, ldb);
            }

            // B3 <- B3 - A2^T X2 - X2^T A2
            fsyr2k (F, Uplo, Uplo==FFLAS::FflasUpper?FFLAS::FflasTrans:FFLAS::FflasNoTrans, N2, N1, F.mOne, A2, lda, B2, ldb, F.one, B3, ldb);

            // Solving  A3^T X3 + X3^T A3 = B3 in B3
            ftrssyr2k (F, Uplo, diagA, N2, A3, lda, B3, ldb, threshold);
        }
    }

} // FFPACK

#endif // __FFLASFFPACK_ffpack_ftrssyr2k_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
