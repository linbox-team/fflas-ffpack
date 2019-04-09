/*
 * Copyright (C) 2014 FFLAS-FFACK group
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

#ifndef __FFLASFFPACK_ffpack_fgetrs_INL
#define __FFLASFFPACK_ffpack_fgetrs_INL


namespace FFPACK {


    template <class Field>
    void
    fgetrs (const Field& F,
            const FFLAS::FFLAS_SIDE Side,
            const size_t M, const size_t N, const size_t R,
            typename Field::Element_ptr A, const size_t lda,
            const size_t *P, const size_t *Q,
            typename Field::Element_ptr B, const size_t ldb,
            int * info)
    {

        *info =0;
        if (Side == FFLAS::FflasLeft) { // Left looking solve A X = B

                // B <- P^T B
            applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N, 0, M, B, ldb, P);
                // B <- L^-1 B
            ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit,
                   R, N, F.one, A, lda , B, ldb);

                // Verifying that the system is consistent @todo: disable this check optionnally
            fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-R, N, R, F.mOne, A+R*lda, lda, B, ldb, F.one, B+R*ldb, ldb);
            if (! FFLAS::fiszero (F, M-R, N, B+R*ldb, ldb)) *info = 1;

                // B <- U^-1 B
            ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
                   R, N, F.one, A, lda , B, ldb);
                // B <- P^T B
            applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans, N, 0, M, B, ldb, Q);

        }
        else { // Right Looking X A = B

                // B <- B Q^T
            applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans, M, 0, N, B, ldb, Q);
                // B <- B U^-1
            ftrsm (F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
                   M, R, F.one, A, lda , B, ldb);
                //  Verifying that the system is consistent @todo: disable this check optionnally
            fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M, N-R, R, F.mOne, B, ldb, A+R, lda, F.one, B+R, ldb);
            if (! FFLAS::fiszero (F, M, N-R, B+R, ldb)) *info = 1;

                // B <- B L^-1
            ftrsm (F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit,
                   M, R, F.one, A, lda , B, ldb);

                // B <- B P^T
            applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, M, 0, N, B, ldb, P);
        }
    }

    template <class Field>
    typename Field::Element_ptr
    fgetrs (const Field& F,
            const FFLAS::FFLAS_SIDE Side,
            const size_t M, const size_t N, const size_t NRHS, const size_t R,
            typename Field::Element_ptr A, const size_t lda,
            const size_t *P, const size_t *Q,
            typename Field::Element_ptr X, const size_t ldx,
            typename Field::ConstElement_ptr B, const size_t ldb,
            int * info)
    {
        *info =0;
        typename Field::Element_ptr W;
        size_t ldw;

        if (Side == FFLAS::FflasLeft) { // Left looking solve A X = B

            if (M > N){ // Cannot copy B into X
                W = FFLAS::fflas_new (F, M, NRHS);
                ldw = NRHS;

                FFLAS::fassign(F,M,NRHS,B,ldb,W,ldw);
                    // B <- P^T B
                applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, NRHS, 0, M, W, ldw, P);
                    // B <- L^-1 B
                ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit,
                       R, NRHS, F.one, A, lda , W, ldw);

                    // Verifying that the system is consistent @todo: disable this check optionnally
                fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-R, NRHS, R, F.mOne, A+R*lda, lda, W, ldw, F.one, W+R*ldw, ldw);
                if (! FFLAS::fiszero (F, M-R, NRHS, W+R*ldw, ldw)) *info = 1;

                    // B <- U^-1 B
                ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
                       R, NRHS, F.one, A, lda , W, ldw);

                FFLAS::fassign (F,R,NRHS,W,ldw,X,ldx);
                FFLAS::fflas_delete (W);

                FFLAS::fzero (F, N-R, NRHS, X+R*ldx, ldx);

                    // B <- Q^T B
                applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans, NRHS, 0, N, X, ldx, Q);

            }
            else { // Copy B to X directly

                FFLAS::fassign (F,M,NRHS,B,ldb,X,ldx);

                    // B <- P^T B
                applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, NRHS, 0, M, X, ldx, P);
                    // B <- L^-1 B
                ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit,
                       R, NRHS, F.one, A, lda , X, ldx);

                    // Verifying that the system is consistent @todo: disable this check optionnally
                fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-R, NRHS, R, F.mOne, A+R*lda, lda, X, ldx, F.one, X+R*ldx, ldx);
                if (! FFLAS::fiszero (F, M-R, NRHS, X+R*ldx, ldx)) *info = 1;

                    // B <- U^-1 B
                ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
                       R, NRHS, F.one, A, lda , X, ldx);

                FFLAS::fzero (F, N-M, NRHS, X+M*ldx, ldx);

                    // B <- Q^T B
                applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans, NRHS, 0, N, X, ldx, Q);
            }
            return X;

        }
        else { // Right Looking X A = B

            if (M < N) {
                W = FFLAS::fflas_new (F, NRHS, N);
                ldw = N;
                FFLAS::fassign (F,NRHS, N, B, ldb, W, ldw);

                    // B <- B Q^T
                applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans, NRHS, 0, N, W, ldw, Q);

                    // B <- B U^-1
                ftrsm (F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
                       NRHS, R, F.one, A, lda , W, ldw);

                    // Verifying that the system is consistent @todo: disable this check optionnally
                fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, NRHS, N-R, R, F.mOne,
                       W, ldw, A+R, lda, F.one, W+R, ldw);

                if (! FFLAS::fiszero (F, NRHS, N-R, W+R, ldw)) *info = 1;

                    // B <- B L^-1
                ftrsm (F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit,
                       NRHS, R, F.one, A, lda , W, ldw);

                    // X <- B
                FFLAS::fassign (F, NRHS, R, W, ldw, X, ldx);
                FFLAS::fflas_delete (W);

                FFLAS::fzero (F, NRHS, M-R, X+R, ldx);

                    // X <- X P^T
                applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, NRHS, 0, M, X, ldx, P);

            } else { // M >=N
                FFLAS::fassign (F, NRHS, N, B, ldb, X, ldx);

                    // B <- B Q^T
                applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans, NRHS, 0, N, X, ldx, Q);

                    // B <- B U^-1
                ftrsm (F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
                       NRHS, R, F.one, A, lda , X, ldx);

                    // Verifying that the system is consistent @todo: disable this check optionnally
                fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, NRHS, N-R, R, F.mOne,
                       X, ldx, A+R, lda, F.one, X+R, ldx);

                if (! FFLAS::fiszero (F, NRHS, N-R, X+R, ldx)) *info = 1;

                // B <- B L^-1
                ftrsm (F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit,
                       NRHS, R, F.one, A, lda , X, ldx);

                FFLAS::fzero (F, NRHS, M-N, X+N, ldx);

                    // X <- X P^T
                applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, NRHS, 0, M, X, ldx, P);
            }
            return X;
        }
    }

} // FFPACK

#endif // __FFLASFFPACK_ffpack_fgetrs_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
