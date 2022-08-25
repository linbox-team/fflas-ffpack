/*
 * Copyright (C) 2020 the FFLAS-FFPACK group
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@univ-grenoble-alpes.fr>
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

#ifndef __FFLASFFPACK_fflas_fsyr2k_INL
#define __FFLASFFPACK_fflas_fsyr2k_INL
#include <fflas-ffpack/utils/fflas_io.h>
namespace FFLAS {
    template<class Field>
    inline typename Field::Element_ptr
    fsyr2k (const Field& F,
            const FFLAS_UPLO UpLo,
            const FFLAS_TRANSPOSE trans,
            const size_t N,
            const size_t K,
            const typename Field::Element alpha,
            typename Field::ConstElement_ptr A, const size_t lda,
            typename Field::ConstElement_ptr B, const size_t ldb,
            const typename Field::Element beta,
            typename Field::Element_ptr C, const size_t ldc){

        if (!N) return C;
        if (!K) fscalin(F, N,N, beta, C, ldc);

        typename Field::Element_ptr diagC(nullptr);

        if (F.characteristic() != 2) {
                // diagonal Cii <- Cii/2
            typename Field::Element itwo; F.init(itwo, 2); F.invin(itwo);
            fscalin(F, N, itwo, C, ldc+1);
        } else {
                // store diagonal diagC <- beta Cii
            diagC = FFLAS::fflas_new(F,N);
            fscal(F, N, beta, C, ldc+1, diagC, 1);
        }

            // Set unused triangular part of C to zero
        if (UpLo == FflasUpper) {
                // Lower part set to zero
            for(size_t i=0; i<N; ++i)
                fzero(F, i, C+i*ldc, 1);
        } else {
                // Upper part set to zero
            for(size_t i=0; i<N; ++i)
                fzero(F, i, C+i, ldc);
        }

        FFLAS_TRANSPOSE oppTrans;
        if (trans==FflasNoTrans) {oppTrans=FflasTrans;}
        else {oppTrans=FflasNoTrans;}

            // alpha(A B^T) + beta(C)
        fgemm (F, trans, oppTrans, N, N, K, alpha, A, lda, B, ldb, beta, C, ldc);

                // adds alpha(B A^T) from the unused part
        if (F.characteristic() != 2) {
                // diagonal included (so that initial div by 2 is compensated)
            if (UpLo == FflasUpper) {
                for(size_t i=0; i<N; ++i)
                    faddin(F, i+1, C+i*ldc, 1, C+i, ldc);
            } else {
                for(size_t i=0; i<N; ++i)
                    faddin(F, N-i, C+i*ldc+i, 1, C+i*ldc+i, ldc);
            }
        } else {
                // diagonal excluded
            if (UpLo == FflasUpper) {
                for(size_t i=0; i<N; ++i)
                    faddin(F, i, C+i*ldc, 1, C+i, ldc);
            } else {
                for(size_t i=0; i<N; ++i)
                    faddin(F, N-i-1, C+i*ldc+(i+1), 1, C+(i+1)*ldc+i, ldc);
            }
                // restore diagonal Cii <- beta Cii
            fassign(F, N, diagC, 1, C, ldc+1);
            FFLAS::fflas_delete(diagC);
        }

        return C;
    }
} // namespace FFLAS
#endif //__FFLASFFPACK_fflas_fsyr2k_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
