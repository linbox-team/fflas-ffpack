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
        if (!K) FFLAS::fscalin(F, N,N, beta, C, ldc);

        typename Field::Element two,itwo; F.init(two); F.init(itwo);
        typename Field::Element_ptr diagC(NULL);

        if (F.characteristic() != 2) {
                // Cii <- Cii/2
            F.init(two, 2);
            F.inv(itwo,two);
            for(size_t i=0; i<N; ++i)
                F.mulin(*(C+i*ldc+i), itwo);
        } else {
                // store diagC <- beta Cii
            diagC = FFLAS::fflas_new(F,N);
            for(size_t i=0; i<N; ++i)
                F.mul(diagC[i],beta,*(C+i*ldc+i));
        }

            // Set unused triangular part of C to zero
        if (UpLo == FflasUpper) {
                // Lower part set to zero
            for(size_t i=0; i<N; ++i)
                for(size_t j=0; j<i; ++j)
                    F.assign(*(C+i*ldc+j), F.zero);
        } else {
                // Upper part set to zero
            for(size_t i=0; i<N; ++i)
                for(size_t j=i+1; j<N; ++j)
                    F.assign(*(C+i*ldc+j), F.zero);
        }

        FFLAS_TRANSPOSE oppTrans;
        if (trans==FflasNoTrans) {oppTrans=FflasTrans;}
        else {oppTrans=FflasNoTrans;}

            // alpha(A B^T) + beta(C)
        fgemm (F, trans, oppTrans, N, N, K, alpha, A, lda, B, ldb, beta, C, ldc);

            // adds alpha(B A^T) from the unused part
        if (UpLo == FflasUpper) {
            for(size_t i=0; i<N; ++i)
                for(size_t j=0; j<i; ++j)
                    F.addin(*(C+j*ldc+i), *(C+i*ldc+j));
        } else {
            for(size_t i=0; i<N; ++i)
                for(size_t j=i+1; j<N; ++j)
                    F.addin(*(C+j*ldc+i), *(C+i*ldc+j));
        }

        if (F.characteristic() != 2) {
                // Cii <- 2*Cii
            for(size_t i=0; i<N; ++i)
                F.addin(*(C+i*ldc+i), *(C+i*ldc+i));
        } else {
                // restore Cii <- beta Cii
            for(size_t i=0; i<N; ++i)
                F.assign(*(C+i*ldc+i),diagC[i]);
            FFLAS::fflas_delete(diagC);
        }

        return C;
    }
} // namespace FFLAS
#endif //__FFLASFFPACK_fflas_fsyr2k_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
