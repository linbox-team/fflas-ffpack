/*
 * Copyright (C) 2019 the FFLAS-FFPACK group
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

#ifndef __FFLASFFPACK_fflas_fsyrk_strassen_INL
#define __FFLASFFPACK_fflas_fsyrk_strassen_INL

namespace FFLAS {
    template<class Field>
    inline typename Field::Element_ptr
    fsyrk_strassen (const Field& F,
                    const FFLAS_UPLO UpLo,
                    const FFLAS_TRANSPOSE trans,
                    const size_t N,
                    const size_t K,
                    const typename Field::Element alpha,
                    typename Field::ConstElement_ptr A, const size_t lda,
                    const typename Field::Element beta,
                    typename Field::Element_ptr C, const size_t ldc){
        

            // S1 = A11 - A21 in C12

            // A22' = A22 Y in T1

            // S2 = A21 + A22' in T2

            // P4 = S1 x S2^T in  T3
        fgemm (beta=0);


            // P1 = A11 x A11^T in C12
        fsyrk_strassen(beta=0);

            // C11 = P1 + beta.C11 in C11
        faxpyin();
        
            // U3 = P2 + P1 = A12 x A12^T + P1 in C11
        fsyrk_strassen(beta=1);

            // S3 = A11 - S2 in T2
        fsubin();
        
            // U1 = P1 - P5 = -S3 x S3^T + P1 in C12
        fsyrk_strassen(beta=1);

            // U2 = U1 - P4 in C12
        fsubin();

            // U5 = U2 - P4^T + C22 in C22 (add)
        fsubin();
        faddin();

            // S4 = S3  + A12 x Y in T2
            // A12xY on the fly

            // - P3 = - A22 Y x S4^T + C21 in C21
        fgemm(beta);

            // U4 = U2 - P3  in C21
        faddin ();

////////////////////////////
            // version without accumulation (WIP)
            // S1 = A11 - A21

            // S2 = A21 + A22 x Y^T

            // S3 = A11 - S2

            // S4 = S3 x Y + A12

            // P1 = A11 x A11^T in C11

            // P2 = A12 x A12^T

            // P3 = A22 Y x S4^T

            // P4 = S1 x S2^T

            // P5 = S3 x S3^T

            // U1 = P1 - P5

            // U2 = U1 - P4

            // U3 = P1 + P2 in C11

            // U4 = U2 - P3 in C21

            // U5 = U2 - P4^T in C22
        
    }
}

#endif // __FFLASFFPACK_fflas_fsyrk_strassen_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
