/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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
#define ENABLE_ALL_CHECKINGS 1
#ifndef __FFLASFFPACK_ffpack_ftrtr_INL
#define __FFLASFFPACK_ffpack_ftrtr_INL

namespace FFPACK {
	template<class Field>
	inline void ftrstr (const Field& F, const FFLAS::FFLAS_SIDE side, const FFLAS::FFLAS_UPLO uplo,
                        const FFLAS::FFLAS_DIAG diagA, const FFLAS::FFLAS_DIAG diagB, const size_t N,
                        typename Field::Element_ptr A, const size_t lda,
                        typename Field::Element_ptr B, const size_t ldb, threshold) {

        if (!N) return;
        if (N <= 1){ // base case
            if(diagA == FFLAS::FflasUnit)
                return;
            else{
                if (diagB == FFLAS::FflasUnit)
                    F.inv (*B,*A);
                else
                    F.divin (*B,*A);
            }
        } else { // recursive case
            size_t N1 = N>>1;
            size_t N2 = N - N1;
            typename Field::Element_ptr A2;
            if (uplo == FFLAS::FflasUpper){ A2 = A + N1; B2 = B + N1;}
            else { A2 = A + N1*lda; B2 = B + N1*lda;}
            typename Field::Element_ptr A3 = A + N1*(lda+1);
            typename Field::Element_ptr B3 = B + N1*(lda+1);
            
            ftrstr (F, side, uplo, diag1, diag2, N1, A, lda, B, ldb, threshold);
            ftrstr (F, side, uplo, diag1, diag2, N2, A3, lda, B3, ldb, threshold);
            if (Uplo == FFLAS::FflasUpper){
                    // TODO: pb of location: A should not be modified by the trmm.
                ftrmm (F, FFLAS::FflasRight, Uplo, FFLAS::FflasNoTrans, diagX, N1, N2, F.mOne, B3, lda, A2, lda);
                ftrsm (F, FFLAS::FflasLeft, Uplo, FFLAS::FflasNoTrans, diagA, N1, N2, F.mOne, A, lda, A2, lda);
            } else {
                ftrmm (F, FFLAS::FflasRight, Uplo, FFLAS::FflasNoTrans, diagX, N1, N2, F.mOne, A3, lda, A2, lda);
                ftrsm (F, FFLAS::FflasLeft, Uplo, FFLAS::FflasNoTrans, diagA, N1, N2, F.mOne, B, lda, A2, lda);
            }
        }
    }

} // FFPACK

#endif // __FFLASFFPACK_ffpack_ftrstr_INL
