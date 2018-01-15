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
#ifndef __FFLASFFPACK_ffpack_ftrstr_INL
#define __FFLASFFPACK_ffpack_ftrstr_INL

#include "fflas-ffpack/utils/fflas_io.h"

namespace FFPACK {
	template<class Field>
	inline void ftrstr (const Field& F, const FFLAS::FFLAS_SIDE side, const FFLAS::FFLAS_UPLO UpLo,
                        const FFLAS::FFLAS_DIAG diagA, const FFLAS::FFLAS_DIAG diagB, const size_t N,
                        typename Field::Element_ptr A, const size_t lda,
                        typename Field::Element_ptr B, const size_t ldb, const size_t threshold) {

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
                // Warning/TODO: only implemented for diagA=diagB for the moment
            FFLAS::FFLAS_SIDE oppSide = (side == FFLAS::FflasLeft) ? FFLAS::FflasRight: FFLAS::FflasLeft;
            size_t N1 = N>>1;
            size_t N2 = N - N1;
            typename Field::Element_ptr A2,B2;
            if (UpLo == FFLAS::FflasUpper){ A2 = A + N1; B2 = B + N1;}
            else { A2 = A + N1*lda; B2 = B + N1*ldb;}
            typename Field::Element_ptr A3 = A + N1*(lda+1);
            typename Field::Element_ptr B3 = B + N1*(ldb+1);
            
                /* Solving [ A1 A2 ] [ X1 X2 ] = [ B1 B2 ]
                 *         [    A3 ] [    X3 ] = [    B3 ]
                 */
                // B1 <- A1^-1 . B1
            // FFLAS::WriteMatrix(std::cerr<<"avant rec1 A = "<<std::endl,F,N,N,A,lda);
            ftrstr (F, side, UpLo, diagA, diagB, N1, A, lda, B, ldb, threshold);
            
            // FFLAS::WriteMatrix(std::cerr<<"after rec1 A = "<<std::endl,F,N,N,A,lda);

                // B3 <- A3^-1 . B3
            ftrstr (F, side, UpLo, diagA, diagB, N2, A3, lda, B3, ldb, threshold);

            // FFLAS::WriteMatrix(std::cerr<<"after rec2 A = "<<std::endl,F,N,N,A,lda);

                            
            if ((UpLo == FFLAS::FflasUpper && side==FFLAS::FflasLeft) ||
                (UpLo == FFLAS::FflasLower && side==FFLAS::FflasRight)){
                    // B2 <- B2 - A2 . B3 
                // FFLAS::WriteMatrix(std::cerr<<"before trmm A = "<<std::endl,F,N,N,A,lda);
                // FFLAS::WriteMatrix(std::cerr<<"before trmm B = "<<std::endl,F,N,N,B,ldb);
                ftrmm (F, oppSide, UpLo, FFLAS::FflasNoTrans, diagB, N1, N2, F.mOne, B3, ldb, A2, lda, F.one, B2, ldb);
                // FFLAS::WriteMatrix(std::cerr<<"after trmm B = "<<std::endl,F,N,N,B,ldb);
                    // B2 <- A1^-1 . B2
                ftrsm (F, side, UpLo, FFLAS::FflasNoTrans, diagA, N1, N2, F.one, A, lda, B2, ldb);
                // FFLAS::WriteMatrix(std::cerr<<"after trsm B = "<<std::endl,F,N,N,B,ldb);
            } else {
                    // B2 <- B2 - A2 . B1
                ftrmm (F, oppSide, UpLo, FFLAS::FflasNoTrans, diagB, N2, N1, F.mOne, B, ldb, A2, lda, F.one, B2, ldb);
                    // B2 <- A3^-1 . B2
                ftrsm (F, side, UpLo, FFLAS::FflasNoTrans, diagA, N2, N1, F.one, A3, lda, B2, ldb);
            }
        }
    }

} // FFPACK

#endif // __FFLASFFPACK_ffpack_ftrstr_INL
