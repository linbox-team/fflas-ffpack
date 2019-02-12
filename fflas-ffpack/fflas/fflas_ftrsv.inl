/* fflas/fflas_ftrsv.inl
 * Copyright (C) 2005 Clement Pernet
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
#ifndef __FFLASFFPACK_ftrsv_INL
#define __FFLASFFPACK_ftrsv_INL

namespace FFLAS {
    //---------------------------------------------------------------------
    // ftrsv: TRiangular System solve with vector
    // Computes  X <- op(A^-1).X
    // size of X is m
    //---------------------------------------------------------------------
    template<class Field>
    inline void
    ftrsv (const Field& F, const FFLAS_UPLO Uplo,
           const FFLAS_TRANSPOSE TransA, const FFLAS_DIAG Diag,
           const size_t N,typename Field::ConstElement_ptr A, size_t lda,
           typename Field::Element_ptr X, int incX)
    {

        typename Field::Element_ptr Xi;
        typename Field::ConstElement_ptr Ai;
        if ( Uplo == FflasLower ){
            if ( TransA == FflasTrans){
                Ai = A+(N-1)*(lda+1); // bottom right entry of A
                Xi = X+(int)(N-1)*incX;
                for(size_t i=0; i<N; Ai-=(lda+1), Xi-=incX, ++i) {
                    if (i)
                        F.subin(*Xi, fdot(F, i, Ai+lda, lda, Xi+incX, incX));
                    if ( Diag==FflasNonUnit )
                        F.divin(*Xi,*Ai);
                }
            } // End FflasTrans
            else{
                Ai = A;
                Xi = X;
                for(size_t i=0 ; i<N; Ai+=lda,Xi+=incX, ++i ){
                    F.subin (*Xi, fdot (F, i, Ai, 1, X, incX));
                    if ( Diag==FflasNonUnit )
                        F.divin(*Xi,*(Ai+i));
                }
            }
        } // End EFflasLower
        else{
            if ( TransA == FflasTrans){
                Ai = A;
                Xi = X;
                for(size_t i=0; i<N; ++Ai,Xi+=incX, ++i) {
                    F.subin(*Xi, fdot(F, i, Ai, lda, X, incX));
                    if ( Diag == FflasNonUnit )
                        F.divin(*Xi, *(Ai+i*lda));
                }

            } // End FflasTrans
            else{
                Ai = A+(lda+1)*(N-1);
                Xi = X+incX*(int)(N-1);
                for(size_t i=0; Xi>=X; Ai-=lda+1,Xi-=incX, ++i ){
                    if (i)
                        F.subin (*Xi, fdot (F, i, Ai+1, 1, Xi+incX, incX));
                    if ( Diag==FflasNonUnit )
                        F.divin(*Xi,*Ai);
                }
            }
        }
    }

}
#endif // __FFLASFFPACK_ftrsv_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
