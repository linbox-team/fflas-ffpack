/* ffpack/ffpack_charpoly_danilevski.inl
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
#ifndef __FFLASFFPACK_ffpack_charpoly_danilveski_INL
#define __FFLASFFPACK_ffpack_charpoly_danilveski_INL

namespace FFPACK {

    //---------------------------------------------------------------------
    // CharPoly: Compute the characteristic polynomial of A using
    // Danilevski's algorithm.
    //---------------------------------------------------------------------

    template <class Field, class Polynomial>
    std::list<Polynomial>&
    Danilevski (const Field& F, std::list<Polynomial>& charp,
                const size_t N, typename Field::Element_ptr A,
                const size_t lda)
    {
        charp.clear();
        size_t dtot=0;
        typename Field::Element_ptr pivot,e,u1;
        typename Field::Element invp;
        for (size_t k=0; k<N; ++k){
            size_t i = k+1;
            e = pivot = A + (k+1) * lda + k; // coef
            while ((i<N) && F.isZero(*e)) { e += lda; i++; }
            if (i < N){
                if (i > k + 1) {
                    FFLAS::fswap (F, N-k, e, 1, pivot, 1);
                    FFLAS::fswap (F, N, A+i, lda, A+k+1, lda);
                }
                F.inv (invp, *pivot);
                FFLAS::fscalin (F, N-k-1, invp, pivot+1, 1);
                FFLAS::fscalin (F, N-dtot, *pivot, A+dtot*lda+k+1, lda);
                // X <- X - uw
                FFLAS::fger (F, k + 1-dtot, N - k -1, F.mOne,
                             A + dtot*lda + k, lda, pivot+1, 1,
                             A+k+1+dtot*lda, lda);
                if (k<N-2){

                    // Y <- Y - vw
                    FFLAS::fger (F, N-k-2, N-k-1,  F.mOne, pivot+lda, lda, pivot+1, 1,
                                 pivot+lda+1,lda);
                    //6
                    fgemv (F, FFLAS::FflasNoTrans, N-dtot, N-k-2,
                           F.one, A+dtot*lda+k+2, lda, pivot+lda, lda,  F.one,
                           A+dtot*lda+k+1,lda);
                }
                //5
                u1 = A+dtot*lda+k;
                for (i = dtot; i <= k; ++i){
                    F.addin( *(u1+lda+1), *u1);
                    u1+=lda;
                }
            }
            if (i==N){// completed one companion block
                size_t d;
                d = k+1-dtot;
                typename Field::Element_ptr Ai = A+k+dtot*lda;
                Polynomial  P(d+1);
                for (i = 0; i < d; ++i){
                    F.neg (P[i], *(Ai+i*lda));
                }
                F.assign( P[d],  F.one);
                charp.push_front(P);
                dtot+=d;
            }
        }
        return charp;
    }

} // FFPACK

#endif // __FFLASFFPACK_ffpack_charpoly_danilveski_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
