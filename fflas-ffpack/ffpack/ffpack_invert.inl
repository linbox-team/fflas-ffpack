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

#ifndef __FFLASFFPACK_ffpack_invert_INL
#define __FFLASFFPACK_ffpack_invert_INL


namespace FFPACK {
    template <class Field>
    typename Field::Element_ptr
    Invert (const Field& F, const size_t M,
            typename Field::Element_ptr A, const size_t lda,
            int& nullity)
    {
        FFLASFFPACK_check(lda >= M);

        Checker_invert<Field> checker(F,M,A,lda);

        if (M == 0) {
            nullity = 0 ;
            return NULL ;
        }
        size_t * P = FFLAS::fflas_new<size_t>(M);
        size_t * Q = FFLAS::fflas_new<size_t>(M);
        size_t R =  ReducedRowEchelonForm (F, M, M, A, lda, P, Q, true, FfpackGaussJordanTile);
        nullity = (int)(M - R);
        applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, M, 0, (int)R, A, lda, P);
        FFLAS::fflas_delete(P);
        FFLAS::fflas_delete(Q);

        checker.check(A,nullity);
        return A;
    }

    template <class Field>
    typename Field::Element_ptr
    Invert (const Field& F, const size_t M,
            typename Field::ConstElement_ptr A, const size_t lda,
            typename Field::Element_ptr X, const size_t ldx,
            int& nullity)
    {
        FFLASFFPACK_check(lda >= M);
        FFLASFFPACK_check(ldx >= M);
        if (M == 0) {
            nullity = 0 ;
            return NULL ;
        }
        FFLAS::fassign(F,M,M,A,lda,X,ldx);
        Invert (F, M, X, ldx, nullity);
        return X;
    }

    template <class Field>
    typename Field::Element_ptr
    Invert2( const Field& F, const size_t M,
             typename Field::Element_ptr A, const size_t lda,
             typename Field::Element_ptr X, const size_t ldx,
             int& nullity)
    {
        FFLASFFPACK_check(lda >= M);
        FFLASFFPACK_check(ldx >= M);

        if (M == 0) {
            nullity = 0 ;
            return NULL ;
        }

        size_t *P = FFLAS::fflas_new<size_t>(M);
        size_t *rowP = FFLAS::fflas_new<size_t>(M);


        nullity = int(M - LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, M, A, lda, P, rowP));

        if (nullity > 0){
            FFLAS::fflas_delete( P);
            FFLAS::fflas_delete( rowP);
            return NULL;
        }
        else {
            // Initializing X to 0
#if 0/*  timer remnants */
            t1.clear();
            t1.start();
#endif
            //! @todo this init is not all necessary (done after ftrtri)
            FFLAS::fzero(F,M,M,X,ldx);

            // X = L^-1 in n^3/3
            ftrtri (F, FFLAS::FflasLower, FFLAS::FflasUnit, M, A, lda);
            for (size_t i=0; i<M; ++i){
                for (size_t j=i; j<M; ++j)
                    F.assign(*(X +i*ldx+j), F.zero);
                F.assign (*(X+i*(ldx+1)), F.one);
            }
            for (size_t i=1; i<M; ++i)
                FFLAS::fassign (F, i, (A+i*lda), 1, (X+i*ldx), 1);

            ftrsm( F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
                   M, M, F.one, A, lda , X, ldx);

            // X = P^-1.X
            applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans,
                    M, 0,(int) M, X, ldx, P );

            FFLAS::fflas_delete( P);
            FFLAS::fflas_delete( rowP);
            return X;
        }
    }

} // FFPACK

#endif // __FFLASFFPACK_ffpack_invert_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
