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

#ifndef __FFLASFFPACK_ffpack_fgesv_INL
#define __FFLASFFPACK_ffpack_fgesv_INL


namespace FFPACK {



    template <class Field>
    size_t
    fgesv (const Field& F,
           const FFLAS::FFLAS_SIDE Side,
           const size_t M, const size_t N,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::Element_ptr B, const size_t ldb,
           int * info)
    {

        size_t Na;
        if (Side == FFLAS::FflasLeft)
            Na = M;
        else
            Na = N;

        size_t* P = FFLAS::fflas_new<size_t>(Na);
        size_t* Q = FFLAS::fflas_new<size_t>(Na);

        size_t R = PLUQ (F, FFLAS::FflasNonUnit, Na, Na, A, lda, P, Q);

        fgetrs (F, Side, M, N, R, A, lda, P, Q, B, ldb, info);

        FFLAS::fflas_delete( P);
        FFLAS::fflas_delete( Q);

        return R;
    }

    template <class Field>
    size_t
    fgesv (const Field& F,
           const FFLAS::FFLAS_SIDE Side,
           const size_t M, const size_t N, const size_t NRHS,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::Element_ptr X, const size_t ldx,
           typename Field::ConstElement_ptr B, const size_t ldb,
           int * info)
    {

        size_t* P = FFLAS::fflas_new<size_t>(M);
        size_t* Q = FFLAS::fflas_new<size_t>(N);

        size_t R = PLUQ (F, FFLAS::FflasNonUnit, M, N, A, lda, P, Q);

        fgetrs (F, Side, M, N, NRHS, R, A, lda, P, Q, X, ldx, B, ldb, info);

        FFLAS::fflas_delete (P);
        FFLAS::fflas_delete (Q);

        return R;
    }

} //FFPACK

#endif // __FFLASFFPACK_ffpack_fgesv_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
