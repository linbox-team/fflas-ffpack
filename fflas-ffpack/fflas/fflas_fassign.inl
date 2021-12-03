/* fflas/fflas_fassign.inl
 * Copyright (C) 2007 Clement Pernet
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

#ifndef __FFLASFFPACK_fassign_INL
#define __FFLASFFPACK_fassign_INL

#include <string.h>
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include <givaro/zring.h>

#include "fflas-ffpack/utils/debug.h"

namespace FFLAS {


    /***************************/
    /*         LEVEL 1         */
    /***************************/


    template<class Field>
    inline void
    fassign (const Field& F, const size_t N,
             typename Field::ConstElement_ptr Y, const size_t incY,
             typename Field::Element_ptr X, const size_t incX)
    {
        typename Field::Element_ptr Xi = X;
        typename Field::ConstElement_ptr Yi=Y;

        if (incX == 1 && incY == 1) {
            for (; Xi < X+N; ++Xi, ++Yi)
                F.assign(*Xi,*Yi);

        }
        else {
            for (; Xi < X+N*incX; Xi+=incX, Yi+=incY )
                F.assign(*Xi,*Yi);
        }
        return;
    }

    template<>
    inline void
    fassign (const Givaro::Modular<float>& F, const size_t N,
             const float * Y, const size_t incY,
             float * X, const size_t incX)
    {

#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_scopy((int)N,Y,(int)incY,X,(int)incX);

        return;
    }

    template<>
    inline void
    fassign (const Givaro::ModularBalanced<float>& F, const size_t N,
             const float * Y, const size_t incY,
             float * X, const size_t incX)
    {

#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_scopy((int)N,Y,(int)incY,X,(int)incX);

        return;
    }

    template<>
    inline void
    fassign (const Givaro::ZRing<float>& F, const size_t N,
             const float * Y, const size_t incY,
             float * X, const size_t incX)
    {

#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_scopy((int)N,Y,(int)incY,X,(int)incX);

        return;
    }

    template<>
    inline void
    fassign (const Givaro::Modular<double>& F, const size_t N,
             const double * Y, const size_t incY,
             double * X, const size_t incX)
    {

#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_dcopy((int)N,Y,(int)incY,X,(int)incX);

        return;
    }

    template<>
    inline void
    fassign (const Givaro::ModularBalanced<double>& F, const size_t N,
             const double * Y, const size_t incY,
             double * X, const size_t incX)
    {

#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_dcopy((int)N,Y,(int)incY,X,(int)incX);

        return;
    }

    template<>
    inline void
    fassign (const Givaro::ZRing<double>& F, const size_t N,
             const double * Y, const size_t incY ,
             double * X, const size_t incX)
    {

#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_dcopy((int)N,Y,(int)incY,X,(int)incX);

        return;
    }


    /***************************/
    /*         LEVEL 2         */
    /***************************/


    template<class Field>
    void fassign (const Field& F, const size_t m, const size_t n,
                  typename Field::ConstElement_ptr B, const size_t ldb ,
                  typename Field::Element_ptr A, const size_t lda)
    {
        if (!m || !n) return;
        FFLASFFPACK_check(n<=std::min(lda,ldb));
        // if possible, copy one big block
        if (lda == n && ldb == n) {
            fassign(F,m*n,B,1,A,1);
            return ;
        }
        // else, copy row after row
        for (size_t i = 0 ; i < m ; ++i) {
            fassign(F,n,B+i*ldb,1,A+i*lda,1);
        }
        return;

    }


} // FFLAS


#endif // __FFLASFFPACK_fassign_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
