/* fflas/fflas_freduce.inl
 * Copyright (C) 2014 FFLAS FFPACK group
 *
 * Written by  Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

#ifndef __FFLASFFPACK_fflas_freduce_H
#define __FFLASFFPACK_fflas_freduce_H

#include "fflas-ffpack/fflas/fflas_simd.h"
#include "fflas-ffpack/field/field-traits.h"
#include "fflas-ffpack/utils/cast.h"

namespace FFLAS {

    template<class T>
    struct support_simd_mod  : public std::false_type {} ;

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
    template<>
    struct support_simd_mod<float> : public std::true_type {} ;
    template<>
    struct support_simd_mod<double> : public std::true_type {} ;
#ifdef __x86_64__
    template<>
    struct support_simd_mod<int64_t> : public std::true_type {} ;
#endif  // __x86_64__

#endif // __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

    /* Note that support_simd_mod => support_fast_mod */
    template<class T>
    struct support_fast_mod  : public std::false_type {} ;
    template<>
    struct support_fast_mod<float> : public std::true_type {} ;
    template<>
    struct support_fast_mod<double> : public std::true_type {} ;
    template<>
    struct support_fast_mod<int64_t> : public std::true_type {} ;

} // FFLAS

#include "fflas-ffpack/fflas/fflas_freduce.inl"

namespace FFLAS {

    /***************************/
    /*         LEVEL 1         */
    /***************************/

    template<class Field>
    void
    freduce (const Field & F, const size_t m,
             typename Field::ConstElement_ptr  B, const size_t incY,
             typename Field::Element_ptr A, const size_t incX)
    {
        return details::freduce (F,m,B,incY,A,incX,typename FieldTraits<Field>::category());
    }

    template<class Field>
    void
    freduce (const Field & F, const size_t m,
             typename Field::Element_ptr A, const size_t incX)
    {
        return details::freduce (F,m,A,incX,typename FieldTraits<Field>::category());
    }

    template<class Field>
    void
    freduce_constoverride(const Field & F, const size_t m,
                          typename Field::ConstElement_ptr A, const size_t incX)
    {
        return freduce(F, m, FFPACK::fflas_const_cast<typename Field::Element_ptr>(A), incX);
    }

    // OOOPS
    // CP: to be moved to a fflas_finit field, if ever needed
    template<class Field, class ConstOtherElement_ptr>
    void
    finit (const Field& F, const size_t n,
           ConstOtherElement_ptr Y, const size_t incY,
           typename Field::Element_ptr X, const size_t incX)
    {
        typename Field::Element_ptr Xi = X ;
        ConstOtherElement_ptr Yi = Y ;

        if (incX == 1 && incY == 1)
            for (; Yi < Y + n ; ++Xi, ++Yi) {
                F.init(*Xi, *Yi);
            }
        else
            for (; Yi < Y+n*incY; Xi+=incX, Yi += incY ) {
                F.init(*Xi, *Yi);
            }
    }


    template<class Field>
    void
    finit (const Field& F, const size_t n,
           typename Field::Element_ptr X, const size_t incX)
    {
        typename Field::Element_ptr Xi = X ;

        if (incX == 1)
            for (; Xi < X + n ; ++Xi) {
                F.init(*Xi);
            }
        else
            for (; Xi < X+n*incX; Xi+=incX ) {
                F.init(*Xi);
            }
    }

    /***************************/
    /*         LEVEL 2         */
    /***************************/


    template<class Field>
    void
    freduce (const Field& F, const size_t m , const size_t n,
             typename Field::Element_ptr A, const size_t lda)
    {
        if (n == lda)
            freduce (F, n*m, A, 1);
        else
            for (size_t i = 0 ; i < m ; ++i)
                freduce (F, n, A+i*lda, 1);
        return;
    }
    template<class Field>
    void
    freduce (const Field& F, const FFLAS_UPLO UpLo, const size_t N,
             typename Field::Element_ptr A, const size_t lda)
    {
        typename Field::Element_ptr Ai = A;
        if (UpLo == FflasUpper){
            for (size_t i = N ; i > 0 ; --i, Ai+=lda+1)
                freduce (F, i, Ai, 1);
        } else { // Lower
            for (size_t i = 1 ; i <= N ; ++i, Ai+=lda)
                freduce (F, i, Ai, 1);
        }
        return;
    }
    template<class Field>
    void
    pfreduce (const Field& F, const size_t m , const size_t n,
              typename Field::Element_ptr A, const size_t lda, const size_t numths)
    {
        SYNCH_GROUP(
                    FORBLOCK1D(iter, m, SPLITTER(numths),
                               size_t rowsize= iter.end()-iter.begin();
                               TASK(MODE(CONSTREFERENCE(F) READWRITE(A[iter.begin()*lda])),
                                    freduce (F, rowsize, n, A+iter.begin()*lda, lda);
                                   );
                              );
                   );
        return;
    }

    template<class Field>
    void
    freduce (const Field& F, const size_t m , const size_t n,
             typename Field::ConstElement_ptr B, const size_t ldb,
             typename Field::Element_ptr A, const size_t lda)
    {
        for (size_t i = 0 ; i < m ; ++i) {
            freduce(F,n,B+i*ldb,1,A+i*lda,1);
        }
    }


    template<class Field>
    void
    freduce_constoverride(const Field & F, const size_t m, const size_t n,
                          typename Field::ConstElement_ptr A, const size_t lda)
    {
        return freduce(F, m, n,
                       FFPACK::fflas_const_cast<typename Field::Element_ptr>(A), lda);
    }

    // CP: to be moved to a fflas_finit field, if ever needed
    template<class Field, class OtherElement_ptr>
    void
    finit (const Field& F, const size_t m , const size_t n,
           const OtherElement_ptr B, const size_t ldb,
           typename Field::Element_ptr A, const size_t lda)
    {
        if (n == lda && n == ldb)
            finit (F, n*m, B, 1, A, 1);
        else
            for (size_t i = 0 ; i < m ; ++i)
                finit (F, n, B + i*ldb, 1, A + i*lda, 1);
        return;
    }

    template<class Field>
    void
    finit (const Field& F, const size_t m , const size_t n,
           typename Field::Element_ptr A, const size_t lda)
    {
        if (n == lda)
            finit (F, n*m, A, 1);
        else
            for (size_t i = 0 ; i < m ; ++i)
                finit (F, n, A + i*lda, 1);
        return;
    }

} // end of namespace FFLAS

//#include "fflas_freduce_mp.inl" moved to fflas.h

#endif // __FFLASFFPACK_fflas_freduce_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
