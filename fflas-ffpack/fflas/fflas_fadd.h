/*
 * Copyright (C) 2014 FFLAS-FFPACK group
 *
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

#ifndef __FFLASFFPACK_fadd_H
#define __FFLASFFPACK_fadd_H

#include "fflas-ffpack/fflas/fflas_simd.h"

namespace FFLAS {

    template<class T>
    struct support_simd_add  : public std::false_type {} ;

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
    template<>
    struct support_simd_add<float> : public std::true_type {} ;
    template<>
    struct support_simd_add<double> : public std::true_type {} ;
#ifdef SIMD_INT // why? PK - 2019/12
    template<>
    struct support_simd_add<int32_t> : public std::true_type {} ;
    template<>
    struct support_simd_add<int64_t> : public std::true_type {} ;
#endif // SIMD_INT
#endif // __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

} // FFLAS

#include "fflas_fadd.inl"

namespace FFLAS {

    /***************************/
    /*         LEVEL 1         */
    /***************************/

    template <class Field>
    void
    fadd (const Field & F,  const size_t N,
          typename Field::ConstElement_ptr A, const size_t inca,
          typename Field::ConstElement_ptr B, const size_t incb,
          typename Field::Element_ptr C, const size_t incc)
    {
        details::fadd<Field, true>(F,N,A,inca,B,incb,C,incc
                                   , typename FieldTraits<Field>::category() );
    }



    template <class Field>
    void
    faddin (const Field& F,  const size_t N,
            typename Field::ConstElement_ptr B, const size_t incb,
            typename Field::Element_ptr C, const size_t incc)
    {
        fadd(F,N,B,incb,C,incc,C,incc);
        return;
    }

    template <class Field>
    void
    fsub(const Field & F,  const size_t N,
         typename Field::ConstElement_ptr A, const size_t inca,
         typename Field::ConstElement_ptr B, const size_t incb,
         typename Field::Element_ptr C, const size_t incc)
    {

        details::fadd<Field, false>(F,N,A,inca,B,incb,C,incc
                                    , typename FieldTraits<Field>::category() );
    }



    template <class Field>
    void
    fsubin (const Field& F,  const size_t N,
            typename Field::ConstElement_ptr B, const size_t incb,
            typename Field::Element_ptr C, const size_t incc)
    {
        fsub(F,N,C,incc,B,incb,C,incc);
        return;
    }

    // C = A + a B
    template <class Field>
    void
    fadd (const Field& F, const size_t N,
          typename Field::ConstElement_ptr A, const size_t inca,
          const typename Field::Element alpha,
          typename Field::ConstElement_ptr B, const size_t incb,
          typename Field::Element_ptr C, const size_t incc)
    {
        if (C == A && inca == incc)
            return faxpy(F,N,alpha,B,incb,C,incc);
        if (F.isOne(alpha))
            return fadd(F,N,A,inca,B,incb,C,incc);
        if (F.isMOne(alpha)){
            return fsub(F,N,A,inca,B,incb,C,incc);
        }
        if (F.isZero(alpha))
            return fassign(F,N,A,inca,C,incc);

        if (inca == 1 && incb == 1 && incc == 1) {
            for (size_t i = 0 ; i < N ; ++i) {
                //!@todo optimise here
                F.mul(C[i],alpha,B[i]);
                F.addin(C[i],A[i]);
            }
            return;
        }

        typename Field::ConstElement_ptr Ai = A, Bi = B;
        typename Field::Element_ptr Ci = C;
        for (; Ai < A+N*inca; Ai+=inca, Bi+=incb, Ci+=incc) {
            F.mul(*Ci,alpha,*Bi);
            F.addin (*Ci, *Ai);
        }
    }


    /***************************/
    /*         LEVEL 2         */
    /***************************/


    template <class Field>
    void
    pfadd (const Field & F,  const size_t M, const size_t N,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr B, const size_t ldb,
           typename Field::Element_ptr C, const size_t ldc, const size_t numths){
        SYNCH_GROUP(
                    FORBLOCK1D(iter, M, SPLITTER(numths),
                               size_t rowsize= iter.end()-iter.begin();
                               TASK(MODE(CONSTREFERENCE(F) READWRITE(C[iter.begin()*ldc]) READ(A[iter.begin()*lda], B[iter.begin()*ldb])),
                                    fadd(F, rowsize, N, A+iter.begin()*lda, lda, B+iter.begin()*ldb, ldb, C+iter.begin()*ldc, ldc);
                                   );
                              );
                   );
    }

    template <class Field>
    void
    pfsub (const Field & F,  const size_t M, const size_t N,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr B, const size_t ldb,
           typename Field::Element_ptr C, const size_t ldc, const size_t numths){
        SYNCH_GROUP(
                    FORBLOCK1D(iter, M, SPLITTER(numths),
                               size_t rowsize= iter.end()-iter.begin();
                               TASK(MODE(CONSTREFERENCE(F) READWRITE(C[iter.begin()*ldc]) READ(A[iter.begin()*lda], B[iter.begin()*ldb])),
                                    fsub(F, rowsize, N, A+iter.begin()*lda, lda, B+iter.begin()*ldb, ldb, C+iter.begin()*ldc, ldc);
                                   );
                              );
                   );
    }


    template <class Field>
    void
    pfaddin (const Field& F, const size_t M, const size_t N,
             typename Field::ConstElement_ptr B, const size_t ldb,
             typename Field::Element_ptr C, const size_t ldc, size_t numths){

        SYNCH_GROUP(
                    FORBLOCK1D(iter, M, SPLITTER(numths),
                               size_t rowsize= iter.end()-iter.begin();
                               TASK(MODE(CONSTREFERENCE(F) READWRITE(C[iter.begin()*ldc]) READ(B[iter.begin()*ldb])),
                                    faddin(F, rowsize, N, B+iter.begin()*ldb, ldb, C+iter.begin()*ldc, ldc);
                                   );
                              );
                   );
    }

    template <class Field>
    void
    pfsubin (const Field& F, const size_t M, const size_t N,
             typename Field::ConstElement_ptr B, const size_t ldb,
             typename Field::Element_ptr C, const size_t ldc, size_t numths){
        SYNCH_GROUP(
                    FORBLOCK1D(iter, M, SPLITTER(numths),
                               size_t rowsize= iter.end()-iter.begin();
                               TASK(MODE(CONSTREFERENCE(F) READWRITE(C[iter.begin()*ldc]) READ(B[iter.begin()*ldb])),
                                    fsubin(F, rowsize, N, B+iter.begin()*ldb, ldb, C+iter.begin()*ldc, ldc);
                                   );
                              );
                   );
    }

    template <class Field>
    void
    fadd (const Field& F, const size_t M, const size_t N,
          typename Field::ConstElement_ptr A, const size_t lda,
          typename Field::ConstElement_ptr B, const size_t ldb,
          typename Field::Element_ptr C, const size_t ldc)
    {
        if (N == lda && N == ldb && N == ldc)
            return fadd(F,M*N,A,1,B,1,C,1);
        typename Field::ConstElement_ptr Ai = A, Bi = B;
        typename Field::Element_ptr Ci = C;
        for (; Ai < A+M*lda; Ai+=lda, Bi+=ldb, Ci+=ldc)
            fadd(F,N,Ai,1,Bi,1,Ci,1);
    }

    template <class Field>
    void
    fsub (const Field& F, const size_t M, const size_t N,
          typename Field::ConstElement_ptr A, const size_t lda,
          typename Field::ConstElement_ptr B, const size_t ldb,
          typename Field::Element_ptr C, const size_t ldc)
    {
        if (N == lda && N == ldb && N == ldc)
            return fsub(F,M*N,A,1,B,1,C,1);
        typename Field::ConstElement_ptr Ai = A, Bi = B;
        typename Field::Element_ptr Ci = C;
        for (; Ai < A+M*lda; Ai+=lda, Bi+=ldb, Ci+=ldc)
            fsub(F,N,Ai,1,Bi,1,Ci,1);
    }

    template <class Field>
    void
    faddin (const Field& F, const size_t M, const size_t N,
            typename Field::ConstElement_ptr B, const size_t ldb,
            typename Field::Element_ptr C, const size_t ldc)
    {
        if (N == ldb && N == ldc)
            return faddin(F,M*N,B,1,C,1);
        const typename Field::Element  *Bi = B;
        typename Field::Element_ptr Ci = C;
        for (; Bi < B+M*ldb;  Bi+=ldb, Ci+=ldc)
            faddin(F,N,Bi,1,Ci,1);
    }

    template <class Field>
    void
    faddin (const Field& F,
            const FFLAS_UPLO uplo,
            const size_t N,
            typename Field::ConstElement_ptr B, const size_t ldb,
            typename Field::Element_ptr C, const size_t ldc)
    {
        const typename Field::Element  *Bi = B;
        typename Field::Element_ptr Ci = C;
        if (uplo == FflasUpper){
            for (size_t i=N; i>0; --i,  Bi+=ldb+1, Ci+=ldc+1)
                faddin(F,i,Bi,1,Ci,1);
        } else {
            for (size_t i=1; i <= N; ++i, Bi+=ldb, Ci+=ldc)
                faddin(F,i,Bi,1,Ci,1);
        }
    }


    template <class Field>
    void
    fsubin (const Field& F, const size_t M, const size_t N,
            typename Field::ConstElement_ptr B, const size_t ldb,
            typename Field::Element_ptr C, const size_t ldc)
    {
        if (N == ldb && N == ldc)
            return fsubin(F,M*N,B,1,C,1);
        typename Field::ConstElement_ptr Bi = B;
        typename Field::Element_ptr Ci = C;
        for (; Bi < B+M*ldb;  Bi+=ldb, Ci+=ldc)
            fsubin(F,N,Bi,1,Ci,1);
    }


    // C = A + a B
    template <class Field>
    void
    fadd (const Field& F, const size_t M, const size_t N,
          typename Field::ConstElement_ptr A, const size_t lda,
          const typename Field::Element alpha,
          typename Field::ConstElement_ptr B, const size_t ldb,
          typename Field::Element_ptr C, const size_t ldc)
    {
        if (C == A && lda == ldc)
            return faxpy(F,M,N,alpha,B,ldb,C,ldc);
        if (F.isOne(alpha))
            return fadd(F,M,N,A,lda,B,ldb,C,ldc);
        if (F.isMOne(alpha))
            return fsub(F,M,N,A,lda,B,ldb,C,ldc);
        if (F.isZero(alpha))
            return fassign(F,M,N,A,lda,C,ldc);

        if (N == lda && N == ldb && N == ldc)
            return fadd(F,M*N,A,1,alpha,B,1,C,1);

        typename Field::ConstElement_ptr Ai = A, Bi = B;
        typename Field::Element_ptr Ci = C;
        for (; Ai < A+M*lda; Ai+=lda, Bi+=ldb, Ci+=ldc)
            for (size_t i=0; i<N; i++) {
                F.mul(Ci[i],alpha,Bi[i]);
                F.addin (Ci[i], Ai[i]);
            }
    }


} // FFLAS


#endif // __FFLASFFPACK_fadd_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
