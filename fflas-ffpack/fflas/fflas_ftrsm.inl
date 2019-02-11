/* fflas/fflas_ftrsm.inl
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

#ifndef __FFLASFFPACK_ftrsm_INL
#define __FFLASFFPACK_ftrsm_INL


namespace FFLAS {

    //---------------------------------------------------------------------
    // ftrsm: TRiangular System solve with matrix
    // Computes  B <- alpha.op(A^-1).B,  B <- alpha.B.op(A^-1)
    // B is M*N, A is M*M if Side==FflasLeft, N*N if Side==FflasRight
    //---------------------------------------------------------------------
    template<class Field>
    inline void
    ftrsm (const Field& F, const FFLAS_SIDE Side,
           const FFLAS_UPLO Uplo,
           const FFLAS_TRANSPOSE TransA,
           const FFLAS_DIAG Diag,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
#ifdef __FFLAS__TRSM_READONLY
           typename Field::ConstElement_ptr
#else
           typename Field::Element_ptr
#endif
           A, const size_t lda,
           typename Field::Element_ptr B, const size_t ldb)
    {
        ParSeqHelper::Sequential PSH;
        TRSMHelper<StructureHelper::Recursive, ParSeqHelper::Sequential> H(PSH);
        Checker_ftrsm<Field> checker(F, M, N, alpha, B, ldb);
        ftrsm(F, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb, H);
        checker.check(Side, Uplo, TransA, Diag, M, N, A, lda, B, ldb);
    }

    template<class Field>
    inline void
    ftrsm (const Field& F, const FFLAS_SIDE Side,
           const FFLAS_UPLO Uplo,
           const FFLAS_TRANSPOSE TransA,
           const FFLAS_DIAG Diag,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
#ifdef __FFLAS__TRSM_READONLY
           typename Field::ConstElement_ptr
#else
           typename Field::Element_ptr
#endif
           A, const size_t lda,
           typename Field::Element_ptr B, const size_t ldb,
           const ParSeqHelper::Sequential& PSH)
    {
        TRSMHelper<StructureHelper::Recursive, ParSeqHelper::Sequential> H(PSH);
        ftrsm(F, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb, H);
    }

    template<class Field, class Cut, class Param>
    inline void
    ftrsm (const Field& F, const FFLAS_SIDE Side,
           const FFLAS_UPLO Uplo,
           const FFLAS_TRANSPOSE TransA,
           const FFLAS_DIAG Diag,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
#ifdef __FFLAS__TRSM_READONLY
           typename Field::ConstElement_ptr
#else
           typename Field::Element_ptr
#endif
           A, const size_t lda,
           typename Field::Element_ptr B, const size_t ldb,
           const ParSeqHelper::Parallel<Cut,Param>& PSH)
    {
        TRSMHelper<StructureHelper::Iterative, ParSeqHelper::Parallel<Cut,Param> > H(PSH);
        ftrsm(F, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb, H);
    }

    template<class Field, class ParSeqTrait=ParSeqHelper::Sequential>
    inline void
    ftrsm (const Field& F, const FFLAS_SIDE Side,
           const FFLAS_UPLO Uplo,
           const FFLAS_TRANSPOSE TransA,
           const FFLAS_DIAG Diag,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
#ifdef __FFLAS__TRSM_READONLY
           typename Field::ConstElement_ptr
#else
           typename Field::Element_ptr
#endif
           A, const size_t lda,
           typename Field::Element_ptr B, const size_t ldb,
           TRSMHelper<StructureHelper::Recursive, ParSeqTrait> & H)
    {
        if (!M || !N ) return;

        if ( Side==FflasLeft ){
            if ( Uplo==FflasUpper){
                if (TransA == FflasNoTrans){
                    if (Diag == FflasUnit)
                        Protected::ftrsmLeftUpperNoTransUnit<typename Field::Element> ()(F,M,N,A,lda,B,ldb,H);
                    else
                        Protected::ftrsmLeftUpperNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,H);
                } else {
                    if (Diag == FflasUnit)
                        Protected::ftrsmLeftUpperTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,H);
                    else
                        Protected::ftrsmLeftUpperTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,H);
                }
            } else {
                if (TransA == FflasNoTrans){
                    if (Diag == FflasUnit)
                        Protected::ftrsmLeftLowerNoTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,H);
                    else
                        Protected::ftrsmLeftLowerNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,H);
                } else {
                    if (Diag == FflasUnit)
                        Protected::ftrsmLeftLowerTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,H);
                    else
                        Protected::ftrsmLeftLowerTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,H);
                }
            }
        } else {
            if ( Uplo == FflasUpper){
                if (TransA == FflasNoTrans){
                    if (Diag == FflasUnit)
                        Protected::ftrsmRightUpperNoTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,H);
                    else
                        Protected::ftrsmRightUpperNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,H);
                } else {
                    if (Diag == FflasUnit)
                        Protected::ftrsmRightUpperTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,H);
                    else
                        Protected::ftrsmRightUpperTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,H);
                }
            } else {
                if (TransA == FflasNoTrans){
                    if (Diag == FflasUnit)
                        Protected::ftrsmRightLowerNoTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,H);
                    else
                        Protected::ftrsmRightLowerNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,H);
                } else {
                    if (Diag == FflasUnit)
                        Protected::ftrsmRightLowerTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,H);
                    else
                        Protected::ftrsmRightLowerTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,H);
                }
            }
        }
        if (!F.isOne(alpha))
            fscalin(F,M,N,alpha,B,ldb);

    }


#ifndef DOXYGEN_SHOULD_SKIP_THIS

    namespace Protected {

#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT



#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT
        //==

#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT



#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT



#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT

    } // Protected

#endif // SKIPPED BY DOXYGEN

} // FFLAS

#endif // __FFLASFFPACK_ftrsm_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
