/* fflas/fflas_ftrmm.inl
 * Copyright (C) 2007 Clement Pernet
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

#ifndef __FFLASFFPACK_ftrmm_INL
#define __FFLASFFPACK_ftrmm_INL

namespace FFLAS { namespace Protected {

    template <class NewField, class Field>
    inline void
    ftrmm_convert (const Field& F,
                   const FFLAS_SIDE Side,
                   const FFLAS_UPLO Uplo,
                   const FFLAS_TRANSPOSE TransA,
                   const FFLAS_DIAG Diag,
                   const size_t M, const size_t N,
                   const typename Field::Element alpha,
                   typename Field::ConstElement_ptr A, const size_t lda,
                   typename Field::Element_ptr B, const size_t ldb)
    {
        typedef typename NewField::Element FloatElement;
        NewField G((FloatElement) F.characteristic());
        FloatElement tmp, alphaf;
        F.convert (tmp, alpha);
        G.init(alphaf, tmp);

        size_t K = (Side == FflasLeft) ? M : N;
        FloatElement* Af = FFLAS::fflas_new(G, K, K);
        FloatElement* Bf = FFLAS::fflas_new(G, M, N);

        fconvert(F, K, K, Af, K, A, lda);
        freduce(G, K, K, Af, K);
        fconvert(F, M, N, Bf, N, B, ldb);
        freduce(G, M, N, Bf, N);

        ftrmm(G, Side, Uplo, TransA, Diag, M, N, alphaf, Af, K, Bf, N);

        finit(F, M, N, Bf, N, B, ldb);

        fflas_delete(Af);
        fflas_delete(Bf);
    }

    template <class NewField, class Field>
    inline void
    ftrmm_convert (const Field& F,
                   const FFLAS_SIDE Side,
                   const FFLAS_UPLO Uplo,
                   const FFLAS_TRANSPOSE TransA,
                   const FFLAS_DIAG Diag,
                   const size_t M, const size_t N,
                   const typename Field::Element alpha,
                   typename Field::ConstElement_ptr A, const size_t lda,
                   typename Field::ConstElement_ptr B, const size_t ldb,
                   const typename Field::Element beta,
                   typename Field::Element_ptr C, const size_t ldc)
    {
        typedef typename NewField::Element FloatElement;
        NewField G((FloatElement) F.characteristic());
        FloatElement tmp, alphaf, betaf;
        F.convert (tmp, alpha);
        G.init(alphaf, tmp);
        F.convert (tmp, beta);
        G.init(betaf, tmp);

        size_t K = (Side == FflasLeft) ? M : N;
        FloatElement* Af = FFLAS::fflas_new(G, K, K);
        FloatElement* Bf = FFLAS::fflas_new(G, M, N);
        FloatElement* Cf = FFLAS::fflas_new(G, M, N);

        fconvert(F, K, K, Af, K, A, lda);
        freduce(G, K, K, Af, K);
        fconvert(F, M, N, Bf, N, B, ldb);
        freduce(G, M, N, Bf, N);
        if (!F.isZero(beta)){
            fconvert(F, M, N, Cf, N, C, ldc);
            freduce(G, M, N, Cf, N);
        }

        ftrmm(G, Side, Uplo, TransA, Diag, M, N, alphaf, Af, K, Bf, N, betaf, Cf, N);

        finit(F, M, N, Cf, N, C, ldc);

        fflas_delete(Af);
        fflas_delete(Bf);
        fflas_delete(Cf);
    }

    // SFINAE helpers for ConvertTo dispatch (in-place)
    template<class Field>
    inline typename std::enable_if<
        std::is_same<typename ModeTraits<Field>::value,
                     ModeCategories::ConvertTo<ElementCategories::MachineFloatTag>>::value, bool
    >::type
    ftrmm_try_convert (const Field& F, const FFLAS_SIDE Side,
                       const FFLAS_UPLO Uplo, const FFLAS_TRANSPOSE TransA,
                       const FFLAS_DIAG Diag, const size_t M, const size_t N,
                       const typename Field::Element alpha,
                       typename Field::ConstElement_ptr A, const size_t lda,
                       typename Field::Element_ptr B, const size_t ldb)
    {
        if (!std::is_same<Field,Givaro::Modular<float> >::value){
            if (F.cardinality() == 2)
                ftrmm_convert<Givaro::Modular<float>,Field>(F,Side,Uplo,TransA,Diag,M,N,alpha,A,lda,B,ldb);
            else if (!std::is_same<Field,Givaro::ModularBalanced<float> >::value){
                if (F.cardinality() < DOUBLE_TO_FLOAT_CROSSOVER)
                    ftrmm_convert<Givaro::ModularBalanced<float>,Field>(F,Side,Uplo,TransA,Diag,M,N,alpha,A,lda,B,ldb);
                else if (!std::is_same<Field,Givaro::ModularBalanced<double> >::value && 16*F.cardinality() < Givaro::ModularBalanced<double>::maxCardinality())
                    ftrmm_convert<Givaro::ModularBalanced<double>,Field>(F,Side,Uplo,TransA,Diag,M,N,alpha,A,lda,B,ldb);
                else
                    FFPACK::failure()(__func__,__LINE__,"Invalid ConvertTo Mode for this field");
            }
        }
        return true;
    }

    template<class Field>
    inline typename std::enable_if<
        !std::is_same<typename ModeTraits<Field>::value,
                      ModeCategories::ConvertTo<ElementCategories::MachineFloatTag>>::value, bool
    >::type
    ftrmm_try_convert (const Field& F, const FFLAS_SIDE Side,
                       const FFLAS_UPLO Uplo, const FFLAS_TRANSPOSE TransA,
                       const FFLAS_DIAG Diag, const size_t M, const size_t N,
                       const typename Field::Element alpha,
                       typename Field::ConstElement_ptr A, const size_t lda,
                       typename Field::Element_ptr B, const size_t ldb)
    {
        return false;
    }

    // SFINAE helpers for ConvertTo dispatch (out-of-place)
    template<class Field>
    inline typename std::enable_if<
        std::is_same<typename ModeTraits<Field>::value,
                     ModeCategories::ConvertTo<ElementCategories::MachineFloatTag>>::value, bool
    >::type
    ftrmm_try_convert (const Field& F, const FFLAS_SIDE Side,
                       const FFLAS_UPLO Uplo, const FFLAS_TRANSPOSE TransA,
                       const FFLAS_DIAG Diag, const size_t M, const size_t N,
                       const typename Field::Element alpha,
                       typename Field::ConstElement_ptr A, const size_t lda,
                       typename Field::ConstElement_ptr B, const size_t ldb,
                       const typename Field::Element beta,
                       typename Field::Element_ptr C, const size_t ldc)
    {
        if (!std::is_same<Field,Givaro::Modular<float> >::value){
            if (F.cardinality() == 2)
                ftrmm_convert<Givaro::Modular<float>,Field>(F,Side,Uplo,TransA,Diag,M,N,alpha,A,lda,B,ldb,beta,C,ldc);
            else if (!std::is_same<Field,Givaro::ModularBalanced<float> >::value){
                if (F.cardinality() < DOUBLE_TO_FLOAT_CROSSOVER)
                    ftrmm_convert<Givaro::ModularBalanced<float>,Field>(F,Side,Uplo,TransA,Diag,M,N,alpha,A,lda,B,ldb,beta,C,ldc);
                else if (!std::is_same<Field,Givaro::ModularBalanced<double> >::value && 16*F.cardinality() < Givaro::ModularBalanced<double>::maxCardinality())
                    ftrmm_convert<Givaro::ModularBalanced<double>,Field>(F,Side,Uplo,TransA,Diag,M,N,alpha,A,lda,B,ldb,beta,C,ldc);
                else
                    FFPACK::failure()(__func__,__LINE__,"Invalid ConvertTo Mode for this field");
            }
        }
        return true;
    }

    template<class Field>
    inline typename std::enable_if<
        !std::is_same<typename ModeTraits<Field>::value,
                      ModeCategories::ConvertTo<ElementCategories::MachineFloatTag>>::value, bool
    >::type
    ftrmm_try_convert (const Field& F, const FFLAS_SIDE Side,
                       const FFLAS_UPLO Uplo, const FFLAS_TRANSPOSE TransA,
                       const FFLAS_DIAG Diag, const size_t M, const size_t N,
                       const typename Field::Element alpha,
                       typename Field::ConstElement_ptr A, const size_t lda,
                       typename Field::ConstElement_ptr B, const size_t ldb,
                       const typename Field::Element beta,
                       typename Field::Element_ptr C, const size_t ldc)
    {
        return false;
    }

}} // Protected, FFLAS

namespace FFLAS {

    //---------------------------------------------------------------------
    // ftrmm: TRiangular Matrix Multiply
    // Computes  B <- alpha.op(A).B,  B <- alpha.B.op(A)
    // B is M*N, A is M*M if Side==FflasLeft, N*N if Side==FflasRight
    // //---------------------------------------------------------------------
    template<class Field>
    inline void
    ftrmm (const Field& F, const FFLAS_SIDE Side,
           const FFLAS_UPLO Uplo,
           const FFLAS_TRANSPOSE TransA,
           const FFLAS_DIAG Diag,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::Element_ptr B, const size_t ldb)
    {
        if (!M || !N ) return;
        if (Protected::ftrmm_try_convert(F, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb))
            return;

        if ( Side==FflasLeft ){
            if ( Uplo==FflasUpper){
                if (TransA == FflasNoTrans){
                    if (Diag == FflasUnit)
                        Protected::ftrmmLeftUpperNoTransUnit<typename Field::Element> ()(F,M,N,A,lda,B,ldb);
                    else
                        Protected::ftrmmLeftUpperNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
                } else {
                    if (Diag == FflasUnit)
                        Protected::ftrmmLeftUpperTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
                    else
                        Protected::ftrmmLeftUpperTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
                }
            } else {
                if (TransA == FflasNoTrans){
                    if (Diag == FflasUnit)
                        Protected::ftrmmLeftLowerNoTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
                    else
                        Protected::ftrmmLeftLowerNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
                } else {
                    if (Diag == FflasUnit)
                        Protected::ftrmmLeftLowerTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
                    else
                        Protected::ftrmmLeftLowerTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
                }
            }
        } else {
            if ( Uplo == FflasUpper){
                if (TransA == FflasNoTrans){
                    if (Diag == FflasUnit)
                        Protected::ftrmmRightUpperNoTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
                    else
                        Protected::ftrmmRightUpperNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
                } else {
                    if (Diag == FflasUnit)
                        Protected::ftrmmRightUpperTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
                    else
                        Protected::ftrmmRightUpperTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
                }
            } else {
                if (TransA == FflasNoTrans){
                    if (Diag == FflasUnit)
                        Protected::ftrmmRightLowerNoTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
                    else
                        Protected::ftrmmRightLowerNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
                } else {
                    if (Diag == FflasUnit)
                        Protected::ftrmmRightLowerTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
                    else
                        Protected::ftrmmRightLowerTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
                }
            }
        }
        if (!F.isOne(alpha))
            fscalin(F,M,N,alpha,B,ldb);

    }
    // //---------------------------------------------------------------------
    // //---------------------------------------------------------------------
    template<class Field>
    inline void
    ftrmm (const Field& F, const FFLAS_SIDE Side,
           const FFLAS_UPLO Uplo,
           const FFLAS_TRANSPOSE TransA,
           const FFLAS_DIAG Diag,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr B, const size_t ldb,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc)
    {
        if (!M || !N ) return;
        if (Protected::ftrmm_try_convert(F, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb, beta, C, ldc))
            return;

        typename Field::Element bet;
        F.init(bet);
        F.assign(bet,beta);
        if (!F.isOne(alpha))
            F.divin(bet,alpha);
        if ( Side==FflasLeft ){
            if ( Uplo==FflasUpper){
                if (TransA == FflasNoTrans){
                    if (Diag == FflasUnit)
                        Protected::ftrmmLeftUpperNoTransUnit<typename Field::Element> ()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                    else
                        Protected::ftrmmLeftUpperNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                } else {
                    if (Diag == FflasUnit)
                        Protected::ftrmmLeftUpperTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                    else
                        Protected::ftrmmLeftUpperTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                }
            } else {
                if (TransA == FflasNoTrans){
                    if (Diag == FflasUnit)
                        Protected::ftrmmLeftLowerNoTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                    else
                        Protected::ftrmmLeftLowerNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                } else {
                    if (Diag == FflasUnit)
                        Protected::ftrmmLeftLowerTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                    else
                        Protected::ftrmmLeftLowerTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                }
            }
        } else {
            if ( Uplo == FflasUpper){
                if (TransA == FflasNoTrans){
                    if (Diag == FflasUnit)
                        Protected::ftrmmRightUpperNoTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                    else
                        Protected::ftrmmRightUpperNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                } else {
                    if (Diag == FflasUnit)
                        Protected::ftrmmRightUpperTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                    else
                        Protected::ftrmmRightUpperTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                }
            } else {
                if (TransA == FflasNoTrans){
                    if (Diag == FflasUnit)
                        Protected::ftrmmRightLowerNoTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                    else
                        Protected::ftrmmRightLowerNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                } else {
                    if (Diag == FflasUnit)
                        Protected::ftrmmRightLowerTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                    else
                        Protected::ftrmmRightLowerTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb,bet,C,ldc);
                }
            }
        }
        if (!F.isOne(alpha))
            fscalin(F,M,N,alpha,C,ldc);

    }
#ifndef DOXYGEN_SHOULD_SKIP_THIS

    namespace Protected {

#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
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
#include "fflas_ftrmm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT

    } // Protected

#endif // SKIPPED BY DOXYGEN

} // FFLAS

#endif // __FFLASFFPACK_ftrmm_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
