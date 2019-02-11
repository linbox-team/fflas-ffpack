/*
 * Copyright (C) 2014 FFLAS-FFPACK group
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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

#ifndef __FFLASFFPACK_fgemv_mp_INL
#define __FFLASFFPACK_fgemv_mp_INL

#include "fflas-ffpack/field/rns-integer-mod.h"

namespace FFLAS {


    // specialization of the fgemv function for the field RNSInteger<rns_double>
    inline FFPACK::rns_double::Element_ptr
    fgemv (const FFPACK::RNSInteger<FFPACK::rns_double>& F, const FFLAS_TRANSPOSE ta,
           const size_t M, const size_t N,
           const FFPACK::rns_double::Element alpha,
           FFPACK::rns_double::ConstElement_ptr A, const size_t lda,
           FFPACK::rns_double::ConstElement_ptr X, const size_t incX,
           const FFPACK::rns_double::Element beta,
           FFPACK::rns_double::Element_ptr Y, const size_t incY,
           MMHelper<FFPACK::RNSInteger<FFPACK::rns_double>, MMHelperAlgo::Classic, ModeCategories::DefaultTag> & H)
    {
        if (M!=0 && N !=0){
            for (size_t i=0;i<F.size();i++)
                fgemv(F.rns()._field_rns[i], ta,
                      M, N,
                      alpha._ptr[i*alpha._stride],
                      A._ptr+i*A._stride, lda,
                      X._ptr+i*X._stride, incX,
                      beta._ptr[i*beta._stride],
                      Y._ptr+i*Y._stride, incY
                     );
        }
        return Y;
    }


    // specialization of the fgemv function for the field RNSIntegerMod<rns_double>
    inline FFPACK::rns_double::Element_ptr
    fgemv (const FFPACK::RNSIntegerMod<FFPACK::rns_double>& F, const FFLAS_TRANSPOSE ta,
           const size_t M, const size_t N,
           const FFPACK::rns_double::Element alpha,
           FFPACK::rns_double::ConstElement_ptr A, const size_t lda,
           FFPACK::rns_double::ConstElement_ptr X, const size_t incX,
           const FFPACK::rns_double::Element beta,
           FFPACK::rns_double::Element_ptr Y, const size_t incY,
           MMHelper<FFPACK::RNSIntegerMod<FFPACK::rns_double>, MMHelperAlgo::Classic, ModeCategories::DefaultTag> & H)
    {
        //std::cout<<"HERE 1"<<std::endl;
        MMHelper<FFPACK::RNSInteger<FFPACK::rns_double>, MMHelperAlgo::Classic, ModeCategories::DefaultTag >  H2;
        //std::cout<<"HERE 2"<<std::endl;
        fgemv(F.delayed(),ta,M,N,alpha,A,lda,X,incX, beta,Y,incY,H2);
        //std::cout<<"HERE 3"<<std::endl;
        size_t Ydim = (ta == FflasNoTrans)?M:N;
        freduce (F, Ydim, Y, incY);
        return Y;
    }


    // BB hack. might not work.
    // Calling fgemm, TODO: really specialize fgemv
    // specialization of the fgemv function for the field Givaro::ZRing<Givaro::Integer>
    inline Givaro::Integer* fgemv (const Givaro::ZRing<Givaro::Integer>& F,
                                   const FFLAS_TRANSPOSE ta,
                                   const size_t m, const size_t n,
                                   const Givaro::Integer alpha,
                                   Givaro::Integer* A, const size_t lda,
                                   Givaro::Integer* X, const size_t ldx,
                                   Givaro::Integer beta,
                                   Givaro::Integer* Y, const size_t ldy,
                                   MMHelper<Givaro::ZRing<Givaro::Integer>, MMHelperAlgo::Classic, ModeCategories::ConvertTo<ElementCategories::RNSElementTag> > & H)
    {
        MMHelper<Givaro::ZRing<Givaro::Integer>, MMHelperAlgo::Classic, ModeCategories::ConvertTo<ElementCategories::RNSElementTag>, ParSeqHelper::Sequential> H2;
        fgemm(F,ta,FFLAS::FflasNoTrans, (ta==FFLAS::FflasNoTrans)?m:n, 1,(ta==FFLAS::FflasNoTrans)?n:m, alpha,A,lda,X,ldx,beta,Y,ldy,H2);
        return Y;
    }

    // specialization of the fgemv function for the field Givaro::Modular<Givaro::Integer>
    // Calling fgemm, TODO: really specialize fgemv
    inline Givaro::Integer* fgemv (const Givaro::Modular<Givaro::Integer>& F,
                                   const FFLAS_TRANSPOSE ta,
                                   const size_t m, const size_t n,
                                   const Givaro::Integer alpha,
                                   Givaro::Integer* A, const size_t lda,
                                   Givaro::Integer* X, const size_t ldx,
                                   Givaro::Integer beta,
                                   Givaro::Integer* Y, const size_t ldy,
                                   MMHelper<Givaro::Modular<Givaro::Integer>, MMHelperAlgo::Classic, ModeCategories::ConvertTo<ElementCategories::RNSElementTag> > & H)
    {
        MMHelper<Givaro::Modular<Givaro::Integer>, MMHelperAlgo::Classic, ModeCategories::ConvertTo<ElementCategories::RNSElementTag>, ParSeqHelper::Sequential> H2;
        fgemm(F,ta,FFLAS::FflasNoTrans,(ta==FFLAS::FflasNoTrans)?m:n,1,(ta==FFLAS::FflasNoTrans)?n:m,alpha,A,lda,X,ldx,beta,Y,ldy,H2);
        return Y;
    }

    // specialization of the fgemv function for the field Givaro::Modular<RecInt::ruint<K>>
    // Calling fgemm, TODO: really specialize fgemv
    template <size_t K1, size_t K2, class ParSeq>
    inline RecInt::ruint<K1>*
    fgemv (const Givaro::Modular<RecInt::ruint<K1>,RecInt::ruint<K2> >& F,
           const FFLAS_TRANSPOSE ta,
           const size_t m, const size_t n,
           const RecInt::ruint<K1> alpha,
           const RecInt::ruint<K1>* A, const size_t lda,
           const RecInt::ruint<K1>* X, const size_t incx,
           RecInt::ruint<K1> beta,
           RecInt::ruint<K1>* Y, const size_t incy,
           MMHelper<Givaro::Modular<RecInt::ruint<K1>,RecInt::ruint<K2> >,
           MMHelperAlgo::Classic,
           ModeCategories::ConvertTo<ElementCategories::RNSElementTag>,
           ParSeq >  & H) {
        MMHelper<Givaro::Modular<RecInt::ruint<K1>,RecInt::ruint<K2> >, MMHelperAlgo::Classic, ModeCategories::ConvertTo<ElementCategories::RNSElementTag>, ParSeqHelper::Sequential> H2;
        fgemm (F,ta,FflasNoTrans,(ta==FFLAS::FflasNoTrans)?m:n,1,(ta==FFLAS::FflasNoTrans)?n:m,alpha,A,lda,X,incx,beta,Y,incy,H2);
        return Y;
    }


} // end namespace FFLAS

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
