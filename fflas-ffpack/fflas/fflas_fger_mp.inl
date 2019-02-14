/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
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
/** @file fflas_fgemm/fgemm_classical_mp.inl
 * @brief matrix multiplication with multiprecision input (either over Z or over Z/pZ)
 */


#ifndef __FFPACK_fger_mp_INL
#define __FFPACK_fger_mp_INL

#include <givaro/modular-integer.h>
#include <givaro/zring.h>

#include "fflas-ffpack/fflas/fflas_helpers.inl"
#include "fflas-ffpack/fflas/fflas_fgemm/fgemm_classical_mp.inl"
#include "fflas-ffpack/field/rns-integer.h"
#include "fflas-ffpack/field/rns-integer-mod.h"

namespace FFLAS{


    inline void
    fger (const Givaro::Modular<Givaro::Integer>& F, const size_t M, const size_t N,
          const typename Givaro::Integer alpha,
          typename Givaro::Integer* x, const size_t incx,
          typename Givaro::Integer* y, const size_t incy,
          typename Givaro::Integer* A, const size_t lda,
          MMHelper<Givaro::Modular<Givaro::Integer>, MMHelperAlgo::Classic, ModeCategories::ConvertTo<ElementCategories::RNSElementTag> > & H)
    {
        MMHelper<Givaro::Modular<Givaro::Integer>, MMHelperAlgo::Classic, ModeCategories::DefaultTag>  H2;
        FFLAS::fger(F,M,N,alpha,x,incx,y,incy,A,lda,H2);
    }

    template<typename RNS>
    inline void
    fger (const FFPACK::RNSInteger<RNS>& F, const size_t M, const size_t N,
          const typename FFPACK::RNSInteger<RNS>::Element alpha,
          typename FFPACK::RNSInteger<RNS>::Element_ptr x, const size_t incx,
          typename FFPACK::RNSInteger<RNS>::Element_ptr y, const size_t incy,
          typename FFPACK::RNSInteger<RNS>::Element_ptr A, const size_t lda,
          MMHelper<FFPACK::RNSInteger<RNS>, MMHelperAlgo::Classic, ModeCategories::DefaultTag> & H)
    {
        for(size_t i=0;i<F.size();i++){
            FFLAS::fger(F.rns()._field_rns[i],M,N,
                        alpha._ptr[i*alpha._stride],
                        x._ptr+i*x._stride,incx,y._ptr+i*y._stride,incy,A._ptr+i*A._stride,lda);
        }
    }

    template<typename RNS>
    inline void
    fger (const FFPACK::RNSIntegerMod<RNS>& F, const size_t M, const size_t N,
          const typename FFPACK::RNSIntegerMod<RNS>::Element alpha,
          typename FFPACK::RNSIntegerMod<RNS>::Element_ptr x, const size_t incx,
          typename FFPACK::RNSIntegerMod<RNS>::Element_ptr y, const size_t incy,
          typename FFPACK::RNSIntegerMod<RNS>::Element_ptr A, const size_t lda,
          MMHelper<FFPACK::RNSIntegerMod<RNS>, MMHelperAlgo::Classic> & H)
    {
        typedef FFPACK::RNSInteger<RNS> RnsDomain;
        MMHelper<RnsDomain, MMHelperAlgo::Classic>  H2;
        RnsDomain Zrns(F.rns());
        FFLAS::fger(Zrns,M,N,alpha,x,incx,y,incy,A,lda,H2);

        // reduce the result mod p
        freduce (F, M, N, A, lda);
    }


} // namespace FFLAS

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
