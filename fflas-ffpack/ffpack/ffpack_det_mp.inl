/*
 * Copyright (C) 2016 the FFLAS-FFPACK group
 *
 * Written by Clément Pernet <clement.pernet@imag.fr>
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

#ifndef __FFPACK_det_mp_INL
#define __FFPACK_det_mp_INL
#include <givaro/zring.h>
#include "givaro/givinteger.h"

#include "fflas-ffpack/field/rns-integer.h"
#include "fflas-ffpack/fflas-ffpack.h"

namespace FFPACK {

    template<class PSHelper>
    inline typename FFPACK::RNSInteger<FFPACK::rns_double>::Element_ptr&
    Det (const FFPACK::RNSInteger<FFPACK::rns_double>& F,
         typename FFPACK::RNSInteger<FFPACK::rns_double>::Element_ptr& det,
         const size_t N,
         typename FFPACK::RNSInteger<FFPACK::rns_double>::Element_ptr A, const size_t lda,
         const PSHelper& psH){

        for(size_t i=0;i<F.size();i++){
            const FFPACK::rns_double::ModField & Fmod =  F.rns()._field_rns[i];
            FFPACK::rns_double::ModField::Element dmf;
            Fmod.assign (*(det._ptr+i*det._stride), FFPACK::Det (Fmod, dmf, N, A._ptr+i*A._stride, lda, psH));
        }
        return det;
    }


    template <class PSHelper>
    inline Givaro::Integer&
    Det (const Givaro::ZRing<Givaro::Integer>& F, Givaro::Integer& det,
         const size_t N,  Givaro::Integer * A, const size_t lda,
         const PSHelper& psH, size_t*P,size_t*Q){

        if (N==0)
            return  F.assign(det,F.one);

        size_t Abs = FFLAS::bitsize(F,N,N,A,lda);
        // Hadamard's bound on the bitsize of the determinant over Z
        int64_t Detbs = (int64_t) ceil (N * (log(double(N))/(log(2.0)*2.0) + Abs));
        Givaro::Integer Detbound = Givaro::Integer(1) << Detbs;
        FFPACK::rns_double RNS(Detbound, 23);
        typedef FFPACK::RNSInteger<FFPACK::rns_double> RnsDomain;
        RnsDomain Zrns(RNS);
        typename RnsDomain::Element_ptr Arns, Detrns;
        Arns = FFLAS::fflas_new(Zrns,N,N);
        Detrns = FFLAS::fflas_new(Zrns,1,1);

        FFLAS::finit_rns(Zrns,N,N,(Abs/16)+((Abs%16)?1:0),A,lda,Arns);
        Det(Zrns, Detrns, N, Arns, N, psH);
        FFLAS::fconvert_rns (Zrns,1,1, Givaro::Integer(1),&det, 1, Detrns);

        FFLAS::fflas_delete(Arns);
        FFLAS::fflas_delete(Detrns);
        return det;
    }

}

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
