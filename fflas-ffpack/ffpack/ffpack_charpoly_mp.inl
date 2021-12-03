/*
 * Copyright (C) 2015 the FFLAS-FFPACK group
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

#ifndef __FFPACK_charpoly_mp_INL
#define __FFPACK_charpoly_mp_INL
#include <givaro/zring.h>
#include "givaro/givinteger.h"
#include "givaro/givpoly1.h"

#include "fflas-ffpack/field/rns-integer.h"
#include "fflas-ffpack/fflas-ffpack.h"

namespace FFPACK {

    typename FFPACK::RNSInteger<FFPACK::rns_double>::Element_ptr
    inline CharPoly (const FFPACK::RNSInteger<FFPACK::rns_double>& F,
                     typename FFPACK::RNSInteger<FFPACK::rns_double>::Element_ptr charp,
                     const size_t N,
                     typename FFPACK::RNSInteger<FFPACK::rns_double>::Element_ptr A, const size_t lda,
                     Givaro::ZRing<Givaro::Integer>::RandIter& G, const FFPACK_CHARPOLY_TAG CharpTag, size_t degree){
        //std::cerr<<"Using "<<F.size()<<" moduli";
        Givaro::Timer t;
        for(size_t i=0;i<F.size();i++){
            t.clear();t.start();
            typedef FFPACK::rns_double::ModField Field;
            typedef Givaro::Poly1Dom<Field> PolRing;
            PolRing::Element cp(N+1);
            Field::RandIter Gp(F.rns()._field_rns[i]); //TODO set the seed from G's seed
            PolRing R(F.rns()._field_rns[i]);
            FFPACK::CharPoly (R, cp, N, A._ptr+i*A._stride, lda, Gp, CharpTag, degree);
            FFLAS::fassign(Givaro::ZRing<double>(), N+1,  &(cp[0]),1, charp._ptr+i*charp._stride, 1);
            t.stop();
            //std::cerr<<"Iteration "<<i<<" --> "<<t.realtime()<<std::endl;
        }

        return charp;
    }
    template <>
    inline Givaro::Poly1Dom<Givaro::ZRing<Givaro::Integer> >::Element&
    CharPoly(const Givaro::Poly1Dom<Givaro::ZRing<Givaro::Integer> >& R,
             Givaro::Poly1Dom<Givaro::ZRing<Givaro::Integer> >::Element& charp,
             const size_t N,  Givaro::Integer * A, const size_t lda,
             Givaro::ZRing<Givaro::Integer>::RandIter& G, const FFPACK_CHARPOLY_TAG CharpTag, size_t degree){

        const Givaro::ZRing<Givaro::Integer>& F = R.getdomain();
        size_t Abs = FFLAS::bitsize(F,N,N,A,lda);
        // See [Dumas Pernet Wang ISSAC'05] for the following bound on the bitsize
        // of the coefficients of the characteristic polynomial
        int64_t CPbs = (int64_t) ceil(N/2.0*(log(double(N))/log(2.0)+2*Abs+0.21163275));
        Givaro::Integer CPbound = Givaro::Integer(1) << CPbs;

        FFPACK::rns_double RNS(CPbound, 23);
        typedef FFPACK::RNSInteger<FFPACK::rns_double> RnsDomain;
        RnsDomain Zrns(RNS);
        typename RnsDomain::Element_ptr Arns, CPrns;
        Arns = FFLAS::fflas_new(Zrns,N,N);
        CPrns = FFLAS::fflas_new(Zrns,1,N+1);
        //std::cerr<<"finit...";
        FFLAS::finit_rns(Zrns,N,N,(Abs/16)+((Abs%16)?1:0),A,lda,Arns);
        //std::cerr<<"...done"<<std::endl;

        //std::cerr<<"charpoly...";
        CharPoly(Zrns, CPrns, N, Arns, N, G, CharpTag, degree);
        //std::cerr<<"...done"<<std::endl;

        //std::cerr<<"fconvert...";
        charp.resize(N+1);
        FFLAS::fconvert_rns (Zrns,1,N+1, Givaro::Integer(1),&(charp[0]), N+1, CPrns);
        //std::cerr<<"...done"<<std::endl;

        FFLAS::fflas_delete(Arns);
        FFLAS::fflas_delete(CPrns);
        return charp;
    }


}

#endif // __FFPACK_charpoly_mp_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
