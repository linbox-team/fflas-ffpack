/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 * This file is Free Software and part of FFLAS-FFPACK.
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

/*! @file tests/test-utils.h
 * @ingroup tests
 * @brief Utilities to create matrices with prescribed shapes, properties,...
 * To be used in the tests
 */

#ifndef __FFLASFFPACK_tests_test_utils_H
#define __FFLASFFPACK_tests_test_utils_H

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/utils/debug.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include <givaro/givinteger.h>
#include <givaro/givintprime.h>
#include <givaro/givranditer.h>
#include <givaro/givtimer.h>

#include <random>
#include <functional>

namespace FFLAS {
    uint64_t getSeed(){
        return Givaro::BaseTimer::seed();
    }
}
namespace FFPACK {

    /*** Field chooser for test according to characteristic q and bitsize b ***/
    /* if q=-1 -> field is chosen randomly with a charateristic of b bits
       if b=0 -> bitsize is chosen randomly according to maxFieldElt
       */
    template<typename Field>
    Field* chooseField(Givaro::Integer q, uint64_t b, uint64_t seed){
        Givaro::Integer maxV = static_cast<Givaro::Integer>(FFLAS::maxCardinality<Field>());
        if (maxV>0 && (q> maxV || b> maxV.bitsize()))
            return nullptr;

        std::mt19937 mt_rand(seed);
        Givaro::Integer::seeding(mt_rand());
        Givaro::IntPrimeDom IPD;
        Givaro::Integer p;

        if (q == -1){ // No modulus specified
            do{
                Givaro::Integer _p;
                if (b <= 1){ // No bitsize specified
                    Givaro::Integer::random_lessthan<true> (_p, maxV);
                } else {
                    Givaro::Integer::random_exact_2exp (_p, b);
                }
                IPD.prevprime(p, _p+1);
            } while (p < Givaro::Integer(FFLAS::minCardinality<Field>()) || (maxV>0 && p > maxV));
        } else { // modulus specified
            p=q;
        }
        return new Field(p);
    }

    template<> Givaro::ZRing<int32_t>* chooseField<Givaro::ZRing<int32_t> >(Givaro::Integer q, uint64_t b, uint64_t seed){return new Givaro::ZRing<int32_t>();}
    template<> Givaro::ZRing<int64_t>* chooseField<Givaro::ZRing<int64_t> >(Givaro::Integer q, uint64_t b, uint64_t seed){return new Givaro::ZRing<int64_t>();}
    template<> Givaro::ZRing<float>* chooseField<Givaro::ZRing<float> >(Givaro::Integer q, uint64_t b, uint64_t seed){return new Givaro::ZRing<float>();}
    template<> Givaro::ZRing<double>* chooseField<Givaro::ZRing<double> >(Givaro::Integer q, uint64_t b, uint64_t seed){return new Givaro::ZRing<double>();}

} // FFPACK
#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
