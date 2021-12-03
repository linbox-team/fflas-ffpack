/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Clément Pernet
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


//--------------------------------------------------------------------------
//                        Test for det
//
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iomanip>
#include <iostream>
#include <givaro/modular-balanced.h>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"

#include "fflas-ffpack/utils/test-utils.h"

using namespace FFLAS;

template<class Field, class RandIter>
bool test_det(Field &F, size_t n, int iter, RandIter& G)
{
    typedef typename Field::Element Element;
    //! @todo test with stride
    Element * A = fflas_new (F, n, n);
    FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::Threads> parH;
    FFLAS::ParSeqHelper::Sequential seqH;

    bool pass = true;
    Element d,dt;
    F.init(d); F.init(dt);
    PAR_BLOCK{
        for(int i = 0;i<iter;++i){
            G.random(dt);
            FFPACK::RandomMatrixWithDet(F, n, dt, A, n, G);
            FFPACK::Det (F, d, n, A, n);
            if (!F.areEqual(dt,d)) {
                pass = false;
                break;
            }

            FFPACK::RandomMatrixWithDet(F, n, dt, A, n, G);
            FFPACK::Det(F,d,n,A,n,seqH);
            if (!F.areEqual(dt,d)) {
                pass = false;
                break;
            }
            FFPACK::RandomMatrixWithDet(F, n, dt, A, n, G);
            parH.set_numthreads(NUM_THREADS);
            FFPACK::Det(F,d,n,A,n,parH);
            if (!F.areEqual(dt,d)) {
                pass = false;
                break;
            }

            ++dt;
        }
    }

    //test the wrapped Det without the need of using a PAR_BLOCK macro to label the desired parallel region
    for(int i = 0;i<iter;++i){
        FFPACK::RandomMatrixWithDet(F, n, dt, A, n, G);
        FFPACK::pDet(F,d,n,A,n);
        if (!F.areEqual(dt,d)) {
            pass = false;
            break;
        }
    }

    fflas_delete( A);
    return pass;
}

int main(int argc, char** argv)
{

    int iters =10 ;
    Givaro::Integer p = 65521 ;
    size_t n = 200 ;
    uint64_t seed = getSeed();
    Argument as[] = {
        { 'p', "-p P", "Set the field characteristic.",         TYPE_INTEGER , &p },
        { 'n', "-n N", "Set the dimension of the matrix.",      TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);

    bool pass = true ;
    typedef Givaro::ModularBalanced<double> Field;
    Field F(p);
    Field::RandIter G(F,seed);
    Givaro::ZRing<Givaro::Integer> ZZ;
    Givaro::ZRing<Givaro::Integer>::RandIter GZZ(ZZ,seed);

    pass = pass && test_det(F,n,iters,G);
        // pass = pass && test_det(ZZ,n,iters,GZZ); @fixme: need a specific random matrix generator over ZZ

    return ((pass==true)?0:1);
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
