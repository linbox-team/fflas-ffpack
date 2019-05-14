/*
 * Copyright (C) the FFLAS-FFPACK group
 * Written by Clément Pernet <clement.pernet@imag.fr>
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

#define ENABLE_ALL_CHECKINGS 1

#define  __FFLASFFPACK_SEQUENTIAL

#include "fflas-ffpack/fflas-ffpack-config.h"

#include <iomanip>
#include <iostream>
#include <random>
#include <chrono>

#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>


using namespace std;
using namespace FFLAS;
using namespace FFPACK;
using Givaro::Modular;
using Givaro::ModularBalanced;

template <class Field>
bool run_with_field (Givaro::Integer q, size_t b, size_t n, size_t iters, uint64_t seed){
    bool ok = true ;
    int nbit=(int)iters;
    while (ok && nbit){
        Field* F= chooseField<Field>(q,b,seed);
        if (F==nullptr)
            return true;

        cout<<"Checking with ";F->write(cout)<<endl;

        size_t lda = n + (rand() % 4);
        size_t ldx = n + (rand() % 4);

        typename Field::Element_ptr A = fflas_new(*F, n, lda);
        typename Field::Element_ptr X = fflas_new(*F, n, ldx);

        typename Field::RandIter R(*F,seed++);

        RandomMatrixWithRankandRandomRPM (*F, A, lda, n, n, n, G);

        int nullity;
        FFPACK::Invert(*F, n, A, lda, X, ldx, nullity);

        if (nullity != 0){
            std::cerr<<"Error: Singular matrix detected"<<std::endl;
            fflas_delete(A);
            fflas_delete(X);
            return ok = false;
        }

        typename Field::Element_ptr Y = fflas_new(*F, n, n);
        fidentity(*F, n, n, Y, n);

        fgemm(*F, FflasNoTrans, FflasNoTrans, n,n,n, F->one, A, lda, X, ldx, F->mOne, Y, n);

        if (! fiszero(*F,n,n,Y,n)){
            //FFLAS::WriteMatrix (std::cerr<<"Y = "<<std::endl,*F,n,n,Y,n);
            std::cerr<<"Error: A * A^{-1} != Id"<<std::endl;
            fflas_delete(A);
            fflas_delete(X);
            fflas_delete(Y);
            return ok = false;
        }

        nbit--;
        fflas_delete(A);
        fflas_delete(X);
        fflas_delete(Y);
    }
    return ok;
}

int main(int argc, char** argv)
{
    cerr<<setprecision(10);
    Givaro::Integer q=-1;
    size_t b=0;
    size_t n=300;
    size_t iters=3;
    bool loop=false;
    uint64_t seed = getSeed();
    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
        { 'n', "-n N", "Set the dimension of the square matrix.", TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    srand(seed);

    bool ok = true;
    do{
        ok = ok && run_with_field<Modular<double> >(q,b,n,iters,seed);
        ok = ok && run_with_field<ModularBalanced<double> >(q,b,n,iters,seed);
        ok = ok && run_with_field<Modular<float> >(q,b,n,iters,seed);
        ok = ok && run_with_field<ModularBalanced<float> >(q,b,n,iters,seed);
        ok = ok && run_with_field<Modular<int32_t> >(q,b,n,iters,seed);
        ok = ok && run_with_field<ModularBalanced<int32_t> >(q,b,n,iters,seed);
        ok = ok && run_with_field<Modular<int64_t> >(q,b,n,iters,seed);
        ok = ok && run_with_field<ModularBalanced<int64_t> >(q,b,n,iters,seed);
        ok = ok && run_with_field<Modular<Givaro::Integer> >(q,(b?b:512),n/4+1,iters,seed);
    } while (loop && ok);

    return !ok ;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
