/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Hongguang ZHU <zhuhongguang2014@gmail.com>
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
//#define  __FFLASFFPACK_SEQUENTIAL

//#define ENABLE_ALL_CHECKINGS 1

//#include "fflas-ffpack/fflas-ffpack-config.h"
#include <givaro/modular-integer.h>

#include <iomanip>
#include <iostream>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>


using namespace std;
using namespace FFLAS;
using namespace FFPACK;
using Givaro::Modular;
using Givaro::ModularBalanced;

template<typename Field, class RandIter>
bool check_solve(const Field &F, size_t m, RandIter& Rand, bool isParallel){

    typename Field::Element_ptr A, A2, B, B2, x;

    size_t lda,incb,incx;
    lda=m;
    incb=1;
    incx=1;
    A  = FFLAS::fflas_new(F,m,lda);
    A2 = FFLAS::fflas_new(F,m,lda);
    B  = FFLAS::fflas_new(F,m,incb);
    B2 = FFLAS::fflas_new(F,m,incb);
    x  = FFLAS::fflas_new(F,m,incx);

    RandomMatrixWithRank (F,  m,  m, m, A, lda, Rand);

    RandomMatrix (F, m, 1, B, incb, Rand);

    FFLAS::fassign (F, m, B, incb, B2, incb);
    FFLAS::fassign (F, m, m, A, lda, A2, lda);
#ifdef DEBUG
    FFLAS::WriteMatrix(std::cout<<"b:="<<std::endl,F,m,1,B,incb)<<std::endl;
#endif

    FFLAS::Timer t; t.clear();
    double time=0.0;
    t.clear();
    t.start();
    if(isParallel){
        FFPACK::pSolve(F, m, A, lda, x, incx, B, incb);
    }else{
        FFPACK::Solve(F, m, A, lda, x, incx, B, incb);
    }

    t.stop();
    time+=t.realtime();

    FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A2, lda, x, incx, F.zero, B, incb);

    bool ok = true;
    if (FFLAS::fequal (F, m, 1, B2, incb, B, incb)){

        cout << " PASSED ("<<time<<")";

    } else{
#ifdef DEBUG
        FFLAS::WriteMatrix(std::cout<<"A*x:="<<std::endl,F,m,1,B2,incb)<<std::endl;
        FFLAS::WriteMatrix(std::cout<<"b:="<<std::endl,F,m,1,B,incb)<<std::endl;
#endif
        cout << " FAILED ("<<time<<")";
        ok=false;

    }

    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(A2);
    FFLAS::fflas_delete(B);
    FFLAS::fflas_delete(B2);
    FFLAS::fflas_delete(x);
    return ok;
}

template <class Field>
bool run_with_field (Givaro::Integer q, size_t b, size_t m, size_t iters, uint64_t seed){
    bool ok = true ;
    int nbit=(int)iters;

    while (ok &&  nbit){
        //typedef typename Field::Element Element ;
        // choose Field
        Field* F= chooseField<Field>(q,b,seed);
        typename Field::RandIter G(*F,seed++);
        if (F==nullptr)
            return true;

        std::ostringstream oss;
        F->write(oss);
        std::cout.fill('.');
        std::cout<<"Checking ";
        std::cout.width(50);
        std::cout<<oss.str();
        std::cout<<" ... ";

            // testing a sequential run
        std::cout<<" seq: ";
        ok = ok && check_solve(*F,m,G,false);

            // testing a parallel run
        std::cout<<" par: ";
        ok = ok && check_solve(*F,m,G,true);

        std::cout<<std::endl;
        nbit--;
        delete F;

    }

    return ok;
}

int main(int argc, char** argv)
{
    cerr<<setprecision(10);
    Givaro::Integer q=-1;
    size_t b=0;
    size_t m=600;

    size_t iters=4;
    bool loop=false;
    uint64_t seed = getSeed();
    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
        { 'm', "-m M", "Set the dimension of unknown square matrix.",      TYPE_INT , &m },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    bool ok = true;

    do{
        ok = ok && run_with_field<Modular<double> >(q,b,m,iters,seed);
        ok = ok && run_with_field<ModularBalanced<double> >(q,b,m,iters,seed);
        ok = ok && run_with_field<Modular<float> >(q,b,m,iters,seed);
        ok = ok && run_with_field<ModularBalanced<float> >(q,b,m,iters,seed);
        ok = ok && run_with_field<Modular<int32_t> >(q,b,m,iters,seed);
        ok = ok && run_with_field<ModularBalanced<int32_t> >(q,b,m,iters,seed);
        ok = ok && run_with_field<Modular<int64_t> >(q,b,m,iters,seed);
        ok = ok && run_with_field<ModularBalanced<int64_t> >(q,b,m,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,5,m/6,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,(b?b:512),m/6,iters,seed);

    } while (loop && ok);

    return !ok ;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
