/*
 * Copyright (C) 2017 FFLAS-FFPACK
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

#define ENABLE_ALL_CHECKINGS 1

#include "fflas-ffpack/fflas-ffpack-config.h"

#include <iomanip>
#include <iostream>

#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/paladin/parallel.h"
#include "fflas-ffpack/paladin/fflas_plevel1.h"
#include <givaro/zring.h>
#include <givaro/modular.h>

#include <random>
#include <chrono>

using namespace std;
using namespace FFLAS;
using Givaro::Modular;
using Givaro::ModularBalanced;

template<typename Field>
bool check_fdot (const Field &F, size_t n,
                 typename Field::ConstElement_ptr a, size_t inca,
                 typename Field::ConstElement_ptr b, size_t incb){


    std::string st = " inca = " +to_string(inca) + " incb = " + to_string(incb);

    cout<<std::left<<"Checking FDOT";
    cout.fill('.');
    cout.width(30);
    cout<<st;
    FFLAS::Timer t; t.clear();
    double time=0.0;
    t.clear(); t.start();

    typename Field::Element d; F.init(d); F.assign(d, F.zero);

    F.assign(d, fdot (F, n, a, inca, b, incb));

    t.stop();
    time+=t.usertime();

    typename Field::Element dcheck; F.init(dcheck); F.assign(dcheck,F.zero);
    for(size_t i=0; i<n; ++i)
        F.axpyin (dcheck, a[i*inca], b[i*incb]);

    F.subin(d, dcheck);

    cout.fill('.');
    cout.width(17);
    st = "Seq("+to_string(time)+")";
    cout<<st;
    if (! F.areEqual(d,F.zero)){
        F.write(std::cout << " FAILED: diff = ",d)<<endl;
        return false;
    }

    t.clear(); t.start();
    F.assign(d, F.zero);

    PAR_BLOCK {
        FFLAS::ParSeqHelper::Parallel<
        FFLAS::CuttingStrategy::Block,
        FFLAS::StrategyParameter::Threads> Par(NUM_THREADS);

        // d <- d + <A,B>
        fdot(F, n, a, inca, b, incb, d, Par);
    }
    t.stop();
    time=t.usertime();

    F.subin(d, dcheck);
    st = "Par(" +to_string(time)+")";

    cout.fill('.');
    cout.width(17);
    cout << st;
    if (! F.areEqual(d,F.zero)){
        F.write(std::cout << " FAILED diff = ",d)<<endl;
        return false;
    }
    cout << "PASSED"<<endl;
    return true;
}

template <class Field>
bool run_with_field (Givaro::Integer q, size_t BS, size_t n, size_t iters, uint64_t seed){
    bool ok = true ;
    int nbit=(int)iters;

    while (ok &&  nbit){
        Field* F= FFPACK::chooseField<Field>(q,BS,seed);
        typename Field::RandIter G(*F,seed++);
        if (F==nullptr)
            return true;

        cout<<"Checking with ";F->write(cout)<<endl;

        size_t inca = 1 + rand() % n;
        size_t incb = 1 + rand() % n;
        typename Field::Element_ptr a = fflas_new (*F, n, inca);
        typename Field::Element_ptr b = fflas_new (*F, n, incb);

        FFPACK::RandomMatrix (*F, n, inca, a, inca, G);
        FFPACK::RandomMatrix (*F, n, incb, b, incb, G);

        ok = ok && check_fdot(*F,n,a,1,b,1);
        ok = ok && check_fdot(*F,n,a,inca,b,incb);
        ok = ok && check_fdot(*F,n,a,1,b,incb);
        ok = ok && check_fdot(*F,n,a,inca,b,1);

        fflas_delete(a);
        fflas_delete(b);
        nbit--;
        delete F;
    }
    return ok;
}

bool run_with_Integer (size_t BS, size_t n, size_t iters, uint64_t seed){
    bool ok = true ;
    int nbit=(int)iters;
    Givaro::GivRandom generator;
    Givaro::IntegerDom IPD;
    typedef Givaro::ZRing<Givaro::Integer> Field;
    Field F;
    Field::RandIter G(F,seed);
    G.setBitsize(BS);

    while (ok &&  nbit){

        F.write(cout<<"Checking with ") << " and bitsize " << BS <<endl;

        size_t inca = 1 + rand() % n;
        size_t incb = 1 + rand() % n;
        typename Field::Element_ptr a = fflas_new (F, n, inca);
        typename Field::Element_ptr b = fflas_new (F, n, incb);

        PAR_BLOCK { pfrand(F, G, n*inca,1, a); pfrand(F, G, n*incb, 1, b); }

        ok = ok && check_fdot(F,n,a,1,b,1);
        ok = ok && check_fdot(F,n,a,inca,b,incb);
        ok = ok && check_fdot(F,n,a,1,b,incb);
        ok = ok && check_fdot(F,n,a,inca,b,1);

        fflas_delete(a);
        fflas_delete(b);
        nbit--;
    }
    return ok;
}


int main(int argc, char** argv)
{
    cerr<<setprecision(10);
    Givaro::Integer q=-1;
    size_t b=0;
    int n=2578;
    size_t iters=3;
    bool loop=false;
    uint64_t seed = getSeed();
    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
        { 'n', "-n N", "Set the  dimension.", TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);

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
        ok = ok && run_with_field<Modular<Givaro::Integer> >(q,5,n/4+1,iters,seed);
        ok = ok && run_with_field<Modular<Givaro::Integer> >(q,(b?b:512),n/4+1,iters,seed);
        ok = ok && run_with_Integer((b?b:512),n/4+1,iters,seed);
    } while (loop && ok);

    if (!ok) std::cerr<<"with seed = "<<seed<<std::endl;

    return !ok ;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
