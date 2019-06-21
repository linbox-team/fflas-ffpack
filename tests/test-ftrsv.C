/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 * and Clement Pernet <clement.pernet@univ-grenoble-alpes.fr>
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
#define  __FFLASFFPACK_SEQUENTIAL

#define ENABLE_ALL_CHECKINGS 1

#include "fflas-ffpack/fflas-ffpack-config.h"

#include <iomanip>
#include <iostream>
#include <random>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"
#include <givaro/modular.h>

using namespace std;
using namespace FFPACK;
using namespace FFLAS;
using Givaro::Modular;
using Givaro::ModularBalanced;


template<typename Field, class RandIter>
bool check_ftrsv (const Field &F, size_t n, FFLAS_UPLO uplo, FFLAS_TRANSPOSE trans, FFLAS_DIAG diag, RandIter& Rand){

    typedef typename Field::Element Element;
    Element * A, *b, *b2, *c;
    size_t lda = n + (rand() % n );
    size_t incb = 1 + (rand() % 2);
    A  = fflas_new(F,n,lda);
    b  = fflas_new(F,n,incb);
    b2 = fflas_new(F,n,incb);
    c  = fflas_new(F,n,incb);

    RandomTriangularMatrix (F, n, n, uplo, diag, true, A, lda, Rand);
    RandomMatrix (F, n, incb, b, incb, Rand);
    fassign (F, n, incb, b, incb, b2, incb);

    string ss=string((uplo == FflasLower)?"Lower_":"Upper_")+string((trans == FflasTrans)?"Trans_":"NoTrans_")+string((diag == FflasUnit)?"Unit":"NonUnit");

    cout<<std::left<<"Checking FTRSV_";
    cout.fill('.');
    cout.width(30);
    cout<<ss;


    Timer t; t.clear();
    double time=0.0;
    t.clear();
    t.start();
    ftrsv (F, uplo, trans, diag, n, A, lda, b, incb);
    t.stop();
    time+=t.usertime();

    fgemv(F, trans, n, n, F.one, A, lda, b, incb, F.zero, c, incb);

    bool ok = true;
    if (fequal (F, n,  b2, incb, c, incb)){
        //cout << "\033[1;32mPASSED\033[0m ("<<time<<")"<<endl;
        cout << "PASSED ("<<time<<")"<<endl;
        //cerr<<"PASSED ("<<time<<")"<<endl;
    } else{
        //cout << "\033[1;31mFAILED\033[0m ("<<time<<")"<<endl;
        cout << "FAILED ("<<time<<")"<<endl;
        ok=false;
        //cerr<<"FAILED ("<<time<<")"<<endl;
    }

    fflas_delete(A);
    fflas_delete(b);
    fflas_delete(b2);
    fflas_delete(c);
    return ok;
}
template <class Field>
bool run_with_field (Givaro::Integer q, size_t b, size_t n, size_t iters, uint64_t seed){
    bool ok = true ;
    int nbit=(int)iters;

    while (ok &&  nbit){
        //typedef typename Field::Element Element ;
        // choose Field
        Field* F= chooseField<Field>(q,b,seed);
        typename Field::RandIter G(*F,seed++);
        if (F==nullptr)
            return true;

        cout<<"Checking with ";F->write(cout)<<endl;

        ok = ok && check_ftrsv(*F,n,FflasLower,FflasNoTrans,FflasUnit,G);
        ok = ok && check_ftrsv(*F,n,FflasUpper,FflasNoTrans,FflasUnit,G);
        ok = ok && check_ftrsv(*F,n,FflasLower,FflasTrans,FflasUnit,G);
        ok = ok && check_ftrsv(*F,n,FflasUpper,FflasTrans,FflasUnit,G);
        ok = ok && check_ftrsv(*F,n,FflasLower,FflasNoTrans,FflasNonUnit,G);
        ok = ok && check_ftrsv(*F,n,FflasUpper,FflasNoTrans,FflasNonUnit,G);
        ok = ok && check_ftrsv(*F,n,FflasLower,FflasTrans,FflasNonUnit,G);
        ok = ok && check_ftrsv(*F,n,FflasUpper,FflasTrans,FflasNonUnit,G);
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
    size_t n=483;
    size_t iters=1;
    bool loop=false;
    uint64_t seed = getSeed();
    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
        { 'n', "-n N", "Set the dimension of the system.", TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);

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
    } while (loop && ok);

    return !ok ;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
