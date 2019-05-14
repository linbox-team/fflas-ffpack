/*
 * Copyright (C) FFLAS-FFPACK 2017
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

// CP: disabling ALL_CHECKINGS as they create
// - conditional jumps on unitialized values (valgrind)
// - segfaults on travis virtual machines
// TODO: fix this (unsuccessful so far) https://github.com/linbox-team/fflas-ffpack/issues/123
//#define ENABLE_ALL_CHECKINGS 1

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <givaro/modular-integer.h>

#include <iomanip>
#include <iostream>
#include <random>

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
bool check_ftrmm (const Field &F, size_t m, size_t n, const typename Field::Element &alpha, FFLAS::FFLAS_SIDE side, FFLAS::FFLAS_UPLO uplo, FFLAS::FFLAS_TRANSPOSE trans, FFLAS::FFLAS_DIAG diag, RandIter& Rand){

    typedef typename Field::Element Element;
    Element * A, *B, *B2, *C;
    size_t k = (side==FFLAS::FflasLeft?m:n);
    size_t lda,ldb,ldc;
    lda=k+13;
    ldb=n+14;
    ldc=n+15;
    A  = FFLAS::fflas_new(F,k,lda);
    B  = FFLAS::fflas_new(F,m,ldb);
    B2 = FFLAS::fflas_new(F,m,ldb);
    C  = FFLAS::fflas_new(F,m,ldc);

    RandomTriangularMatrix (F, k, k, uplo, diag, true, A, lda, Rand);
    RandomMatrix (F, m, n, B, ldb, Rand);
    FFLAS::fassign (F, m, n, B, ldb, B2, ldb);

    string ss=string((uplo == FFLAS::FflasLower)?"Lower_":"Upper_")+string((side == FFLAS::FflasLeft)?"Left_":"Right_")+string((trans == FFLAS::FflasTrans)?"Trans_":"NoTrans_")+string((diag == FFLAS::FflasUnit)?"Unit":"NonUnit");

    cout<<std::left<<"Checking FTRMM_";
    cout.fill('.');
    cout.width(35);
    cout<<ss;


    FFLAS::Timer t; t.clear();
    double time=0.0;
    t.clear();
    t.start();
    FFLAS::ftrmm (F, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
    t.stop();
    time+=t.usertime();


    if (side == FFLAS::FflasLeft)
        FFLAS::fgemm(F, trans, FFLAS::FflasNoTrans, m, n, m, alpha, A, lda, B2, ldb, F.zero, C, ldc);
    else
        FFLAS::fgemm(F, FFLAS::FflasNoTrans, trans, m, n, n, alpha, B2, ldb, A, lda, F.zero, C, ldc);

    bool ok = true;
    if (FFLAS::fequal (F, m, n, B, ldb, C, ldc)){
        //cout << "\033[1;32mPASSED\033[0m ("<<time<<")"<<endl;
        cout << "PASSED ("<<time<<")"<<endl;
        //cerr<<"PASSED ("<<time<<")"<<endl;
    } else{
        //cout << "\033[1;31mFAILED\033[0m ("<<time<<")"<<endl;
        cout << "FAILED ("<<time<<")"<<endl;
        ok=false;
        //cerr<<"FAILED ("<<time<<")"<<endl;
    }

    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);
    FFLAS::fflas_delete(B2);
    FFLAS::fflas_delete(C);
    return ok;
}
template <class Field>
bool run_with_field (Givaro::Integer q, size_t b, size_t m, size_t n, uint64_t a, size_t iters, uint64_t seed){
    bool ok = true ;
    int nbit=(int)iters;

    while (ok &&  nbit){
        //typedef typename Field::Element Element ;
        // choose Field
        Field* F= chooseField<Field>(q,b,seed);
        typename Field::RandIter G(*F,seed++);
        if (F==nullptr)
            return true;

        typename Field::Element alpha;
        F->init (alpha, (typename Field::Element)a);
        cout<<"Checking with ";F->write(cout)<<endl;

        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit,G);
        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit,G);
        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit,G);
        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit,G);
        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit,G);
        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit,G);
        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit,G);
        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit,G);
        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,G);
        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,G);
        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit,G);
        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit,G);
        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,G);
        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,G);
        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit,G);
        ok = ok && check_ftrmm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit,G);
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
    size_t m=128;
    size_t n=128;
    size_t a=1;
    size_t iters=1;
    bool loop=false;
    uint64_t seed =  getSeed();
    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
        { 'm', "-m M", "Set the row dimension of unknown matrix.",      TYPE_INT , &m },
        { 'n', "-n N", "Set the column dimension of the unknown matrix.", TYPE_INT , &n },
        { 'a', "-a A", "Set the scaling of trmm",                         TYPE_INT , &a },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    bool ok = true;
    do{
        ok = ok && run_with_field<Modular<double> >(q,b,m,n,a,iters,seed);
        ok = ok && run_with_field<ModularBalanced<double> >(q,b,m,n,a,iters,seed);
        ok = ok && run_with_field<Modular<float> >(q,b,m,n,a,iters,seed);
        ok = ok && run_with_field<ModularBalanced<float> >(q,b,m,n,a,iters,seed);
        ok = ok && run_with_field<Modular<int32_t> >(q,b,m,n,a,iters,seed);
        ok = ok && run_with_field<ModularBalanced<int32_t> >(q,b,m,n,a,iters,seed);
        ok = ok && run_with_field<Modular<int64_t> >(q,b,m,n,a,iters,seed);
        ok = ok && run_with_field<ModularBalanced<int64_t> >(q,b,m,n,a,iters,seed);
        ok = ok && run_with_field<Modular<Givaro::Integer> >(q,5,m/4+1,n/4+1,a,iters,seed);
        ok = ok && run_with_field<Modular<Givaro::Integer> >(q,(b?b:512),m/4+1,n/4+1,a,iters,seed);
    } while (loop && ok);

    return !ok ;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
