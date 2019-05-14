/*
 * Copyright (C) 2018 FFLAS-FFPACK
 * Written by Clement Pernet <clement.pernet@univ-grenoble-alpes.fr>
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
#include "fflas-ffpack/utils/fflas_io.h"

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <givaro/modular-integer.h>

#include <iomanip>
#include <iostream>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>

#include <random>

using namespace std;
using namespace FFLAS;
using namespace FFPACK;
using Givaro::Modular;
using Givaro::ModularBalanced;

template<typename Field, class RandIter>
bool check_ftrstr (const Field &F, size_t n, FFLAS::FFLAS_SIDE side, FFLAS::FFLAS_UPLO uplo, FFLAS::FFLAS_DIAG diagA, FFLAS::FFLAS_DIAG diagB, RandIter& Rand){

    typedef typename Field::Element Element;
    Element * A, *B, *B2, *C;
    size_t lda,ldb;
    lda=n+(rand() % 13);
    ldb=n+(rand() % 14);

    A  = FFLAS::fflas_new(F,n,lda);
    B  = FFLAS::fflas_new(F,n,ldb);
    B2 = FFLAS::fflas_new(F,n,ldb);
    C  = FFLAS::fflas_new(F,n,n);

    RandomTriangularMatrix (F, n, n, uplo, diagA, true, A, lda, Rand);
    RandomTriangularMatrix (F, n, n, uplo, diagB, true, B, ldb, Rand);
    FFLAS::fassign (F, n, n, B, ldb, B2, ldb);

    string ss=string((side == FFLAS::FflasLeft)?"L":"R")+string((uplo == FFLAS::FflasLower)?"Lo":"Up")+string((diagA == FFLAS::FflasUnit)?"U":"N")+string((diagB == FFLAS::FflasUnit)?"U":"N");

    cout<<std::left<<"Checking FTRSTR_";
    cout.fill('.');
    cout.width(25);
    cout<<ss;


    FFLAS::Timer t; t.clear();
    double time=0.0;
    t.clear();
    t.start();
    FFPACK::ftrstr (F, side, uplo, diagA, diagB, n, A, lda, B, ldb);
    t.stop();
    time+=t.usertime();

    FFLAS::ftrmm(F, side, uplo, FflasNoTrans, diagA, n,n,F.one, A, lda, B,ldb);
    bool ok = true;
    if (FFLAS::fequal (F, n, n, B2, ldb, B, ldb)){
        cout << "PASSED ("<<time<<")"<<endl;
    } else{
        cout << "FAILED ("<<time<<")"<<endl;
        WriteMatrix(std::cerr<<"B2 = "<<std::endl,F,n,n,B2,ldb);
        WriteMatrix(std::cerr<<"B = "<<std::endl,F,n,n,B,ldb);
        ok=false;
    }
    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);
    FFLAS::fflas_delete(B2);
    FFLAS::fflas_delete(C);
    return ok;
}
template <class Field>
bool run_with_field (Givaro::Integer q, size_t b,  size_t n, size_t iters, uint64_t seed){
    bool ok = true ;
    int nbit=(int)iters;

    while (ok &&  nbit){
        // choose Field
        Field* F= chooseField<Field>(q,b,seed);
        typename Field::RandIter G(*F,seed++);
        if (F==nullptr)
            return true;

        cout<<"Checking with ";F->write(cout)<<endl;
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNonUnit,FFLAS::FflasNonUnit,G);
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNonUnit,FFLAS::FflasUnit,G);
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasUnit,FFLAS::FflasNonUnit,G);
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasUnit,FFLAS::FflasUnit,G);
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNonUnit,FFLAS::FflasNonUnit,G);
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNonUnit,FFLAS::FflasUnit,G);
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasUnit,FFLAS::FflasNonUnit,G);
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasUnit,FFLAS::FflasUnit,G);
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNonUnit,FFLAS::FflasNonUnit,G);
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNonUnit,FFLAS::FflasUnit,G);
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasUnit,FFLAS::FflasNonUnit,G);
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasUnit,FFLAS::FflasUnit,G);
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNonUnit,FFLAS::FflasNonUnit,G);
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNonUnit,FFLAS::FflasUnit,G);
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasUnit,FFLAS::FflasNonUnit,G);
        ok = ok && check_ftrstr(*F,n,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasUnit,FFLAS::FflasUnit,G);
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
    size_t n=83;
    size_t iters=1;
    bool loop=false;
    uint64_t seed = getSeed();
    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
        { 'n', "-n N", "Set the column dimension of the unknown matrix.", TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);
    std::cerr<<"seed = "<<seed<<std::endl;
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
    } while (loop && ok);

    return !ok ;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
