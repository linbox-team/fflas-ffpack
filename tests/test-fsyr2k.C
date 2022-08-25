/*
 * Copyright (C) 2018 FFLAS-FFPACK
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
#include <random>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"
#include <givaro/modular.h>


using namespace std;
using namespace FFLAS;
using Givaro::Modular;
using Givaro::ModularBalanced;

template<typename Field, class RandIter>
bool check_fsyr2k (const Field &F, size_t n, size_t k,
                   const typename Field::Element &alpha, const typename Field::Element &beta,
                   FFLAS::FFLAS_UPLO uplo, FFLAS::FFLAS_TRANSPOSE trans, RandIter& Rand){

    typedef typename Field::Element Element;
    Element * A, *B, *C, *C2;
    size_t ldc = n+15;
    size_t Arows = (trans==FFLAS::FflasNoTrans)?n:k;
    size_t Acols = (trans==FFLAS::FflasNoTrans)?k:n;
    size_t lda = Acols+(rand()%13);
    size_t ldb = Acols+(rand()%13);

    A  = FFLAS::fflas_new(F,Arows,lda);
    B  = FFLAS::fflas_new(F,Arows,ldb);
    C  = FFLAS::fflas_new(F,n,ldc);
    C2  = FFLAS::fflas_new(F,n,ldc);

    FFPACK::RandomTriangularMatrix (F, n, n, uplo, FflasNonUnit, true, C, ldc, Rand);
    FFPACK::RandomMatrix (F, Arows, Acols, A, lda, Rand);
    FFPACK::RandomMatrix (F, Arows, Acols, B, ldb, Rand);
    FFLAS::fassign (F, n, n, C, ldc, C2, ldc);

    string ss=string((uplo == FFLAS::FflasLower)?"Lower_":"Upper_")+string((trans == FFLAS::FflasTrans)?"Trans":"NoTrans");

    cout<<std::left<<"Checking FSYR2K_";
    cout.fill('.');
    cout.width(35);
    cout<<ss;

    FFLAS::Timer t; t.clear();
    double time=0.0;
    t.clear(); t.start();

    fsyr2k (F, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

    t.stop();
    time+=t.usertime();

    fgemm (F, trans, (trans==FflasNoTrans)?FflasTrans:FflasNoTrans, n, n, k, alpha, A, lda, B, ldb, beta, C2, ldc);
    fgemm (F, trans, (trans==FflasNoTrans)?FflasTrans:FflasNoTrans, n, n, k, alpha, B, ldb, A, lda, F.one, C2, ldc);

    bool ok = true;
    if (uplo == FflasUpper){
        for (size_t i=0; i<n; i++)
            for (size_t j=i; j<n; j++){
                ok = ok && F.areEqual(C2[i*ldc+j], C[i*ldc+j]);
                if (!ok) std::cerr<<"Fail at i = "<<i<<" j = "<<j<<std::endl;}
    } else {
        for (size_t i=0; i<n; i++)
            for (size_t j=0; j<=i; j++)
                ok = ok && F.areEqual(C2[i*ldc+j], C[i*ldc+j]);
    }
    if (ok)
        cout << "PASSED ("<<time<<")"<<endl;
    else{
        cout << "FAILED ("<<time<<")"<<endl;
        WriteMatrix(std::cerr<<"C = "<<std::endl,F,n,n, C, ldc);
        WriteMatrix(std::cerr<<"C2 = "<<std::endl,F,n,n, C2, ldc);
    }

    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);
    FFLAS::fflas_delete(C2);
    FFLAS::fflas_delete(C);
    return ok;
}

template <class Field>
bool run_with_field (Givaro::Integer q, size_t b, size_t n, size_t k, int a, int c, size_t iters, uint64_t seed){
    bool ok = true ;
    int nbit=(int)iters;

    while (ok &&  nbit){
        //typedef typename Field::Element Element ;
        // choose Field
        Field* F= FFPACK::chooseField<Field>(q,b,seed);
        typename Field::RandIter G(*F,seed++);
        if (F==nullptr)
            return true;

        typename Field::Element alpha, beta;
        F->init (alpha, (typename Field::Element)a);
        F->init (beta, (typename Field::Element)c);
        cout<<"Checking with ";F->write(cout)<<endl;

        ok = ok && check_fsyr2k(*F,n,k,alpha,beta,FflasUpper,FflasNoTrans,G);
        ok = ok && check_fsyr2k(*F,n,k,alpha,beta,FflasUpper,FflasTrans,G);
        ok = ok && check_fsyr2k(*F,n,k,alpha,beta,FflasLower,FflasNoTrans,G);
        ok = ok && check_fsyr2k(*F,n,k,alpha,beta,FflasLower,FflasTrans,G);

        // checking with k > n (=k+n)
        ok = ok && check_fsyr2k(*F,n,k+n,alpha,beta,FflasUpper,FflasNoTrans,G);
        ok = ok && check_fsyr2k(*F,n,k+n,alpha,beta,FflasUpper,FflasTrans,G);
        ok = ok && check_fsyr2k(*F,n,k+n,alpha,beta,FflasLower,FflasNoTrans,G);
        ok = ok && check_fsyr2k(*F,n,k+n,alpha,beta,FflasLower,FflasTrans,G);
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
    int k=55;
    int n=109;
    int a=-1;
    int c=1;
    size_t iters=3;
    bool loop=false;
    uint64_t seed = getSeed();

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
        { 'k', "-k K", "Set the  dimension",      TYPE_INT , &k },
        { 'n', "-n N", "Set the column dimension.", TYPE_INT , &n },
        { 'a', "-a A", "Set the scaling alpha",                         TYPE_INT , &a },
        { 'c', "-c C", "Set the scaling beta",                         TYPE_INT , &c },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);

    bool ok = true;
    do{
        ok = ok && run_with_field<Modular<double> >(q,b,n,k,a,c,iters,seed);
        ok = ok && run_with_field<Modular<double> >(2,b,n,k,a,c,iters,seed); // forcing Z/2Z test
        if (q!=2) ok = ok && run_with_field<ModularBalanced<double> >(q,b,n,k,a,c,iters,seed);
        ok = ok && run_with_field<Modular<float> >(q,b,n,k,a,c,iters,seed);
        if (q!=2) ok = ok && run_with_field<ModularBalanced<float> >(q,b,n,k,a,c,iters,seed);
        ok = ok && run_with_field<Modular<int32_t> >(q,b,n,k,a,c,iters,seed);
        if (q!=2) ok = ok && run_with_field<ModularBalanced<int32_t> >(q,b,n,k,a,c,iters,seed);
        ok = ok && run_with_field<Modular<int64_t> >(q,b,n,k,a,c,iters,seed);
        if (q!=2) ok = ok && run_with_field<ModularBalanced<int64_t> >(q,b,n,k,a,c,iters,seed);
        ok = ok && run_with_field<Modular<Givaro::Integer> >(q,5,n/4+1,k/4+1,a,c,iters,seed);
        ok = ok && run_with_field<Modular<Givaro::Integer> >(q,(b?b:512),n/4+1,k/4+1,a,c,iters,seed);
    } while (loop && ok);

    if (!ok) std::cerr<<"with seed = "<<seed<<std::endl;

    return !ok ;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
