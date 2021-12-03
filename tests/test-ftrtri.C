/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Clément Pernet <clement.pernet@imag.fr>
 * Philippe LEDENT <philippe.ledent@etu.univ-grenoble-alpes
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
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/ffpack/ffpack.h"


using namespace std;
using namespace FFPACK;
using namespace FFLAS;
using Givaro::Modular;
using Givaro::ModularBalanced;


template<typename Field, class RandIter>
bool check_ftrtri (const Field &F, size_t n, FFLAS_UPLO uplo, FFLAS_DIAG diag, RandIter& Rand){
    typedef typename Field::Element Element;
    Element * A, * B;
    size_t lda = n + (rand() % n );
    A  = fflas_new(F,n,lda);
    B  = fflas_new(F,n,lda);

    RandomTriangularMatrix (F, n, n, uplo, FFLAS::FflasNonUnit, true, A, lda, Rand);
    fassign (F, n, n, A, lda, B, lda); // copy of A

    if (diag == FFLAS::FflasUnit) // Making the implicit unit diagonal explicit on B
        for (size_t i=0; i<n; i++)
            F.assign (B[i*(lda+1)],F.one);

    string ss=string((uplo == FflasLower)?"Lower_":"Upper_")+string((diag == FflasUnit)?"Unit":"NonUnit");

    cout<<std::left<<"Checking FTRTRI_";
    cout.fill('.');
    cout.width(30);
    cout<<ss;
    //	<< endl;


    Timer t; t.clear();
    double time=0.0;
    t.clear();
    t.start();
    ftrtri (F, uplo, diag, n, A, lda);
    t.stop();
    time+=t.usertime();

    // B <- A times B
    ftrmm(F, FFLAS::FflasRight, uplo, FFLAS::FflasNoTrans, diag, n, n, F.one, A, lda, B, lda);

    // Is B the identity matrix ?
    bool ok = true;
    for(size_t li = 0; (li < n) && ok; li++){
        for(size_t co = 0; (co < n) && ok; co++){
            ok = ((li == co) && (F.areEqual(B[li*lda+co],F.one))) || (F.areEqual(B[li*lda+co],F.zero));
        }
    }


    if (ok){
        cout << "PASSED ("<<time<<")"<<endl;
    } else{
        cout << "FAILED ("<<time<<")"<<endl;
        WriteMatrix(std::cout << "\nA^-1" << std::endl, F,n,n,A,lda);
    }

    fflas_delete(A);
    fflas_delete(B);
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

        ok = ok && check_ftrtri(*F,n,FflasLower,FflasUnit,G);
        ok = ok && check_ftrtri(*F,n,FflasUpper,FflasUnit,G);
        ok = ok && check_ftrtri(*F,n,FflasLower,FflasNonUnit,G);
        ok = ok && check_ftrtri(*F,n,FflasUpper,FflasNonUnit,G);
        nbit--;
        delete F;
    }
    if (!ok)
        std::cout << "with seed = "<< seed << std::endl;

    return ok;
}

int main(int argc, char** argv)
{
    cerr<<setprecision(10);
    Givaro::Integer q=-1;
    size_t b=0;
    size_t n=157;
    size_t iters=3;
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
        ok = ok && run_with_field<Modular<Givaro::Integer> >(q,5,n/6+1,iters,seed);
        ok = ok && run_with_field<Modular<Givaro::Integer> >(q,(b?b:512),n/6+1,iters,seed);
    } while (loop && ok);
    return !ok ;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
