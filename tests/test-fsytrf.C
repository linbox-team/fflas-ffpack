/*
 * Copyright (C) 2016 FFLAS-FFPACK
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
//          Test for the computations of the LDLT factorization
//--------------------------------------------------------------------------

#include <iostream>
#include <iterator>
#include <vector>
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &x){
    std::ostream_iterator<T> out_it (os,", ");
    std::copy ( x.begin(), x.end(), out_it );
    return os;
}
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"

#include <iomanip>
#include <random>
#include <chrono>

#include <givaro/modular.h>

#include "fflas-ffpack/utils/test-utils.h"

using namespace std;
using namespace FFPACK;
using namespace FFLAS;


template <class Field,  class RandIter>
bool test_RPM_fsytrf (Field& F, FFLAS_UPLO uplo, string file, size_t n, size_t r, RandIter& G, size_t threshold){

    typename Field::Element_ptr A;
    size_t lda;
    if (!file.empty()){
        ReadMatrix (file, F, n, n, A);
        lda = n;
    } else {
        lda = n+(rand() % 13);
        A = fflas_new (F, n,lda);
        RandomSymmetricMatrixWithRankandRandomRPM (F, n, r, A, lda, G);
    }

    typename Field::Element_ptr B = fflas_new(F, n,lda);
    fassign (F,n,n,A,lda, B, lda);

    size_t * P = fflas_new<size_t>(n);
    size_t rank = fsytrf_RPM (F, uplo, n, A, lda, P, threshold);

    typename Field::Element_ptr T = fflas_new(F, n, n);
    typename Field::Element_ptr U = fflas_new(F, n, n);
    getTridiagonal(F,n,rank,A,lda, P, T, n);
    getTriangular(F,FflasUpper, FflasUnit, n,n,rank,A,lda,  U, n, false);
    fgemm(F,FflasTrans,FflasNoTrans, n,n,n,F.one, U,n,T,n,F.zero,A,lda);
    fgemm(F,FflasNoTrans,FflasNoTrans, n,n,n,F.one,A,lda,U,n,F.zero,T,n);
    for (size_t i=0; i<n; i++)
        if ((int)P[i]<0){
            P[i]=-P[i]-1;
            P[i+1]=-P[i+1]-1;
            i++;
        }
    // FFLAS::WritePermutation(std::cerr<<"test LAPACK Perm =",P,n);
    applyP(F,FflasRight, FflasNoTrans, n, 0, n, T,n,P);
    applyP(F,FflasLeft, FflasTrans, n, 0, n, T,n,P);

    if (!fequal(F,n,n,T,n,B,lda)){
        cerr<<"ERROR: A != P^T U^T D U P"<<endl;
        WriteMatrix(std::cerr<<"PTUTDUP ="<<std::endl,F,n,n,T,n);
        WriteMatrix(std::cerr<<"A ="<<std::endl,F,n,n,B,lda);
        WriteMatrix(std::cerr<<"sortie de l'algo ="<<std::endl,F,n,n,A,lda);
        WriteMatrix(std::cerr<<"Tridiagonal ="<<std::endl,F,n,n,T,n);
        WriteMatrix(std::cerr<<"Triangular ="<<std::endl,F,n,n,U,n);

        fflas_delete(P,A,T,B,U);
        return false;
    }
    if (rank != r){
        cerr<<"ERROR: rank incorrect"<<endl;
        fflas_delete(P,A,T,B,U);
        return false;
    }
    fflas_delete(P,A,T,B,U);
    return true;
}

template <class Field,  class RandIter>
bool test_generic_fsytrf (Field& F, FFLAS_UPLO uplo, string file, size_t n, RandIter& G, size_t threshold){

    typename Field::Element_ptr A;
    size_t lda;
    if (!file.empty()){
        ReadMatrix (file, F, n, n, A);
        lda = n;
    } else {
        lda = n+(rand() % 13);
        A = fflas_new (F, n,lda);
        RandomSymmetricMatrix (F, n, true, A, lda, G);
    }

    typename Field::Element_ptr B = fflas_new(F, n,lda);
    fassign (F,n,n,A,lda, B, lda);

    bool success=FFPACK::fsytrf (F, uplo, n, A, lda, threshold);
    if (!success) cerr<<"Non definite matrix"<<endl;

    if (uplo == FflasLower) { // Testing is B ==  L D L^T
        cout<<"Low";
        // copying L on L^T
        for (size_t i=0; i<n; i++)
            fassign (F, n-i-1, A+i*(lda+1)+lda, lda, A+i*(lda+1)+1, 1);
    } else { // Testing is B ==  U^T D U
        cout<<"Up";
        // copying U on U^T
        for (size_t i=0; i<n; i++)
            fassign(F, n-i-1, A+i*(lda+1)+1, 1, A+i*(lda+1)+lda, lda);
    }
    // L^T <- D L^T or U <- D U
    for (size_t i=0; i<n; i++){
        fscalin (F, n-i-1, A[i*(lda+1)], A+i*(lda+1)+1, 1);
    }
    // A <- L x L^T  or A <- U^T x U
    ftrtrm (F, FflasRight, FflasNonUnit, n, A, lda);

    bool ok = fequal(F, n, n, A, lda, B, lda);

    cout<<"...";

    fflas_delete(A);
    fflas_delete(B);
    return ok;
}

template<class Field>
bool run_with_field(Givaro::Integer q, uint64_t b, size_t n, size_t r, size_t iters, string file, size_t threshold, uint64_t& seed){
    bool ok = true ;
    int nbit=(int)iters;

    while (ok &&  nbit){
        // choose Field
        Field* F= chooseField<Field>(q,b,seed);
        if (F==nullptr)
            return true;
        typename Field::RandIter G(*F,seed++);

        ostringstream oss;
        F->write(oss);

        cout.fill('.');
        cout<<"Checking ";
        cout.width(40);
        cout<<oss.str();
        cout<<" ... ";

        cout<<"GRP..";
        ok = ok && test_generic_fsytrf (*F, FflasUpper, file, n, G, threshold);
        ok = ok && test_generic_fsytrf (*F, FflasLower, file, n, G, threshold);
        cout<<"RPM..";
        ok = ok && test_RPM_fsytrf (*F, FflasUpper, file, n, r, G, threshold);
        //ok = ok && test_RPM_fsytrf (*F, FflasLower, file, n, r, G, threshold);
        size_t NN = (rand() % n)+50;
        size_t RR = (rand() % NN);
        size_t THRESHOLD = (rand() % NN/3)+2;
        cout<<"Random dim...";
        cout<<"GRP..";
        ok = ok && test_generic_fsytrf (*F, FflasUpper, file, NN, G, THRESHOLD);
        ok = ok && test_generic_fsytrf (*F, FflasLower, file, NN, G, THRESHOLD);
        cout<<"RPM..";
        ok = ok && test_RPM_fsytrf (*F, FflasUpper, file, NN, RR, G, THRESHOLD);
        //ok = ok && test_RPM_fsytrf (*F, FflasLower, file, NN,RR, G, THRESHOLD);

        delete F;

        nbit--;
        if (!ok)
            cout << "FAILED "<<endl;
        else
            cout << "PASSED "<<endl;
    }
    return ok;
}

int main(int argc, char** argv){
    cerr<<setprecision(20);

    Givaro::Integer q = -1;
    size_t b = 0;
    size_t n = 280;
    size_t r = 107;
    size_t iters = 6 ;
    bool loop=false;
    uint64_t seed = getSeed();
    size_t threshold =64;
    string file = "";
    Argument as[] = {
        { 'q', "-q Q", "Set the field cardinality.",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
        { 'n', "-n N", "Set the order of the matrix.", TYPE_INT , &n },
        { 'r', "-r R", "Set the rank of the matrix.", TYPE_INT , &r },
        { 'i', "-i I", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
        { 't', "-t T", "Set the threshold to the base case.",    TYPE_INT , &threshold },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        { 'f', "-f file", "Set input file", TYPE_STR, &file },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    srand(seed);
    cerr << "with seed = " << seed << endl;
    if (r>n) r=n;
    bool ok=true;
    do{
        ok = ok && run_with_field<Givaro::Modular<float> >           (q,b,n, r,iters,file,threshold,seed);
        ok = ok && run_with_field<Givaro::Modular<double> >          (q,b,n, r,iters,file,threshold,seed);
        ok = ok && run_with_field<Givaro::ModularBalanced<float> >   (q,b,n, r,iters,file,threshold,seed);
        ok = ok && run_with_field<Givaro::ModularBalanced<double> >   (q,b,n, r,iters,file,threshold,seed);
        ok = ok && run_with_field<Givaro::Modular<int32_t> >   (q,b,n, r,iters,file,threshold,seed);
        ok = ok && run_with_field<Givaro::ModularBalanced<int32_t> >   (q,b,n, r,iters,file,threshold,seed);
        ok = ok && run_with_field<Givaro::Modular<int64_t> >   (q,b,n, r,iters,file,threshold,seed);
        ok = ok && run_with_field<Givaro::ModularBalanced<int64_t> >   (q,b,n, r,iters,file,threshold,seed);
        ok = ok && run_with_field<Givaro::Modular<Givaro::Integer> >(q,(b?b:128),n/4+1,r/4+1,iters,file,threshold,seed);
    } while (loop && ok);

    if (!ok) cerr << "with seed = " << seed << endl;

    return !ok;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
