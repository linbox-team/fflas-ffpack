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

#define ENABLE_ALL_CHECKINGS 1

#include "fflas-ffpack/fflas-ffpack-config.h"

#include <iomanip>
#include <iostream>
#include <random>

#include "fflas-ffpack/utils/fflas_io.h"
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
bool check_fsyrk (const Field &F, size_t n, size_t k, size_t w,
                  const typename Field::Element &alpha, const typename Field::Element &beta,
                  FFLAS::FFLAS_UPLO uplo, FFLAS::FFLAS_TRANSPOSE trans, RandIter& Rand){

    typedef typename Field::Element Element;
    Element * A, *C, *C2;
    size_t ldc = n+15;
    size_t Arows = (trans==FFLAS::FflasNoTrans)?n:k;
    size_t Acols = (trans==FFLAS::FflasNoTrans)?k:n;
    size_t lda = Acols+13;

    A  = FFLAS::fflas_new(F,Arows,lda);
    C  = FFLAS::fflas_new(F,n,ldc);
    C2  = FFLAS::fflas_new(F,n,ldc);

    FFPACK::RandomTriangularMatrix (F, n, n, uplo, FflasNonUnit, true, C, ldc, Rand);
    FFPACK::RandomMatrix (F, Arows, Acols, A, lda, Rand);
    FFLAS::fassign (F, n, n, C, ldc, C2, ldc);

    string ss=string((uplo == FFLAS::FflasLower)?"Lower_":"Upper_")+string((trans == FFLAS::FflasTrans)?"Trans":"NoTrans");

    cout<<std::left<<"Checking FSYRK_";
    cout.fill('.');
    cout.width(35);
    cout<<ss;

    FFLAS::Timer t; t.clear();
    double time=0.0;
    t.clear(); t.start();

    //     // find a, b such that a^2 + b^2 = -1 mod p
    // Givaro::Integer a,b;
    // Givaro::IntSqrtModDom<> ISM;
    // ISM.sumofsquaresmodprime (a, b, -1, F.characteristic())h;
    // typename Field::Element y1, y2;
    // F.init (y1, a);
    // F.init (y2, b);
    //std::cerr<<"Launching fsyrk_strassen with alpha = "<<alpha<<" beta = "<<beta<<" w = "<<w
            //<<" and "<<a<<"^2 + "<<b<<"^2 = -1"
        //<<std::endl;
    // WriteMatrix (std::cerr, F, n, k, A, lda);
    // WriteMatrix(std::cerr, F, n, k, A, lda,FflasSageMath );
    // WriteMatrix (std::cerr, F, n, n, C, ldc);
    // WriteMatrix(std::cerr, F, n, n, C, ldc,FflasSageMath );
    if (w == size_t(-1))
            //w= (rand() % 5);
        w=1;
    MMHelper<Field, MMHelperAlgo::Winograd> H(F,w);
    fsyrk (F, uplo, trans, n, k, alpha, A, lda, beta, C, ldc, H);

    t.stop();
    time+=t.usertime();
        //WriteMatrix (std::cerr<<"Result C = "<<std::endl, F, n, n, C, ldc);

    fgemm (F, trans, (trans==FflasNoTrans)?FflasTrans:FflasNoTrans, n, n, k, alpha, A, lda, A, lda, beta, C2, ldc);

    bool ok = true;
        //std::cerr<<std::endl<<std::endl;
    if (uplo == FflasUpper){
        for (size_t i=0; i<n; i++){
            for (size_t j=0;j<i;j++)
                    //std::cerr<<" ";
            for (size_t j=i; j<n; j++){
                ok = ok && F.areEqual(C2[i*ldc+j], C[i*ldc+j]);
                    //if (F.areEqual(C2[i*ldc+j], C[i*ldc+j]))
                        //std::cerr<<".";
                    //else
                        //std::cerr<<"X";
            }
                //std::cerr<<std::endl;
        }
    } else {
        for (size_t i=0; i<n; i++){
            for (size_t j=0; j<=i; j++){
                ok = ok && F.areEqual(C2[i*ldc+j], C[i*ldc+j]);
                    if (F.areEqual(C2[i*ldc+j], C[i*ldc+j]))
                            std::cerr<<".";
                    else
                        std::cerr<<"X";
            }
            std::cerr<<std::endl;
        }
    }
    if (ok)
        //cout << "\033[1;32mPASSED\033[0m ("<<time<<")"<<endl;
        cout << "PASSED ("<<time<<")"<<endl;
    //cerr<<"PASSED ("<<time<<")"<<endl;
    else
        //cout << "\033[1;31mFAILED\033[0m ("<<time<<")"<<endl;
        cout << "FAILED ("<<time<<")"<<endl;
    //cerr<<"FAILED ("<<time<<")"<<endl;

    if (!ok){
    }
    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(C2);
    FFLAS::fflas_delete(C);
    return ok;
}
template<typename Field, class RandIter>
bool check_fsyrk_diag (const Field &F, size_t n, size_t k,
                       const typename Field::Element &alpha, const typename Field::Element &beta,
                       FFLAS::FFLAS_UPLO uplo, FFLAS::FFLAS_TRANSPOSE trans, RandIter& Rand){

    typedef typename Field::Element Element;
    Element * A, *B, *C, *C2, *D;
    size_t ldc = n+(rand()%50);
    size_t Arows = (trans==FFLAS::FflasNoTrans)?n:k;
    size_t Acols = (trans==FFLAS::FflasNoTrans)?k:n;
    size_t lda = Acols+(rand()%50);
    size_t incD = 1+(rand()%100);
    A  = FFLAS::fflas_new(F,Arows,lda);
    B  = FFLAS::fflas_new(F,Arows,lda);
    C  = FFLAS::fflas_new(F,n,ldc);
    C2  = FFLAS::fflas_new(F,n,ldc);
    D  = FFLAS::fflas_new(F,k,incD);
    Givaro::GeneralRingNonZeroRandIter<Field,RandIter> nzRand (Rand);
    for (size_t i=0; i<k; i++)
        nzRand.random(D[i*incD]);
    FFPACK::RandomTriangularMatrix (F, n, n, uplo, FflasNonUnit, true, C, ldc, Rand);
    FFPACK::RandomMatrix (F, Arows, Acols, A, lda, Rand);
    FFLAS::fassign (F, n, n, C, ldc, C2, ldc);
    FFLAS::fassign (F, Arows, Acols, A, lda, B, lda);

    string ss=string((uplo == FFLAS::FflasLower)?"Lower_":"Upper_")+string((trans == FFLAS::FflasTrans)?"Trans":"NoTrans");

    cout<<std::left<<"Checking FSYRK_DIAG_";
    cout.fill('.');
    cout.width(30);
    cout<<ss;


    // FFLAS::WriteMatrix ( std::cerr<<"A = "<<std::endl,F,Arows, Acols, A, lda);
    // FFLAS::WriteMatrix ( std::cerr<<"C = "<<std::endl,F,n,n, C, ldc);
    // FFLAS::WriteMatrix ( std::cerr<<"D = "<<std::endl,F,k,1,D, incD);
    FFLAS::Timer t; t.clear();
    double time=0.0;
    t.clear(); t.start();

    fsyrk (F, uplo, trans, n, k, alpha, A, lda, D, incD, beta, C, ldc, 13);

    t.stop();
    time+=t.usertime();

    // std::cerr<<"After fsyrk_diag"<<std::endl;
    //  FFLAS::WriteMatrix (std::cerr<<"A = "<<std::endl,F,Arows, Acols, A, lda);
    //  FFLAS::WriteMatrix (std::cerr<<"C = "<<std::endl,F,n,n,C,ldc);

    bool ok = true;

    typename Field::Element tmp;
    F.init(tmp);
    if (trans==FflasNoTrans){
        // Checking whether  A = B x D
        for (size_t i=0; i < Arows; i++)
            for (size_t j=0; j < Acols; j++)
                ok = ok && F.areEqual(A[i*lda+j], F.mul (tmp, B[i*lda+j], D[j*incD]));
    } else {
        // Checking whether  A = D x B
        for (size_t i=0; i < Arows; i++)
            for (size_t j=0; j < Acols; j++){
                if(!F.areEqual(A[i*lda+j], F.mul (tmp, B[i*lda+j], D[i*incD]))){
                    std::cerr<<"B["<<i<<","<<j<<"] = "<<B[i*lda+j]<<" != "<<D[i*incD]<<" * "<<A[i*lda+j]<<std::endl;
                    ok=false;
                }
            }
    }
    if (!ok){
        std::cerr<<"Scaling failed"<<std::endl;
        FFLAS::fflas_delete(A, B, C, C2, D);
        return ok;
    }

    fgemm (F, trans, (trans==FflasNoTrans)?FflasTrans:FflasNoTrans, n, n, k, alpha, A, lda, B, lda, beta, C2, ldc);

    if (uplo == FflasUpper){
        for (size_t i=0; i<n; i++)
            for (size_t j=i; j<n; j++){
                if(!F.areEqual(C2[i*ldc+j], C[i*ldc+j])){
                    std::cerr<<"C2["<<i<<","<<j<<"] = "<<C2[i*ldc+j]<<" != C["<<i<<","<<j<<"] = "<< C[i*ldc+j]<<std::endl;
                    ok=false;
                }
            }
    } else {
        for (size_t i=0; i<n; i++)
            for (size_t j=0; j<=i; j++)
                ok = ok && F.areEqual(C2[i*ldc+j], C[i*ldc+j]);
    }
    if (ok)
        cout << "PASSED ("<<time<<")"<<endl;
    else
        cout << "FAILED ("<<time<<")"<<endl;

    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);
    FFLAS::fflas_delete(C2);
    FFLAS::fflas_delete(C);
    FFLAS::fflas_delete(D);
    return ok;
}
template<typename Field, class RandIter>
bool check_fsyrk_bkdiag (const Field &F, size_t n, size_t k,
                         const typename Field::Element &alpha, const typename Field::Element &beta,
                         FFLAS::FFLAS_UPLO uplo, FFLAS::FFLAS_TRANSPOSE trans, RandIter& Rand){

    typedef typename Field::Element Element;
    Element * A, *B, *C, *C2, *D;
    size_t ldc = n+(rand()%50);
    size_t Arows = (trans==FFLAS::FflasNoTrans)?n:k;
    size_t Acols = (trans==FFLAS::FflasNoTrans)?k:n;
    size_t lda = Acols+(rand()%50);
    size_t incD = 1+(rand()%100);
    A  = FFLAS::fflas_new(F,Arows,lda);
    B  = FFLAS::fflas_new(F,Arows,lda);
    C  = FFLAS::fflas_new(F,n,ldc);
    C2  = FFLAS::fflas_new(F,n,ldc);
    D  = FFLAS::fflas_new(F,k,incD);
    Givaro::GeneralRingNonZeroRandIter<Field,RandIter> nzRand (Rand);
    std::vector<bool> tb(k,false);
    for (size_t i=0; i<k; i++){
        if (rand()% 10 != 0 || i==k-1)
            nzRand.random(D[i*incD]);
        else{
            nzRand.random(D[i*incD]);
            F.assign(D[(i+1)*incD],D[i*incD]);
            tb[i]=true;
            i++;
        }
    }

    FFPACK::RandomTriangularMatrix (F, n, n, uplo, FflasNonUnit, true, C, ldc, Rand);
    FFPACK::RandomMatrix (F, Arows, Acols, A, lda, Rand);
    FFLAS::fassign (F, n, n, C, ldc, C2, ldc);
    FFLAS::fassign (F, Arows, Acols, A, lda, B, lda);

    string ss=string((uplo == FFLAS::FflasLower)?"Lower_":"Upper_")+string((trans == FFLAS::FflasTrans)?"Trans":"NoTrans");

    cout<<std::left<<"Checking FSYRK_BK_DIAG_";
    cout.fill('.');
    cout.width(30);
    cout<<ss;


    // FFLAS::WriteMatrix ( std::cerr<<"A = "<<std::endl,F,Arows, Acols, A, lda);
    // FFLAS::WriteMatrix ( std::cerr<<"C = "<<std::endl,F,n,n, C, ldc);
    // FFLAS::WriteMatrix ( std::cerr<<"D = "<<std::endl,F,k,1,D, incD);
    FFLAS::Timer t; t.clear();
    double time=0.0;
    t.clear(); t.start();

    fsyrk (F, uplo, trans, n, k, alpha, A, lda, D, incD, tb, beta, C, ldc, 13);

    t.stop();
    time+=t.usertime();

    // std::cerr<<"After fsyrk_bk_diag"<<std::endl;
    // FFLAS::WriteMatrix (std::cerr<<"A = "<<std::endl,F,Arows, Acols, A, lda);
    // FFLAS::WriteMatrix (std::cerr<<"C = "<<std::endl,F,n,n,C,ldc);
    bool ok = true;

    typename Field::Element tmp;
    F.init(tmp);
    if (trans==FflasNoTrans){
        // Checking whether  A = B x D
        for (size_t i=0; i < Arows; i++)
            for (size_t j=0; j < Acols; j++)
                if (!tb[j])
                    ok = ok && F.areEqual(A[i*lda+j], F.mul (tmp, B[i*lda+j], D[j*incD]));
                else{
                    ok = ok && F.areEqual(A[i*lda+j], F.mul (tmp, B[i*lda+j+1], D[j*incD]));
                    ok = ok && F.areEqual(A[i*lda+j+1], F.mul (tmp, B[i*lda+j], D[j*incD]));
                    j++;
                }

    } else {
        // Checking whether  A = D x B
        for (size_t j=0; j < Acols; j++){
            for (size_t i=0; i < Arows; i++)
                if (!tb[i])
                    ok = ok && F.areEqual(A[i*lda+j], F.mul (tmp, B[i*lda+j], D[i*incD]));
                else{
                    ok = ok && F.areEqual(A[i*lda+j], F.mul (tmp, B[(i+1)*lda+j], D[i*incD]));
                    ok = ok && F.areEqual(A[(i+1)*lda+j], F.mul (tmp, B[i*lda+j], D[i*incD]));
                    i++;
                }
            // {
            // 	std::cerr<<"B["<<i<<","<<j<<"] = "<<B[i*lda+j]<<" != "<<D[i*incD]<<" * "<<A[i*lda+j]<<std::endl;
            // 	ok=false;
            // }
        }
    }
    if (!ok){
        std::cerr<<"Scaling failed"<<std::endl;
        std::cerr<<"alpha = "<<alpha<<" beta="<<beta<<std::endl;
        std::cerr<<"tb = "<<tb<<std::endl;
        FFLAS::fflas_delete(A, B, C, C2, D);
        return ok;
    }

    fgemm (F, trans, (trans==FflasNoTrans)?FflasTrans:FflasNoTrans, n, n, k, alpha, A, lda, B, lda, beta, C2, ldc);

    if (uplo == FflasUpper){
        for (size_t i=0; i<n; i++)
            for (size_t j=i; j<n; j++){
                if(!F.areEqual(C2[i*ldc+j], C[i*ldc+j])){
                    std::cerr<<"C2["<<i<<","<<j<<"] = "<<C2[i*ldc+j]<<" != C["<<i<<","<<j<<"] = "<< C[i*ldc+j]<<std::endl;
                    ok=false;
                }
            }
    } else {
        for (size_t i=0; i<n; i++)
            for (size_t j=0; j<=i; j++)
                ok = ok && F.areEqual(C2[i*ldc+j], C[i*ldc+j]);
    }
    if (ok)
        //cout << "\033[1;32mPASSED\033[0m ("<<time<<")"<<endl;
        cout << "PASSED ("<<time<<")"<<endl;
    //cerr<<"PASSED ("<<time<<")"<<endl;
    else
        //cout << "\033[1;31mFAILED\033[0m ("<<time<<")"<<endl;
        cout << "FAILED ("<<time<<")"<<endl;
    //cerr<<"FAILED ("<<time<<")"<<endl;

    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);
    FFLAS::fflas_delete(C2);
    FFLAS::fflas_delete(C);
    FFLAS::fflas_delete(D);
    return ok;
}

template <class Field>
bool run_with_field (Givaro::Integer q, size_t b, size_t n, size_t k, size_t w, int a, int c, size_t iters, uint64_t seed){
    bool ok = true ;
    int nbit=(int)iters;

    while (ok &&  nbit){
        //typedef typename Field::Element Element ;
        // choose Field
        Field* F= FFPACK::chooseField<Field>(q,b,seed);
        if (F==nullptr)
            return true;
        typename Field::RandIter G(*F,seed++);

        typename Field::Element alpha, beta;
        F->init (alpha, a);
        F->init (beta, c);
        cout<<"Checking with ";F->write(cout)<<endl;

            // ok = ok && check_fsyrk(*F,n,k,w,alpha,beta,FflasUpper,FflasNoTrans,G);
        // ok = ok && check_fsyrk(*F,n,k,w,alpha,beta,FflasUpper,FflasTrans,G);
            // ok = ok && check_fsyrk(*F,n,k,w,alpha,beta,FflasLower,FflasNoTrans,G);
        ok = ok && check_fsyrk(*F,n,k,w,alpha,beta,FflasLower,FflasTrans,G);
        // ok = ok && check_fsyrk_diag(*F,n,k,w,alpha,beta,FflasUpper,FflasNoTrans,G);
        // ok = ok && check_fsyrk_diag(*F,n,k,w,alpha,beta,FflasUpper,FflasTrans,G);
        // ok = ok && check_fsyrk_diag(*F,n,k,w,alpha,beta,FflasLower,FflasNoTrans,G);
        // ok = ok && check_fsyrk_diag(*F,n,k,w,alpha,beta,FflasLower,FflasTrans,G);
        // ok = ok && check_fsyrk_bkdiag(*F,n,k,w,alpha,beta,FflasUpper,FflasNoTrans,G);
        // ok = ok && check_fsyrk_bkdiag(*F,n,k,w,alpha,beta,FflasUpper,FflasTrans,G);
        // ok = ok && check_fsyrk_bkdiag(*F,n,k,w,alpha,beta,FflasLower,FflasNoTrans,G);
        // ok = ok && check_fsyrk_bkdiag(*F,n,k,w,alpha,beta,FflasLower,FflasTrans,G);

        // // checking with k > n (=k+n)
        // ok = ok && check_fsyrk(*F,n,k+n,w,alpha,beta,FflasUpper,FflasNoTrans,G);
        // ok = ok && check_fsyrk(*F,n,k+n,w,alpha,beta,FflasUpper,FflasTrans,G);
        // ok = ok && check_fsyrk(*F,n,k+n,w,alpha,beta,FflasLower,FflasNoTrans,G);
        // ok = ok && check_fsyrk(*F,n,k+n,w,alpha,beta,FflasLower,FflasTrans,G);
        // ok = ok && check_fsyrk_diag(*F,n,k+n,w,alpha,beta,FflasUpper,FflasNoTrans,G);
        // ok = ok && check_fsyrk_diag(*F,n,k+n,w,alpha,beta,FflasUpper,FflasTrans,G);
        // ok = ok && check_fsyrk_diag(*F,n,k+n,w,alpha,beta,FflasLower,FflasNoTrans,G);
        // ok = ok && check_fsyrk_diag(*F,n,k+n,w,alpha,beta,FflasLower,FflasTrans,G);
        // ok = ok && check_fsyrk_bkdiag(*F,n,k+n,w,alpha,beta,FflasUpper,FflasNoTrans,G);
        // ok = ok && check_fsyrk_bkdiag(*F,n,k+n,w,alpha,beta,FflasUpper,FflasTrans,G);
        // ok = ok && check_fsyrk_bkdiag(*F,n,k+n,w,alpha,beta,FflasLower,FflasNoTrans,G);
        // ok = ok && check_fsyrk_bkdiag(*F,n,k+n,w,alpha,beta,FflasLower,FflasTrans,G);
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
    int k=125;
    int n=219;
    int w=-1;
    int a=-1;
    int c=1;
    size_t iters=3;
    bool loop=false;
    uint64_t seed = getSeed();

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
        { 'k', "-k K", "Set the dimension k",      TYPE_INT , &k },
        { 'n', "-n N", "Set the dimension n.", TYPE_INT , &n },
        { 'w', "-w W", "Set the number of recursive levels for Strassen.", TYPE_INT , &w },
        { 'a', "-a A", "Set the scaling alpha",                         TYPE_INT , &a },
        { 'c', "-c C", "Set the scaling beta",                         TYPE_INT , &c },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);
    srand(seed);
    bool ok = true;
    do{
        ok = ok && run_with_field<Modular<double> >(q,b,n,k,w,a,c,iters,seed);
        ok = ok && run_with_field<ModularBalanced<double> >(q,b,n,k,w,a,c,iters,seed);
        ok = ok && run_with_field<Modular<float> >(q,b,n,k,w,a,c,iters,seed);
        ok = ok && run_with_field<ModularBalanced<float> >(q,b,n,k,w,a,c,iters,seed);
        // ok = ok && run_with_field<Modular<int32_t> >(q,b,n,k,a,c,iters,seed);
        // ok = ok && run_with_field<ModularBalanced<int32_t> >(q,b,n,k,a,c,iters,seed);
        // ok = ok && run_with_field<Modular<int64_t> >(q,b,n,k,a,c,iters,seed);
        // ok = ok && run_with_field<ModularBalanced<int64_t> >(q,b,n,k,a,c,iters,seed);

        // conversion to RNS basis not available yet for fsyrk
        // ok = ok && run_with_field<Modular<Givaro::Integer> >(q,5,n/4+1,k/4+1,a,c,iters,seed);
        // ok = ok && run_with_field<Modular<Givaro::Integer> >(q,(b?b:512),n/4+1,k/4+1,a,c,iters,seed);
    } while (loop && ok);

    if (!ok) std::cerr<<"with seed = "<<seed<<std::endl;

    return !ok ;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
