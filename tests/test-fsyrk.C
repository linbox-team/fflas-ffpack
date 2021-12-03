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

//#define DEBUG_TRACE

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

#ifdef DEBUG_TRACE
    WriteMatrix (std::cerr, F, n, k, A, lda);
    WriteMatrix(std::cerr<<"A=", F, n, k, A, lda,FflasSageMath );
    std::cerr<<"A11=A[:"<<Arows/2<<",:"<<Acols/2
             <<"]; A21=A["<<Arows/2<<":,:"<<Acols/2<<"]; A12=A[:"
             <<Arows/2<<","<<Acols/2<<":]; A22=A["<<Arows/2<<":,"<<Acols/2<<":];"<<std::endl;
    WriteMatrix (std::cerr, F, n, n, C, ldc);
    WriteMatrix(std::cerr, F, n, n, C, ldc,FflasSageMath );
#endif
    if (w == size_t(-1))
        w= (rand() % 5);
        
    MMHelper<Field, MMHelperAlgo::Winograd> H(F,w);
    fsyrk (F, uplo, trans, n, k, alpha, A, lda, beta, C, ldc, H);

    t.stop();
    time+=t.usertime();
#ifdef DEBUG_TRACE
    WriteMatrix (std::cerr<<"Result C = "<<std::endl, F, n, n, C, ldc);
#endif
    
    fgemm (F, trans, (trans==FflasNoTrans)?FflasTrans:FflasNoTrans, n, n, k, alpha, A, lda, A, lda, beta, C2, ldc);

    bool ok = true;
#ifdef DEBUG_TRACE
    std::cerr<<std::endl<<std::endl;
#endif
    if (uplo == FflasUpper){
        for (size_t i=0; i<n; i++){
#ifdef DEBUG_TRACE
            for (size_t j=0;j<i;j++)
                std::cerr<<" ";
#endif
            for (size_t j=i; j<n; j++){
                ok = ok && F.areEqual(C2[i*ldc+j], C[i*ldc+j]);
#ifdef DEBUG_TRACE
                if (F.areEqual(C2[i*ldc+j], C[i*ldc+j]))
                    std::cerr<<".";
                else
                    std::cerr<<"X";
#endif
            }
#ifdef DEBUG_TRACE
            std::cerr<<std::endl;
#endif
        }
    } else {
        for (size_t i=0; i<n; i++){
            for (size_t j=0; j<=i; j++){
                ok = ok && F.areEqual(C2[i*ldc+j], C[i*ldc+j]);
#ifdef DEBUG_TRACE
                if (F.areEqual(C2[i*ldc+j], C[i*ldc+j]))
                    std::cerr<<".";
                else
                    std::cerr<<"X";
#endif
            }
#ifdef DEBUG_TRACE
            std::cerr<<std::endl;
#endif
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


#ifdef DEBUG_TRACE
    FFLAS::WriteMatrix ( std::cerr<<"A = "<<std::endl,F,Arows, Acols, A, lda);
    FFLAS::WriteMatrix ( std::cerr<<"C = "<<std::endl,F,n,n, C, ldc);
    FFLAS::WriteMatrix ( std::cerr<<"D = "<<std::endl,F,k,1,D, incD);
#endif
    FFLAS::Timer t; t.clear();
    double time=0.0;
    t.clear(); t.start();

    fsyrk (F, uplo, trans, n, k, alpha, A, lda, D, incD, beta, C, ldc, 13);

    t.stop();
    time+=t.usertime();

#ifdef DEBUG_TRACE
    std::cerr<<"After fsyrk_diag"<<std::endl;
    FFLAS::WriteMatrix (std::cerr<<"A = "<<std::endl,F,Arows, Acols, A, lda);
    FFLAS::WriteMatrix (std::cerr<<"C = "<<std::endl,F,n,n,C,ldc);
#endif
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
        fflas_delete(A, B, C, C2, D);
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

    fflas_delete(A);
    fflas_delete(B);
    fflas_delete(C2);
    fflas_delete(C);
    fflas_delete(D);
    return ok;
}
template<typename Field, class RandIter>
bool check_fsyrk_bkdiag (const Field &F, size_t n, size_t k,
                         const typename Field::Element &alpha, const typename Field::Element &beta,
                         FFLAS_UPLO uplo, FFLAS_TRANSPOSE trans, RandIter& Rand){

    typedef typename Field::Element Element;
    Element * A, *B, *C, *C2, *D;
    size_t ldc = n+(rand()%50);
    size_t Arows = (trans==FflasNoTrans)?n:k;
    size_t Acols = (trans==FflasNoTrans)?k:n;
    size_t lda = Acols+(rand()%50);
    size_t incD = 1+(rand()%100);
    A  = fflas_new(F,Arows,lda);
    B  = fflas_new(F,Arows,lda);
    C  = fflas_new(F,n,ldc);
    C2  = fflas_new(F,n,ldc);
    D  = fflas_new(F,k,incD);
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
    fassign (F, n, n, C, ldc, C2, ldc);
    fassign (F, Arows, Acols, A, lda, B, lda);

    string ss=string((uplo == FflasLower)?"Lower_":"Upper_")+string((trans == FflasTrans)?"Trans":"NoTrans");

    cout<<std::left<<"Checking FSYRK_BK_DIAG_";
    cout.fill('.');
    cout.width(30);
    cout<<ss;


#ifdef DEBUG_TRACE
    WriteMatrix ( std::cerr<<"A = "<<std::endl,F,Arows, Acols, A, lda);
    WriteMatrix ( std::cerr<<"C = "<<std::endl,F,n,n, C, ldc);
    WriteMatrix ( std::cerr<<"D = "<<std::endl,F,k,1,D, incD);
#endif
    Timer t; t.clear();
    double time=0.0;
    t.clear(); t.start();

    fsyrk (F, uplo, trans, n, k, alpha, A, lda, D, incD, tb, beta, C, ldc, 13);

    t.stop();
    time+=t.usertime();

#ifdef DEBUG_TRACE
    std::cerr<<"After fsyrk_bk_diag"<<std::endl;
    WriteMatrix (std::cerr<<"A = "<<std::endl,F,Arows, Acols, A, lda);
    WriteMatrix (std::cerr<<"C = "<<std::endl,F,n,n,C,ldc);
#endif
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
        fflas_delete(A, B, C, C2, D);
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

    fflas_delete(A);
    fflas_delete(B);
    fflas_delete(C2);
    fflas_delete(C);
    fflas_delete(D);
    return ok;
}

template <class Field, class RandIter>
bool check_computeS1S2 (const Field& F, size_t N, size_t K, FFLAS_TRANSPOSE trans, RandIter& G){

    string ss=string((trans == FflasTrans)?"Trans":"NoTrans");

    cout<<std::left<<"Checking ComputeS1S2_";
    cout.fill('.');
    cout.width(30);
    cout<<ss;

    bool ok = true;
    typename Field::Element_ptr A, S, T;
    size_t STrows = 2*((trans==FflasNoTrans)?N:K);
    size_t STcols = 2*((trans==FflasNoTrans)?K:N);
    size_t Arows = 2*STrows;
    size_t Acols = 2*STcols;
    size_t lda = Acols+(rand()%50);
    size_t lds = STcols+(rand()%50);
    size_t ldt = STcols+(rand()%50);
    A  = fflas_new(F,Arows,lda);
    S  = fflas_new(F,STrows,lds);
    T  = fflas_new(F,STrows,ldt);

    FFPACK::RandomMatrix (F, Arows, Acols, A, lda, G);

    MMHelper<Field, MMHelperAlgo::Winograd, ModeCategories::LazyTag> WH(F,1);
    Givaro::Integer a,b;
    Givaro::IntSqrtModDom<> ISM;
    Givaro::ZRing<Givaro::Integer> Z;
    Z.init(a);
    Z.init(b);
    ISM.sumofsquaresmodprime (a, b, -1, F.characteristic());
    typename Field::Element y1, y2,negy2;
    F.init (y1, a);
    F.init (y2, b);
    
    F.neg(negy2,y2);
#ifdef DEBUG_TRACE
    WriteMatrix(std::cerr<<"A = "<<std::endl, F, Arows, Acols, A, lda);
    WriteMatrix(std::cerr<<"A = ", F, Arows, Acols, A, lda, FflasSageMath);
    std::cerr<<"A11=A[:2,:2]; A21=A[2:,:2]; A12=A[:2,2:]; A22=A[2:,2:]; Y=Matrix(GF("
             <<F.cardinality()<<"),"<<STrows<<","<<STcols<<",["<<y1<<","<<y2<<","<<negy2
             <<","<<y1<<"])"<<std::endl
             <<"S = (A21-A11)*Y; T = A22 - A21*Y"<<std::endl;
    std::cerr<<"y1, y2 = "<<y1<<" "<<y2<<std::endl;
#endif
    computeS1S2 (F, trans, 4*N, 4*K, y1, y2, A, lda, S, lds, T, ldt, WH); 

#ifdef DEBUG_TRACE
    WriteMatrix(std::cerr<<"S = "<<std::endl, F, STrows, STcols, S, lds);
    WriteMatrix(std::cerr<<"T = "<<std::endl, F, STrows, STcols, T, ldt);
#endif
    
    typename Field::Element_ptr V, Y, W;
    typename Field::Element_ptr A21, A22;
    size_t Ydim = 2*K;
    size_t ldy = Ydim;
    size_t ldv = STcols;
    size_t ldw = STcols;
    if (trans == FflasNoTrans){
        A21 = A + 2*N*lda;
        A22 = A21 + 2*K;
    } else {
        A21 = A + 2*N;
        A22 = A21 + 2*K*lda;
    }
    V  = fflas_new (F, STrows, ldv);
    Y  = fflas_new (F, Ydim, Ydim);
    W  = fflas_new (F, STrows, ldw);
    for (size_t i=0; i<Ydim; ++i){
        for (size_t j=0; j<Ydim; ++j){
            F.assign(Y[i*ldy+j], F.zero);
        }
        F.assign (Y [i*ldy+i], y1);
        if (i < K)
            F.assign (Y [i*ldy + K + i], y2);
        else
            F.neg (Y [i*ldy + i - K], y2);
    }
        // W <- A21 - A11
    fsub (F, STrows, STcols, A21, lda, A, lda, W, ldw);
        // S < W x Y (NoTrans) or Y x W (Trans)
    if (trans == FflasNoTrans)
        fgemm (F, FflasNoTrans, FflasNoTrans, STrows, STcols, Ydim, F.one, W, ldw, Y, ldy, F.zero, V, ldv);
    else
        fgemm (F, FflasTrans, FflasNoTrans, STrows, STcols, Ydim, F.one, Y, ldy, W, ldw, F.zero, V, ldv);
    
    if (!fequal (F, STrows, STcols, S, lds, V, ldv)){
        std::cerr<< "FAILED: computing S"<<std::endl;
        ok = false;
    }

        // W <- -A21 x Y (NoTrans) or Y x (-A21) (Trans)
    if (trans == FflasNoTrans)
        fgemm (F, FflasNoTrans, FflasNoTrans, STrows, STcols, Ydim, F.mOne, A21, lda, Y, ldy, F.zero, W, ldw);
    else
        fgemm (F, FflasTrans, FflasNoTrans, STrows, STcols, Ydim, F.mOne, Y, ldy, A21, lda, F.zero, W, ldw);
    
        // V <- A22 + W
    fadd (F, STrows, STcols, A22, lda, W, ldw, V, ldv);
    
    if (!fequal (F, STrows, STcols, T, ldt, V, ldv)){
        std::cerr<< "FAILED: computing T"<<std::endl;
        ok = false;
    }
    
    fflas_delete(A,S,T,V,W,Y);
    
    if (ok) cout << "PASSED"<<endl;

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

        ok = ok && check_fsyrk(*F,n,k,w,alpha,beta,FflasUpper,FflasNoTrans,G);
        ok = ok && check_fsyrk(*F,n,k,w,alpha,beta,FflasUpper,FflasTrans,G);
        ok = ok && check_fsyrk(*F,n,k,w,alpha,beta,FflasLower,FflasNoTrans,G);
        ok = ok && check_fsyrk(*F,n,k,w,alpha,beta,FflasLower,FflasTrans,G);
        ok = ok && check_fsyrk_diag(*F,n,k,alpha,beta,FflasUpper,FflasNoTrans,G);
        ok = ok && check_fsyrk_diag(*F,n,k,alpha,beta,FflasUpper,FflasTrans,G);
        ok = ok && check_fsyrk_diag(*F,n,k,alpha,beta,FflasLower,FflasNoTrans,G);
        ok = ok && check_fsyrk_diag(*F,n,k,alpha,beta,FflasLower,FflasTrans,G);
        ok = ok && check_fsyrk_bkdiag(*F,n,k,alpha,beta,FflasUpper,FflasNoTrans,G);
        ok = ok && check_fsyrk_bkdiag(*F,n,k,alpha,beta,FflasUpper,FflasTrans,G);
        ok = ok && check_fsyrk_bkdiag(*F,n,k,alpha,beta,FflasLower,FflasNoTrans,G);
        ok = ok && check_fsyrk_bkdiag(*F,n,k,alpha,beta,FflasLower,FflasTrans,G);

            // checking with k > n (=k+n)
        ok = ok && check_fsyrk(*F,n,k+n,w,alpha,beta,FflasUpper,FflasNoTrans,G);
        ok = ok && check_fsyrk(*F,n,k+n,w,alpha,beta,FflasUpper,FflasTrans,G);
        ok = ok && check_fsyrk(*F,n,k+n,w,alpha,beta,FflasLower,FflasNoTrans,G);
        ok = ok && check_fsyrk(*F,n,k+n,w,alpha,beta,FflasLower,FflasTrans,G);
        ok = ok && check_fsyrk_diag(*F,n,k+n,alpha,beta,FflasUpper,FflasNoTrans,G);
        ok = ok && check_fsyrk_diag(*F,n,k+n,alpha,beta,FflasUpper,FflasTrans,G);
        ok = ok && check_fsyrk_diag(*F,n,k+n,alpha,beta,FflasLower,FflasNoTrans,G);
        ok = ok && check_fsyrk_diag(*F,n,k+n,alpha,beta,FflasLower,FflasTrans,G);
        ok = ok && check_fsyrk_bkdiag(*F,n,k+n,alpha,beta,FflasUpper,FflasNoTrans,G);
        ok = ok && check_fsyrk_bkdiag(*F,n,k+n,alpha,beta,FflasUpper,FflasTrans,G);
        ok = ok && check_fsyrk_bkdiag(*F,n,k+n,alpha,beta,FflasLower,FflasNoTrans,G);
        ok = ok && check_fsyrk_bkdiag(*F,n,k+n,alpha,beta,FflasLower,FflasTrans,G);

            // Checking the preadditions with the skew othogonal matrix
        ok = ok && check_computeS1S2(*F, n, k, FflasNoTrans, G);
        ok = ok && check_computeS1S2(*F, n, k, FflasTrans, G);
        nbit--;
        delete F;
    }
    return ok;
}

int main(int argc, char** argv)
{
    cerr<<setprecision(10);
    cout<<setprecision(10);

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
        ok = ok && run_with_field<Modular<int64_t> >(q,b,n,k,w,a,c,iters,seed);
        ok = ok && run_with_field<ModularBalanced<int64_t> >(q,b,n,k,w,a,c,iters,seed);
        // conversion to RNS basis not available yet for fsyrk
        // ok = ok && run_with_field<Modular<Givaro::Integer> >(q,5,n/4+1,k/4+1,a,c,iters,seed);
        // ok = ok && run_with_field<Modular<Givaro::Integer> >(q,(b?b:512),n/4+1,k/4+1,a,c,iters,seed);
        seed++;
    } while (loop && ok);

    if (!ok) std::cerr<<"with seed = "<<seed-1<<std::endl;

    return !ok ;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
