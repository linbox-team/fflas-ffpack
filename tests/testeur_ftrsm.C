/*
 * Copyright (C) 2007 FFLAS-FFPACK
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
//                        Sanity check for ftrsm and ftrmm
//
//--------------------------------------------------------------------------

#include <iomanip>
#include <iostream>
#include "givaro/modular-balanced.h"
//#include "fflas-ffpack/field/modular-int.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "givaro/givintprime.h"

using namespace std;
using namespace FFPACK;

//typedef Givaro::Modular<int> Field;
//typedef Givaro::Modular<float> Field;
typedef ModularBalanced<double> Field;

int main(int argc, char** argv){


    FFLAS::Timer tim;
    Givaro::IntPrimeDom IPD;
    uint64_t p;
    size_t M, N, K ;
    bool keepon = true;
    Givaro::Integer _p,tmp;
    cerr<<setprecision(10);

    size_t TMAX = 300;
    size_t PRIMESIZE = 23;
    if (argc > 1 )
        TMAX = atoi(argv[1]);
    if (argc > 2 )
        PRIMESIZE = atoi(argv[2]);

    FFLAS::FFLAS_TRANSPOSE trans;
    FFLAS::FFLAS_SIDE side;
    FFLAS::FFLAS_UPLO uplo;
    FFLAS::FFLAS_DIAG diag;
    size_t lda, ldb;

    Field::Element * A, *Abis, *B,* Bbis;
    Field::Element alpha;

    while (keepon){
        srandom(_p);
        do{
            //		max = Integer::random(2);
            _p = random();//max % (2<<30);
            IPD.prevprime( tmp, (_p% (1<<PRIMESIZE)) );
            p =  tmp;
        }while( (p <= 2) );

        Field F (p);
        Field::RandIter RValue (F);

        do{
            M = (size_t)  random() % TMAX;
            N = (size_t)  random() % TMAX;
        } while ((M == 0) || (N == 0));

        ldb = N;

        if (random()%2)
            trans = FFLAS::FflasNoTrans;
        else
            trans = FFLAS::FflasTrans;


        if (random()%2)
            diag = FFLAS::FflasUnit;
        else
            diag = FFLAS::FflasNonUnit;

        if (random()%2){
            side = FFLAS::FflasLeft;
            K = M;
            lda = M;
        } else {
            side = FFLAS::FflasRight;
            K = N;
            lda = N;
        }

        if (random()%2)
            uplo = FFLAS::FflasUpper;
        else
            uplo = FFLAS::FflasLower;

        while (F.isZero(RValue.random (alpha)));

        A = FFLAS::fflas_new<Field::Element>(K*K);
        B = FFLAS::fflas_new<Field::Element>(M*N);
        Abis = FFLAS::fflas_new<Field::Element>(K*K);
        Bbis = FFLAS::fflas_new<Field::Element>(M*N);
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j){
                RValue.random (*(B + i*N + j));
                *(Bbis + i*N + j) = *(B + i*N + j);
            }
        for (size_t i = 0; i < K; ++i)
            for (size_t j = 0; j < K; ++j)
                *(Abis + i*K + j) = RValue.random (*(A + i*K + j));
        for (size_t i = 0; i < K; ++i){
            while (F.isZero(RValue.random (*(A + i*(K+1)))));
            *(Abis + i*(K +1)) = *(A + i*(K+1));
        }

        cout <<"p = "<<(size_t)p
        <<" M = "<<M
        <<" N = "<<N
        <<((side==FFLAS::FflasLeft)?" Left ":" Right ")
        <<((uplo==FFLAS::FflasLower)?" Lower ":" Upper ")
        <<((trans==FFLAS::FflasTrans)?" Trans ":" NoTrans ")
        <<((diag==FFLAS::FflasUnit)?" Unit ":" NonUnit ")
        <<"....";


        tim.clear();
        tim.start();
        FFLAS::ftrsm (F, side, uplo, trans, diag, M, N, alpha,
                      A, lda, B, ldb);
        tim.stop();

        // Verification
        Field::Element invalpha;
        F.inv(invalpha, alpha);
        FFLAS::ftrmm (F, side, uplo, trans, diag, M, N, invalpha,
                      A, K, B, N);
        for (size_t i = 0;i < M;++i)
            for (size_t j = 0;j < N; ++j)
                if ( !F.areEqual (*(Bbis + i*N+ j ), *(B + i*N + j))){
                    cerr<<endl
                    <<"Bbis ["<<i<<", "<<j<<"] = "<<(*(Bbis + i*N + j))
                    <<" ; B ["<<i<<", "<<j<<"] = "<<(*(B + i*N + j));

                    keepon = false;
                }
        for (size_t i = 0;i < K; ++i)
            for (size_t j = 0;j < K; ++j)
                if ( !F.areEqual (*(A + i*K + j), *(Abis + i*K + j))){
                    cerr<<endl
                    <<"A ["<<i<<", "<<j<<"] = "<<(*(A + i*K + j))
                    <<" ; Abis ["<<i<<", "<<j<<"] = "<<(*(Abis + i*K + j));
                    keepon = false;
                }
        if (keepon) {
            cout<<" Passed "
            <<double(M*N)/1000000.0*double(K)/tim.usertime()<<" Mfops"<<endl;

            FFLAS::fflas_delete( B);
            FFLAS::fflas_delete( Bbis);
            FFLAS::fflas_delete( A);
            FFLAS::fflas_delete( Abis);
        } else {

            cerr<<endl;
            FFLAS::WriteMatrix (cerr<<"A = "<<endl, F,  K, K, Abis, K);
            FFLAS::WriteMatrix (cerr<<"B = "<<endl, F,  M, N, Bbis, N);
        }
    }

    cout<<endl;
    cerr<<"FAILED with p = "<<(size_t)p
    <<" M = "<<M
    <<" N = "<<N
    <<" alpha = "<<alpha
    <<((side==FFLAS::FflasLeft)?" Left ":" Right ")
    <<((uplo==FFLAS::FflasLower)?" Lower ":" Upper ")
    <<((trans==FFLAS::FflasTrans)?" Trans ":" NoTrans ")
    <<((diag==FFLAS::FflasUnit)?" Unit ":" NonUnit ")
    <<endl;

    cerr<<"A:"<<endl;
    cerr<<K<<" "<<K<<" M"<<endl;
    for (size_t i=0; i<K; ++i)
        for (size_t j=0; j<K; ++j)
            if ((*(Abis + i*lda + j)))
                cerr<<i+1<<" "<<j+1<<" "
                <<((int) *(Abis+i*lda+j) )
                <<endl;
    cerr<<"0 0 0"<<endl<<endl;

    cerr<<"B:"<<endl;
    cerr<<M<<" "<<N<<" M"<<endl;
    for (size_t i=0; i<M; ++i)
        for (size_t j=0; j<N; ++j)
            if ((*(Bbis + i*ldb + j)))
                cerr<<i+1<<" "<<j+1<<" "
                <<((int) *(Bbis+i*ldb+j) )
                <<endl;
    cerr<<"0 0 0"<<endl<<endl;

    FFLAS::fflas_delete( A);
    FFLAS::fflas_delete( Abis);
    FFLAS::fflas_delete( B);
    FFLAS::fflas_delete( Bbis);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
