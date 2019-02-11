/*
 * Copyright (C) FFLAS-FFPACK
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
//                        Test for the  lqup decomposition
//
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
using namespace std;
//#include "fflas-ffpack/field/modular-int.h"
//#include "fflas-ffpack/field/modular-positive.h"
#include "givaro/modular-balanced.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "givaro/givintprime.h"


using namespace FFPACK;

//typedef Givaro::Modular<double> Field;
typedef ModularBalanced<double> Field;
//typedef Givaro::Modular<float> Field;
//typedef ModularBalanced<float> Field;
//typedef Givaro::Modular<int> Field;
//typedef GivaroZpz<int32_t> Field;
//typedef GivaroGfq Field;

int main(int argc, char** argv){
    FFLAS::Timer tim;
    Givaro::IntPrimeDom IPD;
    uint64_t p;
    size_t M, N ;
    bool keepon = true;
    Givaro::Integer _p,tmp;
    cerr<<setprecision(10);
    size_t TMAX = 100;
    size_t PRIMESIZE = 23;

    if (argc > 1 )
        TMAX = atoi(argv[1]);
    if (argc > 2 )
        PRIMESIZE = atoi(argv[2]);

    FFLAS::FFLAS_TRANSPOSE ta;
    FFLAS::FFLAS_DIAG diag;
    size_t lda;

    Field::Element * A, *Abis, *X,* U, *L;
    size_t *P, *Q;
    while (keepon){
        srandom(_p);
        do{
            //		max = Integer::random(2);
            _p = random();//max % (2<<30);
            IPD.prevprime( tmp, (_p% (1<<PRIMESIZE)) );
            p =  tmp;

        }while( (p <= 2) );

        Field F( p);
        Field::RandIter RValue( F );

        do{
            M = (size_t)  random() % TMAX;
            N = (size_t)  random() % TMAX;
        } while ((M == 0) || (N == 0));
        lda = N;
        if (random()%2)
            diag = FFLAS::FflasUnit;
        else
            diag = FFLAS::FflasNonUnit;


        if (random()%2){
            ta = FFLAS::FflasTrans;
            L = FFLAS::fflas_new<Field::Element>(M*N);
            U = FFLAS::fflas_new<Field::Element>(N*N);
            P = FFLAS::fflas_new<size_t>(M);
            Q = FFLAS::fflas_new<size_t>(N);
            for (size_t i=0; i<M; ++i) P[i] = 0;
            for (size_t i=0; i<N; ++i) Q[i] = 0;
        }
        else{
            ta = FFLAS::FflasNoTrans;
            L = FFLAS::fflas_new<Field::Element>(M*M);
            U = FFLAS::fflas_new<Field::Element>(M*N);
            P = FFLAS::fflas_new<size_t>(N);
            Q = FFLAS::fflas_new<size_t>(M);
            for (size_t i=0; i<N; ++i) P[i] = 0;
            for (size_t i=0; i<M; ++i) Q[i] = 0;
        }

        size_t R=0;
        Field::Element * G = FFLAS::fflas_new<Field::Element>(M*M);
        Field::Element * H = FFLAS::fflas_new<Field::Element>(M*N);
        size_t t;
        do{
            t = (size_t) random() % 10;
        } while ((!t)||(t==1));
        for (size_t i=0; i<M; ++i)
            if (!(random() % t))
                for (size_t j=0; j < M; ++j)
                    RValue.random (*(G+i*M+j));
            else
                for (size_t j=0; j < M; ++j)
                    F.assign(*(G+i*M+j), F.zero);



        for (size_t j=0; j < N; ++j)
            if (!(random() % t))
                for (size_t i=0; i<M; ++i)
                    RValue.random (*(H+i*N+j));
            else
                for (size_t i=0; i<M; ++i)
                    F.assign(*(H+i*N+j), F.zero);

        // 		FFLAS::WriteMatrix (cerr<<"G = "<<endl,F,M,M,G,M);
        // 		FFLAS::WriteMatrix (cerr<<"H = "<<endl,F,M,N,H,N);
        A = FFLAS::fflas_new<Field::Element>(M*N);
        FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M,
                      N, M, F.one, G, M, H, N, F.zero, A, N);
        FFLAS::fflas_delete( G);
        FFLAS::fflas_delete( H);

        Abis = FFLAS::fflas_new<Field::Element>(M*N);
        for (size_t i=0; i<M*N; ++i)
            *(Abis+i) = *(A+i);

        X = FFLAS::fflas_new<Field::Element>(M*N);


        cout <<"p = "<<(size_t)p<<" M = "<<M
        <<" N = "<<N
        <<((diag==FFLAS::FflasUnit)?" Unit ":" Non Unit ")
        <<((ta==FFLAS::FflasNoTrans)?"LQUP ( A ) ":"LQUP ( A^T ) ")
        <<"....";


        tim.clear();
        tim.start();
        R = FFPACK::LUdivine (F, diag, ta, M, N, A, lda, P, Q);
        tim.stop();


        //FFLAS::WriteMatrix (cerr<<"Result = "<<endl,F,M,N,Abis,lda);

        if (ta == FFLAS::FflasNoTrans){

            for (size_t i=0; i<R; ++i){
                for (size_t j=0; j<i; ++j)
                    F.assign ( *(U + i*N + j), F.zero);
                for (size_t j=i+1; j<N; ++j)
                    F.assign (*(U + i*N + j), *(A+ i*N+j));
            }
            for (size_t i=R;i<M; ++i)
                for (size_t j=0; j<N; ++j)
                    F.assign(*(U+i*N+j), F.zero);
            for ( size_t i=0; i<M; ++i ){
                size_t j=0;
                for (; j< ((i<R)?i:R) ; ++j )
                    F.assign( *(L + i*M+j), *(A+i*N+j));
                for (; j<M; ++j )
                    F.assign( *(L+i*M+j), F.zero);
            }

            //FFLAS::WriteMatrix (cerr<<"L = "<<endl,F,M,M,U,M);
            //FFLAS::WriteMatrix (cerr<<"U = "<<endl,F,M,N,U,N);
            FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
                            M,0,(int) R, L, M, Q);
            for ( size_t  i=0; i<M; ++i )
                F.assign(*(L+i*(M+1)), F.one);

            if (diag == FFLAS::FflasNonUnit)
                for ( size_t  i=0; i<R; ++i )
                    F.assign (*(U+i*(N+1)), *(A+i*(lda+1)));

            else{
                for (size_t i=0; i<R; ++i ){
                    *(L+Q[i]*(M+1)) = *(A+Q[i]*lda+i);
                    F.assign (*(U+i*(N+1)), F.one);
                }
            }

            FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
                            M,0,(int) R, U, N, P);
            FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
                            N,0,(int) R, U, N, Q);
            FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M,N,M, 1.0, L,M, U,N, 0.0, X,N);
            //FFLAS::fflas_delete( A);
        } else {

            for (size_t i=0; i<R; ++i){
                for (size_t j=0; j<i; ++j)
                    F.assign ( *(L + i + j*N), F.zero);
                for (size_t j=i+1; j<M; ++j)
                    F.assign (*(L + i + j*N), *(A+ i+j*N));
            }

            for (size_t i=R;i<N; ++i)
                for (size_t j=0; j<M; ++j)
                    F.assign(*(L+i+j*N), F.zero);
            for ( size_t i=0; i<N; ++i ){
                size_t j=0;
                for (;  j< ((i<R)?i:R) ; ++j )
                    F.assign( *(U + i+j*N), *(A+i+j*N));
                for (; j<N; ++j )
                    F.assign( *(U+i+j*N), F.zero);
            }

            FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans,
                            N,0,(int) R, U, N, Q);
            for (size_t i=0; i<N; ++i)
                F.assign (*(U+i*(N+1)), F.one);
            if (diag == FFLAS::FflasNonUnit)
                for ( size_t i=0; i<R; ++i )
                    F.assign (*(L+i*(N+1)), *(A+i*(lda+1)));
            else{
                for ( size_t i=0; i<R; ++i ){
                    *(U+Q[i]*(N+1)) = *(A+Q[i]+i*N);
                    F.assign (*(L+i*(N+1)), F.one);
                }
            }
            // FFLAS::WriteMatrix (cerr<<"L = "<<endl,F,M,N,L,N);
            // 			FFLAS::WriteMatrix (cerr<<"U = "<<endl,F,N,N,U,N);

            FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
                            N,0,(int) R, L, N, P);
            FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
                            M,0,(int) R, L, N, Q);
            FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M,N,N, 1.0, L,N, U,N, 0.0, X,N);
        }
        for (size_t i=0; i<M; ++i)
            for (size_t j=0; j<N; ++j)
                if (!F.areEqual (*(Abis+i*N+j), *(X+i*N+j))){
                    cerr<<"error for i,j="<<i<<" "<<j<<" "<<*(Abis+i*N+j)<<" "<<*(X+i*N+j)<<endl;
                    keepon = false;
                }

        //FFLAS::WriteMatrix (cerr<<"X = "<<endl,F,m,n,X,n);
        //FFLAS::WriteMatrix (cerr<<"B = "<<endl,F,m,n,B,n);

        if (keepon){
            cout<<"R = "<<R
            <<" Passed "
            <<(double(M*M)/1000.0*(double(N)-double(M)/3.0)/tim.usertime()/1000.0)<<"Mfops"<<endl;
            FFLAS::fflas_delete( A);
            FFLAS::fflas_delete( L);
            FFLAS::fflas_delete( U);
            FFLAS::fflas_delete( Abis);
            FFLAS::fflas_delete( X);
            FFLAS::fflas_delete( P);
            FFLAS::fflas_delete( Q);
        }
        else{
            cerr<<"Abis = "<<endl;
            FFLAS::WriteMatrix (cerr, F, M, N, Abis, N);
            cerr<<"X = "<<endl;
            FFLAS::WriteMatrix (cerr, F, M, N, X, N);
        }
    }
    cout<<endl;
    cerr<<"FAILED with p = "<<(size_t)p<<" M = "<<M<<" N = "<<N
    <<" trans = "<<ta<<" diag = "<<diag<<endl;

    cerr<<"A:"<<endl;
    cerr<<M<<" "<<N<<" M"<<endl;
    for (size_t i=0; i<M; ++i)
        for (size_t j=0; j<N; ++j)
            if (*(Abis+i*lda+j))
                cerr<<i+1<<" "<<j+1<<" "<<((int) *(Abis+i*lda+j) )<<endl;
    cerr<<"0 0 0"<<endl<<endl;

    FFLAS::fflas_delete( A);
    FFLAS::fflas_delete( Abis);
    FFLAS::fflas_delete( L);
    FFLAS::fflas_delete( U);
    FFLAS::fflas_delete( X);
    FFLAS::fflas_delete( P);
    FFLAS::fflas_delete( Q);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
