/*
 * Copyright (C) 2015 the FFLAS-FFPACK group
 * Written by
 *
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
 *
 */


/*
 *******************************************************
 Parallel PLUQ quad recurisve with OpenMP
 *******************************************************

 g++ -D__FFLASFFPACK_HAVE_CBLAS -Wall -g -fopenmp -O3 -march=native -mavx -I/home/sultan/soft/fflas-ffpack/ -I/usr/local/soft/givaro-3.7.1/include  test-ppluq.C -L/home/pernet/Logiciels/ATLAS_1TH/lib -lcblas -latlas -L/usr/local/soft/givaro-3.7.1/lib -lgivaro -lm -lrt -Wl,-rpath -Wl,/usr/local/soft/givaro-3.7.1/lib  -o test-ppluq
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
//#include "omp.h"

#define __FFLASFFPACK_USE_OPENMP

#define __FFLAS__TRSM_READONLY
#define __PFTRSM_FOR_PLUQ
#include "fflas-ffpack/utils/fflas_io.h"
#include <givaro/modular-balanced.h>
//#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "sys/time.h"

//#define BASECASE_K 256
#include "fflas-ffpack/utils/test-utils.h"
//#include "fflas-ffpack/ffpack/parallel.h"

using namespace std;
using namespace FFLAS;
using namespace FFPACK;
#ifndef MODULO
#define MODULO 1
#endif

#if(MODULO==1)
typedef Givaro::Modular<double> Field;
#else
typedef Givaro::ZRing<double> Field;
#endif

#ifndef __FFLASFFPACK_DEBUG
#define __FFLASFFPACK_DEBUG 1
#endif

#ifndef SEQ
#define SEQ 1
#endif

void verification_PLUQ(const Field & F, typename Field::Element * B, typename Field::Element * A,
                       size_t * P, size_t * Q, size_t m, size_t n, size_t R)
{

    Field::Element * X = FFLAS::fflas_new<Field::Element>(m*n);
    Field::Element * L, *U;
    L = FFLAS::fflas_new<Field::Element>(m*R);
    U = FFLAS::fflas_new<Field::Element>(R*n);
    ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::Threads> H;

    PARFOR1D (i,m*R, H, { F.init(L[i]); F.assign(L[i], F.zero); } );
    PARFOR1D (i,m*R, H, { F.init(U[i]); F.assign(U[i], F.zero); } );
    PARFOR1D (i,m*n, H, { F.init(X[i]); F.assign(X[i], F.zero); } );


    PARFOR1D (i,R, H,
              for (size_t j=0; j<i; ++j)
              F.assign ( *(U + i*n + j), F.zero);
              for (size_t j=i; j<n; ++j)
              F.assign (*(U + i*n + j), *(A+ i*n+j));
             );

    PARFOR1D (j,R, H,
              for (size_t i=0; i<=j; ++i )
              F.assign( *(L+i*R+j), F.zero);
              F.assign(*(L+j*R+j), F.one);
              for (size_t i=j+1; i<m; i++)
              F.assign( *(L + i*R+j), *(A+i*n+j));
             );


    FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, R,0,m, L, R, P);

    FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, R,0,n, U, n, Q);
    FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,R,
                  1.0, L,R, U,n, 0.0, X,n);
    bool fail = false;
    PARFOR1D (i,m, H,
              for (size_t j=0; j<n; ++j)
              if (!F.areEqual (*(B+i*n+j), *(X+i*n+j))){
              std::stringstream errs;
              errs << " B["<<i<<","<<j<<"] = " << (*(B+i*n+j))
              << " X["<<i<<","<<j<<"] = " << (*(X+i*n+j))
              << std::endl;
//              std::cerr << errs;
              fail=true;
              }
             );

    if (fail)
        std::cerr<<"FAIL"<<std::endl;
    else
        std::cerr<<"PASS"<<std::endl;

    FFLAS::fflas_delete( U);
    FFLAS::fflas_delete( L);
    FFLAS::fflas_delete( X);
}

int main(int argc, char** argv)
{

    size_t p, n, m;
    int nbf;

    if (argc > 6){
        std::cerr<<"usage : PLUQ-rec-omp <p> <m> <n> <i> <file>"<<std::endl
        //		std::cerr<<"usage : PLUQ-rec-omp <m> <n> <p> <r> <i>"<<std::endl
        <<std::endl;
        exit(-1);
    }

    p = (argc>1 ? atoi( argv[1] ) : 1009);

    m = (argc>2 ? atoi( argv[2] ) : 1024);
    n = (argc>3 ? atoi( argv[3] ) : 1024);
    // r = atoi( argv[4] );
    nbf = (argc>4 ? atoi( argv[4] ) : 1);

    //	size_t lda = n;

    // random seed
    // ifstream f("/dev/urandom");
    // size_t seed1, seed2, seed3,seed4;
    // f.read(reinterpret_cast<char*>(&seed1), sizeof(seed1));
    // f.read(reinterpret_cast<char*>(&seed2), sizeof(seed2));
    // f.read(reinterpret_cast<char*>(&seed3), sizeof(seed3));
    // f.read(reinterpret_cast<char*>(&seed4), sizeof(seed4));

    //     seed1=10;seed2=12;
    //     seed3=13;seed4=14;

    enum FFLAS::FFLAS_DIAG diag = FFLAS::FflasNonUnit;
    size_t R;
    size_t RR;

    const Field F((double)p);
    // Field::RandIter G(F, seed1);

    // Field::Element * U = FFLAS::fflas_new<Field::Element>(n*n);

    ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::Threads> H;

    typename Field::Element* Acop;
    if (argc > 5) {
        FFLAS::ReadMatrix (argv[5],F,m,n,Acop,FFLAS::FflasAuto);
    } else {
        Field::RandIter G(F);
        Acop = FFLAS::fflas_new<Field::Element>(m*n);
        PARFOR1D(i,(size_t)m, H,
                 for (size_t j=0; j<(size_t)n; ++j)
                 G.random (*(Acop+i*n+j));
                );

    }

    // FFLAS::fflas_new<Field::Element>(n*m);
    Field::Element* A = FFLAS::fflas_new<Field::Element>(n*m);
    Field::Element* A2 = FFLAS::fflas_new<Field::Element>(n*m);
#if __FFLASFFPACK_DEBUG
    Field::Element* Adebug = FFLAS::fflas_new<Field::Element>(n*m);
    Field::Element* Adebug2 = FFLAS::fflas_new<Field::Element>(n*m);
#endif
    // std::vector<size_t> Index_P(r);

    // U = construct_U(F,G, n, r, Index_P, seed4, seed3);
    // A = construct_L(F,G, m, r, Index_P, seed2);
    // M_randgen(F, A, U, r, m, n);
    // size_t taille=m*n;
    // for(size_t i=0; i<taille;++i) U[i]=A[i];

    struct timespec t0, t1;// tt0, tt1;
    double delay, avrg;//, avrgg;
    double t_total=0;

    size_t maxP, maxQ;
    maxP = m;
    maxQ = n;

    size_t *P = FFLAS::fflas_new<size_t>(maxP);
    size_t *Q = FFLAS::fflas_new<size_t>(maxQ);

    size_t *P2 = FFLAS::fflas_new<size_t>(maxP);
    size_t *Q2 = FFLAS::fflas_new<size_t>(maxQ);

    PARFOR1D(i, (size_t)m, H,
             for (size_t j=0; j<(size_t)n; ++j) {
             *(A+i*n+j) = *(Acop+i*n+j) ;
             *(A2+i*n+j) = *(Acop+i*n+j) ;
#if __FFLASFFPACK_DEBUG
             *(Adebug+i*n+j) = *(Acop+i*n+j) ;
             *(Adebug2+i*n+j) = *(Acop+i*n+j) ;
#endif
             }
            );



    for ( int i=0;i<nbf+1;i++){
        for (size_t j=0;j<maxP;j++){
            P[j]=0;
            P2[j]=0;
        }
        for (size_t j=0;j<maxQ;j++){
            Q[j]=0;
            Q2[j]=0;
        }

        PARFOR1D(i, (size_t)m, H,
                 for (size_t j=0; j<(size_t)n; ++j){
                 *(A+i*n+j) = *(Acop+i*n+j) ;
                 *(A2+i*n+j) = *(Acop+i*n+j) ;
                 }
                );



        clock_gettime(CLOCK_REALTIME, &t0);
        PAR_BLOCK{
            H.set_numthreads(NUM_THREADS);
            R = PLUQ(F, diag, (size_t)m, (size_t)n, A, (size_t)n, P, Q, H);// Parallel PLUQ
        }
        clock_gettime(CLOCK_REALTIME, &t1);
        delay = (double)(t1.tv_sec-t0.tv_sec)+(double)(t1.tv_nsec-t0.tv_nsec)/1000000000;

        RR = pPLUQ(F, diag, (size_t)m, (size_t)n, A2, (size_t)n, P2, Q2);// Parallel pPLUQ

        if(i)
            t_total +=delay;

    }
    avrg = t_total/nbf;
    std::cerr << "MODULO: " << (MODULO?p:0) << std::endl;

    PAR_BLOCK{
        std::cerr<<"Parallel --- m: "<<m<<" , n: " << n << " , r: " <<R<<" "
        <<avrg<<" "<<(2.0*n*n*n)/(double(3.0*(1000000000)*avrg))<<" "
        //#ifdef  __FFLASFFPACK_USE_OPENMP
        <<NUM_THREADS<<endl;
        //#else
    }
    //<<endl;
    //#endi

    //	std::cout<<typeid(A).name()<<endl;
#if __FFLASFFPACK_DEBUG
    cout<<"check equality A == PLUQ ?"<<endl;
    verification_PLUQ(F,Adebug,A,P,Q,m,n,R);
    verification_PLUQ(F,Adebug2,A,P2,Q2,m,n,RR);
    FFLAS::fflas_delete( Adebug);
    FFLAS::fflas_delete( Adebug2);
#endif
#if(SEQ==1)
    struct timespec  tt0, tt1;
    double avrgg;
    //call sequential PLUQ
    size_t * PP = FFLAS::fflas_new<size_t>(maxP);
    size_t * QQ = FFLAS::fflas_new<size_t>(maxQ);
    for (size_t j=0;j<maxP;j++)
        PP[j]=0;
    for (size_t j=0;j<maxQ;j++)
        QQ[j]=0;

    clock_gettime(CLOCK_REALTIME, &tt0);
    size_t R2 = PLUQ(F, diag, m, n, Acop, n, PP, QQ);
    clock_gettime(CLOCK_REALTIME, &tt1);
    FFLAS::fflas_delete( Acop);
    avrgg = (double)(tt1.tv_sec-tt0.tv_sec)+(double)(tt1.tv_nsec-tt0.tv_nsec)/1000000000;
    //verification
    std::cerr<<"Sequential : "<<m<<" "<<R2<<" "
    <<avrgg<<" "<<(2.0*n*n*n)/(double(3.0*(1000000000)*avrgg))<<endl;
#endif

    FFLAS::fflas_delete( A);
    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
