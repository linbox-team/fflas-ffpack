/* Copyright (c) FFLAS-FFPACK
 * Written by Ziad Sultan <ziad.sultan@imag.fr>
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
 */

//#include "goto-def.h"

//#define __FFLASFFPACK_USE_OPENMP
//#define __FFLASFFPACK_USE_TBB

//#define __FFLASFFPACK_USE_DATAFLOW
//#define  __FFLASFFPACK_FORCE_SEQ
//#define WINOPAR_KERNEL
//#define CLASSIC_SEQ
// #define PROFILE_PLUQ
// #define MONOTONIC_CYCLES
// #define MONOTONIC_MOREPIVOTS
// #define MONOTONIC_FEWPIVOTS

#ifdef MONOTONIC_CYCLES
#define MONOTONIC_APPLYP
#endif
#ifdef MONOTONIC_MOREPIVOTS
#define MONOTONIC_APPLYP
#endif
#ifdef MONOTONIC_FEWPIVOTS
#define MONOTONIC_APPLYP
#endif

// declare that the call to openblas_set_numthread will be made here, hence don't do it
// everywhere in the call stack
#define __FFLASFFPACK_OPENBLAS_NT_ALREADY_SET 1

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <givaro/modular.h>
#include <givaro/givranditer.h>
#include <iostream>

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/ffpack/ffpack.h"

#ifdef __FFLASFFPACK_USE_KAAPI
#include "libkomp.h"
#endif

using namespace std;

//typedef Givaro::ModularBalanced<double> Field;
typedef Givaro::Modular<double> Field;
//typedef Givaro::ModularBalanced<int64_t> Field;
//typedef Givaro::Modular<int64_t> Field;
//typedef Givaro::Modular<float> Field;
//typedef Givaro::ModularBalanced<float> Field;
//typedef Givaro::ZRing<double> Field;
//typedef Givaro::UnparametricZRing<double> Field;

void verification_PLUQ(const Field & F, typename Field::Element * B, typename Field::Element * A,
                       size_t * P, size_t * Q, size_t m, size_t n, size_t R)
{

    FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads> H;

    Field::Element_ptr X = FFLAS::fflas_new (F, m,n);
    Field::Element_ptr L, U;
    L = FFLAS::fflas_new(F, m,R);
    U = FFLAS::fflas_new(F, R,n);

    PARFOR1D (i, m*R,H, F.init(L[i], 0.0); );

    PARFOR1D (i,n*R,H, F.init(U[i], 0.0); );

    PARFOR1D (i,m*n,H, F.init(X[i], 0.0); );

    PARFOR1D (i,R,H,
              for (size_t j=0; j<i; ++j)
              F.assign ( *(U + i*n + j), F.zero);
              for (size_t j=i; j<n; ++j)
              F.assign (*(U + i*n + j), *(A+ i*n+j));
             );

    PARFOR1D (j,R,H,
              for (size_t i=0; i<=j; ++i )
              F.assign( *(L+i*R+j), F.zero);
              F.assign(*(L+j*R+j), F.one);
              for (size_t i=j+1; i<m; i++)
              F.assign( *(L + i*R+j), *(A+i*n+j));
             );

    PAR_BLOCK{
        SYNCH_GROUP(

                    TASK(MODE(CONSTREFERENCE(F,P,L)),
                         FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, R,0,m, L, R, P););
                    TASK(MODE(CONSTREFERENCE(F,Q,U)),
                         FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, R,0,n, U, n, Q););
                    WAIT;
                    typename FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads> pWH (MAX_THREADS);

                    TASK(MODE(CONSTREFERENCE(F,U,L,X)),
                         FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,R,
                                       F.one, L,R, U,n, F.zero, X,n, pWH););
                   );
    }
    bool fail = false;
    for(size_t i=0; i<m; ++i)
        for (size_t j=0; j<n; ++j)
            if (!F.areEqual (*(B+i*n+j), *(X+i*n+j))){
                std::cout << " Initial["<<i<<","<<j<<"] = " << (*(B+i*n+j))
                << " Result["<<i<<","<<j<<"] = " << (*(X+i*n+j))
                << std::endl;
                fail=true;
            }

    if (fail)
        std::cout<<"FAIL"<<std::endl;
    else
        std::cout<<"PASS"<<std::endl;

    FFLAS::fflas_delete( U);
    FFLAS::fflas_delete( L);
    FFLAS::fflas_delete( X);
}

void Rec_Initialize(Field &F, Field::Element * C, size_t m, size_t n, size_t ldc)
{
    if(std::min(m,n) <= ldc/NUM_THREADS){
        for(size_t i=0; i<m; i++)
            FFLAS::fzero(F, 1, n, C+i*n, n);
    }
    else{
        size_t M2 = m >> 1;
        size_t N2 = n >> 1;
        typename Field::Element * C2 = C + N2;
        typename Field::Element * C3 = C + M2*ldc;
        typename Field::Element * C4 = C3 + N2;

        SYNCH_GROUP(
                    TASK(MODE(CONSTREFERENCE(F)), Rec_Initialize(F,C,M2,N2, ldc););
                    TASK(MODE(CONSTREFERENCE(F)), Rec_Initialize(F,C2,M2,n-N2, ldc););
                    TASK(MODE(CONSTREFERENCE(F)), Rec_Initialize(F,C3,m-M2,N2, ldc););
                    TASK(MODE(CONSTREFERENCE(F)), Rec_Initialize(F,C4,m-M2,n-N2, ldc););
                   );
    }
}

int main(int argc, char** argv) {

#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t iter = 3 ;
    bool slab=false;
    Givaro::Integer q = 131071 ;
    int m = 2000 ;
    int n = 2000 ;
    int r = 2000 ;
    int v = 0;
    int t=MAX_THREADS;
    int NBK = -1;
    bool par=false;
    bool grp =true;
    Argument as[] = {
        { 's', "-s S", "Use the Slab recursive algorithm (LUdivine)instead of the tile recursive algorithm (PLUQ).", TYPE_BOOL , &slab },
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'm', "-m M", "Set the row dimension of A.",      TYPE_INT , &m },
        { 'n', "-n N", "Set the col dimension of A.",      TYPE_INT , &n },
        { 'r', "-r R", "Set the rank of matrix A.",            TYPE_INT , &r },
        { 'g', "-g yes/no", "Generic rank profile (yes) or random rank profile (no).", TYPE_BOOL , &grp },
        { 'i', "-i I", "Set number of repetitions.",            TYPE_INT , &iter },
        { 'v', "-v V", "Set 1 if need verification of result else 0.",            TYPE_INT , &v },
        { 't', "-t T", "number of virtual threads to drive the partition.", TYPE_INT , &t },
        { 'b', "-b B", "number of numa blocks per dimension for the numa placement", TYPE_INT , &NBK },
        { 'p', "-p P", "whether to run or not the parallel PLUQ", TYPE_BOOL , &par },
        END_OF_ARGUMENTS
    };
    FFLAS::parseArguments(argc,argv,as);
    Field F(q);
    if (q>F.maxCardinality())
        return 0;
    if (r > std::min(m,n)){
        std::cerr<<"Warning: rank can not be greater than min (m,n). It has been forced to min (m,n)"<<std::endl;
        r=std::min(m,n);
    }
    if (!par) { t=1;NBK=1;}
    if (NBK==-1) NBK = t;

    Field::Element_ptr A,  Acop;
    A = FFLAS::fflas_new(F,m,n);

    PAR_BLOCK{
        Rec_Initialize(F, A, m, n, n);
        //              FFLAS::pfzero(F,m,n,A,m/NBK);
        if (grp){
            size_t * cols = FFLAS::fflas_new<size_t>(n);
            size_t * rows = FFLAS::fflas_new<size_t>(m);
            for (int i=0; i<n; ++i)
                cols[i] = i;
            for (int i=0; i<m; ++i)
                rows[i] = i;
            FFPACK::RandomMatrixWithRankandRPM (F, m, n ,r, A, n, rows, cols);
            FFLAS::fflas_delete(cols);
            FFLAS::fflas_delete(rows);
        } else
            FFPACK::RandomMatrixWithRankandRandomRPM (F, m, n ,r, A, n);
    }
    size_t R;
    FFLAS::Timer chrono;
    double *time=new double[iter];

    enum FFLAS::FFLAS_DIAG diag = FFLAS::FflasNonUnit;
    size_t maxP, maxQ;
    maxP = m;
    maxQ = n;

    size_t *P = FFLAS::fflas_new<size_t>(maxP);
    size_t *Q = FFLAS::fflas_new<size_t>(maxQ);

    Acop = FFLAS::fflas_new(F,m,n);
    FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,
                                  FFLAS::StrategyParameter::Threads> parH;
    PARFOR1D(i,(size_t)m,parH,
             FFLAS::fassign(F, n, A + i*n, 1, Acop + i*n, 1);
             // for (size_t j=0; j<(size_t)n; ++j)
             //     Acop[i*n+j]= A[i*n+j];
            );
    for (size_t i=0;i<=iter;++i){

        PARFOR1D(j,maxP,parH, P[j]=0; );
        PARFOR1D(j,maxQ,parH, Q[j]=0; );
        PARFOR1D(k,(size_t)m,parH,
                 FFLAS::fassign(F, n, Acop + k*n, 1, A + k*n, 1);
                 // for (size_t j=0; j<(size_t)n; ++j)
                 //     F.assign( A[k*n+j] , Acop[k*n+j]) ;
                );
        chrono.clear();

        if (i) chrono.start();
        if (par){
/*
            PAR_BLOCK{
                parH.set_numthreads(t);
                R = FFPACK::PLUQ(F, diag, m, n, A, n, P, Q, parH);
            }
*/
            R = FFPACK::pPLUQ(F, diag, m, n, A, n, P, Q);
        }
        else{
            if (slab)
                R = FFPACK::LUdivine (F, diag, FFLAS::FflasNoTrans, m, n, A, n, P, Q);
            else
                R = FFPACK::PLUQ(F, diag, m, n, A, n, P, Q);
        }
        if (i) {chrono.stop(); time[i-1]=chrono.realtime();}

    }
    std::sort(time, time+iter);
    double mediantime = time[iter/2];
    delete[] time;
    // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14
#define CUBE(x) ((x)*(x)*(x))
    double gflop =  2.0/3.0*CUBE(double(r)/1000.0) +2*m/1000.0*n/1000.0*double(r)/1000.0  - double(r)/1000.0*double(r)/1000.0*(m+n)/1000;
    std::cout << "Time: " << mediantime
    << " Gfops: " << gflop / mediantime;
    FFLAS::writeCommandString(std::cout, as) << std::endl;

    //verification
    if(v)
        verification_PLUQ(F,Acop,A,P,Q,m,n,R);

    FFLAS::fflas_delete (P);
    FFLAS::fflas_delete (Q);
    FFLAS::fflas_delete (A);
    FFLAS::fflas_delete (Acop);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
