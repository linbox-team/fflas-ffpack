/* Copyright (c) FFLAS-FFPACK
 * Written by Clement Pernet <clement.pernet@imag.fr>, from benchmark-pluq by Ziad Sultan
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

// declare that the call to openblas_set_numthread will be made here, hence don't do it
// everywhere in the call stack
#define __FFLASFFPACK_OPENBLAS_NT_ALREADY_SET

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas/fflas.h"
#include <givaro/modular-balanced.h>
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/ffpack/ffpack.h"
using namespace std;

typedef Givaro::ModularBalanced<double> Field;

// random generator function:                                    
ptrdiff_t myrandom (ptrdiff_t i) { return rand()%i;}

// pointer object to it:                                                                                               
ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;

typename Field::Element* construct_U(const Field& F, Field::RandIter& G, size_t n, size_t r, std::vector<size_t>& P, size_t commonseed, size_t seed, int nt)
{
    size_t lda = n;
    Field::Element *U=new Field::Element[n*n];

    FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::Threads>  H(nt);

    std::vector<size_t> E(r);
    PARFOR1D(i,r,H, E[i]=i; );


    srand48(commonseed);
    std::vector<size_t> Z(n);
    PARFOR1D(i,n,H, Z[i]=i;);

    P.resize(r);
    for(size_t i=0; i<r; ++i) {
        size_t index=lrand48() % Z.size();
        P[i] = Z[ index ];
        Z.erase(Z.begin()+index);
    }

    PARFOR1D(i,r,H,
             while( F.isZero( G.random(U[ P[i]*lda+P[i] ]) ) ) {}
             for(size_t j=P[i]+1;j<n;++j)
             G.random(U[ P[i]*lda+j]);
            );

    return U;
}

typename Field::Element* construct_L(const Field& F, Field::RandIter& G, size_t m, size_t r, const std::vector<size_t>& P, size_t seed, int nt)
{
    FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::Threads>  H(nt);

    size_t lda = m;
    size_t taille=m*m;
    Field::Element * L= new Field::Element[taille];
    PARFOR1D(i,taille,H, { F.init(L[i]); F.assign(L[i],F.zero); } );

    std::vector<size_t> E(r);
    PARFOR1D(i,r,H, E[i]=i;);

    srand48(seed);
    std::vector<size_t> Z(m);
    PARFOR1D(i,m,H, Z[i]=i; );

    std::vector<size_t> Q(r);
    for(size_t i=0; i<r; ++i) {
        size_t index=lrand48() % Z.size();
        Q[i] = Z[ index ];
        Z.erase(Z.begin()+index);
    }

    for(size_t i=0; i<r; ++i) {
        size_t index=lrand48() % E.size();
        size_t perm = E[ index ];

        E.erase(E.begin()+index);
        F.init(L[Q[perm]*lda+P[perm]]);
        F.assign(L[Q[perm]*lda+P[perm]], F.one);
        for(size_t j=Q[perm]+1;j<m;++j)
            G.random(L[j*lda+P[perm]]);
    }
    return L;
}


typename Field::Element* M_randgen(const Field& F, typename Field::Element* L,typename Field::Element* U, size_t r, size_t m, size_t n, int nt)
{
    Field::Element alpha, beta;
    F.init(alpha);
    F.assign(alpha, F.one);
    F.init(beta);
    F.assign(beta, F.zero);
    size_t lda = n;
    Field::Element * A = new Field::Element[m*n];

    // Computing produit L * U (ideally should use parallel ftrmm

    /*
       FFLAS::ftrmm(F,  FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, m,n,1.0, U, lda, L, lda);
       */


    FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads>  pWH(nt);
    PAR_BLOCK{
        FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                      m,n,r, alpha, L,r, U,
                      lda,beta,A,lda,pWH);
    }
    return L;
}

void verification_PLUQ(const Field & F, typename Field::Element * B, typename Field::Element * A,
                       size_t * P, size_t * Q, size_t m, size_t n, size_t R, int nt)
{

    FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::Threads>  H(nt);

    Field::Element * X = FFLAS::fflas_new<Field::Element>(m*n);
    Field::Element * L, *U;
    L = FFLAS::fflas_new<Field::Element>(m*R);
    U = FFLAS::fflas_new<Field::Element>(R*n);

    PARFOR1D (i,m*R,H, { F.init(L[i]); F.assign(L[i], 0.0); } );
    PARFOR1D (i,n*R,H, { F.init(U[i]); F.assign(U[i], 0.0); } );
    PARFOR1D (i,m*n,H, { F.init(X[i]); F.assign(X[i], 0.0); } );


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
#pragma omp task shared(F, P, L)
        FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, R,0,m, L, R, P);
#pragma omp task shared(F, Q, U)
        FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, R,0,n, U, n, Q);
#pragma omp taskwait
        FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads>  pWH (nt);
#pragma omp task shared(F, L, U, X)
        FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,R,
                      F.one, L,R, U,n, F.zero, X,n, pWH);

    }
    bool fail = false;
    //  PARFOR1D (size_t i=0; i<m; ++i)
    for(size_t i=0; i<m; ++i)
        for (size_t j=0; j<n; ++j)
            if (!F.areEqual (*(B+i*n+j), *(X+i*n+j))){
                std::cout << " Initial["<<i<<","<<j<<"] = " << (*(B+i*n+j))
                << " Result"<<i<<","<<j<<"] = " << (*(X+i*n+j))
                << std::endl;

                std::stringstream errs;
                errs << " B["<<i<<","<<j<<"] = " << (*(B+i*n+j))
                << " X["<<i<<","<<j<<"] = " << (*(X+i*n+j))
                << std::endl;
                std::cout << errs.str();
                fail=true;
                std::cout<<" end verification"<<std::endl;
                exit(1);
            }

    if (fail)
        std::cout<<"FAIL"<<std::endl;
    else
        std::cout<<"PASS"<<std::endl;

    FFLAS::fflas_delete( U);
    FFLAS::fflas_delete( L);
    FFLAS::fflas_delete( X);
}



int main(int argc, char** argv) {

#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t iter = 3 ;
    int q = 131071 ;
    Field F(q);
    int m = 2000 ;
    int n = 2000 ;
    size_t r = 2000 ;
    int v = 0;

    int nt = -1;

    int NBK = -1;
    int  a=1;
    bool transform = false;
    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INT , &q },
        { 'm', "-m M", "Set the row dimension of A.",      TYPE_INT , &m },
        { 'n', "-n N", "Set the col dimension of A.",      TYPE_INT , &n },
        { 'r', "-r R", "Set the rank of matrix A.",            TYPE_INT , &r },
        { 'i', "-i I", "Set number of repetitions.",            TYPE_INT , &iter },
        { 'v', "-v V", "Set 1 if need verification of result else 0.",            TYPE_INT , &v },
        { 't', "-t T", "whether to compute the transformation matrix.", TYPE_BOOL , &transform },
        { 'a', "-a A", "Algorithm for PLUQ decomposition: 0=LUdivine 1=PLUQ.", TYPE_INT , &a },
        { 'b', "-b B", "number of numa blocks per dimension for the numa placement", TYPE_INT , &NBK },
        { 'p', "-p P", "number of threads to use", TYPE_INT , &nt },
        END_OF_ARGUMENTS
    };
    //		{ 'p', "-p P", "0 for sequential, 1 for 2D iterative,
    //2 for 2D rec, 3 for 2D rec adaptive, 4 for 3D rc in-place, 5 for 3D rec, 6 for 3D rec adaptive.", TYPE_INT , &p },
    FFLAS::parseArguments(argc,argv,as);
    std::cerr<<"nt:="<<nt<<std::endl; if(nt==-1) PAR_BLOCK{ nt=NUM_THREADS; }  std::cerr<<"nt:="<<nt<<std::endl;
    FFPACK::FFPACK_LU_TAG LuTag = a?FFPACK::FfpackTileRecursive:FFPACK::FfpackSlabRecursive;
    if (NBK==-1) NBK = nt;
    typedef Field::Element Element;
    Element * A, * Acop;
    A = FFLAS::fflas_new(F,m,n,Alignment::CACHE_PAGESIZE);
    Acop = FFLAS::fflas_new(F,m,n,Alignment::CACHE_PAGESIZE);

    Field::Element * U = new Field::Element[n*n];
    // random seed                                                                                         
    ifstream f("/dev/urandom");
    size_t seed1, seed2, seed3,seed4;
    f.read(reinterpret_cast<char*>(&seed1), sizeof(seed1));
    f.read(reinterpret_cast<char*>(&seed2), sizeof(seed2));
    f.read(reinterpret_cast<char*>(&seed3), sizeof(seed3));
    f.read(reinterpret_cast<char*>(&seed4), sizeof(seed4));
    std::vector<size_t> Index_P(r);
    Field::RandIter GG(F, seed1);

    PAR_BLOCK{ FFLAS::pfrand(F,GG,m,n,A,m/NBK); }


    //       std::cout<<"Construct U"<<endl;
    U = construct_U(F,GG, n, r, Index_P, seed4, seed3,nt);
    //       std::cout<<"Construct L"<<endl;
    A = construct_L(F,GG, m, r, Index_P, seed2, nt);
    //       std::cout<<"randgen"<<endl;
    A = M_randgen(F, A, U, r, m, n, nt);
    size_t R=0;
    FFLAS::Timer chrono;
    double time=0.0;

    //       enum FFLAS::FFLAS_DIAG diag = FFLAS::FflasNonUnit;
    size_t maxP, maxQ;
    maxP = m;
    maxQ = n;

    size_t *P = FFLAS::fflas_new<size_t>(maxP);
    size_t *Q = FFLAS::fflas_new<size_t>(maxQ);

    FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::Threads>  H(nt);

    PARFOR1D(i,(size_t)m,H,
             for (size_t j=0; j<(size_t)n; ++j)
             Acop[i*n+j]= (*(A+i*n+j));
            );

    for (size_t i=0;i<=iter;++i){

        PARFOR1D(j,maxP,H, P[j]=0; );

        PARFOR1D(j,maxQ,H, Q[j]=0; );

        PARFOR1D(k,(size_t)m,H,
                 for (size_t j=0; j<(size_t)n; ++j)
                 *(A+k*n+j) = *(Acop+k*n+j) ;  
                );

    //Additional tests for templated helper
    FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::Threads> parH(nt);
    

        chrono.clear();

        if (i) chrono.start();
        // Added by AB 2014-12-15
        //#ifdef __FFLASFFPACK_USE_OPENMP
        PAR_BLOCK{
            r = FFPACK::RowEchelonForm(F,m,n,A,n,P,Q,transform,LuTag,parH);
        }
        if (i) {chrono.stop(); time+=chrono.realtime();}

    }

    // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14
#define CUBE(x) ((x)*(x)*(x))
    double gflop =  2.0/3.0*CUBE(double(r)/1000.0) +2*m/1000.0*n/1000.0*double(r)/1000.0  - double(r)/1000.0*double(r)/1000.0*(m+n)/1000;
    if (transform)
        gflop += CUBE(double(r)/1000.0)/3.0 + double(r)/1000.0*double(r)/1000.0*double(n-r)/1000.0;
    std::cout << "Time: " << time / double(iter)  << " Gfops: " << gflop / time * double(iter-1);
    FFLAS::writeCommandString(std::cout, as) << std::endl;

    //verification
    if(v)
        verification_PLUQ(F,Acop,A,P,Q,m,n,R,nt);
    FFLAS::fflas_delete( A);
    FFLAS::fflas_delete( Acop);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

