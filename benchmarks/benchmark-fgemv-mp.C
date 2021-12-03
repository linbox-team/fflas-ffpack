/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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
 *.
 */

// declare that the call to openblas_set_numthread will be made here, hence don't do it
// everywhere in the call stack
#define __FFLASFFPACK_OPENBLAS_NT_ALREADY_SET 1

#if not defined(MG_DEFAULT)
#define MG_DEFAULT MG_ACTIVE
#endif
#if not defined(STD_RECINT_SIZE)
#define STD_RECINT_SIZE 8
#endif


#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <typeinfo>
#include <vector>
#include <string>
using namespace std;

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "givaro/modular-integer.h"
#include "givaro/givcaster.h"
#include "fflas-ffpack/paladin/parallel.h"
#ifdef BENCH_RECINT
#include "recint/recint.h"
#endif


template<typename T>
std::ostream& write_matrix(std::ostream& out, Givaro::Integer p, size_t m, size_t n, T* C, size_t ldc){

    size_t www(size_t((double(p.bitsize())*log(2.))/log(10.)));
    out<<"Matrix("<<m<<','<<n<<",[[";
    out.width(www+1);
    out<<std::right<<C[0];
    for (size_t j=1;j<n;++j){
        out<<',';
        out.width(www);
        out<<std::right<<C[j];
    }
    out<<']';
    for (size_t i=1;i<m;++i){
        out<<endl<<",[";
        out.width(www+1);
        out<<std::right<<C[i*ldc];
        for (size_t j=1;j<n;++j){
            out<<',';
            out.width(www);
            out<<std::right<<C[i*ldc+j];
        }
        out<<']';
    }
    return out<<"])";
}


static size_t iters = 3 ;
static Givaro::Integer q = -1 ;
static unsigned long b = 512 ;
static size_t m = 512 ;
static size_t k = 512 ;
static int nbw = -1 ;
static size_t seed= time(NULL);
static Argument as[] = {
    { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
    { 'b', "-b B", "Set the bitsize of the random characteristic.",         TYPE_INT , &b },
    { 'm', "-m M", "Set the dimension m of the matrix.",                    TYPE_INT , &m },
    { 'k', "-k K", "Set the dimension k of the matrix.",                    TYPE_INT , &k },
    { 'w', "-w N", "Set the number of winograd levels (-1 for random).",    TYPE_INT , &nbw },
    { 'i', "-i R", "Set number of repetitions.",                            TYPE_INT , &iters },
    { 's', "-s S", "Sets seed.",                            				TYPE_INT , &seed },
    END_OF_ARGUMENTS
};

template<typename Ints>
int tmain(){
    srand( (int)seed);
    srand48(seed);
    Givaro::Integer::seeding(seed);

    typedef Givaro::Modular<Ints> Field;
    Givaro::Integer p;
    FFLAS::Timer chrono, TimFreivalds;
    double time=0.;
    for (size_t loop=0;loop<iters;loop++){
        Givaro::Integer::random_exact_2exp(p, b);
        Givaro::IntPrimeDom IPD;
        IPD.nextprimein(p);
        Ints ip; Givaro::Caster<Ints,Givaro::Integer>(ip,p);
        Givaro::Caster<Givaro::Integer,Ints>(p,ip); // to check consistency

        Field F(ip);
        size_t lda,ldb,ldc;
        lda=k;
        ldb=1;
        ldc=1;

        typename Field::RandIter Rand(F,seed);
        typename Field::Element_ptr A,B,C;
        A= FFLAS::fflas_new(F,m,lda);
        B= FFLAS::fflas_new(F,k,ldb);
        C= FFLAS::fflas_new(F,m,ldc);

        // 		for (size_t i=0;i<m;++i)
        // 			for (size_t j=0;j<k;++j)
        // 				Rand.random(A[i*lda+j]);
        // 		for (size_t i=0;i<k;++i)
        // 			for (size_t j=0;j<n;++j)
        // 				Rand.random(B[i*ldb+j]);
        // 		for (size_t i=0;i<m;++i)
        // 			for (size_t j=0;j<n;++j)
        // 				Rand.random(C[i*ldc+j]);

        PAR_BLOCK { FFLAS::pfrand(F,Rand, m,k,A,m/size_t(MAX_THREADS)); }
        PAR_BLOCK { FFLAS::pfrand(F,Rand, k,1,B,k/MAX_THREADS); }
        PAR_BLOCK { FFLAS::pfzero(F, m,1,C,m/MAX_THREADS); }


        Ints alpha,beta;
        alpha=F.one;
        beta=F.zero;


        using  FFLAS::CuttingStrategy::Recursive;
        using  FFLAS::StrategyParameter::TwoDAdaptive;
        // RNS MUL_LA
        chrono.clear();chrono.start();
        {
            FFLAS::ParSeqHelper::Sequential seqH;
            FFLAS::fgemv(F,FFLAS::FflasNoTrans,m,k,alpha,A,lda,B,ldb,beta,C,ldc,seqH);
        }
        chrono.stop();
        time+=chrono.realtime();

        FFLAS::fflas_delete(A);
        FFLAS::fflas_delete(B);
        FFLAS::fflas_delete(C);

    }

    double Mflops=((2.*double(m)-1)/1000.*double(k)/1000.0) /time * double(iters);
    // 	Mflops*=p.bitsize()/16.;
    cout << "Time: "<< (time/double(iters))  <<" Gfops: "<<Mflops*1.0/1000.0
    << " (total:" << time <<") "
    <<typeid(Ints).name()
    <<" perword: "<< (Mflops*double(p.bitsize()))/64. ;
    FFLAS::writeCommandString(std::cout << " | " << p << " (" << p.bitsize()<<")|", as)  << std::endl;
    return 0;
}



int main(int argc, char** argv){

#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    FFLAS::parseArguments(argc,argv,as);

    int r1 = tmain<Givaro::Integer>();

#ifdef BENCH_RECINT
    r1 += tmain<RecInt::rint<STD_RECINT_SIZE>>();
#endif
    return r1;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
