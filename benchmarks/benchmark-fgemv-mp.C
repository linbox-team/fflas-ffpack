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


static size_t iters = 3 ;
static Givaro::Integer q = -1 ;
static unsigned long b = 512 ;
static size_t m = 512 ;
static size_t k = 512 ;

static size_t seed= time(NULL);
static int par = 0;
int t = 1;

size_t GrainSize = 64;
static Argument as[] = {
    { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
    { 'b', "-b B", "Set the bitsize of the random characteristic.",         TYPE_INT , &b },
    { 'm', "-m M", "Set the dimension m of the matrix.",                    TYPE_INT , &m },
    { 'k', "-k K", "Set the dimension k of the matrix.",                    TYPE_INT , &k },

    { 'i', "-i R", "Set number of repetitions.",                            TYPE_INT , &iters },
    { 's', "-s S", "Sets seed.",                            				TYPE_INT , &seed },
    { 'p', "-p P", "0 for sequential, 1 for <Recursive,Thread>, 2 for <Row,Thread>, 3 for <Row,Grain>.",
                                                                                    TYPE_INT , &par },
    { 'g', "-g G", "Sets GrainSize.",                            			        TYPE_INT , &GrainSize },
    { 't', "-t T", "number of virtual threads to drive the partition.",             TYPE_INT , &t },
    END_OF_ARGUMENTS
};

template <class Field, class Matrix, class Vector>
bool check_result(Field& F, size_t m, size_t lda, Matrix& A, Vector& X, size_t incX, Vector& Y, size_t incY){
  //Naive result checking by comparing result from pfgemv against the one from fgemv
  typename Field::Element_ptr Y2 = FFLAS::fflas_new(F,m,1);
  FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y2,  incY);

  for(size_t j=0; j<m; ++j){
    if(!F.areEqual(Y2[j],Y[j])){
      FFLAS::fflas_delete(Y2);
      return false;
    }
  }
  FFLAS::fflas_delete(Y2);
  return true;
}

template <class Field, class arg>
void benchmark_disp(Field& F, double& time, size_t iters, int p,  size_t m, size_t k, arg& as){

  std::cout << "Time: " << time / double(iters)
     << " Gflops: " << (2.*double(m)/1000.*double(k)/1000.0/1000.0) / time * double(iters);
        FFLAS::writeCommandString(std::cout, as) << std::endl;
}

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
        
        //Fill with random value datum
        PAR_BLOCK { FFLAS::pfrand(F,Rand, m,k,A,m/size_t(MAX_THREADS)); }
        PAR_BLOCK { FFLAS::pfrand(F,Rand, k,1,B,k/MAX_THREADS); }
        PAR_BLOCK { FFLAS::pfzero(F, m,1,C,m/MAX_THREADS); }


        Ints alpha,beta;
        alpha=F.one;
        beta=F.zero;


        using  FFLAS::CuttingStrategy::Recursive;
        using  FFLAS::StrategyParameter::TwoDAdaptive;
        // RNS MUL_LA
        chrono.clear();

        //@TODO: Still need to use PAR_BLOCK to label the parallel region, impl as pDet to wrap this into one function
        PAR_BLOCK {
            if (par){
              typedef FFLAS::CuttingStrategy::Row row;
              typedef FFLAS::CuttingStrategy::Recursive rec;
              typedef FFLAS::StrategyParameter::Threads threads;
              typedef FFLAS::StrategyParameter::Grain grain;

              if (loop) { chrono.start(); }

              switch (par){

              case 1:{
	            FFLAS::ParSeqHelper::Parallel<rec, threads>  H(t);
	            FFLAS::fgemv(F,FFLAS::FflasNoTrans,m,k,alpha,A,lda,B,ldb,beta,C,ldc, H);
	            break;
	            }
              case 2:{
	            FFLAS::ParSeqHelper::Parallel<row, threads>  H(t);
	            FFLAS::fgemv(F,FFLAS::FflasNoTrans,m,k,alpha,A,lda,B,ldb,beta,C,ldc, H);
	            break;
              }
              case 3:{
	            FFLAS::ParSeqHelper::Parallel<row, grain>  H(GrainSize);
	            FFLAS::fgemv(F,FFLAS::FflasNoTrans,m,k,alpha,A,lda,B,ldb,beta,C,ldc, H);
	            break;
              }
              default:{
	            FFLAS::ParSeqHelper::Sequential  H;
	            FFLAS::fgemv(F,FFLAS::FflasNoTrans,m,k,alpha,A,lda,B,ldb,beta,C,ldc, H);
	            break;
              }
              }

              if (loop) {chrono.stop(); time+=chrono.realtime();}
            }else{
              if (loop) chrono.start();
              FFLAS::fgemv(F,FFLAS::FflasNoTrans,m,k,alpha,A,lda,B,ldb,beta,C,ldc);
              if (loop) {chrono.stop(); time+=chrono.realtime();}
            }

                time+=chrono.realtime();

                FFLAS::fflas_delete(A);
                FFLAS::fflas_delete(B);
                FFLAS::fflas_delete(C);

            }
/*
            if(!check_result(F, m, lda,  A,  B, ldb,  C, ldc)){
              std::cerr<<"Computation failed with wrong result"<<std::endl;
              break;
            }
*/
        }

    Field F;
    benchmark_disp(F, time, iters, p, m, k, as);
    return 0;
}

int main(int argc, char** argv){
    //Set the defaut value to the number of all available threads
    //PAR_BLOCK { t = NUM_THREADS; };

    FFLAS::parseArguments(argc,argv,as);

    int r1 = tmain<Givaro::Integer>();

#ifdef BENCH_RECINT
    r1 += tmain<RecInt::rint<STD_RECINT_SIZE>>();
#endif
    return r1;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
