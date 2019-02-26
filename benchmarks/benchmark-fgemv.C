/* Copyright (c) FFLAS-FFPACK
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

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular-balanced.h>

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/args-parser.h"

#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/test-utils.h"

#include "fflas-ffpack/utils/timer.h"
#include "givaro/modular-integer.h"
#include "givaro/givcaster.h"

using namespace FFPACK;

using namespace std;
using namespace FFLAS;

template <typename Field>
struct need_field_characteristic { static constexpr bool value = false; };
template <typename Field>
struct need_field_characteristic<Givaro::Modular<Field>>{ static constexpr bool value = true; };
template <typename Field>
struct need_field_characteristic<Givaro::ModularBalanced<Field>>{ static constexpr bool value = true; };


template <class Field, class RandIter, class Matrix, class Vector>
void fill_value(Field& F, RandIter& Rand, 
    Matrix& A, Vector& X, Vector& Y, 
    size_t m, size_t k, size_t incX, size_t incY, size_t lda, int NBK){

        // TODO: replace by a 1D pfrand
        PAR_BLOCK {

            SYNCH_GROUP(
                        FORBLOCK1D(iter, m, SPLITTER(NBK, 
                        CuttingStrategy::Row, StrategyParameter::Threads
                        ),
                                   TASK(MODE(CONSTREFERENCE(F,Rand)),
                                        {
                                        frand(F, Rand,
                                              iter.end()-iter.begin(),
                                              k,
                                              A+iter.begin()*lda,
                                              lda);
                                        }
                                       );
                                  );
                       );

            //FFLAS::pfrand(F,Rand, m,k,A,m/NBK);
        }

    FFLAS::frand(F,Rand, k,1,X,incX);
    FFLAS::fzero(F, m,1,Y,incY);

}


template <class Field, class Matrix, class Vector>
void genData(Field& F, 
    Matrix& A, Vector& X, Vector& Y, 
    size_t m, size_t k, size_t incX, size_t incY, size_t lda, int NBK,
    int bitsize, uint64_t seed){

    typename Field::RandIter Rand(F,bitsize,seed);
    fill_value(F, Rand, A, X, Y, m, k, incX, incY, lda, NBK);
}


#if 0

template <class Field, class Matrix, class Vector, class T>
bool benchmark_with_timer(Field& F, int p,
    Matrix& A, Vector& X, Vector& Y, 
    size_t m, size_t k, size_t incX, size_t incY, size_t lda, size_t iters, int t,
    T& chrono, double& time){

    bool pass = true;
    for (size_t i=0;i<=iters;++i){

        chrono.clear();

        if (p){

            typedef CuttingStrategy::Row row;
            typedef CuttingStrategy::Recursive rec;
            typedef StrategyParameter::Threads threads;
            typedef StrategyParameter::Grain grain;


            PAR_BLOCK{
                if (i) { chrono.start(); }

                switch (p){
                case 1:{
                           ParSeqHelper::Parallel<rec, threads>  H(t);
                           FFLAS::pfgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y,  incY, H);
                           break;}
                case 2:{
                           ParSeqHelper::Parallel<row, threads>  H(t);
                           FFLAS::pfgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y,  incY, H);
                           break;
                       }
                case 3:{
                           size_t BS = 64;
                           ParSeqHelper::Parallel<row, grain>  H(BS);
                           FFLAS::pfgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y,  incY, H);
                           break;
                       }
                default:{
                            FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y,  incY);
                            break;
                        }
                }
            }//PAR_BLOCK
            if (i) {chrono.stop(); time+=chrono.realtime();}
        }else{
            if (i) chrono.start();
            FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y,  incY);
            if (i) {chrono.stop(); time+=chrono.realtime();}
        }
        //Naive result checking by comparing result from pfgemv against the one from fgemv
        typename Field::Element_ptr Y2 = FFLAS::fflas_new(F,m,1);

        FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y2,  incY);



        for(size_t j=0; j<m; ++j){
            pass &= F.areEqual(Y2[j],Y[j]);
            if(!F.areEqual(Y2[j],Y[j])) std::cerr<<std::setprecision(500)<<" : Y2["<<j<<"]<"<<Y2[j]<<">=!=Y["<<j<<"]<"<<Y[j]<<">"<<std::endl;

        }
    }
    return pass;
}

#else

enum Index { One, Two, Three };

template <template <Index> class Switcher, class Field, class Matrix, class Vector> 
void test_case_switch(int k, Field& F, int t, size_t m, size_t lda, Matrix& A, Vector& X, size_t incX, Vector& Y, size_t incY)
{
            typedef CuttingStrategy::Row row;
            typedef CuttingStrategy::Recursive rec;
            typedef StrategyParameter::Threads threads;
            typedef StrategyParameter::Grain grain;
            size_t BS = 64;
            ParSeqHelper::Parallel<rec, threads>  H(t); 
            ParSeqHelper::Parallel<row, threads>  H2(t);
            ParSeqHelper::Parallel<row, grain>  H3(BS);
    switch (k)
    {
        case 1: Switcher<One>::template test(F, m, lda, A, X, incX, Y,  incY, H); break;
        case 2: Switcher<Two>::template test(F, m, lda, A, X, incX, Y,  incY, H2); break;
        case 3:  Switcher<Three>::template test(F, m, lda, A, X, incX, Y,  incY, H3); break;
        default: FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y,  incY); break;
    }
}

template <Index I> struct CaseSwitch
{
    template <class Field, class Matrix, class Vector, class Helper>
    static void test(Field& F, size_t m, size_t lda, Matrix& A, Vector& X, size_t incX, Vector& Y, size_t incY, Helper& H){
        FFLAS::pfgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y,  incY, H);
    };
};


template <class Field, class Matrix, class Vector, class T>
bool benchmark_with_timer(Field& F, int p,
    Matrix& A, Vector& X, Vector& Y, 
    size_t m, size_t k, size_t incX, size_t incY, size_t lda, size_t iters, int t,
    T& chrono, double& time){
    
    bool pass = true;
    for (size_t i=0;i<=iters;++i){

        chrono.clear();

        if (p){

            PAR_BLOCK{
                if (i) { chrono.start(); }

                test_case_switch<CaseSwitch>(p, F, t,  m, lda, A, X, incX, Y,  incY);

            }//PAR_BLOCK
            if (i) {chrono.stop(); time+=chrono.realtime();}
        }else{
            if (i) chrono.start();
            FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y,  incY);
            if (i) {chrono.stop(); time+=chrono.realtime();}
        }
        //Naive result checking by comparing result from pfgemv against the one from fgemv
        typename Field::Element_ptr Y2 = FFLAS::fflas_new(F,m,1);
        FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y2,  incY);

        for(size_t j=0; j<m; ++j){
            pass &= F.areEqual(Y2[j],Y[j]);
            if(!F.areEqual(Y2[j],Y[j])) std::cerr<<std::setprecision(500)<<" : Y2["<<j<<"]<"<<Y2[j]<<">=!=Y["<<j<<"]<"<<Y[j]<<">"<<std::endl;

        }
    }
    return pass;
}

#endif


template <class Field, class arg>
void benchmark_in_ZZ(Field& F, int p,  size_t m, size_t k,
    int NBK, int bitsize, uint64_t seed, size_t iters, int t, 
    arg& as
    ){

    Timer chrono, TimFreivalds;
    double time=0.0;
    size_t lda,incX,incY;
    lda=k;
    incX=1;
    incY=1;
    typename Field::Element_ptr A,X,Y; 

    A = FFLAS::fflas_new(F,m,lda);
    X = FFLAS::fflas_new(F,k,incX);
    Y = FFLAS::fflas_new(F,m,incY);

    genData(F, A, X, Y, m, k, incX, incY, lda, NBK, bitsize, seed);
/*
            FFLAS::WriteMatrix (std::cout<<"A:=",F,m,k,A,lda)<<';'<<std::endl;
            FFLAS::WriteMatrix(std::cout<<"X:=",F,k,1,X,1)<<';'<<std::endl;
            FFLAS::WriteMatrix(std::cout<<"Y:=",F,k,1,Y,1)<<';'<<std::endl;
*/
    bool pass=benchmark_with_timer( F, p, A, X, Y, m, k, incX, incY, lda, iters, t, chrono, time);
    
    if(pass){
        std::cout << "Time: " << time / double(iters)
        << " Gflops: " << (2.*double(m)/1000.*double(k)/1000.0/1000.0) / time * double(iters);
        writeCommandString(std::cout, as) << std::endl;
    }else{
            std::cout<<"FAILED for "<<typeid(Field).name()<<std::endl;
            std::cout << "p:=" << p << ';'<<std::endl;

    }
    
    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(X);
    FFLAS::fflas_delete(Y);    
    
}


template <class Field,  class arg >
void benchmark_with_field(int p,  size_t m, size_t k,
    int NBK, int bitsize, uint64_t seed, size_t iters, int t, 
    arg& as
    ){
    
    Field F;
    //static assert to raise compile time error for Non ZRing without providing a characteristic
    static_assert(!need_field_characteristic<Field>::value, "A field characteristic should be provided for Non ZRing data type !");

    benchmark_in_ZZ(F, p,  m, k, NBK, bitsize, seed, iters, t, as);
    
}


template <class Field, class arg>
void benchmark_with_field(const Givaro::Integer& q, int p,  size_t m, size_t k,
    int NBK, int bitsize, uint64_t seed, size_t iters, int t,
    arg& as
    ){
    Field  F(q);
    benchmark_in_ZZ(F, p,  m, k, NBK, bitsize, seed, iters, t, as);    
}


int main(int argc, char** argv) {

    int p=0;    
    
    size_t iters = 10;
    Givaro::Integer q = 131071; 
    size_t m = 800;
    size_t k = 800;
    //static size_t n = 512 ;
    uint64_t seed= getSeed();
    int t=NUM_THREADS;
    int NBK = -1;
    int b=0;
    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of input.",         TYPE_INT , &b },
        { 'p', "-p P", "0 for sequential, 1 for <Recursive,Thread>, 2 for <Row,Thread>, 3 for <Row, Grain>.", TYPE_INT , &p },
        { 'm', "-m M", "Set the dimension m of the matrix.",                    TYPE_INT , &m },
        { 'k', "-k K", "Set the dimension k of the matrix.",                    TYPE_INT , &k },
        { 't', "-t T", "number of virtual threads to drive the partition.", TYPE_INT , &t },
        { 'N', "-n N", "number of numa blocks per dimension for the numa placement", TYPE_INT , &NBK },
        { 'i', "-i R", "Set number of repetitions.",                            TYPE_INT , &iters },
        { 's', "-s S", "Sets seed.",                            				TYPE_INT , &seed },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);

    if (NBK==-1) NBK = t;
    //@Fixme runtime error which appears to be due to precision
    //benchmark_with_field<Givaro::ZRing<float>>( p,  m, k, NBK, b, seed, iters, t, as);
    //benchmark_with_field<Givaro::ZRing<double>>( p,  m, k, NBK, b, seed, iters, t, as);
    
    benchmark_with_field<Givaro::ZRing<int32_t>>( p,  m, k, NBK, b, seed, iters, t, as);
    benchmark_with_field<Givaro::ZRing<Givaro::Integer>>( p,  m, k, NBK, b, seed, iters, t, as);

    benchmark_with_field<Givaro::Modular<float>>(q, p,  m, k, NBK, b, seed, iters, t, as);
    benchmark_with_field<Givaro::Modular<double>>(q, p,  m, k, NBK, b, seed, iters, t, as);
    benchmark_with_field<Givaro::Modular<int32_t>>(q, p,  m, k, NBK, b, seed, iters, t, as);
    benchmark_with_field<Givaro::Modular<Givaro::Integer>>(q, p,  m, k, NBK, b, seed, iters, t, as);

    benchmark_with_field<Givaro::ModularBalanced<float>>(q, p,  m, k, NBK, b, seed, iters, t, as);
    benchmark_with_field<Givaro::ModularBalanced<double>>(q, p,  m, k, NBK, b, seed, iters, t, as);
    benchmark_with_field<Givaro::ModularBalanced<int32_t>>(q, p,  m, k, NBK, b, seed, iters, t, as);

    //@Fixme compile error : ‘Givaro::ModularBalanced<Givaro::Integer> F’ has incomplete type
    //benchmark_with_field<Givaro::ModularBalanced<Givaro::Integer>>(q, p,  m, k, NBK, b, seed, iters, t, as);
    
    
    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
