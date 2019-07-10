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

// declare that the call to openblas_set_numthread will be made here, hence don't do it
// everywhere in the call stack
#define __FFLASFFPACK_OPENBLAS_NT_ALREADY_SET 1

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

template <typename Field>
struct compatible_data_type { static constexpr bool value = true; };
template <>
struct compatible_data_type<Givaro::ZRing<float>>{ static constexpr bool value = false; };
template <>
struct compatible_data_type<Givaro::ZRing<double>>{ static constexpr bool value = false; };


template <class Field, class RandIter, class Matrix, class Vector>
void fill_value(Field& F, RandIter& Rand,
		Matrix& A, Vector& X, Vector& Y,
		size_t m, size_t k, size_t incX, size_t incY, size_t lda, int NBK){
  // TODO: replace by a 1D pfrand
  SYNCH_GROUP(
	      FORBLOCK1D(iter, m, SPLITTER(NBK, CuttingStrategy::Row, StrategyParameter::Threads),
			 TASK(MODE(CONSTREFERENCE(F,Rand,A)),
			      {
                    frand(F, Rand, iter.end()-iter.begin(), k, A+iter.begin()*lda, lda);
			      }
			      );
			 );
	      );
  //FFLAS::pfrand(F,Rand, m,k,A,m/NBK);
  FFLAS::frand(F,Rand, k,1,X,incX);
  FFLAS::fzero(F, m,1,Y,incY);
}

template <class Field, class Matrix, class Vector>
void genData(Field& F,
	     Matrix& A, Vector& X, Vector& Y,
	     size_t m, size_t k, size_t incX, size_t incY, size_t lda, int NBK,
	     uint64_t bitsize, uint64_t seed){
    Givaro::Integer samplesize(1); samplesize <<= bitsize;
    typename Field::RandIter Rand(F,seed,samplesize);
    fill_value(F, Rand, A, X, Y, m, k, incX, incY, lda, NBK);
}

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


template <class Field, class Matrix, class Vector>
bool benchmark_with_timer(Field& F, int p, Matrix& A, Vector& X, Vector& Y, size_t m, size_t k, size_t incX,
			  size_t incY, size_t lda, size_t iters, int t, double& time, size_t GrainSize){
  Timer chrono;
  bool pass = true;
  for (size_t i=0;i<=iters;++i){

    chrono.clear();

    if (p){

      typedef CuttingStrategy::Row row;
      typedef CuttingStrategy::Recursive rec;
      typedef StrategyParameter::Threads threads;
      typedef StrategyParameter::Grain grain;

      if (i) { chrono.start(); }

      switch (p){
      case 1:{
	ParSeqHelper::Parallel<rec, threads>  H(t);
	FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y,  incY, H);
	break;}
      case 2:{
	ParSeqHelper::Parallel<row, threads>  H(t);
	FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y,  incY, H);
	break;
      }
      case 3:{
	ParSeqHelper::Parallel<row, grain>  H(GrainSize);
	FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y,  incY, H);
	break;
      }
      default:{
	FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y,  incY);
	break;
      }
      }

      if (i) {chrono.stop(); time+=chrono.realtime();}
    }else{
      if (i) chrono.start();
      FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, lda, F.one, A, lda, X, incX, F.zero, Y,  incY);
      if (i) {chrono.stop(); time+=chrono.realtime();}
    }

    if(!check_result(F, m, lda,  A,  X, incX,  Y, incY)){
      pass = false;
      break;
    }

  }
  return pass;
}

template <class Field, class arg>
void benchmark_disp(Field& F, bool pass, double& time, size_t iters, int p,  size_t m, size_t k, arg& as){
  if(pass){
    std::cout << "Time: " << time / double(iters)
	      << " Gflops: " << (2.*double(m)/1000.*double(k)/1000.0/1000.0) / time * double(iters);
    writeCommandString(std::cout, as) << std::endl;
  }else{
    std::cout<<"FAILED for "<<typeid(Field).name()<<std::endl;
    std::cout << "p:=" << p << ';'<<std::endl;
  }
}


template <class Field, class arg>
void benchmark_in_Field(Field& F, int p,  size_t m, size_t k, int NBK, uint64_t bitsize, uint64_t seed, size_t iters,
			int t, arg& as, size_t GrainSize){
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

    bool pass=benchmark_with_timer( F, p, A, X, Y, m, k, incX, incY, lda, iters, t, time, GrainSize);

    benchmark_disp(F, pass, time, iters, p, m, k, as);

    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(X);
    FFLAS::fflas_delete(Y);

}

template <class Field,  class arg >
void benchmark_with_field(int p,  size_t m, size_t k, int NBK, uint64_t bitsize, uint64_t seed, size_t iters,
			  int t, arg& as, size_t GrainSize){
  Field F;
  //static assert to raise compile time error for Non ZRing without providing a characteristic
  static_assert(!need_field_characteristic<Field>::value,
		"A field characteristic should be provided for Non ZRing data type !");
  //static assert to raise compile time error for ZRing with either float or double that could lead to inconsistent result
  static_assert(compatible_data_type<Field>::value,
		"The provided data type for ZRing is not compatible for the desired operation and could lead to inconsistent result !");

  benchmark_in_Field(F, p,  m, k, NBK, bitsize, seed, iters, t, as, GrainSize);

}

template <class Field, class arg>
void benchmark_with_field(const Givaro::Integer& q, int p,  size_t m, size_t k,
			  int NBK, uint64_t bitsize, uint64_t seed, size_t iters, int t,
			  arg& as, size_t GrainSize){
    Field  F(q);
    benchmark_in_Field(F, p,  m, k, NBK, bitsize, seed, iters, t, as, GrainSize);
}

int main(int argc, char** argv) {

#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    int p=0;

  size_t iters = 3;
  Givaro::Integer q = 131071;
  size_t m = 4000;
  size_t k = 4000;

  uint64_t seed = getSeed();
  int t;
  PAR_BLOCK { t = NUM_THREADS; }
  int NBK = -1;
  uint64_t b=0;
  size_t GrainSize = 64;

  Argument as[] = {
    { 'q', "-q Q", "Set the field characteristic (-1 for random).",                 TYPE_INTEGER , &q },
    { 'b', "-b B", "Set the bitsize of input.",                                     TYPE_INT , &b },
    { 'p', "-p P", "0 for sequential, 1 for <Recursive,Thread>, 2 for <Row,Thread>, 3 for <Row,Grain>.",
                                                                                    TYPE_INT , &p },
    { 'm', "-m M", "Set the dimension m of the matrix.",                            TYPE_INT , &m },
    { 'k', "-k K", "Set the dimension k of the matrix.",                            TYPE_INT , &k },
    { 't', "-t T", "number of virtual threads to drive the partition.",             TYPE_INT , &t },
    { 'N', "-n N", "number of numa blocks per dimension for the numa placement",    TYPE_INT , &NBK },
    { 'i', "-i R", "Set number of repetitions.",                                    TYPE_INT , &iters },
    { 's', "-s S", "Sets seed.",				        TYPE_INT , &seed },
    { 'g', "-g G", "Sets GrainSize.",			        TYPE_INT , &GrainSize },
    END_OF_ARGUMENTS
  };

  parseArguments(argc,argv,as);

  if (NBK==-1) NBK = t;
  if(q==0){
    PAR_BLOCK {
      //benchmark_with_field<Givaro::ZRing<int32_t>>( p,  m, k, NBK, b, seed, iters, t, as);
      benchmark_with_field<Givaro::ZRing<Givaro::Integer>>( p,  m, k, NBK, b, seed, iters, t, as, GrainSize);
    }
  }else{
    PAR_BLOCK {
      //benchmark_with_field<Givaro::Modular<float>>(q, p,  m, k, NBK, b, seed, iters, t, as);
      //benchmark_with_field<Givaro::Modular<double>>(q, p,  m, k, NBK, b, seed, iters, t, as);
      //benchmark_with_field<Givaro::Modular<int32_t>>(q, p,  m, k, NBK, b, seed, iters, t, as);

      //benchmark_with_field<Givaro::Modular<Givaro::Integer>>(q, p,  m, k, NBK, b, seed, iters, t, as);

      //benchmark_with_field<Givaro::ModularBalanced<float>>(q, p,  m, k, NBK, b, seed, iters, t, as);
      benchmark_with_field<Givaro::ModularBalanced<double>>(q, p,  m, k, NBK, b, seed, iters, t, as, GrainSize);
      //benchmark_with_field<Givaro::ModularBalanced<int32_t>>(q, p,  m, k, NBK, b, seed, iters, t, as);
    }
  }

  return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
