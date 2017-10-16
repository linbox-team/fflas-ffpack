/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
//#include "goto-def.h"

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

// Please do not commit with any of these defines on - AB 2015-01-12
//#define __FFLASFFPACK_USE_TBB
//#define __FFLASFFPACK_USE_OPENMP 
//#define __FFLASFFPACK_USE_DATAFLOW
//#define WINO_PARALLEL_TMPS
//#define __FFLASFFPACK_FORCE_SEQ
//#define PFGEMM_WINO_SEQ 32
//#define CLASSIC_SEQ
#define CLASSIC_HYBRID
//#define WINO_SEQ
//#define FFT_PROFILER
//#define PROFILE_FGEMM_MP
#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular-balanced.h>

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/args-parser.h"

#ifdef __FFLASFFPACK_USE_KAAPI
#include "libkomp.h"
#endif


#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/timer.h"
#include "givaro/modular-integer.h"
#include "givaro/givcaster.h"
#include "fflas-ffpack/paladin/parallel.h"
#include "fflas-ffpack/paladin/pfgemv.inl"
using namespace FFPACK;

using namespace std;
using namespace FFLAS;


int main(int argc, char** argv) {


	int p=0;
	size_t iters = 100;
	Givaro::Integer q = 131071; 
	size_t m = 8000;
	size_t k = 8000;
	//static size_t n = 512 ;
    size_t seed= time(NULL);
	int t=MAX_THREADS;
	int NBK = -1;

	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
		{ 'p', "-p P", "0 for sequential, 1 for Recursive CuttingStrategy using StrategyParameter Thread, 2 for Row CuttingStrategy using StrategyParameter Thread, 3 for Row CuttingStrategy using StrategyParameter Grain.", TYPE_INT , &p },
		{ 'm', "-m M", "Set the dimension m of the matrix.",                    TYPE_INT , &m },
		{ 'k', "-k K", "Set the dimension k of the matrix.",                    TYPE_INT , &k },
		{ 't', "-t T", "number of virtual threads to drive the partition.", TYPE_INT , &t },
		{ 'b', "-b B", "number of numa blocks per dimension for the numa placement", TYPE_INT , &NBK },
		{ 'i', "-i R", "Set number of repetitions.",                            TYPE_INT , &iters },
		{ 's', "-s S", "Sets seed.",                            				TYPE_INT , &seed },
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,as);

	if (NBK==-1) NBK = t;
//  typedef Givaro::Modular<double> Field;
//  typedef Givaro::Modular<Givaro::Integer> Field;
//  typedef Givaro::ModularBalanced<int32_t> Field;
//	typedef Givaro::ModularBalanced<float> Field;
 	typedef Givaro::ModularBalanced<double> Field;

//	typedef Field::Element Element;

  Field F(q);
  Timer chrono, TimFreivalds;
  double time=0.0;

		size_t lda,incX,incY;
		lda=k;
		incX=1;
		incY=1;
  Field::Element_ptr A,X,Y; //, Y2;

  Field::RandIter G(F);
		A = FFLAS::fflas_new(F,m,lda);
		X = FFLAS::fflas_new(F,k,incX);
		Y = FFLAS::fflas_new(F,m,incY);
		//Y2= FFLAS::fflas_new(F,m,incY);
	Field::RandIter Rand(F,seed);


		// TODO: replace by a 1D pfrand
	PAR_BLOCK {
			SYNCH_GROUP(
            FORBLOCK1D(iter, m, SPLITTER(NBK, Row, Threads),
                       TASK(MODE(CONSTREFERENCE(F,Rand)),
                       {
                           frand(F, Rand,
                                 iter.iend()-iter.ibegin(),
                                 k,
                                 A+iter.ibegin()*lda,
                                 lda);
                       }
                            );
                       );
			);
//	FFLAS::pfrand(F,Rand, m,k,A,m/NBK);
	}
		FFLAS::frand(F,Rand, k,1,X,incX);
		FFLAS::fzero(F, m,1,Y,incY);
		//FFLAS::fzero(F, m,1,Y2,incY);
  

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
					FFLAS::pfgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y,  incY, H);
					  break;}
				  case 2:{
					ParSeqHelper::Parallel<row, threads>  H(t);
					FFLAS::pfgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y,  incY, H);
					  break;
				  }
				  case 3:{
					size_t BS = 64;
					ParSeqHelper::Parallel<row, grain>  H(BS);
					FFLAS::pfgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y,  incY, H);
					  break;
				  }
				  default:{
					  FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y,  incY);
					  break;
				  }
			  }
		  }//PAR_BLOCK
	      if (i) {chrono.stop(); time+=chrono.realtime();}
      }else{
		      if (i) chrono.start();
		      FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y,  incY);
		      if (i) {chrono.stop(); time+=chrono.realtime();}
	  }

}
  FFLAS::fflas_delete(A);
  FFLAS::fflas_delete(X);
  FFLAS::fflas_delete(Y);
  

	std::cout << "Time: " << time / double(iters)
			  << " Gflops: " << (2.*double(m)/1000.*double(k)/1000000.0) / time * double(iters);
	writeCommandString(std::cout, as) << std::endl;

  return 0;
}
