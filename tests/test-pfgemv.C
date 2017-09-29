/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/*
 * Copyright (C) the FFLAS-FFPACK group
 * Written by ZHU Hongguang
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

#include <iomanip>
#include <iostream>
using namespace std;

#define  __FFLASFFPACK_USE_OPENMP

#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "givaro/modular.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/fflas/fflas.h"

/*
#ifdef __FFLASFFPACK_USE_KAAPI
#include <kaapi++>
#endif
*/
#ifdef __FFLASFFPACK_USE_OPENMP
#include <omp.h>
#endif

#include "fflas-ffpack/paladin/pfgemv.inl"


using namespace FFLAS;

#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
using namespace std;
using namespace FFPACK;
using Givaro::Modular;
using Givaro::ModularBalanced;


template<typename Field, class RandIter>
bool check_solve(const Field &F, size_t m, RandIter& Rand){

  typename Field::Element_ptr A, X, Y, Y2;

  size_t lda,incX,incY;
  lda=m;  
  incX=1;  
  incY=1;
  A  = FFLAS::fflas_new(F,m,lda);
  
  Y  = FFLAS::fflas_new(F,m,incY);
  Y2 = FFLAS::fflas_new(F,m,incY);
  X  = FFLAS::fflas_new(F,m,incX);
  
  RandomMatrix (F, m, m, A, lda);
  RandomMatrix (F, m, incX, X,  incX);
  RandomMatrix (F, m, incY, Y,  incY);
  RandomMatrix (F, m, incY, Y2, incY);

  FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y,  incY);
  


	FFLAS::Timer t; t.clear();
	double time=0.0;
	bool ok = true;

	t.clear();
	t.start();
 {
    MMHelper<Field, MMHelperAlgo::Classic, ModeTraits<Field>, ParSeqHelper::Parallel<CuttingStrategy::Recursive,StrategyParameter::Threads> >  H;
    FFLAS::pfgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y2,  incY, H);
 }
	t.stop();
	time+=t.usertime();
  if (FFLAS::fequal (F, m, 1, Y2, incY, Y, incY)){
    cout << "PASSED ("<<time<<")"<<endl;
    
  } else{
    ok=false;
  }

 	time=0.0;
	t.clear();
	t.start();
  {
    MMHelper<Field, MMHelperAlgo::Classic, ModeTraits<Field>, ParSeqHelper::Parallel<CuttingStrategy::Row,StrategyParameter::Threads> >  H;
    FFLAS::pfgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y2,  incY, H);
  }
  	t.stop();
	time+=t.usertime();
  if (FFLAS::fequal (F, m, 1, Y2, incY, Y, incY)){
    cout << "PASSED ("<<time<<")"<<endl;
    
  } else{
	ok=false;
  }

if(m>2){
  for(size_t GS=2; GS<m; GS++){
 	time=0.0;
	t.clear();
	t.start();
  {
    MMHelper<Field, MMHelperAlgo::Classic, ModeTraits<Field>, 						 					ParSeqHelper::Parallel<CuttingStrategy::Row,StrategyParameter::Grain> >  H;
    FFLAS::pfgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y2,  incY, GS, H); 
  }
 //FFLAS::WriteMatrix (std::cout << "A:"<< std::endl, F, m, m, A, lda) << std::endl;
   // FFLAS::WriteMatrix (std::cout << "X:"<< std::endl, F, m, incX, X, incX) << std::endl;
	t.stop();
	time+=t.usertime();
  if (FFLAS::fequal (F, m, 1, Y2, incY, Y, incY)){    
    cout << "PASSED ("<<time<<")"<<endl;   
  } else{
	ok=false;
	//cout << "failed	with GS = "<<GS<<endl;
	//break;
  }
 }
}else{
	size_t GS=2;
 	time=0.0;
	t.clear();
	t.start();
  {
    MMHelper<Field, MMHelperAlgo::Classic, ModeTraits<Field>, 						 					ParSeqHelper::Parallel<CuttingStrategy::Row,StrategyParameter::Grain> >  H;
    FFLAS::pfgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y2,  incY, GS, H); 
  }
	t.stop();
	time+=t.usertime();
  if (FFLAS::fequal (F, m, 1, Y2, incY, Y, incY)){    
    cout << "PASSED ("<<time<<")"<<endl;   
  } else{
	ok=false;
	//cout << "failed	with GS = "<<GS<<endl;
  }
}
	FFLAS::fflas_delete(A);
	FFLAS::fflas_delete(X);
	FFLAS::fflas_delete(Y);
	FFLAS::fflas_delete(Y2);


	return ok;
}
template <class Field>
bool run_with_field (Givaro::Integer q, size_t b, size_t m, size_t iters, uint64_t seed){
	bool ok = true ;
	int nbit=(int)iters;

	while (ok &&  nbit){
		//typedef typename Field::Element Element ;
		// choose Field
		//Field* F= chooseField<Field>(3,b);
		Field* F= chooseField<Field>(q,b);
		typename Field::RandIter G(*F,0,seed);
		if (F==nullptr)
			return true;

		cout<<"Checking with ";F->write(cout)<<endl;
		ok = ok && check_solve(*F,m,G);
		
		nbit--;
		delete F;

	}

	return ok;
}







BEGIN_PARALLEL_MAIN(int argc, char** argv)
{

	cerr<<setprecision(10);
	Givaro::Integer q=-1;
	size_t b=0;
	size_t m=256;

	size_t iters=20;
	bool loop=false;
	uint64_t seed = time(NULL);
	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
		{ 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
		{ 'm', "-m M", "Set the dimension of unknown square matrix.",      TYPE_INT , &m },

		{ 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
		{ 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
		{ 's', "-s seed", "Set seed for the random generator", TYPE_INT, &seed },
                END_OF_ARGUMENTS
        };

	FFLAS::parseArguments(argc,argv,as);

	bool ok = true;

	do{
		ok = ok && run_with_field<Modular<double> >(q,b,m,iters,seed);
		ok = ok && run_with_field<ModularBalanced<double> >(q,b,m,iters,seed);
		ok = ok && run_with_field<Modular<float> >(q,b,m,iters,seed); 
		ok = ok && run_with_field<ModularBalanced<float> >(q,b,m,iters,seed); 
		ok = ok && run_with_field<Modular<int32_t> >(q,b,m,iters,seed); 
		ok = ok && run_with_field<ModularBalanced<int32_t> >(q,b,m,iters,seed); 
		ok = ok && run_with_field<Modular<int64_t> >(q,b,m,iters,seed);
		ok = ok && run_with_field<ModularBalanced<int64_t> >(q,b,m,iters,seed); 
		ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,5,m/6,iters,seed);
		ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,(b?b:512),m/6,iters,seed); 

	} while (loop && ok);


  
}
END_PARALLEL_MAIN()


