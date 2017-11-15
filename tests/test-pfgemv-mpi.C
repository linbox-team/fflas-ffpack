/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) the FFLAS-FFPACK group
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

//#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "givaro/modular.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/fflas/fflas.h"


#include "fflas-ffpack/utils/test-utils.h"
/*
#ifdef __FFLASFFPACK_USE_KAAPI
#include <kaapi++>
#endif
*/
#ifdef __FFLASFFPACK_USE_OPENMP
#include <omp.h>
#endif

//#include "fflas-ffpack/paladin/pfgemv.inl"
#include <mpi.h>
#include "fflas-ffpack/paladin/pfgemv-mpi.inl"


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

double time=0.0;
bool ok = true;	
FFLAS::Timer t; t.clear();

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  A  = FFLAS::fflas_new(F,m,lda);
  
  Y  = FFLAS::fflas_new(F,m,incY);
  Y2 = FFLAS::fflas_new(F,m,incY);
  X  = FFLAS::fflas_new(F,m,incX);




if(rank==0){
	RandomMatrix (F, m, m, A, lda);
	RandomMatrix (F, m, incX, X,  incX);
	RandomMatrix (F, m, incY, Y,  incY);
	RandomMatrix (F, m, incY, Y2, incY);  
	
	//FFLAS::WriteMatrix (std::cout << "init A:="<<std::endl, F, m, m, A, lda) << std::endl;
	//FFLAS::WriteMatrix (std::cout << "init X:="<<std::endl, F, m, incX, X, incX) << std::endl;
	//FFLAS::WriteMatrix (std::cout << "init Y:="<<std::endl, F, m, incY, Y, incY) << std::endl;
	
	
	t.clear();
	t.start();
	
	{cout <<"rank("<<rank<< ") begin to execute"<<endl;
		MMHelper<Field, MMHelperAlgo::Classic, ModeTraits<Field>, ParSeqHelper::Parallel<CuttingStrategy::Recursive,StrategyParameter::Threads> >  H;
		
		pfgemv_mpi(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y2,  incY, H);
cout <<"rank("<<rank<< ") finished"<<endl;
	}
	
	t.stop();
	time+=t.usertime();
	
	FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y,  incY);

	if (FFLAS::fequal (F, m, 1, Y2, incY, Y, incY)){
		cout << "PASSED ("<<time<<")"<<endl; 
	} else{
		cout << "*********************************failed1*********************************"<<endl;
		FFLAS::WriteMatrix (std::cout <<"Found Y2:"<<std::endl, F, m, incY, Y2, incY) << std::endl;
		FFLAS::WriteMatrix (std::cout <<"Found Y:"<<std::endl, F, m, incX, Y, incY) << std::endl;
		ok=false; 
	}
 }else{

	{cout <<"rank("<<rank<< ") begin to execute"<<endl;
		MMHelper<Field, MMHelperAlgo::Classic, ModeTraits<Field>, ParSeqHelper::Parallel<CuttingStrategy::Recursive,StrategyParameter::Threads> >  H;
		
		pfgemv_mpi(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y2,  incY, H);
cout <<"rank("<<rank<< ") finished"<<endl;
	}

 }


cout << "process("<<rank<<") continues to NEXT "<<endl;
	

if(rank==0){

	RandomMatrix (F, m, incY, Y,  incY);
	RandomMatrix (F, m, incY, Y2, incY); 	
  
	t.clear();
	t.start();
		
	{
cout <<"0rank("<<rank<< ") begin to execute"<<endl;
		MMHelper<Field, MMHelperAlgo::Classic, ModeTraits<Field>, ParSeqHelper::Parallel<CuttingStrategy::Row,StrategyParameter::Threads> >  H;		
		pfgemv_mpi(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y2,  incY, H);
cout <<"0rank("<<rank<< ") finished"<<endl;
	}
	
	t.stop();
	time+=t.usertime();

	FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y,  incY);
	if (FFLAS::fequal (F, m, 1, Y2, incY, Y, incY)){
		cout << "PASSED ("<<time<<")"<<endl; 
	} else{
		cout << "*********************************failed2********************************"<<endl;
		//FFLAS::WriteMatrix (std::cout <<"Found Y2:"<<std::endl, F, m, incY, Y2, incY) << std::endl;
		//FFLAS::WriteMatrix (std::cout <<"Found Y:"<<std::endl, F, m, incX, Y, incY) << std::endl;
		ok=false; 
	}
	

	
 }else{

	{cout <<"0rank("<<rank<< ") begin to execute"<<endl;
		MMHelper<Field, MMHelperAlgo::Classic, ModeTraits<Field>, ParSeqHelper::Parallel<CuttingStrategy::Row,StrategyParameter::Threads> >  H;
		
		pfgemv_mpi(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y2,  incY, H);
cout <<"0rank("<<rank<< ") finished"<<endl;
	}

 }
 
 


if(rank==0){

	FFLAS::fflas_delete(A);
	FFLAS::fflas_delete(X);
	FFLAS::fflas_delete(Y);
	FFLAS::fflas_delete(Y2);

}

return ok;

}


template <class Field>
bool run_with_field (Givaro::Integer q, size_t b, size_t m, size_t iters, uint64_t seed){
	bool ok = true ;
	int nbit=(int)iters;
//MPI_Init(NULL, NULL);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	while (ok &&  nbit){
		//typedef typename Field::Element Element ;
		// choose Field
		//Field* F= chooseField<Field>(3,b);

		Field* F= FFPACK::chooseField<Field>(q,b,seed);
		typename Field::RandIter G(*F,0,seed);
		if (F==nullptr)
			return true;

		if(rank==0)cout<<"Checking with ";F->write(cout)<<endl;
MPI_Barrier(MPI_COMM_WORLD);
		ok = ok && check_solve(*F,m,G);
MPI_Barrier(MPI_COMM_WORLD);		
		nbit--;
		delete F;

	}
//MPI_Finalize();
	return ok;
}







BEGIN_PARALLEL_MAIN(int argc, char** argv)
{

	cerr<<setprecision(10);
	Givaro::Integer q=-1;
	size_t b=0;
	size_t m=1741;

	size_t iters=3;
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

//MPI_Init(&argc, &argv); 
int provided;
MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
/*if (provided < MPI_THREAD_MULTIPLE) {
	std::cout<<" Error - MPI does not provide needed threading level"<<std::endl;

}else{
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
	if (provided < MPI_THREAD_FUNNELED) {
		std::cout<<" Error - MPI does not provide needed thread safety "<<std::endl;
		return 1;
	}
}*/

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
		//ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,5,m/6,iters,seed);
		//ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,(b?b:512),m/6,iters,seed); 

	} while (loop && ok);

MPI_Finalize();
  
}
END_PARALLEL_MAIN()


