/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
//#include "goto-def.h"

/* Copyright (c) FFLAS-FFPACK
* Written by Clément Pernet <clement.pernet@imag.fr>
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

//#define PFGEMM_WINO_SEQ 32
//#define CLASSIC_SEQ
//#define WINO_SEQ

#include <iostream>
#include <givaro/modular-balanced.h>

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/args-parser.h"

#include "tests/test-utils.h"

#ifdef __FFLASFFPACK_USE_KAAPI
#include "libkomp.h"
#endif


using namespace std;

#ifdef __FFLASFFPACK_USE_DATAFLOW
template<class Element>
void Initialize(Element * C, size_t BS, size_t m, size_t n)
{
//#pragma omp parallel for collapse(2) schedule(runtime) 
	BS=std::max(BS, __FFLASFFPACK_WINOTHRESHOLD_BAL );
	PAR_INSTR{
	for(size_t p=0; p<m; p+=BS) ///row
		for(size_t pp=0; pp<n; pp+=BS) //column
		{
			size_t M=BS, MM=BS;
			if(!(p+BS<m))
				M=m-p;
			if(!(pp+BS<n))
				MM=n-pp;
#pragma omp task 
			{
			for(size_t j=0; j<M; j++)
				for(size_t jj=0; jj<MM; jj++)
					C[(p+j)*n+pp+jj]=0;
			}
		}
	#pragma omp taskwait
	}
	// printf("A = \n");
	// for (size_t i=0; i<m; i+=128)
	//  {
	//  	for (size_t j=0; j<n; j+=128)
	//  	{
	//  		int ld = komp_get_locality_domain_num_for( &C[i*n+j] );
	//  		printf("%i ", ld);
	//  	}
	//  	printf("\n");
	//  }

}
#else
template<class Element>
void Initialize(Element * C, size_t BS, size_t m, size_t n)
{}
#endif
int main(int argc, char** argv) {

	size_t iter = 3 ;
	int q = 131071 ;
	size_t m = 2000 ;
	size_t k = 2000 ;
	size_t n = 2000 ;
	int nbw = -1 ;
	int p=3;
	int t=MAX_THREADS;
	int NBK = -1;

	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INT , &q },
		{ 'm', "-m M", "Set the row dimension of A.",      TYPE_INT , &m },
		{ 'k', "-k K", "Set the col dimension of A.",      TYPE_INT , &k },
		{ 'n', "-n N", "Set the col dimension of B.",      TYPE_INT , &n },
		{ 'w', "-w N", "Set the number of winograd levels (-1 for random).",    TYPE_INT , &nbw },
		{ 'i', "-i R", "Set number of repetitions.",       TYPE_INT , &iter },
		{ 'p', "-p P", "0 for sequential, 1 for 2D iterative, 2 for 2D rec, 3 for 2D rec adaptive, 4 for 3D rc in-place, 5 for 3D rec, 6 for 3D rec adaptive.", TYPE_INT , &p },
		{ 't', "-t T", "number of virtual threads to drive the partition.", TYPE_INT , &t },
		{ 'b', "-b B", "number of numa blocks per dimension for the numa placement", TYPE_INT , &NBK },
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(argc,argv,as);

	if (NBK==-1) NBK = t;
  typedef Givaro::ModularBalanced<double> Field;
//  typedef Givaro::ModularBalanced<float> Field;
  typedef Field::Element Element;

  Field F(q);

  FFLAS::Timer chrono, freivalds;
  double time=0.0, timev=0.0;

  Element * A, * B, * C;

  Field::RandIter G(F); 
  A = FFLAS::fflas_new(F,m,k,Alignment::CACHE_PAGESIZE);
//#pragma omp parallel for collapse(2) schedule(runtime) 
  Initialize(A,m/size_t(NBK),m,k);

  FFLAS::ParSeqHelper::Parallel H;
  size_t i;
//#pragma omp for
  PARFOR1D (i,0, m,H,
	    for (size_t j=0; j<(size_t)k; ++j)
		    G.random (*(A+i*k+j));
	    );
  
  B = FFLAS::fflas_new(F,k,n,Alignment::CACHE_PAGESIZE);
//#pragma omp parallel for collapse(2) schedule(runtime) 
  Initialize(B,k/NBK,k,n);
//#pragma omp parallel for
  PARFOR1D (i, 0, k,H,
            for (size_t j=0; j<(size_t)n; ++j)
            	G.random(*(B+i*n+j));
            );
  
  C = FFLAS::fflas_new(F,m,n,Alignment::CACHE_PAGESIZE);
//#pragma omp parallel for collapse(2) schedule(runtime) 
  Initialize(C,m/NBK,m,n);
  for (i=0;i<=iter;++i){

	  // if (argc > 4){
	  // 	  A = read_field (F, argv[4], &n, &n);
	  // }
	  // else{

      chrono.clear();
      if (p && p!=7){
	      FFLAS::CuttingStrategy meth;
	      switch (p){
		  case 1: meth = FFLAS::BLOCK_THREADS;break;
		  case 2: meth = FFLAS::TWO_D;break;
		  case 3: meth = FFLAS::TWO_D_ADAPT;break;
		  case 4: meth = FFLAS::THREE_D_INPLACE;break;
		  case 5: meth = FFLAS::THREE_D;break;
		  case 6: meth = FFLAS::THREE_D_ADAPT;break;
		  default: meth = FFLAS::BLOCK_THREADS;break;
	      }
	      FFLAS::MMHelper<Field,FFLAS::MMHelperAlgo::Winograd,
			      typename FFLAS::ModeTraits<Field>::value,
			      FFLAS::ParSeqHelper::Parallel> 
		      WH (F, nbw, FFLAS::ParSeqHelper::Parallel(t, meth));	
	      if (i) chrono.start();
	      PAR_INSTR{
		      FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,k, F.one, A, k, B, n, F.zero, C,n,WH);
	      }
	      if (i) {chrono.stop(); time+=chrono.realtime();}

	      
      }else{
	      if(p==7){

		      FFLAS::MMHelper<Field, FFLAS::MMHelperAlgo::WinogradPar>
			      WH (F, nbw, FFLAS::ParSeqHelper::Sequential());
		      //		      cout<<"wino parallel"<<endl;
		      if (i) chrono.start();
		      PAR_INSTR
		      {
			      FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,k, F.one, A, k, B, n, F.zero, C,n,WH);
		      }
		      if (i) {chrono.stop(); time+=chrono.realtime();}
	      }
	      else{

		      FFLAS::MMHelper<Field,FFLAS::MMHelperAlgo::Winograd>//,
			      //typename FFLAS::FieldTraits<Field>::value,
			      //FFLAS::ParSeqHelper::Parallel>
			      WH (F, nbw, FFLAS::ParSeqHelper::Sequential());
		      if (i) chrono.start();
		      FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,k, F.one, A, k, B, n, F.zero, C,n,WH);
		      if (i) {chrono.stop(); time+=chrono.realtime();}
	      }
      }


      
      freivalds.clear();
      freivalds.start();
      bool pass = FFLAS::freivalds(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,k, F.one, A, k, B, n, C,n);
      freivalds.stop();
      timev+=freivalds.usertime();
      if (!pass) 
	      std::cout<<"FAILED"<<std::endl;
	  //std::cout << *A << ' ' << *B << ' ' << *C << ' '<< pass << std::endl;
  }
  FFLAS::fflas_delete( A);
  FFLAS::fflas_delete( B);
  FFLAS::fflas_delete( C);
  
	// -----------
	// Standard output for benchmark - Alexis Breust 2014/11/14
	std::cout << "Time: " << time / double(iter)
			  << " Gflops: " << (2.*double(m)/1000.*double(n)/1000.*double(k)/1000.0) / time * double(iter);
	FFLAS::writeCommandString(std::cout, as) << std::endl;
  
#if DEBUG
        std::cout<<"Freivalds vtime: "<<timev/(double)iter<<std::endl;
#endif

  return 0;
}

