/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

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

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular-balanced.h>

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/args-parser.h"

#include "tests/test-utils.h"

//#define __FFLASFFPACK_USE_DATAFLOW

using namespace std;

//#ifdef __FFLASFFPACK_USE_DATAFLOW
template<class Element>
void Initialize(Element * C, size_t BS, size_t m, size_t n)
{
//#pragma omp parallel for collapse(2) schedule(runtime) 
	BS=std::max(BS, (size_t) __FFLASFFPACK_WINOTHRESHOLD_BAL );
	PAR_BLOCK{
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
//#endif // __FFLASFFPACK_USE_DATAFLOW

int main(int argc, char** argv) {

  size_t iter = 3 ;
  int q = 131071 ;
  size_t m = 2000 ;
  size_t n = 2000 ;
  int p=1; // 0 for sequential 1 for pIter-sRec ; 2 for pRec; 3 for hybrid
  int t=MAX_THREADS;
  int NBK = -1;

  Argument as[] = {
	  { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INT , &q },
	  { 'm', "-m M", "Set the row dimension of A.",      TYPE_INT , &m },
	  { 'n', "-n N", "Set the col dimension of B.",      TYPE_INT , &n },
	  { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iter },
	  { 't', "-t T", "number of virtual threads to drive the partition.", TYPE_INT , &t },
	  { 'b', "-b B", "number of numa blocks per dimension for the numa placement", TYPE_INT , &NBK },
	  { 'p', "-p P", "0 for sequential, 1 for Iterative, 2 for Recursive, 3 for Hybrid.", TYPE_INT , &p },
	  END_OF_ARGUMENTS
	};
  FFLAS::parseArguments(argc,argv,as);
  
  if (NBK==-1) NBK = t;

  typedef Givaro::ModularBalanced<double> Field;
  typedef Field::Element Element;

  Field F(q);
  Element * A;
  Element * B;

   FFLAS::OMPTimer chrono;
  double time=0.0;

  Field::RandIter G(F);
      // if (argc > 5){
      // 	  A = read_field (F, argv[5], &n, &n);
      // }
      // else{
  A = FFLAS::fflas_new (F,m,m,Alignment::CACHE_PAGESIZE);
  Initialize(A,m/NBK,m,m);

  FFLAS::ParSeqHelper::Parallel H;
  size_t i;
  PARFOR1D (i,0,(size_t)m, H,
            for (size_t j = 0; j< (size_t)m; ++j)
            	G.random(*(A+i*m+j));
            );
  
      //}
  
      // if (argc == 7){
      // 	  B = read_field (F, argv[6], &n, &n);
      // }
	  // else{
  B = FFLAS::fflas_new(F,m,n,Alignment::CACHE_PAGESIZE);
  Initialize(B,m/NBK,m,n);
  PARFOR1D (i,0,(size_t)m,H,
            for (size_t j=0 ; j< (size_t)n; ++j)
            	G.random(*(A+i*m+j));
            );
  
      //}
  for (size_t k=0;k<(size_t)m;++k)
	  while (F.isZero( G.random(*(A+k*(m+1)))));
  for (i=0;i<=iter;++i){
      
	  chrono.clear();
	  if (i) chrono.start();

	  if (!p)
		  FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower, 
				FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
				m,n, F.one, A, m, B, n);
	  else{
	  PAR_BLOCK{
	      switch (p) {
		  case 1: {
			  FFLAS::TRSMHelper<FFLAS::StructureHelper::Iterative,
					    FFLAS::ParSeqHelper::Parallel> 
				  PH (FFLAS::ParSeqHelper::Parallel(t,FFLAS::BLOCK_THREADS));
			  FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower, 
					FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
					m,n, F.one, A, m, B, n, PH);
			  break;}
		  case 2: {FFLAS::TRSMHelper<FFLAS::StructureHelper::Recursive, 
					     FFLAS::ParSeqHelper::Parallel> 
				      PH (FFLAS::ParSeqHelper::Parallel(t,FFLAS::ROW_THREADS));
			  FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower, 
					FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
					m,n, F.one, A, m, B, n, PH); 
			  break;}
		  case 3: 
			  FFLAS::TRSMHelper<FFLAS::StructureHelper::Hybrid, FFLAS::ParSeqHelper::Parallel> 
				  PH (FFLAS::ParSeqHelper::Parallel(t,FFLAS::ROW_THREADS));
			  FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower, 
					FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
					m,n, F.one, A, m, B, n, PH);
			  break;
	      }
      
	  }
	  BARRIER;
	  }
	  if (i) {chrono.stop(); time+=chrono.realtime();}

  }
  FFLAS::fflas_delete( A);
  FFLAS::fflas_delete( B);
  
	// -----------
	// Standard output for benchmark - Alexis Breust 2014/11/14
	std::cout << "Time: " << time / double(iter)
			  << " Gflops: " << double(m)/1000. * double(m)/1000. * double(n)/1000. / time * double(iter);
	FFLAS::writeCommandString(std::cout, as) << std::endl;

  return 0;
}
