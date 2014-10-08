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

#include <iostream>

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "tests/test-utils.h"


using namespace std;

int main(int argc, char** argv) {

  // parameter: p, n, iteration, file1, file2

  // int    p    = argc>1 ? atoi(argv[1]) : 1009;
  // int    n    = argc>2 ? atoi(argv[2]) : 2000;
  // size_t iter = argc>3 ? atoi(argv[3]) :    1;
  // size_t strat= argc>4 ? atoi(argv[4]) :    1;
  size_t iter = 3 ;
  int q = 131071 ;
  size_t m = 2000 ;
  size_t n = 2000 ;
  int p=1; // 0 for sequential 1 for pIter-sRec ; 2 for pRec; 3 for hybrid
  Argument as[] = {
	  { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INT , &q },
	  { 'm', "-m M", "Set the row dimension of A.",      TYPE_INT , &m },
	  { 'n', "-n N", "Set the col dimension of B.",      TYPE_INT , &n },
	  { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iter },
	  { 'p', "-p P", "0 for sequential, 1 for Iterative, 2 for Recursive, 3 for Hybrid.", TYPE_INT , &p },
	  END_OF_ARGUMENTS
	};
  FFLAS::parseArguments(argc,argv,as);

  typedef FFPACK::ModularBalanced<double> Field;
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
  A = FFLAS::fflas_new<Element>(m*m);
  for (size_t j = 0; j< (size_t)m*m; ++j)
	  G.random(*(A+j));
      //}
  
      // if (argc == 7){
      // 	  B = read_field (F, argv[6], &n, &n);
      // }
	  // else{
  B = FFLAS::fflas_new<Element>(m*n);
  for (size_t j=0 ; j< (size_t)m*n; ++j)
	  G.random(*(A+j));
      //}
  for (size_t k=0;k<(size_t)m;++k)
	  while (F.isZero( G.random(*(A+k*(m+1)))));
  for (size_t i=0;i<iter;++i){
      
	  chrono.clear();
	  chrono.start();

	  if (!p)
		  FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower, 
				FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
				m,n, F.one, A, n, B, n);
	  else{
	  PAR_REGION{
	      switch (p) {
		  case 1: {
			  FFLAS::TRSMHelper<FFLAS::StructureHelper::Iterative,
					    FFLAS::ParSeqHelper::Parallel> PH 
			      (FFLAS::ParSeqHelper::Parallel(NUM_THREADS,FFLAS::BLOCK_THREADS));
			  FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower, 
					FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
					m,n, F.one, A, n, B, n, PH);
			  break;}
		  case 2: {FFLAS::TRSMHelper<FFLAS::StructureHelper::Recursive, FFLAS::ParSeqHelper::Parallel> PH (FFLAS::ParSeqHelper::Parallel(NUM_THREADS,FFLAS::BLOCK_THREADS));
			  FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower, 
					FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
					m,n, F.one, A, n, B, n, PH); 
			  break;}
		  case 3: 
			  FFLAS::TRSMHelper<FFLAS::StructureHelper::Hybrid, FFLAS::ParSeqHelper::Parallel> PH (FFLAS::ParSeqHelper::Parallel(NUM_THREADS,FFLAS::BLOCK_THREADS));
			  FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower, 
					FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
					m,n, F.one, A, n, B, n, PH);
			  break;
	      }
      
	  }
	  BARRIER;
	  }
	  chrono.stop();
	  time+=chrono.realtime();

  }
  FFLAS::fflas_delete( A);
  FFLAS::fflas_delete( B);

  cerr<<"n: "<<n<<" p: "<<p<<" time: "<<time/(double)iter<<" n^3/time/10^9: "<<(double(n)/1000.*double(n)/1000.*double(n)/1000./time*double(iter))<<endl;


  return 0;
}
