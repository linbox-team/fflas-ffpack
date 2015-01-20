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

#include <iostream>
#include <omp.h>
#include <givaro/modular-balanced.h>

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/args-parser.h"

using namespace std;

int main(int argc, char** argv) {

#ifdef __FFLASFFPACK_USE_OPENMP
  
	size_t iter = 1;
	int    q    = 1009;
	int    n    = 2000;
	int    w    = -1;
	std::string file1 = "";
	std::string file2 = "";
  
	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
		{ 'n', "-n N", "Set the dimension of the matrix.",               TYPE_INT , &n },
		{ 'w', "-w W", "-------.",               TYPE_INT , &w },
		{ 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &iter },
		{ 'f', "-f FILE", "Set the first input file (empty for random).",   TYPE_STR , &file1 },
		{ 'g', "-g FILE", "Set the second input file (empty for random).",  TYPE_STR , &file2 },
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(argc,argv,as);

  // typedef Givaro::Modular<double> Field;
  // typedef Givaro::Modular<float> Field;
  typedef Givaro::ModularBalanced<double> Field;
  // typedef Givaro::ModularBalanced<float> Field;
  typedef Field::Element Element;

  Field F(q);

   FFLAS::OMPTimer chrono;
  double time=0.0;// time2=0.0;

  Element * A, * B, * C;

  for (size_t i=0;i<iter;++i){

	  if (argc > 5){
		  A = read_field (F, argv[5], &n, &n);
	  }
	  else{
		  Field::RandIter G(F);
		  A = FFLAS::fflas_new<Element>(n*n);
          PAR_FOR(size_t i=0; i<(size_t)n; ++i)
              for (size_t j=0; j<(size_t)n; ++j)
                  G.random (*(A+i*n+j));
	  }

	  if (argc == 7){
		  B = read_field (F, argv[6], &n, &n);
	  }
	  else{
		  Field::RandIter G(F);
		  B = FFLAS::fflas_new<Element>(n*n);
          PAR_FOR(size_t i=0; i<(size_t)n; ++i)
              for (size_t j=0; j<(size_t)n; ++j)
                  G.random (*(B+i*n+j));
	  }

	  C = FFLAS::fflas_new<Element>(n*n);

          const FFLAS::CuttingStrategy Strategy = FFLAS::BLOCK_THREADS;
          enum FFLAS::FFLAS_TRANSPOSE ta = FFLAS::FflasNoTrans;
          enum FFLAS::FFLAS_TRANSPOSE tb = FFLAS::FflasNoTrans;
          Field::Element alpha,beta;
          F.init(alpha); F.init(beta);

          F.assign( alpha, F.one);
          F.assign( beta, F.zero);


	  chrono.clear();
	  chrono.start();
      FFLAS::MMHelper<Field,
          FFLAS::MMHelperAlgo::Winograd,
          FFLAS::FieldTraits<Field>::value,
          FFLAS::ParSeqHelper::Parallel>
          pWH (F, w, FFLAS::ParSeqHelper::Parallel(MAX_THREADS,Strategy));
      
      PAR_REGION{
          FFLAS::fgemm(F, ta, tb,n,n,n,alpha, A,n, B,n, beta,C,n, pWH);
          
      }
      BARRIER;
      
	  chrono.stop();
	  time+=chrono.usertime();

      // std::cout << *A << ' ' << *B << ' ' << *C << std::endl;
	  FFLAS::fflas_delete( A);
	  FFLAS::fflas_delete( B);
	  FFLAS::fflas_delete( C);
  }
  
	// -----------
	// Standard output for benchmark - Alexis Breust 2014/11/14
	std::cout << "Time: " << time / double(iter)
			  << " Gflops: " << 2. * double(n)/1000. * double(n)/1000. * double(n)/1000. / time * double(iter);
	FFLAS::writeCommandString(std::cout, as) << std::endl;
  
#else
  
	// -----------
	// Standard output for benchmark - Alexis Breust 2014/11/14
	std::cout << "Error: You need OpenMP to execute this benchmark.";
	
  std::cout << "you need openmp here"
#endif

  return 0;
}

