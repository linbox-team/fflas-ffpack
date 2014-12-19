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
#include <vector>


#ifndef __FFLASFFPACK_HAVE_DGETRF
#define __FFLASFFPACK_HAVE_DGETRF 1
#endif

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/args-parser.h"


using namespace std;

#ifdef __FFLASFFPACK_USE_OPENMP
typedef FFLAS::OMPTimer TTimer;
#else
typedef FFLAS::Timer TTimer;
#endif

#ifdef __FFLASFFPACK_USE_DATAFLOW
template<class Element>
void Initialize(Element * C, int BS, size_t m, size_t n)
{
//#pragma omp parallel for collapse(2) schedule(runtime)
    std::cout << "Initialize PAR_REGION " << BS << ", " << m << 'x' << n << std::endl;
    
        BS=std::max(BS, __FFLASFFPACK_WINOTHRESHOLD_BAL );
        PAR_REGION{
        for(size_t q=0; q<m; q+=BS) ///row
                for(size_t pp=0; pp<n; pp+=BS) //column
                {
                        size_t M=BS, MM=BS;
                        if(!(q+BS<m))
                                M=m-q;
                        if(!(pp+BS<n))
                                MM=n-pp;
#pragma omp task
                        {
                        for(size_t j=0; j<M; j++)
                                for(size_t jj=0; jj<MM; jj++)
                                        C[(q+j)*n+pp+jj]=0;
                        }
                }
        #pragma omp taskwait
        }
        // printf("A = \n");
        // for (size_t i=0; i<m; i+=128)
        //  {
        //      for (size_t j=0; j<n; j+=128)
        //      {
        //              int ld = komp_get_locality_domain_num_for( &C[i*n+j] );
        //              printf("%i ", ld);
        //      }
        //      printf("\n");
        //  }

}
#endif




int main(int argc, char** argv) {
  
	size_t iter = 1;
	int    q    = 1009;
	int    n    = 2000;
	std::string file = "";
	
#ifdef __FFLASFFPACK_USE_DATAFLOW
	size_t NBK = MAX_THREADS;
#endif
  
	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
		{ 'n', "-n N", "Set the dimension of the matrix.",               TYPE_INT , &n },
		{ 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &iter },
		{ 'f', "-f FILE", "Set the input file (empty for random).",  TYPE_STR , &file },
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(argc,argv,as);

  typedef FFPACK::Modular<double> Field;
  typedef Field::Element Element;

  Field F(q);
  Field::Element * A;

  TTimer chrono;
  double time=0.0;

  std::vector<int> Piv(n,0);
  if (iter>1) {
	  if (!file.empty()){
		  A = read_field(F, file.c_str(), &n, &n);
	  }
	  else {
		  A = FFLAS::fflas_new<Element>(n*n);
		  Field::RandIter G(F);
#ifdef __FFLASFFPACK_USE_DATAFLOW
          Initialize(A,n/NBK,n,n);
#endif
#pragma omp parallel for
		  for (int i=0; i<n; ++i)
              for (int j=0; j<n; ++j)
                  G.random(*(A+i*n+j));
	  }
	  clapack_dgetrf(CblasRowMajor,n,n,A,n,&Piv[0]);
	  FFLAS::fflas_delete( A);
  }

  for (size_t i=0;i<iter;++i){
	  if (!file.empty()){
		  A = read_field(F, file.c_str(), &n, &n);
	  }
	  else {
		  A = FFLAS::fflas_new<Element>(n*n);
		  Field::RandIter G(F);
#ifdef __FFLASFFPACK_USE_DATAFLOW
          Initialize(A,n/NBK,n,n);
#endif
#pragma omp parallel for
		  for (int i=0; i<n; ++i)
              for (int j=0; j<n; ++j)
                  G.random(*(A+i*n+j));
	  }

	  chrono.clear();
	  chrono.start();
	  clapack_dgetrf(CblasRowMajor,n,n,A,n,&Piv[0]);
	  chrono.stop();

	  time+=chrono.usertime();
	  FFLAS::fflas_delete( A);
  }
  
	// -----------
	// Standard output for benchmark - Alexis Breust 2014/11/14
	std::cout << "Time: " << time / double(iter)
			  << " Gflops: " << (2.*double(n)/1000.*double(n)/1000.*double(n)/1000.0) / time * double(iter) / 3.;
	FFLAS::writeCommandString(std::cout, as) << std::endl;

  return 0;
}




