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


using namespace std;

#ifdef __FFLASFFPACK_USE_OPENMP
typedef FFLAS::OMPTimer TTimer;
#else
typedef FFLAS::Timer TTimer;
#endif

int main(int argc, char** argv) {

  // parameter: p, n, iteration, file

  int    p    = argc>1 ? atoi(argv[1]) : 1009;
  int    n    = argc>2 ? atoi(argv[2]) : 2000;
  size_t iter = argc>3 ? atoi(argv[3]) :    1;


  typedef FFPACK::Modular<double> Field;
  typedef Field::Element Element;

  Field F(p);
  Field::Element * A;

  TTimer chrono;
  double time=0.0;

  std::vector<int> Piv(n,0);
  if (iter>1) {
	  if (argc > 4){
		  A = read_field(F, argv[4], &n, &n);
	  }
	  else {
		  A = FFLAS::fflas_new<Element>(n*n);
		  Field::RandIter G(F);
#pragma omp parallel for
		  for (size_t i=0; i<n; ++i)
              for (size_t j=0; j<n; ++j)
                  G.random(*(A+i*n+j));
	  }
	  clapack_dgetrf(CblasRowMajor,n,n,A,n,&Piv[0]);
	  FFLAS::fflas_delete( A);
  }

  for (size_t i=0;i<iter;++i){
	  if (argc > 4){
		  A = read_field(F, argv[4], &n, &n);
	  }
	  else {
		  A = FFLAS::fflas_new<Element>(n*n);
		  Field::RandIter G(F);
#pragma omp parallel for
		  for (size_t i=0; i<n; ++i)
              for (size_t j=0; j<n; ++j)
                  G.random(*(A+i*n+j));
	  }

	  chrono.clear();
	  chrono.start();
	  clapack_dgetrf(CblasRowMajor,n,n,A,n,&Piv[0]);
	  chrono.stop();

	  time+=chrono.usertime();
	  FFLAS::fflas_delete( A);
  }

  cerr<<"n: "<<n<<" p: "<<p<<" time: "<<time/(double)iter<<" 2n^3/3/time/10^9: "<<(2.*double(n)/1000.*double(n)/1000.*double(n)/1000./time*double(iter)/3.)<<endl;


  return 0;
}




