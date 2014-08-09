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

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"

using namespace std;

int main(int argc, char** argv) {

  // parameter: p, n, iteration, file1, file2

  int    p    = argc>1 ? atoi(argv[1]) : 1009;
  int    n    = argc>2 ? atoi(argv[2]) : 2000;
  size_t iter = argc>3 ? atoi(argv[3]) :    1;

  typedef FFPACK::ModularBalanced<double> Field;
  typedef Field::Element Element;

  Field F(p);

  Timer chrono, freidvals;
  double time=0.0, timev=0.0;

  Element * A, * B, * C;
  Element *v, *w, *y;

  for (size_t i=0;i<iter;++i){

	  if (argc > 4){
		  A = read_field (F, argv[4], &n, &n);
	  }
	  else{
		  Field::RandIter G(F);
		  A = new Element[n*n];
		  for (size_t j=0; j<(size_t)n*n; ++j)
			  G.random (*(A+j));
	  }

	  if (argc == 6){
		  B = read_field (F, argv[5], &n, &n);
	  }
	  else{
		  Field::RandIter G(F);
		  B = new Element[n*n];
		  for (size_t j=0; j<(size_t)n*n; ++j)
			  G.random(*(B+j));
	  }

	  C = new Element[n*n];

      v = new Element[n];
      {
          Field::RandIter G(F);
          for(size_t j=0; j<(size_t)n; ++j)
              G.random(*(v+j));
      }
      
      w = new Element[n];
      y = new Element[n];

	  chrono.clear();
	  chrono.start();
	  FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n,n,n, F.one,
                    A, n, B, n, F.zero, C,n);
	  chrono.stop();
	  time+=chrono.usertime();

      freidvals.clear();
      freidvals.start();
      FFLAS::fgemv(F, FFLAS::FflasNoTrans, n,n, F.one, 
                   C, n, v, 1, F.zero, w, 1);
      FFLAS::fgemv(F, FFLAS::FflasNoTrans, n,n, F.one, 
                   B, n, v, 1, F.zero, y, 1);
      FFLAS::fgemv(F, FFLAS::FflasNoTrans, n,n, F.one, 
                   A, n, y, 1, F.zero, v, 1);
      bool pass=true;
      for(size_t j=0; j<(size_t)n; ++j) pass &= ( *(w+j) == *(v+j) );
      freidvals.stop();
      timev+=freidvals.usertime();

      std::cerr << *A << ' ' << *B << ' ' << *C << ' '<< pass << std::endl;
	  delete[] A;
	  delete[] B;
	  delete[] C;
  }

  std::cerr<<"n: "<<n<<" p: "<<p<<" time: "<<time/(double)iter<<" 2n^3/time/10^9: "<<(2.*double(n)/1000.*double(n)/1000.*double(n)/1000./time*double(iter))<<std::endl;  
  std::cerr<<"Freidvals vtime: "<<timev/(double)iter<<std::endl;

  return 0;
}

