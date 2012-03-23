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

  int    p    = atoi(argv[1]);
  int n    = atoi(argv[2]);
  size_t iter = atoi(argv[3]);


  // typedef FFPACK::Modular<double> Field;
  // typedef FFPACK::Modular<float> Field;
  // typedef FFPACK::ModularBalanced<double> Field;
  typedef FFPACK::ModularBalanced<float> Field;
  typedef Field::Element Element;

  Field F(p);

  Timer chrono;
  double time=0.0;// time2=0.0;

  Element * A, * B, * C;

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

	  chrono.clear();
	  chrono.start();
	  FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n,n,n, F.one,
			A, n, B, n, F.zero, C,n);
	  chrono.stop();
	  time+=chrono.usertime();

	  delete[] A;
	  delete[] B;
	  delete[] C;
  }

  std::cerr<<"n: "<<n<<" p: "<<p<<" time: "<<time/(double)iter<<std::endl;

  return 0;
}

