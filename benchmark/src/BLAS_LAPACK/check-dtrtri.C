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

#define __FFLASFFPACK_HAVE_DTRTRI 1

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"


using namespace std;

int main(int argc, char** argv) {

  // parameter: p, n, iteration, file

  int    p    = atoi(argv[1]);
  int n    = atoi(argv[2]);
  size_t iter = atoi(argv[3]);


  typedef FFPACK::Modular<double> Field;
  typedef Field::Element Element;

  Field F(p);
  Element * A;

  Timer chrono;
  double time=0.0;

  Field::RandIter G(F);
  for (size_t i=0;i<iter;++i){
    if (argc > 4){
      A = read_field (F, argv[4], &n, &n);
    } else {
      A = new Element[n*n];
      for (size_t j=0; j<(size_t) n*n; ++j)
	G.random(*(A+j));
    }
    for (size_t k=0;k<(size_t)n;++k)
      while (F.isZero( G.random(*(A+k*(n+1)))));

    chrono.clear();
    chrono.start();
    clapack_dtrtri(CblasRowMajor,CblasUpper, CblasNonUnit,n,A,n);
    chrono.stop();

    time+=chrono.usertime();
    delete[] A;

  }

  cerr<<"n: "<<n<<" p: "<<p<<" time: "<<time/(double)iter<<endl;


  return 0;
}
