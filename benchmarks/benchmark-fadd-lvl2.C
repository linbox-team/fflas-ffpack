/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


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
#include <givaro/modular.h>

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/args-parser.h"


using namespace std;
using namespace FFLAS;
using namespace FFPACK;
int main(int argc, char** argv) {
  
  size_t iter = 3;
  int    q    = 131071;
  size_t    n    = 1000;
  size_t    m    = 1000;
  double    a    = 1.0;
  size_t threshold = 64;
  std::string file = "";
  
  Argument as[] = {
    { 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
    { 'n', "-n N", "Set the dimension of the matrix C.",             TYPE_INT , &n },
    { 'm', "-m M", "Set the other dimension of the matrix C.",       TYPE_INT , &m },
    { 'a', "-a A", "Set the value of alpha.",                        TYPE_DOUBLE , &a},
    { 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &iter },
    { 't', "-t T", "Set the threshold to the base case.",            TYPE_INT , &threshold },
    END_OF_ARGUMENTS
  };

  FFLAS::parseArguments(argc,argv,as);
	
  typedef Givaro::ModularBalanced<double> Field;
  typedef Field::Element Element;
	
  Field F(q);
  Field::ConstElement_ptr A, B;
  Field::Element_ptr C;
  Element alpha;
  F.init(alpha, a);

	
  FFLAS::Timer chrono;
  double time=0.0;
	
  for (size_t i=0;i<=iter;++i){
    A = fflas_new(F,n,m);
    size_t lda=m;
    B = fflas_new(F,n,m);
    size_t ldb=m;
    C = fflas_new(F,n,m);
    size_t ldc=m;
    Field::RandIter G(F);
    RandomMatrix (F, n, m, A, lda, G);
    RandomMatrix (F, n, m, B, ldb, G);
    chrono.clear();
    if (i) chrono.start();
    fadd(F,n,m,A,lda,B,ldb,C,ldc);
    if (i) chrono.stop();
		
    time+=chrono.usertime();
    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);
    FFLAS::fflas_delete(C);
  }
  
  // -----------
  // Standard output for benchmark - Alexis Breust 2014/11/14
#define SQUARE(x) ((x)*(x))
  std::cout << "Time: " << time / double(iter)
	    << " Gflops: " << SQUARE(double(n)/1000.)/ time * double(iter);
  FFLAS::writeCommandString(std::cout, as) << std::endl;
  return 0;
}
