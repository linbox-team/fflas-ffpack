/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2016 FFLAS-FFPACK group.
 *
 * Written by Clément Pernet <clement.pernet@imag.fr>
 * Philippe LEDENT <philippe.ledent@etu.univ-grenoble-alpes.fr>
 *
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
 *
 */



#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular-balanced.h>
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"


#ifdef __GIVARO_USE_OPENMP
typedef Givaro::OMPTimer TTimer;
#else
typedef Givaro::Timer TTimer;
#endif

#include <ctime>
#define CUBE(x) ((x)*(x)*(x))
#define GFOPS(n,t) (CUBE(double(n)/1000.0)/(3.0*t))

int main () {
  using namespace std;

  typedef Givaro::ModularBalanced<double> Field;
  Field F(131071);
  size_t n=1000, nmax=1500, k=1000, kmax=1500, prec=1000, nbest=0, count=0;
  TTimer chrono,tim;
  bool bound=false;
  time_t result = std::time(NULL);
  
  // Let C be a random symmetric n by n matrix
  Field::Element_ptr C = FFLAS::fflas_new (F, nmax, nmax);
  size_t ldc = nmax;
  typename Field::RandIter G(F);
  FFPACK::RandomSymmetricMatrix (F, n, true,
				 C, ldc,G);

  // Let B be a copy of C
  Field::Element_ptr B = FFLAS::fflas_new (F, nmax, nmax);
  size_t ldb = ldc;
  FFLAS::fassign (F, n, n, C, ldc, B, ldb);
  
  // Let A be a random n by k matrix
  Field::Element_ptr A = FFLAS::fflas_new (F, nmax, kmax);
  size_t lda = kmax;
  FFPACK::RandomMatrix(F,n,k,A,lda);
  
  // let alpha and beta be scalars in F
  Field::Element alpha = F.mOne, beta = F.one;
  
  cout << std::endl 
       << "---------------------------------------------------------------------"
       << std::endl << std::asctime(std::localtime(&result))
       << std::endl
       << "Threshold for fsyrk base case" ;
  F.write(cout << " (using ") << ')' << endl << endl;

  cout << "fsyrk:  n                   Base case                        Recursive 1 level" << std::endl;
  cout << "                    seconds            Gfops          seconds            Gfops" << std::endl;
  double BCTime, RecTime;
  int iter;
  do{
    iter=300;

    //warm up computation
    FFLAS::fsyrk(F,FFLAS::FflasUpper,FFLAS::FflasNoTrans,n,k,alpha,A,lda,beta,C,ldc,n);
    FFLAS::fassign (F, n, n, B, ldb, C, ldc);
		
    // base case
    chrono.clear();tim.clear();
    for (int i=0;i<iter;i++){
      chrono.start();
      FFLAS::fsyrk(F,FFLAS::FflasUpper,FFLAS::FflasNoTrans,n,k,alpha,A,lda,beta,C,ldc,n);
      chrono.stop();
      tim+=chrono;
      FFLAS::fassign (F, n, n, B, ldb, C, ldc);
    }
    BCTime = tim.usertime()/iter;
		
    tim.clear();chrono.clear();
    for (int i=0;i<iter;i++){
      chrono.start();
      FFLAS::fsyrk(F,FFLAS::FflasUpper,FFLAS::FflasNoTrans,n,k,alpha,A,lda,beta,C,ldc,n-1);
      chrono.stop();
      tim+=chrono;
      FFLAS::fassign (F, n, n, B, ldb, C, ldc);
    }
    RecTime = tim.realtime()/iter;

    cout << "      ";
    cout.width(4);
    cout << n;
    cout << "  ";
    cout.width(15);
    cout << BCTime;
    cout << "  ";
    cout.width(15);
    cout << GFOPS(n, BCTime) << "  ";
    cout.width(15);
    cout << RecTime;
    cout << "  ";
    cout.width(15);
    cout << GFOPS(n, RecTime) << endl;

    if (BCTime > RecTime){
      count++;
      if (count > 2){
	nbest = n;
	bound = true;
	prec = prec >> 1;
	n -= prec;
      }
    }
    else{
      count=0;
      if (bound)
	prec=prec>>1;
      n+=prec;
    }
  } while ((prec > 1 ) && (n < nmax));

  cout<<endl;
  if (nbest != 0 ) {
    cerr << "#ifndef __FFLASFFPACK_FSYRK_THRESHOLD"  << endl;
    cerr << "#define __FFLASFFPACK_FSYRK_THRESHOLD" << ' ' <<  nbest << endl;
    cout << "defined __FFLASFFPACK_FSYRK_THRESHOLD to " << nbest << "" << std::endl;
    std::cerr << "#endif" << endl  << endl;
  }
  FFLAS::fflas_delete(A);
  FFLAS::fflas_delete(B);
  FFLAS::fflas_delete(C);
	
  return 0;
}
