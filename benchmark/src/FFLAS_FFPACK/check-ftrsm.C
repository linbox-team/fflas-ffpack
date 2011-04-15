/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
#include <iostream>

#include "fflas-ffpack/fflas.h"
#include "fflas-ffpack/modular-balanced.h"
#include "timer.h"
#include "Matio.h"

using namespace std;

int main(int argc, char** argv) {

  // parameter: p, n, iteration, file1, file2

  int    p    = atoi(argv[1]);
  int n    = atoi(argv[2]);
  size_t iter = atoi(argv[3]);
  

  typedef Modular<double> Field;
  typedef Field::Element Element;
      
  Field F(p);
  Element one;
  F.init(one, 1.0);
  Element * A;
  Element * B;
  
  Timer chrono;
  double time=0.0;
  int singular;
  
  for (size_t i=0;i<iter;++i){
    Field::RandIter G(F);
    if (argc > 4){
        A = read_field (F, argv[4], &n, &n);    
    }
    else{
      A = new Element[n*n];
      for (size_t i = 0; i< n*n; ++i)
	G.random(*(A+i));      
    }

    if (argc == 6){
      B = read_field (F, argv[5], &n, &n);    
    }
    else{
      B = new Element[n*n];
      Field::RandIter G(F);
      for (size_t i=0 ; i< n*n; ++i)
	G.random(*(A+i));
    }
        
    for (size_t k=0;k<n;++k)
      while (F.isZero( G.random(*(A+k*(n+1)))));
    
    chrono.clear();
    chrono.start();
    FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans,
		  FFLAS::FflasNonUnit, n,n, one, A, n, B, n);
    
    chrono.stop();
    time+=chrono.usertime();
    delete[] A;
    delete[] B;
    
  }

  cerr<<"n: "<<n<<" p: "<<p<<" time: "<<time/iter<<endl;


  return 0;
}
