/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
//#include "goto-def.h"

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
  Element one,zero;
  F.init(one, 1.0);
  F.init(zero,0.0);

  Timer chrono;
  double time=0.0, time2=0.0;
  int singular;
  
  Element * A, * B, * C;

  for (size_t i=0;i<iter;++i){
    
    if (argc > 4){
      A = read_field (F, argv[4], &n, &n);    
    }
    else{
      Field::RandIter G(F);
      A = new Element[n*n];
      for (size_t i=0; i<n*n; ++i)
	G.random (*(A+i));      
    }

    if (argc == 6){
      B = read_field (F, argv[5], &n, &n);    
    }
    else{
      Field::RandIter G(F);
      B = new Element[n*n];
      for (size_t i=0; i<n*n; ++i)
	G.random(*(B+i));
    }
    C = new Element[n*n];

    chrono.clear();
    chrono.start();
    FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n,n,n, one,
		  A, n, B, n, zero, C,n);
    chrono.stop();
    time+=chrono.usertime();

    delete[] A;
    delete[] B;
    delete[] C;
  }
  
  std::cerr<<"n: "<<n<<" p: "<<p<<" time: "<<time/iter<<std::endl;

  return 0;
}

