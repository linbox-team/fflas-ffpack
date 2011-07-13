/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
#include <iostream>
#include <vector>


#ifndef __FFLASFFPACK_HAVE_DGETRF 
#define __FFLASFFPACK_HAVE_DGETRF 1
#endif
#include "fflas-ffpack/fflas.h"
#include "fflas-ffpack/modular-balanced.h"
#include "timer.h"
#include "Matio.h"

using namespace std;

int main(int argc, char** argv) {

  // parameter: p, n, iteration, file

  int    p    = atoi(argv[1]);
  int n    = atoi(argv[2]);
  size_t iter = atoi(argv[3]);
  

  typedef Modular<double> Field;
  typedef Field::Element Element;
      
  Field F(p);
  Element * A;

  Timer chrono;
  double time=0.0;
  int singular;
  std::vector<int> Piv(n,0);
  for (size_t i=0;i<iter;++i){
    if (argc > 4){
      A = read_field (F, argv[4], &n,&n);
    }
    else{
      A = new Element[n*n];
      Field::RandIter G(F);
      for (size_t i=0; i<n*n; ++i)
	G.random(*(A+i));
    }
 
    chrono.clear();
    chrono.start();
    clapack_dgetrf(CblasRowMajor,n,n,A,n,&Piv[0]);
    chrono.stop();
 
    time+=chrono.usertime();
 
  }

  cerr<<"n: "<<n<<" p: "<<p<<" time: "<<time/iter<<endl;


  return 0;
}


