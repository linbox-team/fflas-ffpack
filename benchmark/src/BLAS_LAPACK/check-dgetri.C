#include <iostream>



#include "linbox-config.h"
#undef __LINBOX_HAVE_NTL
#undef __LINBOX_HAVE_GIVARO
#ifndef __LINBOX_HAVE_DGETRI 
#define __LINBOX_HAVE_DGETRI 
#endif
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/field/modular.h"
#include "linbox/util/timer.h"

using namespace LinBox;
using namespace std;

int main(int argc, char** argv) {

  // parameter: p, n, iteration, file

  int    p    = atoi(argv[1]);
  size_t n    = atoi(argv[2]);
  size_t iter = atoi(argv[3]);
  

  typedef Modular<double> Field;
  typedef Field::Element Element;
      
  Field F(p);
  
  LinBox::Timer chrono;
  double time=0.0;
  int singular;
  std::vector<int> Piv(n,0);
  BlasMatrix<Element> A (n,n);

  for (size_t i=0;i<iter;++i){
     if (argc > 4){
      fstream osA(argv[4]);    
      A.read(osA,F);
      osA.close();
    }
    else{
      Field::RandIter G(F);
      BlasMatrix<Element>::RawIterator it = A.rawBegin();
      for (; it != A.rawEnd();++it)
	G.random(*it);      
    }
 
    chrono.clear();
    chrono.start();
    clapack_dgetrf(CblasRowMajor,n,n,A.getPointer(),n,&Piv[0]);
    clapack_dgetri(CblasRowMajor,n,A.getPointer(),n,&Piv[0]);
    chrono.stop();
 
    time+=chrono.usertime();
 
  }
  
  cerr<<"n: "<<n<<" p: "<<p<<" time: "<<time/iter<<endl;

  
  return 0;
}


