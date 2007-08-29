#include <iostream>

#include "linbox-config.h"
#undef __LINBOX_HAVE_NTL
#undef __LINBOX_HAVE_GIVARO
#define __LINBOX_HAVE_DTRTRI 1
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/field/modular.h"
#include "linbox/util/timer.h"

//extern "C" {
//#include "atlas_lapack.h"
//#include "atlas_enum.h"
//#include "clapack.h"
//  int ATL_dtrtriRL(const enum ATLAS_DIAG Diag, const int N, double *A, const int lda);
//}

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
 
    for (size_t k=0;k<n;++k)    
      A.setEntry(k,k,1);
 
    chrono.clear();
    chrono.start();       
    //ATL_dtrtriRL(CblasUnit,n,A.getWritePointer(),n);
    clapack_dtrtri(CblasRowMajor,CblasLower, CblasUnit,n,A.getWritePointer(),n);
    chrono.stop();
 
    time+=chrono.usertime();
 
  }

  cerr<<"n: "<<n<<" p: "<<p<<" time: "<<time/iter<<endl;


  return 0;
}
