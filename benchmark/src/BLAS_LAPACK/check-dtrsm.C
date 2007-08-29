#include <iostream>

//#define __LINBOX_CONFIGURATION
//#include <linbox/config-blas.h>

#include "linbox-config.h"
#undef __LINBOX_HAVE_NTL
#undef __LINBOX_HAVE_GIVARO
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/field/modular.h"
#include "linbox/util/timer.h"

using namespace LinBox;
using namespace std;

int main(int argc, char** argv) {

  // parameter: p, n, iteration, file1, file2

  if (argc != 4)
    std::cout<<"usage: <p> <dim> <iter>\n";

  int    p    = atoi(argv[1]);
  size_t n    = atoi(argv[2]);
  size_t iter = atoi(argv[3]);
  

  typedef Modular<double> Field;
  typedef Field::Element Element;
      
  Field F(p);
  BlasMatrixDomain<Field> BMD(F);
  BlasMatrix<Element> A(n,n), B(n,n);

  LinBox::Timer chrono;
  double time=0.0;
  int singular;
  
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


    if (argc == 6){
      fstream osB(argv[5]);      
      B.read(osB,F);
      osB.close();
    }
    else{
      Field::RandIter G(F);
      BlasMatrix<Element>::RawIterator it = B.rawBegin();
      for (; it != B.rawEnd();++it)
	G.random(*it);
    }
  
    chrono.clear();
    chrono.start();
    cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans,
		CblasUnit, n, n, 1.0, A.getPointer(), n, B.getPointer(), n);
    chrono.stop();
    time+=chrono.usertime();
    
  }

  cerr<<"n: "<<n<<" p: "<<p<<" time: "<<time/iter<<endl;


  return 0;
}
