#include <iostream>

#include "fflas-ffpack/ffpack.h"
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
  Field::Element * A;
  
  Timer chrono;
  double time=0.0;
  int singular;
  
  for (size_t i=0;i<iter;++i){
    if (argc > 4){
      A = read_field(F, argv[4],  &n, &n);
    } else {
      A = new Element[n*n];      
      Field::RandIter G(F);
      for (size_t i=0; i<n*n; ++i)
	G.random(*(A+i));

    }
    int nullity=0;
    chrono.clear();
    chrono.start();
    FFPACK::Invert (F, n, A, n, nullity);
    chrono.stop();
    
    time+=chrono.usertime();
    delete[] A;
  }

  cerr<<"n: "<<n<<" p: "<<p<<" time: "<<time/iter<<endl;


  return 0;
}




