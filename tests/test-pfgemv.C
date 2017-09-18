#define NEWWINO
#ifndef TIME
#define TIME 1
#endif

#define __FFLASFFPACK_DEBUG 1
#include <iomanip>
#include <iostream>
using namespace std;

#define  __FFLASFFPACK_USE_OPENMP
//#define  __FFLASFFPACK_USE_KAAPI

//#define __FFLASFFPACK_FORCE_SEQ

#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "givaro/modular.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "time.h"

/*
#ifdef __FFLASFFPACK_USE_KAAPI
#include <kaapi++>
#endif

#ifdef __FFLASFFPACK_USE_OPENMP
#include <omp.h>
#endif
*/
#include "fflas-ffpack/paladin/pfgemv.inl"

using namespace FFPACK;
using namespace FFLAS;

typedef Givaro::Modular<double> Field;
//typedef Givaro::Modular<float> Field;
//typedef ModularBalanced<double> Field;
//typedef ModularBalanced<float> Field;
//typedef Givaro::Modular<int> Field;

int main(int argc, char** argv) 
//BEGIN_PARALLEL_MAIN(int argc, char** argv)
{
  typename Field::Element_ptr A, X, Y, Y2;
  Field F(17);
  size_t lda,incX,incY,m = 7;
  lda=m;  
  incX=1;  
  incY=1;
  A  = FFLAS::fflas_new(F,m,lda);
  
  Y  = FFLAS::fflas_new(F,m,incY);
  Y2 = FFLAS::fflas_new(F,m,incY);
  X  = FFLAS::fflas_new(F,m,incX);
  
  RandomMatrix (F, m, m, A, lda);
  RandomMatrix (F, m, incX, X,  incX);
  RandomMatrix (F, m, incY, Y,  incY);
  RandomMatrix (F, m, incY, Y2, incY);
  FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y,  incY);
  
  
  PAR_BLOCK
    {
      MMHelper<Field, MMHelperAlgo::Classic, ModeTraits<Field>::value, ParSeqHelper::Parallel<CuttingStrategy::Recursive,StrategyParameter::Threads> >  H;
      FFLAS::pfgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y2,  incY, H);
    }
  
  if (FFLAS::fequal (F, m, 1, Y2, incY, Y, incY)){
    
    cout << "The results are identical"<<endl;
    
  } else{
    
    cout << "The results are identical"<<endl;
    
  }

  PAR_BLOCK
    {
      MMHelper<Field, MMHelperAlgo::Classic, ModeTraits<Field>::value, ParSeqHelper::Parallel<CuttingStrategy::Row,StrategyParameter::Threads> >  H;
      FFLAS::pfgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y2,  incY, H);
    }
  
  if (FFLAS::fequal (F, m, 1, Y2, incY, Y, incY)){
    
    cout << "The results are identical"<<endl;
    
  } else{
    
    cout << "The results are identical"<<endl;
    
  }
size_t BS = 10;
  PAR_BLOCK
    {
      MMHelper<Field, MMHelperAlgo::Classic, ModeTraits<Field>::value, ParSeqHelper::Parallel<CuttingStrategy::Row,StrategyParameter::Grain> >  H;
      FFLAS::pfgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, X, incX, F.zero, Y2,  incY, BS, H);
    }
  
  if (FFLAS::fequal (F, m, 1, Y2, incY, Y, incY)){
    
    cout << "The results are identical"<<endl;
    
  } else{
    
    cout << "The results are identical"<<endl;
    
  }

  
}
//END_PARALLEL_MAIN()


