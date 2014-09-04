/*   
*******************************************************
        Parallel PLUQ quad recurisve with OpenMP
*******************************************************

g++ -D__FFLASFFPACK_HAVE_CBLAS -Wall -g -fopenmp -O3 -march=native -mavx -I/home/sultan/soft/fflas-ffpack/ -I/usr/local/soft/givaro-3.7.1/include  test-ppluq.C -L/home/pernet/Logiciels/ATLAS_1TH/lib -lcblas -latlas -L/usr/local/soft/givaro-3.7.1/lib -lgivaro -lm -lrt -Wl,-rpath -Wl,/usr/local/soft/givaro-3.7.1/lib  -o test-ppluq
*/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
//#include "omp.h"

#define __FFLASFFPACK_USE_OPENMP

#define __FFLAS__TRSM_READONLY
#define __PFTRSM_FOR_PLUQ
#include "fflas-ffpack/utils/Matio.h"
//#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "sys/time.h"

//#define BASECASE_K 256

//#include "fflas-ffpack/ffpack/parallel.h"

using namespace std;
using namespace FFLAS;
using namespace FFPACK;
#define MODULO 1

#if(MODULO==1)
typedef FFPACK::Modular<double> Field;
#else
typedef FFPACK::UnparametricField<double> Field;
#endif

#define DEBUG 
#define SEQ 1


typedef FFPACK::Modular<double> Field;

void verification_PLUQ(const Field & F, typename Field::Element * B, typename Field::Element * A,
		       size_t * P, size_t * Q, size_t m, size_t n, size_t R)
{

  Field::Element * X = FFLAS::fflas_new<Field::Element>(m*n);
  Field::Element * L, *U;
  L = FFLAS::fflas_new<Field::Element>(m*R);
  U = FFLAS::fflas_new<Field::Element>(R*n);

  for (size_t i=0; i<m*R; ++i)
    F.init(L[i], 0.0);

  for (size_t i=0; i<m*R; ++i)
    F.init(U[i], 0.0);

  for (size_t i=0; i<m*n; ++i)
    F.init(X[i], 0.0);


  Field::Element zero,one;
  F.init(zero,0.0);
  F.init(one,1.0);
  for (size_t i=0; i<R; ++i){
    for (size_t j=0; j<i; ++j)
      F.assign ( *(U + i*n + j), zero);
    for (size_t j=i; j<n; ++j)
      F.assign (*(U + i*n + j), *(A+ i*n+j));
  }
  for ( size_t j=0; j<R; ++j ){
    for (size_t i=0; i<=j; ++i )
      F.assign( *(L+i*R+j), zero);
    F.assign(*(L+j*R+j), one);
    for (size_t i=j+1; i<m; i++)
      F.assign( *(L + i*R+j), *(A+i*n+j));
  }

  FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, R,0,m, L, R, P);

  FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, R,0,n, U, n, Q);
  FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,R,
		1.0, L,R, U,n, 0.0, X,n);
  bool fail = false;
  for (size_t i=0; i<m; ++i)
    for (size_t j=0; j<n; ++j)
      if (!F.areEqual (*(B+i*n+j), *(X+i*n+j))){
	std::cerr << " B["<<i<<","<<j<<"] = " << (*(B+i*n+j))
		  << " X["<<i<<","<<j<<"] = " << (*(X+i*n+j))
		  << std::endl;
	fail=true;
      }

  if (fail)
    std::cerr<<"FAIL"<<std::endl;


  else
    std::cerr<<"PASS"<<std::endl;
  delete[] U;
  delete[] L;
  delete[] X;
}

int main(int argc, char** argv)
{

    int p, n, m, nbf;

	if (argc!=4){
		std::cerr<<"usage : PLUQ-rec-omp <p> <file>  <i>"<<std::endl
//		std::cerr<<"usage : PLUQ-rec-omp <m> <n> <p> <r> <i>"<<std::endl
		    <<std::endl;
		exit(-1);
	}
   
	p = atoi( argv[1] );

	// m = atoi( argv[1] );
	// n = atoi( argv[2] );
	// r = atoi( argv[4] );
	nbf = atoi( argv[3] );
	
	//	size_t lda = n;

		// random seed
	// ifstream f("/dev/urandom");
	// size_t seed1, seed2, seed3,seed4;
	// f.read(reinterpret_cast<char*>(&seed1), sizeof(seed1));
	// f.read(reinterpret_cast<char*>(&seed2), sizeof(seed2));
	// f.read(reinterpret_cast<char*>(&seed3), sizeof(seed3));
	// f.read(reinterpret_cast<char*>(&seed4), sizeof(seed4));
    
//     seed1=10;seed2=12;
//     seed3=13;seed4=14;
    
        enum FFLAS::FFLAS_DIAG diag = FFLAS::FflasNonUnit;
        size_t R;
#if(MODULO==1)        
	typedef FFPACK::Modular<double> Field;
#else
	typedef FFPACK::UnparametricField<double> Field;
#endif

	const Field F((double)p);
	// Field::RandIter G(F, seed1);
    
	Field::Element alpha, beta;
	F.init(alpha,1.0);
	F.init(beta,0.0);
	// Field::Element * U = FFLAS::fflas_new<Field::Element>(n*n);

	typename Field::Element* A = read_field(F,argv[2],&m,&n);
// FFLAS::fflas_new<Field::Element>(n*m);
	Field::Element* Acop = FFLAS::fflas_new<Field::Element>(n*m);
#ifdef DEBUG 
	Field::Element* Adebug = FFLAS::fflas_new<Field::Element>(n*m);
#endif
	// std::vector<size_t> Index_P(r);

	// U = construct_U(F,G, n, r, Index_P, seed4, seed3);
	// A = construct_L(F,G, m, r, Index_P, seed2);
	// M_randgen(F, A, U, r, m, n);
	// size_t taille=m*n;
	// for(size_t i=0; i<taille;++i) U[i]=A[i];

	struct timespec t0, t1;// tt0, tt1;
	double delay, avrg;//, avrgg;
	double t_total=0;

        size_t maxP, maxQ;
                maxP = m;
                maxQ = n;

        size_t *P = FFLAS::fflas_new<size_t>(maxP);
        size_t *Q = FFLAS::fflas_new<size_t>(maxQ);

	for(int i=0; i<n*m; i++){
	  Acop[i]=A[i];
#ifdef DEBUG 
          Adebug[i]=A[i];
#endif
        }

        for ( int i=0;i<nbf+1;i++){
            for (size_t j=0;j<maxP;j++)
                P[j]=0;
            for (size_t j=0;j<maxQ;j++)
                Q[j]=0;
	    //            #pragma omp parallel for shared(A, Acop)    
	    for(int k=0; k<n*m; k++)
	      {
		
		A[k]=Acop[k];
	      }
	    
	    clock_gettime(CLOCK_REALTIME, &t0);
	    PAR_REGION{
		R = pPLUQ(F, diag, m, n, A, n, P, Q);// Parallel PLUQ
	    }
	    clock_gettime(CLOCK_REALTIME, &t1);
	    delay = (double)(t1.tv_sec-t0.tv_sec)+(double)(t1.tv_nsec-t0.tv_nsec)/1000000000;

	    if(i)
	      t_total +=delay;
	    
        }
        avrg = t_total/nbf;
	std::cerr<<"Parallel : "<<m<<" "<<R<<" "
                 <<avrg<<" "<<(2.0*n*n*n)/(double(3.0*(1000000000)*avrg))<<" "
	  //#ifdef  __FFLASFFPACK_USE_OPENMP
		 <<NUM_THREADS<<endl;
	//#else
	//<<endl;
	//#endi

	//	std::cout<<typeid(A).name()<<endl;
#ifdef DEBUG
	cout<<"check equality A == PLUQ ?"<<endl;
        verification_PLUQ(F,Adebug,A,P,Q,m,n,R);
        delete[] Adebug;
#endif
#if SEQ
	struct timespec  tt0, tt1;
	double avrgg;
	//call sequential PLUQ
	size_t * PP = FFLAS::fflas_new<size_t>(maxP);
	size_t * QQ = FFLAS::fflas_new<size_t>(maxQ);
	for (size_t j=0;j<maxP;j++)
	  PP[j]=0;
	for (size_t j=0;j<maxQ;j++)
	  QQ[j]=0;
	clock_gettime(CLOCK_REALTIME, &tt0);
	size_t R2 = PLUQ(F, diag, m, n, Acop, n, PP, QQ);
	clock_gettime(CLOCK_REALTIME, &tt1);
        delete[] Acop;
	avrgg = (double)(tt1.tv_sec-tt0.tv_sec)+(double)(tt1.tv_nsec-tt0.tv_nsec)/1000000000;
	//verification
	std::cerr<<"Sequential : "<<m<<" "<<R2<<" "
                 <<avrgg<<" "<<(2.0*n*n*n)/(double(3.0*(1000000000)*avrgg))<<endl;
#endif

        delete[] A;
	return 0;
}
