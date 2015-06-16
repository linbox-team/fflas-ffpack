/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
//#include "goto-def.h"

/* Copyright (c) FFLAS-FFPACK
* Written by Ziad Sultan <ziad.sultan@imag.fr>
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
*/

//#define __FFLASFFPACK_USE_OPENMP
//#define __FFLASFFPACK_USE_TBB

//#define __FFLASFFPACK_USE_DATAFLOW

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/field/nonzero-randiter.h"

#include "tests/test-utils.h"

#ifdef __FFLASFFPACK_USE_KAAPI
#include "libkomp.h"
#endif

using namespace std;


// random generator function:                                                                                          
ptrdiff_t myrandom (ptrdiff_t i) { return rand()%i;}

// pointer object to it:                                                                                               
ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;

template<class Field>
void matrixWithRandRPM (const Field& F, typename Field::Element_ptr A, size_t lda, size_t M, size_t N, size_t R, size_t seed){
	    // generate the r pivots in the rank profile matrix E
	size_t curr = 0;
	std::vector<bool> rows(M,false);
	std::vector<bool> cols(N,false);
	int pivot_r[R];
	int pivot_c[R];
	typedef typename Field::RandIter Randiter ;
	Randiter RI(F);
	FFPACK::NonzeroRandIter<Field,Randiter> nzR(F,RI);
	while (curr<R){
		int i,j;
		while (rows [i = rand() % M]);
		while (cols [j = rand() % N]);
		rows[i] = true;
		cols[i] = true;
		pivot_r[curr] = i;
		pivot_c[curr] = j;
		curr++;
	}
	typename Field::Element_ptr L= FFLAS::fflas_new(F,M,N);
	for (size_t k = 0; k < R; ++k){
		size_t i = pivot_r[k];
		size_t j = pivot_c[k];
		if (!cols [j])
			FFLAS::fzero(F, M, L+j, N);
		else{
			nzR.random (L [i*N+j]);
			for (size_t l=i+1; l < M; ++l)
				RI.random (L [l*N+j]);
		}
	}
	typename Field::Element_ptr U= FFLAS::fflas_new(F,N,N);
	for (size_t i = 0; i < N; ++i){
		nzR.random (U [i*N+i]);
		for (size_t j=i+1; j < N; ++j)
			RI.random (U [i*N+j]);
	}
	
	
	typename Field::Element alpha, beta;
	F.init(alpha,1.0);
	F.init(beta,0.0);
	PAR_BLOCK{
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M,N,N, alpha, L, N, U, N, beta, A, lda, FFLAS::ParSeqHelper::Parallel());
	}
	FFLAS::fflas_delete(L);
	FFLAS::fflas_delete(U);
	
}

//typedef Givaro::ModularBalanced<double> Field;
//typedef Givaro::ZRing<double> Field;
typedef Givaro::ZRing<double> Field;
//typedef Givaro::UnparametricZRing<double> Field;

void verification_PLUQ(const Field & F, typename Field::Element * B, typename Field::Element * A,
		       size_t * P, size_t * Q, size_t m, size_t n, size_t R)
{

	FFLAS::ParSeqHelper::Parallel H;

	Field::Element * X = FFLAS::fflas_new (F, m,n);
	Field::Element * L, *U;
	L = FFLAS::fflas_new(F, m,R);
	U = FFLAS::fflas_new(F, R,n);
	size_t i;
	
	PARFOR1D (i,0, m*R,H, F.init(L[i], 0.0); );
	
	PARFOR1D (i,0,n*R,H, F.init(U[i], 0.0); );
	
	PARFOR1D (i,0,m*n,H, F.init(X[i], 0.0); );	
	
	Field::Element zero,one;
	F.init(zero,0.0);
	F.init(one,1.0);
	PARFOR1D (i,0,R,H,
              for (size_t j=0; j<i; ++j)
              	F.assign ( *(U + i*n + j), zero);
              for (size_t j=i; j<n; ++j)
              	F.assign (*(U + i*n + j), *(A+ i*n+j));
              );
	size_t j;
	
	PARFOR1D (j,0,R,H, 
		  for (i=0; i<=j; ++i )
			  F.assign( *(L+i*R+j), zero);
		  F.assign(*(L+j*R+j), one);
		  for (i=j+1; i<m; i++)
			  F.assign( *(L + i*R+j), *(A+i*n+j));
		  );
	
	PAR_BLOCK{
		SYNCH_GROUP(MAX_THREADS,
		
		//#pragma omp task shared(F, P, L)
		TASK(MODE(CONSTREFERENCE(F,P,L)),
		     FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, R,0,m, L, R, P););
		//#pragma omp task shared(F, Q, U)
		TASK(MODE(CONSTREFERENCE(F,Q,U)),
		     FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, R,0,n, U, n, Q););
		WAIT;
		//#pragma omp taskwait
		const FFLAS::CuttingStrategy method = FFLAS::THREE_D;
		typename FFLAS::ParSeqHelper::Parallel pWH (MAX_THREADS, method);
		//#pragma omp task shared(F, L, U, X)
		TASK(MODE(CONSTREFERENCE(F,U,L,X)),
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,R,
			      F.one, L,R, U,n, F.zero, X,n, pWH););
			    );
	}
	bool fail = false;
	//  PAR_FOR (size_t i=0; i<m; ++i)
	for(i=0; i<m; ++i)
		for (j=0; j<n; ++j)
			if (!F.areEqual (*(B+i*n+j), *(X+i*n+j))){
				std::cout << " Initial["<<i<<","<<j<<"] = " << (*(B+i*n+j))
					  << " Result["<<i<<","<<j<<"] = " << (*(X+i*n+j))
					  << std::endl;

				std::stringstream errs;
				errs << " B["<<i<<","<<j<<"] = " << (*(B+i*n+j))
				     << " X["<<i<<","<<j<<"] = " << (*(X+i*n+j))
				     << std::endl;
				std::cout << errs;
				fail=true;
			}
	
	if (fail)
		std::cout<<"FAIL"<<std::endl;
	else
		std::cout<<"PASS"<<std::endl;
	
	FFLAS::fflas_delete( U);
	FFLAS::fflas_delete( L);
	FFLAS::fflas_delete( X);
}




template<class Element>
void Initialize(Field &F, Element * C, int BS, size_t m, size_t n)
{

	Field::RandIter G(F); 
//#pragma omp parallel for collapse(2) schedule(runtime) 
	SYNCH_GROUP(MAX_THREADS, {
			for(size_t p=0; p<m; p+=BS) ///row
				for(size_t pp=0; pp<n; pp+=BS) //column
			{
				size_t M=BS, MM=BS;
				if(!(p+BS<m))
					M=m-p;
				if(!(pp+BS<n))
					MM=n-pp;
				//#pragma omp task
				TASK(MODE(CONSTREFERENCE(G)),
				{
					for(size_t j=0; j<M; j++)
						for(size_t jj=0; jj<MM; jj++)
							//		C[(p+j)*n+pp+jj]=0;
							G.random (*(C+(p+j)*n+pp+jj));
				});
			}
		    });
	
	//		#pragma omp taskwait
	//	}
	// printf("A = \n");
	// for (size_t i=0; i<m; i+=128)
	//  {
	//  	for (size_t j=0; j<n; j+=128)
	//  	{
	//  		int ld = omp_get_locality_domain_num_for( &C[i*n+j] );
	//  		printf("%i ", ld);
	//  	}
	//  	printf("\n");
	//  }
	
}


int main(int argc, char** argv) {
	
	size_t iter = 3 ;
	int q = 131071 ;
	Field F(q);
	int m = 2000 ;
	int n = 2000 ;
	int r = 2000 ;
	int v = 0;
	//	int p=0;
	int t=MAX_THREADS;
	int NBK = -1;
	bool par=false;
	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INT , &q },
		{ 'm', "-m M", "Set the row dimension of A.",      TYPE_INT , &m },
		{ 'n', "-n N", "Set the col dimension of A.",      TYPE_INT , &n },
		{ 'r', "-r R", "Set the rank of matrix A.",            TYPE_INT , &r },
		{ 'i', "-i I", "Set number of repetitions.",            TYPE_INT , &iter },
		{ 'v', "-v V", "Set 1 if need verification of result else 0.",            TYPE_INT , &v },
		{ 't', "-t T", "number of virtual threads to drive the partition.", TYPE_INT , &t },
		{ 'b', "-b B", "number of numa blocks per dimension for the numa placement", TYPE_INT , &NBK },
		{ 'p', "-p P", "whether to run or not the parallel PLUQ", TYPE_BOOL , &par },
		END_OF_ARGUMENTS
	};
	//		{ 'p', "-p P", "0 for sequential, 1 for 2D iterative,
//2 for 2D rec, 3 for 2D rec adaptive, 4 for 3D rc in-place, 5 for 3D rec, 6 for 3D rec adaptive.", TYPE_INT , &p },
        FFLAS::parseArguments(argc,argv,as);

	if (r > std::min(m,n)){
		std::cerr<<"Warning: rank can not be greater than min (m,n). It has been forced to min (m,n)"<<std::endl;
		r=std::min(m,n);
	}
	if (!par){
		t=1;NBK=1;
	}
	if (NBK==-1) NBK = t;
	typedef Field::Element Element;
	Element * A, * Acop;
	// A = FFLAS::fflas_new(F,m,n,Alignment::CACHE_PAGESIZE);
	// Acop = FFLAS::fflas_new(F,m,n,Alignment::CACHE_PAGESIZE);
	A = FFLAS::fflas_new(F,m,n);

	// random seed                                                                                         
        ifstream f("/dev/urandom");
        size_t seed1, seed2, seed3,seed4;
        f.read(reinterpret_cast<char*>(&seed1), sizeof(seed1));
        f.read(reinterpret_cast<char*>(&seed2), sizeof(seed2));
        f.read(reinterpret_cast<char*>(&seed3), sizeof(seed3));
        f.read(reinterpret_cast<char*>(&seed4), sizeof(seed4));
	std::vector<size_t> Index_P(r);
	Field::RandIter GG(F, seed1);

       Initialize(F,A,m/NBK,m,n);
       
       matrixWithRandRPM(F, A, n, m, n, r, seed1);
       // //       std::cout<<"Construct U"<<endl;
       // construct_U(F, U, GG, n, r, Index_P, seed4, seed3);
       // //       std::cout<<"Construct L"<<endl;
       // A = construct_L(F,GG, m, r, Index_P, seed2);
       // //       std::cout<<"randgen"<<endl;
       // A = M_randgen(F, A, U, r, m, n);
       
       size_t R;
       FFLAS::Timer chrono;
       double *time=new double[iter];
       
       enum FFLAS::FFLAS_DIAG diag = FFLAS::FflasNonUnit;
       size_t maxP, maxQ;
       maxP = m;
       maxQ = n;
       
       size_t *P = FFLAS::fflas_new<size_t>(maxP);
       size_t *Q = FFLAS::fflas_new<size_t>(maxQ);
       
       FFLAS::ParSeqHelper::Parallel H;
       size_t i;
       
       Acop = FFLAS::fflas_new(F,m,n);
       PARFOR1D(i,0,(size_t)m,H, 
                for (size_t j=0; j<(size_t)n; ++j)
                	Acop[i*n+j]= A[i*n+j];
                );
       size_t j;
       size_t k;
              
       for (i=0;i<=iter;++i){
	       	       
	       PARFOR1D(j,0,maxP,H, P[j]=0; );
	       PARFOR1D(j,0,maxQ,H, Q[j]=0; );
	       PARFOR1D(k,0,(size_t)m,H,
                    for (j=0; j<(size_t)n; ++j)
			    F.assign( A[k*n+j] , Acop[k*n+j]) ;  
                    );
	       chrono.clear();
	       
	       if (i) chrono.start();
	       if (par){
		       
		       PAR_BLOCK{
			       R = FFPACK::pPLUQ(F, diag, m, n, A, n, P, Q, t);
		       }
	       }
	       else
		       R = FFPACK::PLUQ(F, diag, m, n, A, n, P, Q);
	       if (i) {chrono.stop(); time[i-1]=chrono.realtime();}
	       
       }
       std::sort(time, time+iter);
       double meantime = time[iter/2];
       delete[] time;
	// -----------
	// Standard output for benchmark - Alexis Breust 2014/11/14
	#define CUBE(x) ((x)*(x)*(x))
       double gflop =  2.0/3.0*CUBE(double(r)/1000.0) +2*m/1000.0*n/1000.0*double(r)/1000.0  - double(r)/1000.0*double(r)/1000.0*(m+n)/1000;
	std::cout << "Time: " << meantime
		  << " Gflops: " << gflop / meantime;
	FFLAS::writeCommandString(std::cout, as) << std::endl;
       
       //verification
       if(v)
	       verification_PLUQ(F,Acop,A,P,Q,m,n,R);
       
       FFLAS::fflas_delete (P);
       FFLAS::fflas_delete (Q);
       FFLAS::fflas_delete (A);
       FFLAS::fflas_delete (Acop);
       
       return 0;
	}

