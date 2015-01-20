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
//#define __FFLASFFPACK_USE_DATAFLOW
//#define FICTIF
#include <iostream>
#include <givaro/modular.h>

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/args-parser.h"

#include "tests/test-utils.h"

#ifdef __FFLASFFPACK_USE_KAAPI
#include "libkomp.h"
#endif

using namespace std;

  typedef Givaro::ModularBalanced<double> Field;
//typedef Givaro::UnparametricRing<double> Field;

// random generator function:                                                                                          
ptrdiff_t myrandom (ptrdiff_t i) { return rand()%i;}

// pointer object to it:                                                                                               
ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;

typename Field::Element* construct_U(const Field& F, Field::RandIter& G, size_t n, size_t r, std::vector<size_t>& P, size_t commonseed, size_t seed)
{
	size_t lda = n;
	Field::Element *U=new Field::Element[n*n];

	std::vector<size_t> E(r);
	PAR_FOR(size_t i=0; i<r; ++i) E[i]=i;
    
	srand48(commonseed);
	std::vector<size_t> Z(n);
	PAR_FOR(size_t i=0; i<n; ++i) Z[i]=i;
	P.resize(r);
	for(size_t i=0; i<r; ++i) {
		size_t index=lrand48() % Z.size();
		P[i] = Z[ index ];
		Z.erase(Z.begin()+index);
	}

	for (size_t i=0; i<r; ++i) {
      		while( F.isZero( G.random(U[ P[i]*lda+P[i] ]) ) ) {}
		for(size_t j=P[i]+1;j<n;++j)
			G.random(U[ P[i]*lda+j]);
	}
    
	return U;
}

typename Field::Element* construct_L(const Field& F, Field::RandIter& G, size_t m, size_t r, const std::vector<size_t>& P, size_t seed)
{
	size_t lda = m;
	size_t taille=m*m;
	Field::Element * L= new Field::Element[taille];
	PAR_FOR(size_t i=0;i<taille; ++i) F.init(L[i],F.zero);
	
	std::vector<size_t> E(r);
	PAR_FOR(size_t i=0; i<r; ++i) E[i]=i;
	
	srand48(seed);
	std::vector<size_t> Z(m);
	PAR_FOR(size_t i=0; i<m; ++i) Z[i]=i;
	std::vector<size_t> Q(r);
	for(size_t i=0; i<r; ++i) {
		size_t index=lrand48() % Z.size();
		Q[i] = Z[ index ];
		Z.erase(Z.begin()+index);
	}
	
	for(size_t i=0; i<r; ++i) {
		size_t index=lrand48() % E.size();
		size_t perm = E[ index ];
		
		E.erase(E.begin()+index);
		F.init(L[Q[perm]*lda+P[perm]],F.one);
		for(size_t j=Q[perm]+1;j<m;++j)
			G.random(L[j*lda+P[perm]]);
	}
	return L;
}


typename Field::Element* M_randgen(const Field& F, typename Field::Element* L,typename Field::Element* U, size_t r, size_t m, size_t n)
{
	Field::Element alpha, beta;
	F.init(alpha,1.0);
	F.init(beta,0.0);
	size_t lda = n;
	Field::Element * A = new Field::Element[m*n];

	// Computing produit L * U (ideally should use parallel ftrmm

	/*
	  FFLAS::ftrmm(F,  FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, m,n,1.0, U, lda, L, lda);
	*/

	const FFLAS::CuttingStrategy method = FFLAS::THREE_D;
	typename FFLAS::ParSeqHelper::Parallel pWH (MAX_THREADS, method);
	PAR_REGION{
	FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				m,n,r, alpha, L,r, U,
		      lda,beta,A,lda,pWH);
	}
	return L;
}

void verification_PLUQ(const Field & F, typename Field::Element * B, typename Field::Element * A,
		       size_t * P, size_t * Q, size_t m, size_t n, size_t R)
{

	Field::Element * X = FFLAS::fflas_new<Field::Element>(m*n);
	Field::Element * L, *U;
	L = FFLAS::fflas_new<Field::Element>(m*R);
	U = FFLAS::fflas_new<Field::Element>(R*n);
	
	PAR_FOR (size_t i=0; i<m*R; ++i)
		F.init(L[i], 0.0);
	
	PAR_FOR (size_t i=0; i<n*R; ++i)
		F.init(U[i], 0.0);
	
	PAR_FOR (size_t i=0; i<m*n; ++i)
		F.init(X[i], 0.0);
	
	
	Field::Element zero,one;
	F.init(zero,0.0);
	F.init(one,1.0);
	PAR_FOR (size_t i=0; i<R; ++i){
		for (size_t j=0; j<i; ++j)
			F.assign ( *(U + i*n + j), zero);
		for (size_t j=i; j<n; ++j)
			F.assign (*(U + i*n + j), *(A+ i*n+j));
	}
	PAR_FOR ( size_t j=0; j<R; ++j ){
		for (size_t i=0; i<=j; ++i )
			F.assign( *(L+i*R+j), zero);
		F.assign(*(L+j*R+j), one);
    for (size_t i=j+1; i<m; i++)
	    F.assign( *(L + i*R+j), *(A+i*n+j));
	}
	
	PAR_REGION{
#pragma omp task shared(F, P, L)
		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, R,0,m, L, R, P);
#pragma omp task shared(F, Q, U)
		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, R,0,n, U, n, Q);
#pragma omp taskwait
		const FFLAS::CuttingStrategy method = FFLAS::THREE_D;
		typename FFLAS::ParSeqHelper::Parallel pWH (MAX_THREADS, method);
#pragma omp task shared(F, L, U, X)
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,R,
			      F.one, L,R, U,n, F.zero, X,n, pWH);

	}
	bool fail = false;
	//  PAR_FOR (size_t i=0; i<m; ++i)
	for(size_t i=0; i<m; ++i)
		for (size_t j=0; j<n; ++j)
			if (!F.areEqual (*(B+i*n+j), *(X+i*n+j))){
				std::cout << " Initial["<<i<<","<<j<<"] = " << (*(B+i*n+j))
					  << " Result"<<i<<","<<j<<"] = " << (*(X+i*n+j))
					  << std::endl;

				std::stringstream errs;
				errs << " B["<<i<<","<<j<<"] = " << (*(B+i*n+j))
				     << " X["<<i<<","<<j<<"] = " << (*(X+i*n+j))
				     << std::endl;
				std::cout << errs;
				fail=true;
				std::cout<<" end verification"<<std::endl;
				exit(1);
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
	
	for(size_t p=0; p<m; p+=BS) ///row
		for(size_t pp=0; pp<n; pp+=BS) //column
			{
				size_t M=BS, MM=BS;
				if(!(p+BS<m))
					M=m-p;
				if(!(pp+BS<n))
					MM=n-pp;
				//#pragma omp task 
				{
					for(size_t j=0; j<M; j++)
						for(size_t jj=0; jj<MM; jj++)
							//		C[(p+j)*n+pp+jj]=0;
							G.random (*(C+(p+j)*n+pp+jj));
				}
			}
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
	if (!par){
		t=1;NBK=1;
	}
	if (NBK==-1) NBK = t;
	typedef Field::Element Element;
	Element * A, * Acop;
	A = FFLAS::fflas_new(F,m,n,Alignment::CACHE_PAGESIZE);
	Acop = FFLAS::fflas_new(F,m,n,Alignment::CACHE_PAGESIZE);

	Field::Element * U = new Field::Element[n*n];
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
       
       //       std::cout<<"Construct U"<<endl;
       U = construct_U(F,GG, n, r, Index_P, seed4, seed3);
       //       std::cout<<"Construct L"<<endl;
       A = construct_L(F,GG, m, r, Index_P, seed2);
       //       std::cout<<"randgen"<<endl;
       A = M_randgen(F, A, U, r, m, n);
       size_t R;
       FFLAS::Timer chrono;
       double *time=new double[iter];
       
       enum FFLAS::FFLAS_DIAG diag = FFLAS::FflasNonUnit;
       size_t maxP, maxQ;
       maxP = m;
       maxQ = n;
       
       size_t *P = FFLAS::fflas_new<size_t>(maxP);
       size_t *Q = FFLAS::fflas_new<size_t>(maxQ);
       
       for/*PAR_FOR*/(size_t i=0; i<(size_t)m; ++i)
	       for (size_t j=0; j<(size_t)n; ++j)
		       Acop[i*n+j]= (*(A+i*n+j));
       
       for (size_t i=0;i<=iter;++i){
	       	       
	       /*PAR_FOR*/for(size_t j=0;j<maxP;j++)
		       P[j]=0;
	       /*PAR_FOR*/for(size_t j=0;j<maxQ;j++)
		       Q[j]=0;
	       
	       /*PAR_FOR*/for(size_t k=0; k<(size_t)m; ++k)
		       for (size_t j=0; j<(size_t)n; ++j)
			       *(A+k*n+j) = *(Acop+k*n+j) ;  
	       chrono.clear();
	       
	       if (i) chrono.start();
	       if (par)
		       PAR_REGION{
			       R = FFPACK::pPLUQ(F, diag, m, n, A, n, P, Q, t);
		       }
	       else
		       R = FFPACK::PLUQ(F, diag, m, n, A, n, P, Q);
	       if (i) {chrono.stop(); time[i-1]=chrono.usertime();}
	       
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
       FFLAS::fflas_delete( A);
       FFLAS::fflas_delete( Acop);
       
       return 0;
	}

