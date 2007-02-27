/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
//--------------------------------------------------------------------------
//          Test for the lsp factorisation
//--------------------------------------------------------------------------
// usage: test-lsp p A n, for n lsp factorization  
// of A over Z/pZ
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
#define DEBUG 1
// Debug option  0: no debug
//               1: check A = LQUP 
//-------------------------------------------------------------------------
using namespace std;


//#define __LUDIVINE_CUTOFF 1
#include <iostream>
#include <iomanip>
#include "Matio.h"
#include "timer.h"
#include "fflas-ffpack/modular-balanced.h"
#include "fflas-ffpack/ffpack.h"

typedef Modular<double> Field;

int main(int argc, char** argv){
	cerr<<setprecision(20);
	int i,j,m,n,nbf,R;

	if (argc!=5){
		cerr<<"usage : test-lqup <p> <A> <c> <i>"<<endl
		    <<"        to do i LQUP factorisation of A with cutoff criterium set to c"
		    <<endl;
		exit(-1);
	}
	Field F(atoi(argv[1]));
	Field::Element * A;
	
	A = read_field(F,argv[2],&m,&n);
	
	size_t *P = new size_t[n];
	size_t *Q = new size_t[m];
		
	size_t cutoff = atoi(argv[3]);
	nbf = atoi(argv[4]);
	
	Timer tim,timc;
	timc.clear();

	enum FFLAS::FFLAS_DIAG diag = FFLAS::FflasNonUnit;
	
	for ( i=0;i<nbf;i++){
		if (i) {		
			delete[] A;
			A = read_field(F,argv[2],&m,&n);
		}
		for (j=0;j<n;j++)
			P[j]=0;
		for (j=0;j<m;j++)
			Q[j]=0;
		tim.clear();      
		tim.start(); 	
		R = FFPACK::LUdivine (F, diag, m, n, A, n, P, Q,
					      FFPACK::FfpackLQUP, cutoff);
		tim.stop();
		timc+=tim;
	}
	//write_field (F,cerr<<"Result = "<<endl, A, m,n,n);
	
#if DEBUG
	Field::Element * L = new Field::Element[m*m];
	Field::Element * U = new Field::Element[m*n];
	Field::Element * X = new Field::Element[m*n];
	
	Field::Element zero,one;
	F.init(zero,0.0);
	F.init(one,1.0);
	for (size_t i=0; i<R; ++i){
		for (size_t j=0; j<i; ++j)
			F.assign ( *(U + i*n + j), zero);
		for (size_t j=i; j<n; ++j)
			F.assign (*(U + i*n + j), *(A+ i*n+j));
	}
	for ( size_t i=0; i<m; ++i ){
		size_t j=0;
		for (; j< ((i<n)?i:n) ; ++j )
			F.assign( *(L + i*m+j), *(A+i*n+j));
		for (; j<m; ++j )
			F.assign( *(L+i*m+j), zero);
	}

	FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans, m,0,m, L, m, Q);
	if (diag == FFLAS::FflasNonUnit)
		for ( size_t i=0; i<m; ++i )
			F.assign (*(L+i*(m+1)),one);
	else{
		size_t i=0;
		while (Q[i]==i){*(L+i*(m+1)) = * (U+i*(n+1)); ++i;}
		for ( size_t i=0; i<R; ++i )
			F.assign (*(U+i*(n+1)),one);
	}
	FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, m,0,R, U, n, P);
	FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans, n,0,R, U, n, Q);
	FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,m, 1.0, L,m, U,n, 0.0, X,n);
	//delete[] A;

	Field::Element * B =  read_field(F,argv[2],&m,&n);
	bool fail = false;
	for (size_t i=0; i<m; ++i)
		for (size_t j=0; j<n; ++j)
			if (!F.areEqual (*(B+i*n+j), *(X+i*n+j)))
				fail=true;
	
	if (fail)
		cerr<<"FAIL"<<endl;


	else
		cerr<<"PASS"<<endl;
		
// 	cout<<m<<" "<<n<<" M"<<endl;
// 	for (size_t i=0; i<m; ++i)
// 		for (size_t j=0; j<n; ++j)
// 			if (!F.isZero(*(A+i*n+j)))
// 				cout<<i+1<<" "<<j+1<<" "<<(*(A+i*n+j))<<endl;
// 	cout<<"0 0 0"<<endl;

#endif
	delete[] A;
	delete[] P;
	delete[] Q;
	delete[] U;
	delete[] L;
	delete[] X;
	double t = timc.usertime();
	double numops = m*m/1000.0*(n-m/3.0);
	
	cerr<<m<<"x"<< n 
	    << " : rank = " << R << "  ["
	    << ((double)nbf/1000.0*(double)numops / t) 
	    << " MFops "
	    << " in "
	    << t/nbf<<"s"
	    <<"]"<< endl;
// 	cout<<m
// 	    <<" "<<((double)nbf/1000.0*(double)numops / t) 
// 	    <<" "<<t/nbf
// 	    <<endl;

	return 0;
}
