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
#include "fflas-ffpack/modular-positive.h"
#include "fflas-ffpack/ffpack.h"

typedef Modular<double> Field;

int main(int argc, char** argv){
	cerr<<setprecision(20);
	int i,j,nbf,m,n;
	int R=0;

	if (argc!=4){
		cerr<<"usage : test-lqup <p> <A>  <i>"<<endl
		    <<"        to do i LQUP factorisation of A"
		    <<endl;
		exit(-1);
	}
	Field F((unsigned long)atoi(argv[1]));
	Field::Element * A;
	
	A = read_field(F,argv[2],&m,&n);
	
	size_t maxP, maxQ;
			
	//	size_t cutoff = atoi(argv[3]);
	nbf = atoi(argv[3]);
	
	Timer tim,timc;
	timc.clear();

	enum FFLAS::FFLAS_DIAG diag = FFLAS::FflasUnit;
	enum FFLAS::FFLAS_TRANSPOSE trans = FFLAS::FflasTrans;
	if (trans == FFLAS::FflasTrans){
		maxP = m;
		maxQ = n;
	} else{
		maxP = n;
		maxQ = m;
	}
	size_t *P = new size_t[maxP];
	size_t *Q = new size_t[maxQ];
		
	//write_field (F,cerr<<"A = "<<endl, A, m,n,n);
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
		R = FFPACK::LUdivine (F, diag, trans, m, n, A, n, P, Q,
				      FFPACK::FfpackLQUP);
		tim.stop();
		timc+=tim;
	}
	write_field (F,cerr<<"Result = "<<endl, A, m,n,n);

	cerr<<"P = [";
	for (size_t i=0; i<maxP; ++i)
		cerr<<P[i]<<" ";
	cerr<<"]"<<endl;
	cerr<<"Q = [";
	for (size_t i=0; i<maxQ; ++i)
		cerr<<Q[i]<<" ";
	cerr<<"]"<<endl;
#if DEBUG
	Field::Element * X = new Field::Element[m*n];
	Field::Element * L, *U;
	if (trans == FFLAS::FflasNoTrans){
		L = new Field::Element[m*m];
		U = new Field::Element[m*n];
				
		Field::Element zero,one;
		F.init(zero,0.0);
		F.init(one,1.0);
		for (int i=0; i<R; ++i){
			for (int j=0; j<i; ++j)
				F.assign ( *(U + i*n + j), zero);
			for (int j=i+1; j<n; ++j)
				F.assign (*(U + i*n + j), *(A+ i*n+j));
		}
		for (int i=R;i<m; ++i)
			for (int j=0; j<n; ++j)
				F.assign(*(U+i*n+j), zero);
		for ( int i=0; i<m; ++i ){
			int j=0;
			for (; j< ((i<R)?i:R) ; ++j )
				F.assign( *(L + i*m+j), *(A+i*n+j));
			for (; j<m; ++j )
				F.assign( *(L+i*m+j), zero);
		}
		
		write_field(F,cerr<<"L = "<<endl,L,m,m,m);
		write_field(F,cerr<<"U = "<<endl,U,m,n,n);
		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans, m,0,R, L, m, Q);
		for ( int i=0; i<m; ++i )
			F.assign(*(L+i*(m+1)), one);

		write_field(F,cerr<<"L = "<<endl,L,m,m,m);
		write_field(F,cerr<<"U = "<<endl,U,m,n,n);
		if (diag == FFLAS::FflasNonUnit)
			for ( int i=0; i<R; ++i )
				F.assign (*(U+i*(n+1)), *(A+i*(n+1)));
			
		else{
			for ( int i=0; i<R; ++i ){
				*(L+Q[i]*(m+1)) = *(A+Q[i]*n+i);
				F.assign (*(U+i*(n+1)),one);
			}
		}
		write_field(F,cerr<<"L = "<<endl,L,m,m,m);
		write_field(F,cerr<<"U = "<<endl,U,m,n,n);

		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, m,0,R, U, n, P);
		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans, n,0,R, U, n, Q);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,m, 1.0, L,m, U,n, 0.0, X,n);
		//delete[] A;
	} else {

		L = new Field::Element[m*n];
		U = new Field::Element[n*n];

		
		Field::Element zero,one;
		F.init(zero,0.0);
		F.init(one,1.0);
		for (int i=0; i<R; ++i){
			for (int j=0; j<i; ++j)
				F.assign ( *(L + i + j*n), zero);
			for (int j=i+1; j<m; ++j)
				F.assign (*(L + i + j*n), *(A+ i+j*n));
		}
		
		for (int i=R;i<n; ++i)
			for (int j=0; j<m; ++j)
				F.assign(*(L+i+j*n), zero);
		for ( int i=0; i<n; ++i ){
			int j=0;
			for (;  j< ((i<R)?i:R) ; ++j )
				F.assign( *(U + i+j*n), *(A+i+j*n));
			for (; j<n; ++j )
				F.assign( *(U+i+j*n), zero);
		}
		write_field(F,cerr<<"L = "<<endl,L,m,n,n);
		write_field(F,cerr<<"U = "<<endl,U,n,n,n);

		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, n,0,R, U, n, Q);


		for (int i=0; i<n; ++i)
			F.assign (*(U+i*(n+1)),one);

		if (diag == FFLAS::FflasNonUnit)
			for ( int i=0; i<R; ++i )
				F.assign (*(L+i*(n+1)), *(A+i*(n+1)));
		else{
			for ( int i=0; i<R; ++i ){
				*(U+Q[i]*(n+1)) = *(A+Q[i]+i*n);
				F.assign (*(L+i*(n+1)),one);
			}
		}
		write_field(F,cerr<<"L = "<<endl,L,m,n,n);
		write_field(F,cerr<<"U = "<<endl,U,n,n,n);

		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans, n,0,R, L, n, P);
		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, m,0,R, L, n, Q);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,n, 1.0, L,n, U,n, 0.0, X,n);
	}
	Field::Element * B =  read_field(F,argv[2],&m,&n);
	bool fail = false;
	for (int i=0; i<m; ++i)
		for (int j=0; j<n; ++j)
			if (!F.areEqual (*(B+i*n+j), *(X+i*n+j)))
				fail=true;
	
 	write_field(F,cerr<<"X = "<<endl,X,m,n,n);
 	write_field(F,cerr<<"B = "<<endl,B,m,n,n);
	delete[] B;
	if (fail)
		cerr<<"FAIL"<<endl;


	else
		cerr<<"PASS"<<endl;
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
