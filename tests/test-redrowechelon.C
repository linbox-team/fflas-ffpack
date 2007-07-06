/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
//--------------------------------------------------------------------------
//          Test for the reduced row echelon factorisation
//--------------------------------------------------------------------------
// usage: test-redrowechelon p A n, for n reduced row echelon computations
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
//#include "fflas-ffpack/modular-balanced.h"
#include "fflas-ffpack/modular-positive.h"
#include "fflas-ffpack/ffpack.h"

typedef Modular<double> Field;

int main(int argc, char** argv){
	//cerr<<setprecision(20);
	int i,j,nbf,m,n;
	int R=0;

	if (argc!=4){
		cerr<<"usage : test-redrowechelon <p> <A> <i>"<<endl
		    <<"        to do i reduced row  echelon computations of A"
		    <<endl;
		exit(-1);
	}
	Field F((unsigned long)atoi(argv[1]));
	Field::Element * A;
	
	A = read_field(F,argv[2],&m,&n);
	
	size_t *P = new size_t[n];
	size_t *Q = new size_t[m];
		
	//	size_t cutoff = atoi(argv[3]);
	nbf = atoi(argv[3]);
	
	Timer tim,timc;
	timc.clear();


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
		R = FFPACK::ReducedRowEchelonForm (F, m, n, A, n, P, Q);
		tim.stop();
		timc+=tim;
	}
	//write_field (F,cerr<<"Result = "<<endl, A, m,n,n);

// 	cerr<<"P = [";
// 	for (size_t i=0; i<n; ++i)
// 		cerr<<P[i]<<" ";
// 	cerr<<"]"<<endl;
// 	cerr<<"Q = [";
// 	for (size_t i=0; i<m; ++i)
// 		cerr<<Q[i]<<" ";
// 	cerr<<"]"<<endl;
#if DEBUG
	Field::Element * L = new Field::Element[m*m];
	Field::Element * U = new Field::Element[m*n];
	Field::Element * X = new Field::Element[m*n];
	
	Field::Element zero,one;
	F.init(zero,0.0);
	F.init(one,1.0);
	
	for (int i=0; i<R; ++i){
		for (int j=0; j<m; ++j)
			F.assign (*(L + i + j*m), *(A+ i+j*n));
	}
	for (int i=R;i<m; ++i){
		for (int j=0; j<m; ++j)
			F.assign(*(L+i+j*m), zero);
		F.init(*(L+i*(m+1)),one);
	}
	FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans, m, 0, R, L, m, P);	

	for ( int i=0; i<R; ++i ){
		for (int j=0; j < m ; ++j)
			F.assign( *(U + i+j*n),zero);
		F.assign(*(U+i*(n+1)), one);
	}
	for ( int i=R; i<n; ++i ){
		for (int j=0; j<R; ++j )
			F.assign (*(U+i+j*n), *(A+i+j*n));
		for (int j=R; j<m; ++j)
			F.assign (*(U+i+j*n), zero);
	}
  	//write_field(F,cerr<<"U = "<<endl,U,m,n,n);
	FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans, m, 0, R, U, n, Q);
  	//write_field(F,cerr<<"U = "<<endl,U,m,n,n);
// 	cerr<<"P = ";
// 	for (size_t i=0; i<n;++i)
// 		cerr<<" "<<P[i];
// 	cerr<<endl;
// 	cerr<<"Q = ";
// 	for (size_t i=0; i<m;++i)
// 		cerr<<" "<<Q[i];
// 	cerr<<endl;
	
	
	//  	write_field(F,cerr<<"R = "<<endl,L,m,n,n);


	Field::Element * B =  read_field(F,argv[2],&m,&n);
	//write_field(F,cerr<<"A = "<<endl,B,m,n,n);
	
	FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,m, 1.0,
		      L, m, B, n, 0.0, X,n);
	//delete[] A;

	bool fail = false;
	for (int i=0; i<m; ++i)
		for (int j=0; j<n; ++j)
			if (!F.areEqual (*(U+i*n+j), *(X+i*n+j)))
				fail=true;
	
//  	write_field(F,cerr<<"X = "<<endl,X,m,n,n);
//    	write_field(F,cerr<<"R = "<<endl,U,m,n,n);
	
	delete[] B;
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

	delete[] U;
	delete[] L;
	delete[] X;
#endif
	delete[] A;
	delete[] P;
	delete[] Q;

	double t = timc.usertime();
	double numops = 2*m*m/1000.0*n;
	
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
