/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
//--------------------------------------------------------------------------
//          Test for the lqup factorisation
//--------------------------------------------------------------------------
// usage: test-lqup p A n, for n lqup factorization  
// of A over Z/pZ
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
#define DEBUG 1
#define __FFLAS__TRSM_READONLY
// Debug option  0: no debug
//               1: check A = LQUP 
//-------------------------------------------------------------------------


#define __FFPACK_LUDIVINE_CUTOFF 60
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "utils/Matio.h"
#include "utils/timer.h"
#include "fflas-ffpack/field/modular-positive.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "test-utils.h"

using namespace std;
using namespace FFPACK;

typedef Modular<double> Field;
//typedef Modular<float> Field;

int main(int argc, char** argv){
	    //cerr<<setprecision(20);
	int m,n;
	size_t R;

	if (argc!=4){
		cerr<<"usage : test-lqup <p> <A>  <i>"<<endl
		    <<"        to do i LQUP factorisation of A"
		    <<endl;
		exit(-1);
	}
	Field F(atof(argv[1]));
	Field::Element * A;
	
	A = read_field(F,argv[2],&m,&n);
	
	size_t maxP, maxQ;
			
	//	size_t cutoff = atoi(argv[3]);
	size_t nbf = atoi(argv[3]);
	
	Timer tim,timc;
	timc.clear();

	enum FFLAS::FFLAS_DIAG diag = FFLAS::FflasNonUnit;
	enum FFLAS::FFLAS_TRANSPOSE trans = FFLAS::FflasNoTrans;
	if (trans == FFLAS::FflasNoTrans){
		maxP = m;
		maxQ = n;
	} else{
		maxP = n;
		maxQ = m;
	}
	size_t *P = new size_t[maxP];
	size_t *Q = new size_t[maxQ];
	
	//write_field (F,cerr<<"A = "<<endl, A, m,n,n);
	size_t * RRP, *CRP;
	for ( size_t i=0;i<nbf;i++){
		if (i) {		
			delete[] A;
			delete[] RRP;
			delete[] CRP;
			A = read_field(F,argv[2],&m,&n);
		}

		for (size_t j=0;j<maxP;j++)
			P[j]=0;
		for (size_t j=0;j<maxQ;j++)
			Q[j]=0;
		tim.clear();      
		tim.start(); 	

		R = FFPACK::PLUQ/*_basecaseCrout*/ (F, diag, m, n, A, n, P, Q);
//		delete[] A;
//		A = read_field(F,argv[2],&m,&n);
		    //	R = FFPACK::LUdivine (F, diag, FFLAS::FflasNoTrans, m, n, A, n, P, Q, FFPACK::FfpackLQUP);
//		std::cerr<<"Fini LUdivine"<<std::endl;
		tim.stop();
		timc+=tim;
		RRP = new size_t[R];
		CRP = new size_t[R];
		RankProfilesFromPLUQ(RRP, CRP, P, Q, m, n, R);
	}
	    // cerr<<"Row Rank Profile = ";
	// for (size_t i=0;i<R;++i)
	// 	cerr<<RRP[i]<<" ";
	// cerr<<endl;
	// cerr<<"Column Rank Profile = ";
	// for (size_t i=0;i<R;++i)
	// 	cerr<<CRP[i]<<" ";
	// cerr<<endl;
	// std::sort(CRP,(CRP+R));
	// std::sort(RRP,(RRP+R));
	// cerr<<"Sorted Row Rank Profile = ";
	// for (size_t i=0;i<R;++i)
	// 	cerr<<RRP[i]<<" ";
	// cerr<<endl;
	// cerr<<"Sorted Column Rank Profile = ";
	// for (size_t i=0;i<R;++i)
	// 	cerr<<CRP[i]<<" ";
	// cerr<<endl;
	
	if (nbf){
		delete[] RRP;
		delete[] CRP;
	}
//	write_field (F,cerr<<"Result = "<<endl, A, m,n,n);

// 	cerr<<"P = [";
// 	for (size_t i=0; i<maxP; ++i)
// 		cerr<<P[i]<<" ";
// 	cerr<<"]"<<endl;
	// cerr<<"Q = [";
	// for (size_t i=0; i<maxQ; ++i)
	// 	cerr<<Q[i]<<" ";
	// cerr<<"]"<<endl;
#if DEBUG
	Field::Element * X = new Field::Element[m*n];
	Field::Element * L, *U;
	L = new Field::Element[m*R];
	U = new Field::Element[R*n];
	
	Field::Element zero,one;
	F.init(zero,0.0);
	F.init(one,1.0);
	for (size_t  i=0; i<R; ++i){
		for (size_t j=0; j<i; ++j)
			F.assign ( *(U + i*n + j), zero);
		for (int j=i; j<n; ++j)
			F.assign (*(U + i*n + j), *(A+ i*n+j));
	}
	for ( size_t j=0; j<R; ++j ){
		for (size_t i=0; i<=j; ++i )
			F.assign( *(L+i*R+j), zero);
		F.assign(*(L+j*R+j), one);
		for (size_t i=j+1; i<(size_t)m; i++)
			F.assign( *(L + i*R+j), *(A+i*n+j));
	}
	
	    //write_field(F,cerr<<"L = "<<endl,L,m,R,R);
	    //write_field(F,cerr<<"U = "<<endl,U,R,n,n);
	// cerr<<endl;
	FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, R,0,m, L, R, P);
	
	    //write_field(F,cerr<<"L = "<<endl,L,m,m,m);
	    //write_field(F,cerr<<"U = "<<endl,U,m,n,n);
// 		write_field(F,cerr<<"L = "<<endl,L,m,m,m);
// 		write_field(F,cerr<<"U = "<<endl,U,m,n,n);
	
	FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, R,0,n, U, n, Q);
	FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,R, 
		      1.0, L,R, U,n, 0.0, X,n);
	    //delete[] A;
	
//////
	    //write_field(F,cerr<<"L = "<<endl,L,m ,n,n);
	    //write_field(F,cerr<<"U = "<<endl,U,n,n,n);
	
	 // cerr<<"P = ";
	 // for (int i=0; i<m; ++i)
	 // 	cerr<<P[i]<<" ";
	 // cerr<<endl;
	 // cerr<<"Q = ";
	 // for (int i=0; i<n; ++i)
	 // 	cerr<<Q[i]<<" ";
	 // cerr<<endl;
	
	Field::Element * B =  read_field(F,argv[2],&m,&n);

	bool fail = false;
	for (size_t i=0; i<(size_t)m; ++i)
		for (size_t j=0; j<(size_t)n; ++j)
			if (!F.areEqual (*(B+i*n+j), *(X+i*n+j))){
				std::cerr << " B["<<i<<","<<j<<"] = " << (*(B+i*n+j))
					  << " X["<<i<<","<<j<<"] = " << (*(X+i*n+j))
					  << endl;
				fail=true;
			}
	// write_field(F,cerr<<"X = "<<endl,X,m,n,n);
	// write_field(F,cerr<<"B = "<<endl,B,m,n,n);
	delete[] B;
	if (fail)
		cerr<<"FAIL"<<endl;


	else
		cerr<<"PASS"<<endl;
	delete[] U;
	delete[] L;
	delete[] X;
#endif
	delete[] A;
	delete[] P;
	delete[] Q;
	
	double t = timc.realtime();
    const int sm = MIN(m,n);
    const int sn = MAX(m,n);
    
	double numops = sm*sm/1000.0*(sn-sm/3.0);
	
	// cerr<<m<<"x"<< n
	//     << " Trans = "<<trans
	//     << " Diag = "<<diag
	//     << " : rank = " << R << "  ["
	//     << ((double)nbf/1000.0*(double)numops / t) 
	//     << " MFops "
	//     << " in "
	//     << t/nbf<<"s"
	//     <<"]"<< endl;
	cerr<<m
	    <<" "<<((double)nbf/1000.0*(double)numops / t) 
	    <<" "<<t/nbf
	    <<" "<<R
	    <<endl;

	return 0;
}
