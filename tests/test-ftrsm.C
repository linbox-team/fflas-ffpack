/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

//--------------------------------------------------------------------------
//                        Test for ftrsm : 1 computation
//
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------

//#define DEBUG 1
#define TIME 1

#include <iomanip>
#include <iostream>
#include "fflas-ffpack/field/modular-balanced.h"
#include "timer.h"
#include "Matio.h"
#include "fflas-ffpack/fflas/fflas.h"


using namespace std;
using namespace FFPACK;

//typedef ModularBalanced<double> Field;
typedef ModularBalanced<float> Field;

int main(int argc, char** argv)
{

	int k,n,m;
	cerr<<setprecision(10);
	Field::Element zero, one;

	if (argc != 10)	{
		cerr<<"Usage : test-ftrsm <p> <A> <B> <iter> <alpha> <left/right> <Up/Low> <NoTrans/Trans> <NonUnit/Unit>"
		    <<endl;
		exit(-1);
	}
	int nbit=atoi(argv[4]); // number of times the product is performed
	Field F(atof(argv[1]));
	F.init(zero,0.0);
	F.init(one,1.0);
	Field::Element * A, *B, *B2;
	A = read_field(F,argv[2],&k,&k);
	B = read_field(F,argv[3],&m,&n);
	B2 = new Field::Element[m*n];


	for (int i=0; i<m;++i){
		for(int j=0; j<n; ++j)
			F.assign(*(B2+i*n+j),*(B+i*n+j));
	}

	Field::Element alpha;
	F.init (alpha, atof(argv[5]));

	FFLAS::FFLAS_SIDE side = (atoi(argv[6])) ? FFLAS::FflasRight :  FFLAS::FflasLeft;
	FFLAS::FFLAS_UPLO uplo = (atoi(argv[7])) ? FFLAS::FflasLower :  FFLAS::FflasUpper;
	FFLAS::FFLAS_TRANSPOSE trans = (atoi(argv[8])) ? FFLAS::FflasTrans :  FFLAS::FflasNoTrans;
	FFLAS::FFLAS_DIAG diag = (atoi(argv[9])) ? FFLAS::FflasUnit :  FFLAS::FflasNonUnit;

	if (   ((side == FFLAS::FflasRight) &&(k != n))
	    || ((side == FFLAS::FflasLeft)&&(k != m))) {
		cerr<<"Error in the dimensions of the input matrices"<<endl;
		exit(-1);
	}

	Timer t; t.clear();
	double time=0.0;
	//write_field(F, cerr<<"A="<<endl, A, k,k,k);

	for(int i = 0;i<nbit;++i){
		t.clear();
		t.start();
		FFLAS::ftrsm (F, side, uplo, trans, diag, m, n, alpha, A, k, B, n);
		t.stop();
		time+=t.usertime();
		if (i+1<nbit)
			for (int i=0; i<m*n;++i)
				F.assign(*(B+i),*(B2+i));
	}

#if DEBUG
	Field::Element invalpha;
	F.inv(invalpha, alpha);

	FFLAS::ftrmm (F, side, uplo, trans, diag, m, n, invalpha, A, k, B, n);
	bool wrong = false;

	for (int i=0;i<m;++i)
		for (int j=0;j<n;++j)
			if ( !F.areEqual(*(B2+i*n+j), *(B+i*n+j))){
				cerr<<"B2 ["<<i<<", "<<j<<"] = "<<(*(B2+i*n+j))
				    <<" ; B ["<<i<<", "<<j<<"] = "<<(*(B+i*n+j))
				    <<endl;
				wrong = true;
			}

	if ( wrong ){
		cerr<<"FAIL"<<endl;
		//write_field (F,cerr<<"B2="<<endl,B2,m,n,n);
		//write_field (F,cerr<<"B="<<endl,B,m,n,n);
	} else
		cerr<<"PASS"<<endl;
#endif

	delete[] A;
	delete[] B;
	delete[] B2;

#if TIME
	double mflops = m*n/1000000.0*nbit*n/time;
	cerr<<"m,n = "<<m<<" "<<n<<". ftrsm "
	    <<((side == FFLAS::FflasLeft)?" Left ":" Right ")
	    <<((uplo == FFLAS::FflasLower)?" Lower ":" Upper ")
	    <<((diag == FFLAS::FflasUnit)?" Unit ":" NonUnit ")
	    <<((trans == FFLAS::FflasTrans)?" Trans ":" NoTrans ")
	    <<"over Z/"<<atoi(argv[1])<<"Z :"
	    <<endl
	    <<"t= "
	    << time/nbit
	    << " s, Mffops = "<<mflops
	    << endl;

	cout<<m<<" "<<n<<" "<<mflops<<" "<<time/nbit<<endl;
#endif
}
