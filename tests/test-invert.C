/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
//--------------------------------------------------------------------------
//                        Test for invert : 1 computation
//                  
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------

#define DEBUG 1
#define TIME 1

#include <iomanip>
#include <iostream>
#include "fflas-ffpack/modular-balanced.h"
#include "timer.h"
#include "Matio.h"
#include "fflas-ffpack/ffpack.h"



using namespace std;

typedef Modular<float> Field;

int main(int argc, char** argv){

	int n;
	int nbit=atoi(argv[3]); // number of times the product is performed
	cerr<<setprecision(10);
	Field::Element zero, one;

	if (argc != 4)	{
		cerr<<"Usage : test-invert <p> <A> <<i>"
		    <<endl
		    <<"         to invert A mod p (i computations)"
		    <<endl;
		exit(-1);
	}
	Field F(atoi(argv[1]));
	F.init(zero,0.0);
	F.init(one,1.0);
	Field::Element * A,*Ab;
	A = read_field(F,argv[2],&n,&n);
	Ab = new Field::Element[n*n];
	for (int i=0; i<n*n;++i)
		F.assign(*(Ab+i),*(A+i));
	Field::Element * X = new Field::Element[n*n];

	Timer tim,t; t.clear();tim.clear(); 
	int nullity=0;
	
	for(int i = 0;i<nbit;++i){
		t.clear();
		t.start();
		FFPACK::Invert2 (F, n, A, n, X, n, nullity);
		t.stop();
		tim+=t;
		if (i+1<nbit)
			for (int i=0; i<n*n;++i)
				F.assign(*(A+i),*(Ab+i));
	}

#if DEBUG
	Field::Element *I = new Field::Element[n*n];
	FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n,
		      1.0, Ab, n, X, n, 0.0, I, n); 
	bool wrong = false;

	for (int i=0;i<n;++i)
		for (int j=0;j<n;++j)
			if ( ((i!=j) && !F.areEqual(*(I+i*n+j),zero))
			     ||((i==j) &&!F.areEqual(*(I+i*n+j),one)))
				wrong = true;
	
	if ( wrong ){
		if (nullity > 0)
			cerr<<"Matrix is singular over Z/"<<argv[1]<<"Z"<<endl;
		else{
			cerr<<"FAIL"<<endl;
			write_field (F,cerr<<"A="<<endl,A,n,n,n);
			write_field (F,cerr<<"X="<<endl,X,n,n,n);
			write_field (F,cerr<<"I="<<endl,I,n,n,n);
		}
	} else {
		cerr<<"PASS"<<endl;
	}
#endif

#if TIME
	double mflops = 2*(n*n/1000000.0)*nbit*n/tim.usertime();
	cerr<<"n = "<<n<<" Inversion over Z/"<<atoi(argv[1])<<"Z : t= "
	     << tim.usertime()/nbit 
	     << " s, Mffops = "<<mflops
	     << endl;
	
	cout<<n<<" "<<mflops<<" "<<tim.usertime()/nbit<<endl;
#endif
}
