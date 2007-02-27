/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
//--------------------------------------------------------------------------
//                        Test for fgemv : 1 computation
//                  
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------

#define DEBUG 0
#define TIME 1

#include <iomanip>
#include <iostream>
#include "fflas-ffpack/modular-balanced.h"
#include "timer.h"
#include "Matio.h"
#include "fflas-ffpack/fflas.h"



using namespace std;

typedef Modular<double> Field;

int main(int argc, char** argv){

	int m,n,k;
	int nbit=atoi(argv[4]); // number of times the product is performed
	cerr<<setprecision(10);
	Field::Element alpha,beta;


	if (argc != 8)	{
		cerr<<"Usage : test-fgemm <p> <A> <b> <i>"
		    <<" <alpha> <beta> <c>"<<endl
		    <<"         to do i computations of c <- alpha AB + beta C"
		    <<endl;
		exit(-1);
	}
	Field F(atoi(argv[1]));

	F.init( alpha, double(atoi(argv[5])));
	F.init( beta, double(atoi(argv[6])));

	Field::Element * A;
	Field::Element * b;
	size_t lda;
	size_t ldb;
	
	b = read_field(F,argv[3],&n,&k);
	A = read_field(F,argv[2],&m,&n);
	
	Field::Element * c;
	c = new Field::Element[n];

	

	Timer tim,t; t.clear();tim.clear(); 
	for(int i = 0;i<nbit;++i){
		c = read_field(F,argv[7],&m,&k);
		t.clear();
		t.start();
		FFLAS::fgemv (F, FFLAS::FflasNoTrans,m,n,alpha, A,n, b,1,
			      beta, c, 1);
		t.stop();
		tim+=t;
	}

#if TIME
	double mflops = (2.0*(m*n/1000000.0)*nbit/tim.usertime());
	cerr << m <<"x" <<n <<" : fgemv over Z/"
	     <<atoi(argv[1])<<"Z : [ "
	     <<mflops<<" MFops in "<<tim.usertime()/nbit<<"s]"
	     << endl;
	
	cerr<<"alpha, beta = "<<alpha <<", "<<beta <<endl;

	cout<<m<<" "<<n<<" "<<mflops<<" "<<tim.usertime()/nbit<<endl;
#endif
}  


