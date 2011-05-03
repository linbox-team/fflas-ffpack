/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) Fflas-Ffpack
 * Written by Clément Pernet
 * This file is Free Software and part of Fflas-Ffpack.
 * See COPYING for license information.
 */


//--------------------------------------------------------------------------
//                        Test for ftrtri : 1 computation
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
#include "fflas-ffpack/ffpack/ffpack.h"



using namespace std;
using namespace FFPACK;

typedef Modular<float> Field;

int main(int argc, char** argv)
{

	int n;
	int nbit=atoi(argv[3]); // number of times the product is performed
	cerr<<setprecision(10);
	Field::Element zero, one;

	if (argc != 4)	{
		cerr<<"Usage : test-ftrtri <p> <A> <<i>"
		    <<endl
		    <<"         to invert a triangular matrix  A mod p (i computations)"
		    <<endl;
		exit(-1);
	}
	Field F(atoi(argv[1]));
	F.init(zero,0.0);
	F.init(one,1.0);
	Field::Element * A,*Ab;
	A = read_field(F,argv[2],&n,&n);
	Ab = new Field::Element[n*n];

	for (int i=0; i<n;++i){
		for(int j=0; j<i; ++j)
			F.assign(*(Ab+i*n+j),*(A+i*n+j));
		F.assign(*(Ab+i*(n+1)), 1.0);
		for(int j=i+1; j<n; ++j)
			F.assign(*(Ab+i*n+j),0.0);
	}


	Field::Element * X = new Field::Element[n*n];

	Timer tim,t; t.clear();tim.clear();

	for(int i = 0;i<nbit;++i){
		t.clear();
		t.start();
				FFPACK::trinv_left (F, n, A, n, X, n);
		//FFPACK::ftrtri(F, FFLAS::FflasLower, FFLAS::FflasUnit, n, A, n);
		t.stop();
		tim+=t;
		if (i+1<nbit)
			for (int i=0; i<n*n;++i)
				F.assign(*(A+i),*(Ab+i));
	}

#if DEBUG
	FFLAS::ftrmm (F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit,
		      n, n, 1.0,
		      X,
		      //A,
		      n, Ab, n);
	bool wrong = false;

	for (int i=0;i<n;++i)
		for (int j=0;j<n;++j)
			if ( ((i!=j) && !F.areEqual(*(Ab+i*n+j),zero))
			     ||((i==j) &&!F.areEqual(*(Ab+i*n+j),one)))
				wrong = true;

	if ( wrong ){
		cerr<<"FAIL"<<endl;
		write_field (F,cerr<<"Ab="<<endl,A,n,n,n);
		//write_field (F,cerr<<"X="<<endl,X,n,n,n);
	}else{

		cerr<<"PASS"<<endl;
	}
#endif

#if TIME
	double mflops = 1.0/3.0*(n*n/1000000.0)*nbit*n/tim.usertime();
	cerr<<"n = "<<n<<" Inversion over Z/"<<atoi(argv[1])<<"Z : t= "
	     << tim.usertime()/nbit
	     << " s, Mffops = "<<mflops
	     << endl;

	cout<<n<<" "<<mflops<<" "<<tim.usertime()/nbit<<endl;
#endif
}
