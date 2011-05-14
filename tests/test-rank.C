/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) Fflas-Ffpack
 * Written by Clément Pernet
 * This file is Free Software and part of Fflas-Ffpack.
 * See COPYING for license information.
 */


//--------------------------------------------------------------------------
//                        Test for rank
//
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------

#include <iomanip>
#include <iostream>
#include "fflas-ffpack/field/modular-balanced.h"
#include "timer.h"
#include "Matio.h"
#include "fflas-ffpack/ffpack/ffpack.h"



using namespace std;
using namespace FFPACK;

typedef ModularBalanced<double> Field;

int main(int argc, char** argv){

	int n,m;
	int nbit=atoi(argv[3]); // number of times the product is performed
	cerr<<setprecision(10);
	if (argc !=  4)	{
		cerr<<"Usage : test-rank <p> <A> <<i>"
		    <<endl
		    <<"         to compute the rank of A mod p (i computations)"
		    <<endl;
		exit(-1);
	}
	Field F(atof(argv[1]));
	Field::Element * A;
	A = read_field(F,argv[2],&m ,&n);

	Timer tim,t;
	t.clear();
	tim.clear();
	size_t r=0;
	for(int i = 0;i<nbit;++i){
		t.clear();
		t.start();
		r = FFPACK::Rank (F, m, n, A, n);
		t.stop();
		tim+=t;
		if (i+1<nbit){
			delete[] A;
			A = read_field(F,argv[2],&m,&n);
		}
	}

	double mflops = 2.0/3.0*(n*(double)r/1000000.0)*nbit*n/tim.usertime();
	cerr<<"m,n = "<<m<<", "<<n<<" Rank (A) = " << r
		     << " mod "<<atoi(argv[1])<<" : t= "
		     << tim.usertime()/nbit
		     << " s, Mffops = "<<mflops
		     << endl;

	cout<<m<<" "<<n<<" "<<r<<" "<<mflops<<" "<<tim.usertime()/nbit<<endl;
}
