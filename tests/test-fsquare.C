/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Clément Pernet
 * This file is Free Software and part of FFLAS-FFPACK.
 * 
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
 *.
 */


//--------------------------------------------------------------------------
//                        Test for fsquare : 1 computation
//
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------

//#define DEBUG 0
#define TIME 1

#include <iomanip>
#include <iostream>
#include "fflas-ffpack/field/modular-balanced.h"
#include "timer.h"
#include "Matio.h"
#include "fflas-ffpack/fflas/fflas.h"


using namespace FFPACK;
using namespace std;

typedef Modular<double> Field;

int main(int argc, char** argv){

	int n;

	cerr<<setprecision(10);
	if (argc != 6)	{
		cerr<<"Usage : test-fsquare <p> <A> <i>"
		    <<"<alpha> <beta>"
		    <<"         to do i computations of C <- AA"
		    <<endl;
		exit(-1);
	}
	Field F(atoi(argv[1]));

	Field::Element * A;
	Field::Element * C;
	// size_t lda;
	// size_t ldb;

	A = read_field(F,argv[2],&n,&n);
	int nbit=atoi(argv[3]); // number of times the product is performed

	Field::Element alpha,beta;
	F.init (alpha, (double)atoi(argv[4]));
	F.init (beta, (double)atoi(argv[5]));

	C = new Field::Element[n*n];
	Timer tim,t; t.clear();tim.clear();
	for(int i = 0;i<nbit;++i){
		t.clear();
		t.start();
		FFLAS::fsquare (F, FFLAS::FflasNoTrans,n, alpha, A,n, beta, C, n);
		t.stop();
		tim+=t;
	}

#if TIME
	double mflops = (2.0*(n*n/1000000.0)*nbit*n/tim.usertime());
	cerr << n <<"x" <<n <<" : fsquare over Z/"
	     <<atoi(argv[1])<<"Z : [ "
	     <<mflops<<" MFops in "<<tim.usertime()/nbit<<"s]"
	     << endl;

	cerr<<"alpha, beta = "<<alpha <<", "<<beta <<endl;

	cout<<n<<" "<<n<<" "<<mflops<<" "<<tim.usertime()/nbit<<endl;
#endif
}


