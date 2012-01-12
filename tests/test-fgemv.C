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
//                        Test for fgemv : 1 computation
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

#if DEBUG
	Field::Element *d = new Field::Element[n];
	for (int i=0; i<m; ++i)
		F.mul (d[i], beta, b[i]);
	for (int i=0; i<m; ++i)
		F.mulin (b[i], alpha);
	for (int i=0; i<m; ++i)
		for (int j=0; j<n; ++j)
			F.axpyin (d[i], *(A+i*m+j), b[j]);
	bool fail = false;
	for (int i=0; i<m; ++i)
		if (!F.areEqual(d[i], c[i]))
			fail = true;

	if (fail)
		cerr<<"FAIL"<<endl;
	else
		cerr<<"PASS"<<endl;
	delete[] d;
#endif
	delete[] A;
	delete[] b;
	delete[] c;
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


