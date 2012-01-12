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
//                        Test for det
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

typedef FFPACK:: ModularBalanced<double> Field;

int main(int argc, char** argv){

	int n;
	int nbit=atoi(argv[3]); // number of times the product is performed
	cerr<<setprecision(10);
	if (argc != 4)	{
		cerr<<"Usage : test-det <p> <A> <<i>"
		    <<endl
		    <<"         to compute the determinant of A mod p (i computations)"
		    <<endl;
		exit(-1);
	}
	Field F(atof(argv[1]));
	Field::Element * A;
	A = read_field(F,argv[2],&n,&n);

	Timer tim,t; t.clear();tim.clear();
	Field::Element d=0;
	for(int i = 0;i<nbit;++i){
		t.clear();
		t.start();
		d = FFPACK::Det (F, n, n, A, n);
		t.stop();
		tim+=t;
		if (i+1<nbit){
			delete[] A;
			A = read_field(F,argv[2],&n,&n);
		}
	}

	double mflops = 2.0/3.0*(n*n/1000000.0)*nbit*n/tim.usertime();
	F.write (cerr<<"n = "<<n<<" Det (A) = ",d)
		     << " mod "<<atoi(argv[1])<<" : t= "
		     << tim.usertime()/nbit
		     << " s, Mffops = "<<mflops
		     << endl;

	cout<<n<<" "<<mflops<<" "<<tim.usertime()/nbit<<endl;
}
