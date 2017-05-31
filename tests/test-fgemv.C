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

#define TIME 1

#if not defined(STD_RECINT_SIZE)
#define STD_RECINT_SIZE 8
#endif

#include <givaro/modular.h>
#include "recint/recint.h"
#include <iomanip>
#include <iostream>
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/fflas/fflas.h"



using namespace std;
using namespace FFPACK;

// typedef Givaro::Modular<double> Field;

typedef RecInt::ruint<STD_RECINT_SIZE> Ints;	
typedef Givaro::Modular<Ints> Field;	

int main(int argc, char** argv){

	size_t m,n,k;
	int nbit=atoi(argv[4]); // number of times the product is performed
	cerr<<setprecision(10);
	Field::Element alpha,beta;


	if (argc != 8)	{
		cerr<<"Usage : test-fgemv <p> <A> <b> <i> <alpha> <beta> <c>"
		    <<" <alpha> <beta> <c>"<<endl
		    <<"         to do i computations of c <- alpha Ab + beta c"
		    <<endl;
		exit(-1);
	}
	Field F(atoi(argv[1]));

	F.init( alpha, double(atoi(argv[5])));
	F.init( beta, double(atoi(argv[6])));

	Field::Element * A;
	Field::Element * b;

	FFLAS::ReadMatrix (argv[2],F,m,n,A);
	FFLAS::ReadMatrix (argv[3],F,n,k,b);

	Field::Element * c;


 FFLAS::Timer tim,t; t.clear();tim.clear();
	for(int i = 0;i<nbit;++i){
		FFLAS::ReadMatrix (argv[7],F,m,k,c);
		t.clear();
		t.start();
		FFLAS::fgemv (F, FFLAS::FflasNoTrans,m,n,alpha, A,n, b,1,
			      beta, c, 1);
		t.stop();
		tim+=t;
	}

#ifdef __FFLASFFPACK_DEBUG
	Field::Element *d;
	FFLAS::ReadMatrix (argv[7],F,m,k,d);
	for (int i=0; i<m; ++i)
		F.mulin (d[i], beta); 					// d <- beta c
	for (int i=0; i<m; ++i)
		F.mulin (b[i], alpha);					// b <- alpha b
	for (int i=0; i<m; ++i)
		for (int j=0; j<n; ++j)
			F.axpyin (d[i], *(A+i*m+j), b[j]);	// d <- Ad+b = alpha Ab + beta c
	bool fail = false;
	for (int i=0; i<m; ++i)
		if (!F.areEqual(d[i], c[i]))
			fail = true;

	if (fail) {
		cerr<<"FAIL"<<endl;
		WriteMatrix(std::cerr<<"i:=",F,n,k,b,k,true)<<std::endl;
		WriteMatrix(std::cerr<<"r:=",F,m,k,c,k,true)<<std::endl;
		WriteMatrix(std::cerr<<"d:=",F,m,k,d,k,true)<<std::endl;
		
		F.write(std::cerr<<"alpha:=", alpha) << ';' << std::endl;
		WriteMatrix(std::cerr<<"A:=",F,m,n,A,n,true)<<std::endl;
		FFLAS::ReadMatrix (argv[3],F,n,k,b);
		WriteMatrix(std::cerr<<"b:=",F,n,k,b,k,true)<<std::endl;
		F.write(std::cerr<<"beta:=", beta) << ';' << std::endl;
		FFLAS::ReadMatrix (argv[7],F,m,k,c);
		WriteMatrix(std::cerr<<"c:=",F,m,k,c,k,true)<<std::endl;
		std::cerr<<"p:=" << F.characteristic() << ';' << std::endl;
	
    } else
		cerr<<"PASS"<<endl;

	FFLAS::fflas_delete( d);
#endif
	FFLAS::fflas_delete( A);
	FFLAS::fflas_delete( b);
	FFLAS::fflas_delete( c);
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


