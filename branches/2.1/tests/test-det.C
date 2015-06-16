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

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iomanip>
#include <iostream>
#include <givaro/modular-balanced.h>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"

#include "test-utils.h"
#include "Matio.h"


// using namespace std;
template<class Field>
bool test_det(Field &F, size_t n, int iter)
{
	typedef typename Field::Element Element;
	//! @todo test with stride
	Element * A = FFLAS::fflas_new<Element>(n*n);
	// A = read_field(F,argv[2],&n,&n);

	bool pass = true;
#ifdef TIME_IT
 FFLAS::Timer tim,t; t.clear();tim.clear();
#endif
	Element d=0;
	Element dt=-4;
	for(int i = 0;i<iter;++i){
		F.init(dt,dt);
		// std::cout << dt << std::endl;
		FFPACK::RandomMatrixWithDet(F,A,dt,n,n);
#ifdef TIME_IT
		t.clear();
		t.start();
#endif
		d = FFPACK::Det (F, n, n, A, n);
		// std::cout << d << std::endl;
#ifdef TIME_IT
		t.stop();
		tim+=t;
#endif
		// if (i+1<iter){
		// FFLAS::fflas_delete( A);
		// A = read_field(F,argv[2],&n,&n);
		if (dt != d) {
			pass = false;
			break;
		}
		++dt;
	}

#ifdef TIME_IT
	double mflops = 2.0/3.0*(n*n/1000000.0)*iter*n/tim.usertime();
	F.write (std::cerr<<"n = "<<n<<" Det (A) = ",d)
	<< " mod "<<atoi(argv[1])<<" : t= "
	<< tim.usertime()/(double)iter
	<< " s, Mffops = "<<mflops
	<< std::endl;

	std::cout<<n<<" "<<mflops<<" "<<tim.usertime()/(double)iter<<std::endl;
#endif
	FFLAS::fflas_delete( A);
	return pass;
	}

int main(int argc, char** argv)
{

	static int iters =10 ;
	static uint64_t p = 65521 ;
	static size_t n = 200 ;

	static Argument as[] = {
		{ 'p', "-p P", "Set the field characteristic.",         TYPE_INT , &p },
		{ 'n', "-n N", "Set the dimension of the matrix.",      TYPE_INT , &n },
		{ 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(argc,argv,as);

	// int n;
	// int iter=atoi(argv[3]); // number of times the product is performed
	std::cerr<<std::setprecision(10);
#if 0 /*  don't know how to do this in parseArguments ; grosse flemme */
	if (argc != 4)	{
		std::cerr<<"Usage : test-det <p> <A> <<i>"
		    <<std::endl
		    <<"         to compute the determinant of A mod p (i computations)"
		    <<std::endl;
		exit(-1);
	}
#endif
	bool pass = true ;
	typedef Givaro::ModularBalanced<double> Field;
	Field F(p);
	pass &= test_det(F,n,iters);
	// pass &= test_det(F,0,iters);

	return ((pass==true)?0:1);
}
