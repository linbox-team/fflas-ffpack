/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2015 the FFLAS-FFPACK group
 * Written by Ashley Lesdalons <Ashley.Lesdalons@e.ujf-grenoble.fr>
 *
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
 *
 */

//--------------------------------------------------------------------------
//          Test for PLUQ check
//--------------------------------------------------------------------------

//-------------------------------------------------------------------------
//#define DEBUG   OR   #define PLUQ_check
// PLUQ_check option  0: no verification
//               	  1: enable PLUQ check in PLUQ and LUdivine
//-------------------------------------------------------------------------

#define DEBUG 1

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"

int exec();

int main(int argc, char** argv) {
	//std::cout << exec() << "\n";
	//exit(0);

	size_t iter = 3 ;
	Givaro::Integer q = 11;//131071;
	size_t m = 0;
	size_t n = 0;
	size_t r;
	bool random_dim = false;

	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).", TYPE_INTEGER , &q },
		{ 'm', "-m M", "Set the row dimension of A.", TYPE_INT , &m },
		{ 'n', "-n N", "Set the col dimension of A.", TYPE_INT , &n },
		{ 'i', "-i R", "Set number of repetitions.", TYPE_INT , &iter },
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(argc,argv,as);
	if (m == 0 || n == 0) random_dim = true;

	srand (time(NULL));

	typedef Givaro::Modular<double> Field;
	Field F(q);

	Field::RandIter RValue(F);

	Field::Element_ptr A;
	size_t *P, *Q;
	size_t pass = 0;	// number of tests that have successfully passed

	for(size_t it=0; it<iter; ++it) {
		if (random_dim) {
			m = rand() % 10000 + 1;
			n = rand() % 10000 + 1;
			std::cout << "m= " << m << "    n= " << n << "\n";
		}

		A = FFLAS::fflas_new(F,m,n);

		P = FFLAS::fflas_new<size_t>(m);
		Q = FFLAS::fflas_new<size_t>(n);

		// generate a random matrix A
		for( size_t i = 0; i < m*n; ++i )
			RValue.random( *(A+i) );

		Checker_PLUQ<Field> checker (F,A,m,n);

		//r = FFPACK::LUdivine_small(F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, m, n, A, n, P, Q);
		//r = FFPACK::LUdivine(F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, m, n, A, n, P, Q, FFPACK::FfpackSingular,60);
		r = FFPACK::PLUQ(F, FFLAS::FflasNonUnit, m, n, A, n, P, Q);

		pass += checker.check(A,r,P,Q) ? 1:0;
	}

	std::cout << pass << "/" << iter << " tests were successful.\n";

	FFLAS::fflas_delete(A,P,Q);

	return 0;
}




int exec() {
	Givaro::Integer q = 131071;
	size_t r,m=0,n=0;

	
	//m = 300;
	//n = 200;
	if (m == 0) m = rand() % 10000 + 1;
	if (n == 0) n = rand() % 10000 + 1;
	//std::cout << "m= " << m << "    n= " << n << "\n";

	typedef Givaro::Modular<double> Field;
	Field F(q);

	Field::RandIter RValue(F);

	Field::Element_ptr A;
	A = FFLAS::fflas_new(F,m,n);
	for( size_t i = 0; i < m*n; ++i )
			RValue.random( *(A+i) );
	//A[0] = 0; A[1] = 3; A[2] = 2; A[3] = 6; A[4] = 6; A[5] = 1;

	//write_field(F,std::cout<<"A = "<<std::endl,A,m,n,n);

	Field::Element_ptr v,av;
	v = FFLAS::fflas_new(F,n,1);
	av = FFLAS::fflas_new(F,m,1);
	typename Field::RandIter G(F);
	FFLAS::frand(F,G,n,v,1);
	//v[0] = 0; v[1] = 8;
	//write_field(F,std::cout<<"v = "<<std::endl,v,n,1,1);
	
	FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, n, F.one, A, n, v, 1, F.zero, av, 1);
	//write_field(F,std::cout<<"av = "<<std::endl,av,m,1,1);

	size_t *P = FFLAS::fflas_new<size_t>(m);
	size_t *Q = FFLAS::fflas_new<size_t>(n);

	r = FFPACK::PLUQ(F, FFLAS::FflasNonUnit, m, n, A, n, P, Q);
	//std::cout << "r = " << r << std::endl;
	//write_field(F,std::cout<<"LU = "<<std::endl,A,m,n,n);

	// v <-- Q.v
	FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans, 1, 0, n, v, 1, Q);
	//std::cout << "B1\n";

 	typename Field::Element_ptr R = FFLAS::fflas_new(F,r,1);

 	// R <- V1
 	FFLAS::fassign(F, r, 1, v, 1, R, 1);
 	//std::cout << "B2\n";

 	// R <- U1.R
 	FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, r, 1, F.one, A, n, R, 1);
 	//std::cout << "B3\n";

 	// R <- U2.V2 + R
 	if (r < n)
 		FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, r, 1, n-r, F.one, A+r, n, v+r, 1, F.one, R, 1);
 	//std::cout << "B4\n";

	typename Field::Element_ptr w = FFLAS::fflas_new(F,m,1);
	for (size_t i=0; i<m; ++i)
		F.assign(*(w+i),F.zero);
	//std::cout << "B5\n";

	// w2 <- L2.R
 	if (r < m)
 		FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m-r, 1, r, F.one, A+r*n, n, R, 1, F.zero, w+r, 1);
 	//std::cout << "B6\n";

 	// R <- L1.R
 	FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, r, 1, F.one, A, n, R, 1);
 	//std::cout << "B7\n";

 	// w1 <- R
 	FFLAS::fassign(F, r, 1, R, 1, w, 1);
 	
 	FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans, 1, 0, m, w, 1, P);

 	//write_field(F,std::cout<<"R = "<<std::endl,R,r,1,1);
 	//write_field(F,std::cout<<"w = "<<std::endl,w,m,1,1);
 	
 	// is av == w ?
	FFLAS::fsub(F, m, 1, w, 1, av, 1, av, 1);
	bool pass = FFLAS::fiszero(F,m,1,av,1);

	FFLAS::fflas_delete(A,av,v,R,P,Q);
	
	//std::cout << pass << std::endl;
	return pass;
}