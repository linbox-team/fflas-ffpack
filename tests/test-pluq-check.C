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

int main(int argc, char** argv) {
	size_t iter = 3 ;
	Givaro::Integer q = 131071;
	size_t m = 0;
	size_t n = 0;
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

	size_t pass = 0;	// number of tests that have successfully passed

	for(size_t it=0; it<iter; ++it) {
		if (random_dim) {
			m = rand() % 10000 + 1;
			n = rand() % 10000 + 1;
		}
			
		std::cout << "m= " << m << "    n= " << n << "\n";

		Field::Element_ptr A = FFLAS::fflas_new(F,m,n);
		size_t *P = FFLAS::fflas_new<size_t>(m);
		size_t *Q = FFLAS::fflas_new<size_t>(n);

		// generate a random matrix A
		for( size_t i = 0; i < m*n; ++i )
			RValue.random( *(A+i) );
		//write_field(F,std::cerr<<"A:=",A,m,n,n,true) <<std::endl;
  
		try {
			//FFPACK::LUdivine_small(F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, m, n, A, n, P, Q);
			//FFPACK::LUdivine(F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, m, n, A, n, P, Q, FFPACK::FfpackSingular,60);
			std::cout << "Computing PLUQ...\n";
			FFPACK::PLUQ(F, FFLAS::FflasNonUnit, m, n, A, n, P, Q);
			std::cout << "Verification successful\n";
			pass++;
		} catch(FailurePLUQcheck &e) {
			std::cout << "Verification failed!\n";
		}

		FFLAS::fflas_delete(A,P,Q);
	}

	std::cout << pass << "/" << iter << " tests were successful.\n";

	return 0;
}
