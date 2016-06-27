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
//          Test for Checker_invert
//--------------------------------------------------------------------------
#define ENABLE_ALL_CHECKINGS 1

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"

int main(int argc, char** argv) {
	srand (time(NULL));
	typedef Givaro::Modular<double> Field;
	Givaro::Integer q = 131071;
	size_t iter = 3;
	size_t MAXM = 1000;
    size_t seed(0);
    
	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).", TYPE_INTEGER , &q },
		{ 'n', "-n N", "Set the maximal size of the matrix.", TYPE_INT , &MAXM },
		{ 'i', "-i R", "Set number of repetitions.", TYPE_INT , &iter },
        { 's', "-s N", "Set the seed.", TYPE_INT , &seed },
		END_OF_ARGUMENTS
	};
	FFLAS::parseArguments(argc,argv,as);
    FFLAS::writeCommandString(std::cout, as) << std::endl;

	Field F(q);

	Field::RandIter Rand(F,0,seed);
	Field::NonZeroRandIter NZRand(Rand);
    srandom(seed);

	int nullity;
	size_t m = MAXM, pass = 0;
	for (size_t i=0; i<iter; ++i) {
		m = random() % MAXM + 1;
		std::cout << "m= " << m << "\n";

		Field::Element_ptr A = FFLAS::fflas_new(F,m,m);

	for (size_t i=0;i<m;++i){
		for (size_t j=0;j<i;++j)
			Rand.random(A[i*m+j]);
        for(size_t j=i+1;j<m;++j)
            F.assign(A[i*m+j],F.zero);
        NZRand.random(A[i*m+i]);
    }

		Checker_invert<Field> checker(Rand,m,A,m);
		FFPACK::Invert(F,m,A,m,nullity);
		try {
			checker.check(A,nullity);
			std::cout << "Verification successful\n";
			pass++;
		} catch (FailureInvertCheck &e) {
			std::cout << "Verification failed!\n";
		}

		FFLAS::fflas_delete(A);
	}

	std::cout << pass << "/" << iter << " tests were successful.\n";


	return 0;
}
