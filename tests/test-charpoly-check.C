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
//          Test for Checker_charpoly
//--------------------------------------------------------------------------

#define ENABLE_ALL_CHECKINGS 1


#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"


template <class Field, class Polynomial>
void printPolynomial (const Field &F, Polynomial &v)
{
	for (int i = v.size() - 1; i >= 0; i--) {
		F.write (std::cout, v[i]);
		if (i > 0)
			std::cout << " x^" << i << " + ";
	}
	std::cout << std::endl;
}

int main(int argc, char** argv) {
	srand (time(NULL));
	typedef Givaro::Modular<double> Field;
	Givaro::Integer q = 131071;
	size_t iter = 3;
    size_t MAXN = 100;
	
	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).", TYPE_INTEGER , &q },
		{ 'i', "-i R", "Set number of repetitions.", TYPE_INT , &iter },
		{ 'n', "-n N", "Set the size of the matrix.", TYPE_INT , &MAXN },
		END_OF_ARGUMENTS
	};
	FFLAS::parseArguments(argc,argv,as);

	Field F(q);
	typedef std::vector<Field::Element> Polynomial;

	Field::RandIter Rand(F);
	Field::Element_ptr A = FFLAS::fflas_new(F,MAXN,MAXN);

	size_t pass = 0;
	for (size_t i=0; i<iter; ++i) {

		size_t n = rand() % MAXN + 1;
// 		std::cout << "n= " << n << "\n";

		Polynomial g(n);

		for( size_t i = 0; i < n*n; ++i )
			Rand.random( *(A+i) );

		try {
			//write_field(F,std::cerr<<"A=",A,n,n,n,true) <<std::endl;
			FFPACK::Checker_charpoly<Field,Polynomial> checker(F,n,A);
			FFPACK::CharPoly(F,g,n,A,n,FFPACK::FfpackLUK);
			//printPolynomial(F,g);
			checker.check(g);
			std::cout << "Verification successful\n";
			pass++;
		} catch(FailureCharpolyCheck &e) {
			std::cout << "Verification failed!\n";
		}
	}

	std::cout << pass << "/" << iter << " tests were successful.\n";	

	return 0;
}
