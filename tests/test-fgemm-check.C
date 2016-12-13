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
//          Test for Checker_fgemm
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
	size_t iter = 10;
	
	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).", TYPE_INTEGER , &q },
		{ 'i', "-i R", "Set number of repetitions.", TYPE_INT , &iter },
		END_OF_ARGUMENTS
	};
	FFLAS::parseArguments(argc,argv,as);

	Field F(q);
	Field::RandIter Rand(F);
	FFLAS::FFLAS_TRANSPOSE ta,tb;

	size_t pass = 0;
	for (size_t i=0; i<iter; ++i) {

		size_t m = rand() % 1000 + 1;
		size_t n = rand() % 1000 + 1;
		size_t k = rand() % 1000 + 1;
		std::cout << "m= " << m << "    n= " << n << "    k= " << k << "   ";

		typename Field::Element alpha,beta;
		F.init(alpha); Rand.random(alpha);
		F.init(beta);  Rand.random(beta);
		
		ta = rand()%2 ? FFLAS::FflasNoTrans : FFLAS::FflasTrans;
		tb = rand()%2 ? FFLAS::FflasNoTrans : FFLAS::FflasTrans;

		std::cout << "m= " << m << "    n= " << n << "    k= " << k
			  <<" ta = ";
		if (ta==FFLAS::FflasNoTrans) std::cout<<"NoTrans";
		else std::cout<<"Trans";
		std::cout<< " tb =  ";
		if (tb==FFLAS::FflasNoTrans) std::cout<<"NoTrans";
		else std::cout<<"Trans";
		std::cout<<" : ";

		size_t lda = ta == FFLAS::FflasNoTrans ? k : m,
			ldb = tb == FFLAS::FflasNoTrans ? n : k,
			ldc = n;

		Field::Element_ptr A = FFLAS::fflas_new(F,m,k);
		Field::Element_ptr B = FFLAS::fflas_new(F,k,n);
		Field::Element_ptr C = FFLAS::fflas_new(F,m,n);

		FFLAS::frand(F,Rand, m,k,A,k);
		FFLAS::frand(F,Rand, k,n,B,n);
		FFLAS::frand(F,Rand, m,n,C,n);

		FFLAS::Checker_fgemm<Field> checker(F,m,n,k,beta,C,ldc);
		FFLAS::fgemm(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
		try {
			checker.check(ta,tb,alpha,A,lda,B,ldb,C);
			std::cout << "PASSED\n";
			pass++;
		} catch (FailureFgemmCheck &e) {
			std::cout << "FAILED\n";
			FFLAS::fflas_delete(A,B,C);
			return -1;
		}

		FFLAS::fflas_delete(A,B,C);
	}

	std::cout << pass << "/" << iter << " tests were successful.\n";

	return 0;
}
