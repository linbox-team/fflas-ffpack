/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

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

#define ENABLE_ALL_CHECKINGS 1 // DO NOT CHANGE
#define _NR_TESTS 10
#define _MAX_SIZE_MATRICES 2000

#include "fflas-ffpack/config-blas.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/checkers/checkers.h"
#include <fstream>

using namespace std;

int main(int argc, char** argv) {
	size_t NR_TESTS = _NR_TESTS;
	int    q    = 131071;
	size_t    MAX_SIZE_MATRICES    = _MAX_SIZE_MATRICES;

	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
		{ 'n', "-n N", "Set the dimension of the matrix.",               TYPE_INT , &MAX_SIZE_MATRICES },
		{ 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &NR_TESTS },
		END_OF_ARGUMENTS
	};

        FFLAS::parseArguments(argc,argv,as);

	srand (time(NULL));
	typedef Givaro::Modular<double> Field;
	typedef std::vector<Field::Element> Polynomial;

	Field F(q);
	Field::RandIter Rand(F);
	Field::NonZeroRandIter NZRand(Rand);

	size_t pass;
	FFLAS::Timer chrono;
	double time1, time2;

	Field::Element_ptr A = FFLAS::fflas_new(F,MAX_SIZE_MATRICES,MAX_SIZE_MATRICES);
	Field::Element_ptr B = FFLAS::fflas_new(F,MAX_SIZE_MATRICES,MAX_SIZE_MATRICES);
	Field::Element_ptr C = FFLAS::fflas_new(F,MAX_SIZE_MATRICES,MAX_SIZE_MATRICES);
	typename Field::Element alpha,beta,tmp;
	F.init(alpha, rand()%1000+1);
	F.init(beta,  rand()%1000+1);
	size_t m,n,k,lda,ldb,ldc;
	FFLAS::FFLAS_TRANSPOSE ta,tb;

        std::cout << "     Matrix size\tSuccess rate\t\tTime comput.\t\tTime checker\n\n";

	// #####   FGEMM   #####
	std::cout << "FGEMM:\n";
	for (size_t i=0; i<MAX_SIZE_MATRICES; i+=500) {
		pass = 0; time1 = 0.0; time2 = 0.0;
		for (size_t j=0; j<NR_TESTS; ++j) {
			m = rand() % 500 + i;
			n = rand() % 500 + i;
			k = rand() % 500 + i;
			ta = FFLAS::FflasNoTrans;//rand()%2 ? FFLAS::FflasNoTrans : FFLAS::FflasTrans,
			tb = FFLAS::FflasNoTrans;//rand()%2 ? FFLAS::FflasNoTrans : FFLAS::FflasTrans;
			lda = ta == FFLAS::FflasNoTrans ? k : m,
			ldb = tb == FFLAS::FflasNoTrans ? n : k,
			ldc = n;

			PAR_BLOCK { FFLAS::pfrand(F,Rand, m,k,A,m/MAX_THREADS); }
			PAR_BLOCK { FFLAS::pfrand(F,Rand, k,n,B,k/MAX_THREADS); }
			PAR_BLOCK { FFLAS::pfrand(F,Rand, m,n,C,n/MAX_THREADS); }

			chrono.clear(); chrono.start();
			FFLAS::Checker_fgemm<Field> checker1(Rand,m,n,k,beta,C,ldc);
			chrono.stop(); time1 += chrono.usertime();

			chrono.clear(); chrono.start();
			FFLAS::fgemm(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
			chrono.stop(); time2 += chrono.usertime();

			chrono.clear(); chrono.start();
			pass += checker1.check(ta,tb,alpha,A,lda,B,ldb,C) ? 1 : 0;
			chrono.stop(); time1 += chrono.usertime();
		}
		time1 /= NR_TESTS;
		time2 /= NR_TESTS;
		std::cout << "     " << i << "-" << i+500 << "\t\t" << pass << "/" << NR_TESTS << "\t\t\t" << time2 
				<< "\t\t" << time1 << endl;
	}
	cout << endl;



	// #####   FTRSM   #####
	std::cout << "FTRSM:\n";
	for (size_t i=0; i<MAX_SIZE_MATRICES; i+=500) {
		pass = 0; time1 = 0.0; time2 = 0.0;
		for (size_t j=0; j<NR_TESTS; ++j) {
			m = rand() % 500 + i;
			n = rand() % 500 + i;
			FFLAS::FFLAS_SIDE side = rand()%2?FFLAS::FflasLeft:FFLAS::FflasRight;
			FFLAS::FFLAS_UPLO uplo = rand()%2?FFLAS::FflasLower:FFLAS::FflasUpper;
			FFLAS::FFLAS_TRANSPOSE trans = rand()%2?FFLAS::FflasNoTrans:FFLAS::FflasTrans;
			FFLAS::FFLAS_DIAG diag = rand()%2?FFLAS::FflasNonUnit:FFLAS::FflasUnit;
			k = (side==FFLAS::FflasLeft?m:n);

			for( size_t i = 0; i < m*n; ++i ) Rand.random( *(B+i) );
			for (size_t i=0;i<k;++i) {
				for (size_t j=0;j<i;++j)
					A[i*k+j]= (uplo == FFLAS::FflasLower)? Rand.random(tmp) : F.zero;
				A[i*k+i]= (diag == FFLAS::FflasNonUnit)? NZRand.random(tmp) : F.one;
				for (size_t j=i+1;j<k;++j)
					A[i*k+j]= (uplo == FFLAS::FflasUpper)? Rand.random(tmp) : F.zero;
			}

			chrono.clear(); chrono.start();
                        FFLAS::Checker_ftrsm<Field> checker2(Rand, m, n, alpha, B, n);
			chrono.stop(); time1 += chrono.usertime();

			chrono.clear(); chrono.start();
			FFLAS::ftrsm(F, side, uplo, trans, diag, m, n, alpha, A, k, B, n);
			chrono.stop(); time2 += chrono.usertime();

			chrono.clear(); chrono.start();
			pass += checker2.check(side, uplo, trans, diag, m, n, A, k, B, n);
			chrono.stop(); time1 += chrono.usertime();
		}
		time1 /= NR_TESTS;
		time2 /= NR_TESTS;
		std::cout << "     " << i << "-" << i+500 << "\t\t" << pass << "/" << NR_TESTS << "\t\t\t" << time2 
				<< "\t\t" << time1 << endl;
	}
	cout << endl;



	// #####   INVERT   #####
	std::cout << "INVERT:\n";
	int nullity;
	for (size_t i=0; i<MAX_SIZE_MATRICES; i+=500) {
		pass = 0; time1 = 0.0; time2 = 0.0;
		for (size_t j=0; j<NR_TESTS; ++j) {
			m = rand() % 500 + i;

			PAR_BLOCK { FFLAS::pfrand(F,Rand, m,m,A,m/MAX_THREADS); }

			try {
				chrono.clear(); chrono.start();
				FFPACK::Checker_invert<Field> checker3(Rand,m,A,m);
				chrono.stop(); time1 += chrono.usertime();
				
				chrono.clear(); chrono.start();
				FFPACK::Invert(F,m,A,m,nullity);
				chrono.stop(); time2 += chrono.usertime();
				
				chrono.clear(); chrono.start();
				pass += checker3.check(A,nullity);
				chrono.stop(); time1 += chrono.usertime();
			} catch(FailureInvertCheck &e) {
				std::cout << " invert verification failed! " << nullity << std::endl;
			} catch(FailurePLUQCheck &e) {
				std::cout << " internal PLUQ verification failed! " << std::endl;
			}
		}
		time1 /= NR_TESTS;
		time2 /= NR_TESTS;
		std::cout << "     " << i << "-" << i+500 << "\t\t" << pass << "/" << NR_TESTS << "\t\t\t" << time2 
				<< "\t\t" << time1 << endl;
	}
	cout << endl;
	



	// #####   PLUQ   #####
	std::cout << "PLUQ:\n";
	for (size_t i=0; i<MAX_SIZE_MATRICES; i+=500) {
		pass = 0; time1 = 0.0; time2 = 0.0;
		for (size_t j=0; j<NR_TESTS; ++j) {
			m = rand() % 500 + i;
			n = rand() % 500 + i;

			PAR_BLOCK { FFLAS::pfrand(F,Rand, m,n,A,m/MAX_THREADS); }

			size_t *P = FFLAS::fflas_new<size_t>(m);
			size_t *Q = FFLAS::fflas_new<size_t>(n);

			chrono.clear(); chrono.start();
                        FFPACK::Checker_PLUQ<Field> checker4 (Rand,m,n,A,n);
			chrono.stop(); time1 += chrono.usertime();

			chrono.clear(); chrono.start();
			k = FFPACK::PLUQ(F, FFLAS::FflasNonUnit, m, n, A, n, P, Q);
			chrono.stop(); time2 += chrono.usertime();

			chrono.clear(); chrono.start();
			pass += checker4.check(A,n,k,P,Q);
			chrono.stop(); time1 += chrono.usertime();

			FFLAS::fflas_delete(P,Q);
		}
		time1 /= NR_TESTS;
		time2 /= NR_TESTS;
		std::cout << "     " << i << "-" << i+500 << "\t\t" << pass << "/" << NR_TESTS << "\t\t\t" << time2 
				<< "\t\t" << time1 << endl;
	}
	cout << endl;




	// #####   CharPoly   #####
	std::cout << "CharPoly:\n";
	for (size_t i=0; i<MAX_SIZE_MATRICES; i+=500) {
		pass = 0; time1 = 0.0; time2 = 0.0;
		for (size_t j=0; j<NR_TESTS; ++j) {
			n = rand() % 500 + i;

			PAR_BLOCK { FFLAS::pfrand(F,Rand, n,n,A,n/MAX_THREADS); }

			try {
			Polynomial g(n);

			chrono.clear(); chrono.start();
                        FFPACK::Checker_charpoly<Field,Polynomial> checker5(Rand,n,A,n);
			chrono.stop(); time1 += chrono.usertime();

			chrono.clear(); chrono.start();
			FFPACK::CharPoly(F,g,n,A,n,FFPACK::FfpackLUK);
			chrono.stop(); time2 += chrono.usertime();

			chrono.clear(); chrono.start();
			pass += checker5.check(g);
			chrono.stop(); time1 += chrono.usertime();
			} catch(FailureCharpolyCheck &e) {
				std::cout << " charpoly verification failed! " << std::endl;
			} catch(FailurePLUQCheck &e) {
				std::cout << " internal PLUQ verification failed! " << std::endl;
			}
		}
		time1 /= NR_TESTS;
		time2 /= NR_TESTS;
		std::cout << "     " << i << "-" << i+500 << "\t\t" << pass << "/" << NR_TESTS << "\t\t\t" << time2 
				<< "\t\t" << time1 << endl;
	}

	return 0;
}
