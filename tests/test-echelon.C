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
//          Test for the echelon factorisation
//--------------------------------------------------------------------------
// usage: test-echelon p A n, for n lsp factorization
// of A over Z/pZ
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
#define DEBUG 1
// Debug option  0: no debug
//               1: check A = LQUP
//-------------------------------------------------------------------------
// using namespace std;



//#define __LUDIVINE_CUTOFF 1
#include <iostream>
#include <iomanip>
#include "Matio.h"
#include "utils/timer.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/ffpack/ffpack.h"

#include "test-utils.h"

#include "fflas-ffpack/utils/args-parser.h"


using namespace FFPACK;

//!@bug does not check that the form is actually correct, just that the product is ok.
template<class Field>
bool
test_echelon(Field &F, size_t m, size_t n, size_t r, size_t iters)
{
	typedef typename Field::Element Element ;
	Element * A = new Element[m*n];
	Element * B = new Element[m*n];

	Element * L = new Element[m*n];
	Element * U = new Element[n*n];
	Element * X = new Element[m*n];


	// A = read_field(F,argv[2],&m,&n);
	size_t lda = n; //!@todo check lda

	size_t *P = new size_t[n];
	size_t *Q = new size_t[m];
	size_t R = (size_t)-1;

	//	size_t cutoff = atoi(argv[3]);
	// iters = atoi(argv[3]);

#ifdef TIME_IT
	Timer tim,timc;
	timc.clear();
#endif

	bool pass=true;


	for (size_t  l=0;l<iters;l++){
		// if (i) {
		// delete[] A;
		// A = read_field(F,argv[2],&m,&n);
		// }
		R = (size_t)-1;
		RandomMatrixWithRank(F,A,r,m,n,lda);
		FFLAS::fcopy(F,m,n,B,lda,A,lda);
		for (size_t j=0;j<n;j++)
			P[j]=0;
		for (size_t j=0;j<m;j++)
			Q[j]=0;
#ifdef TIME_IT
		tim.clear();
		tim.start();
#endif
		R = FFPACK::ColumnEchelonForm (F, m, n, A, n, P, Q);
		if (R != r) {
			pass = false;
			break;
		}
#ifdef TIME_IT
		tim.stop();
		timc+=tim;
#endif
		//write_field (F,std::cerr<<"Result = "<<std::endl, A, m,n,n);

		// 	std::cerr<<"P = [";
		// 	for (size_t i=0; i<n; ++i)
		// 		std::cerr<<P[i]<<" ";
		// 	std::cerr<<"]"<<std::endl;
		// 	std::cerr<<"Q = [";
		// 	for (size_t i=0; i<m; ++i)
		// 		std::cerr<<Q[i]<<" ";
		// 	std::cerr<<"]"<<std::endl;
		// #if DEBUG

		Element one=F.one;
		// F.init(zero,0.0);
		// F.init(one,1.0);
		// for (size_t i=0; i<R; ++i){
		// 	for (size_t j=0; j<=i; ++j)
		// 		F.assign ( *(U + i*n + j), zero);
		// 	F.init (*(U+i*(n+1)),one);
		// 	for (size_t j=i+1; j<n; ++j)
		// 		F.assign (*(U + i*n + j), *(A+ i*n+j));
		// }
		// for (size_t i=R;i<n; ++i){
		// 	for (size_t j=0; j<n; ++j)
		// 		F.assign(*(U+i*n+j), zero);
		// 	F.assign(*(U+i*(n+1)),one);
		// }
		FFPACK::TriangularFromLU (F, FFLAS::FflasUpper, FFLAS::FflasNonUnit, n, n,
					  R, U, n, A, n);
		    // Adding I_{n-R} on the bottom right corner
		for (size_t i=R;i<n; ++i)
			F.assign (*(U+i*(n+1)),one);

		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, n, 0, (int)R, U, n, P);

		FFPACK::EchelonFromLU (F, FFLAS::FflasLower, FFLAS::FflasUnit, m,n,R,Q,L,n,A,n);

		// for ( size_t i=0; i<m; ++i ){
		// 	size_t j=0;
		// 	for (; j <= ((i<R)?i:R) ; ++j )
		// 		F.assign( *(L + i*n+j), *(A+i*n+j));
		// 	for (; j<m; ++j )
		// 		F.assign( *(L+i*n+j), zero);
		// }
		// 	std::cerr<<"P = ";
		// 	for (size_t i=0; i<n;++i)
		// 		std::cerr<<" "<<P[i];
		// 	std::cerr<<std::endl;
		// 	std::cerr<<"Q = ";
		// 	for (size_t i=0; i<m;++i)
		// 		std::cerr<<" "<<Q[i];
		// 	std::cerr<<std::endl;

		// 	write_field(F,std::cerr<<"A = "<<std::endl,A,m,n,n);
		//  	write_field(F,std::cerr<<"L = "<<std::endl,L,m,n,n);
		//  	write_field(F,std::cerr<<"U = "<<std::endl,U,m,n,n);

		// Element * B =  read_field(F,argv[2],&m,&n);

		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,n, 1.0,
			      B, n, U, n, 0.0, X,n);
		//delete[] A;

		bool fail = false;
		for (size_t i=0; i<m; ++i)
			for (size_t j=0; j<n; ++j)
				if (!F.areEqual (*(L+i*n+j), *(X+i*n+j)))
					fail=true;

		// write_field(F,std::cerr<<"X = "<<std::endl,X,m,n,n);
		// write_field(F,std::cerr<<"L = "<<std::endl,L,m,n,n);
		// write_field(F,std::cerr<<"A = "<<std::endl,A,m,n,n);

		if (fail) {
			std::cerr<<"FAIL"<<std::endl;
			pass = false;
			break;
		}


		// else
		// std::cerr<<"PASS"<<std::endl;

		// 	std::cout<<m<<" "<<n<<" M"<<std::endl;
		// 	for (size_t i=0; i<m; ++i)
		// 		for (size_t j=0; j<n; ++j)
		// 			if (!F.isZero(*(A+i*n+j)))
		// 				std::cout<<i+1<<" "<<j+1<<" "<<(*(A+i*n+j))<<std::endl;
		// 	std::cout<<"0 0 0"<<std::endl;

		// delete[] U;
		// delete[] L;
		// delete[] X;
		// #endif
	}

	delete[] U;
	delete[] L;
	delete[] X;
	delete[] B;
	delete[] A;
	delete[] P;
	delete[] Q;

#ifdef TIME_IT
	double t = timc.usertime();
	double numops = m*m/1000.0*(n-m/3.0);

	std::cerr<<m<<"x"<< n
	    << " : rank = " << R << "  ["
	    << ((double)iters/1000.0*(double)numops / t)
	    << " MFops "
	    << " in "
	    << t/iters<<"s"
	    <<"]"<< std::endl;
// 	std::cout<<m
// 	    <<" "<<((double)iters/1000.0*(double)numops / t)
// 	    <<" "<<t/iters
// 	    <<std::endl;
#endif

	return pass;

}

int main(int argc, char** argv){
	std::cerr<<setprecision(20);

	int    p = 101;
	size_t m = 50;
	size_t n = 50;
	size_t r = 20;
	size_t iters = 2 ;

	static Argument as[] = {
		{ 'p', "-p P", "Set the field characteristic.",         TYPE_INT , &p },
		{ 'n', "-n N", "Set the number of cols in the matrix.", TYPE_INT , &n },
		{ 'm', "-m N", "Set the number of rows in the matrix.", TYPE_INT , &m },
		{ 'r', "-r r", "Set the rank of the matrix."          , TYPE_INT , &r },
		{ 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
		// { 'f', "-f file", "Set input file", TYPE_STR, &file },
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(argc,argv,as);


	// if (argc!=4){
		// std::cerr<<"usage : test-lqup <p> <A> <i>"<<std::endl
		    // <<"        to do i Echelon factorisation of A with "
		    // <<std::endl;
		// exit(-1);
	// }
	bool pass = true ;
	typedef Modular<double> Field;
	Field F(p);
	pass &= test_echelon(F,m,n,r,iters);
	return !pass ;
}
