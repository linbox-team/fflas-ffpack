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

//#define __LUDIVINE_CUTOFF 1
#include <iostream>
#include <iomanip>
#include "Matio.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/ffpack/ffpack.h"

#include "test-utils.h"

#include "fflas-ffpack/utils/args-parser.h"


using namespace FFPACK;

//!@bug does not check that the form is actually correct, just that the product is ok.
template<class Field>
bool
test_colechelon(Field &F, size_t m, size_t n, size_t r, size_t iters)
{
	typedef typename Field::Element Element ;
	Element * A = FFLAS::fflas_new<Element>(m*n);
	Element * B = FFLAS::fflas_new<Element>(m*n);
	Element * L = FFLAS::fflas_new<Element>(m*n);
	Element * U = FFLAS::fflas_new<Element>(n*n);
	Element * X = FFLAS::fflas_new<Element>(m*n);     
	size_t lda = n; //!@todo check lda

	size_t *P = FFLAS::fflas_new<size_t>(n);
	size_t *Q = FFLAS::fflas_new<size_t>(m);
	size_t R = (size_t)-1;

	bool pass=true;

	for (size_t  l=0;l<iters;l++){
		R = (size_t)-1;
		RandomMatrixWithRank(F,A,r,m,n,lda);
		FFLAS::fassign(F,m,n,A,lda,B,lda);
		for (size_t j=0;j<n;j++) P[j]=0;
		for (size_t j=0;j<m;j++) Q[j]=0;

		R = FFPACK::ColumnEchelonForm (F, m, n, A, n, P, Q);

		if (R != r) {pass = false; break;}

		FFPACK::getEchelonTransform (F, FFLAS::FflasLower, FFLAS::FflasUnit, m,n,R,P,A,lda,U,n);

		FFPACK::getEchelonForm (F, FFLAS::FflasLower, FFLAS::FflasUnit, m,n,R,Q,A,n,L,n);

		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,n, 1.0, B, n, U, n, 0.0, X,n);

		bool fail = false;
		for (size_t i=0; i<m; ++i)
			for (size_t j=0; j<n; ++j)
				if (!F.areEqual (*(L+i*n+j), *(X+i*n+j)))
					fail=true;

		if (fail) {
			std::cerr<<"FAIL"<<std::endl;
			pass = false;
			break;
		}
	}

	FFLAS::fflas_delete( U);
	FFLAS::fflas_delete( L);
	FFLAS::fflas_delete( X);
	FFLAS::fflas_delete( B);
	FFLAS::fflas_delete( A);
	FFLAS::fflas_delete( P);
	FFLAS::fflas_delete( Q);
	return pass;
}
template<class Field>
bool
test_redcolechelon(Field &F, size_t m, size_t n, size_t r, size_t iters)
{
	typedef typename Field::Element Element ;
	Element * A = FFLAS::fflas_new<Element>(m*n);
	Element * B = FFLAS::fflas_new<Element>(m*n);
	Element * L = FFLAS::fflas_new<Element>(m*n);
	Element * U = FFLAS::fflas_new<Element>(n*n);
	Element * X = FFLAS::fflas_new<Element>(m*n);     
	size_t lda = n; //!@todo check lda

	size_t *P = FFLAS::fflas_new<size_t>(n);
	size_t *Q = FFLAS::fflas_new<size_t>(m);
	size_t R = (size_t)-1;

	bool pass=true;

	for (size_t  l=0;l<iters;l++){
		R = (size_t)-1;
		RandomMatrixWithRank(F,A,r,m,n,lda);
		FFLAS::fassign(F,m,n,A,lda,B,lda);
		for (size_t j=0;j<n;j++) P[j]=0;
		for (size_t j=0;j<m;j++) Q[j]=0;

		R = FFPACK::ReducedColumnEchelonForm (F, m, n, A, n, P, Q);

		if (R != r) {pass = false; break;}

		FFPACK::getReducedEchelonTransform (F, FFLAS::FflasLower, m,n,R,P,A,lda,U,n);

		FFPACK::getReducedEchelonForm (F, FFLAS::FflasLower, m,n,R,Q,A,n,L,n);

		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,n, 1.0, B, n, U, n, 0.0, X,n);

		bool fail = false;
		for (size_t i=0; i<m; ++i)
			for (size_t j=0; j<n; ++j)
				if (!F.areEqual (*(L+i*n+j), *(X+i*n+j)))
					fail=true;

		if (fail) {
			std::cerr<<"FAIL"<<std::endl;
			pass = false;
			break;
		}
	}

	FFLAS::fflas_delete( U);
	FFLAS::fflas_delete( L);
	FFLAS::fflas_delete( X);
	FFLAS::fflas_delete( B);
	FFLAS::fflas_delete( A);
	FFLAS::fflas_delete( P);
	FFLAS::fflas_delete( Q);
	return pass;
}
int main(int argc, char** argv){
	std::cerr<<std::setprecision(20);

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

	bool pass = true ;
	typedef Modular<double> Field;
	Field F(p);
	pass &= test_colechelon(F,m,n,r,iters);
	pass &= test_redcolechelon(F,m,n,r,iters);
	return !pass ;
}
