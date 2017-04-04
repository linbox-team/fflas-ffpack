/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) FFLAS-FFPACK
 * Written by David Lucas
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
//                        Test for minpoly
//
//--------------------------------------------------------------------------
// David Lucas
//-------------------------------------------------------------------------

#include <iomanip>
#include <iostream>

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "test-utils.h"
#include <givaro/modular-integer.h>
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include <givaro/givpoly1.h>

using namespace std;
using namespace FFPACK;
using Givaro::Modular;
using Givaro::ModularBalanced;

typedef Givaro::ModularBalanced<double> Field;
typedef vector<Field::Element> Polynomial;

template<typename Field, class RandIter>
bool check_minpoly(const Field &F, size_t n, RandIter& G)
{
	cout<<"Entering check_minpoly";
	cout<<endl;
	typedef typename Field::Element Element;
	size_t lda, ldv;
	Element *A, *V;

	//Default
	lda = n;
	ldv = n; 

	/*Create variables used for testing (matrices, vectors and polynomials) */

    A = FFLAS::fflas_new(F, n, n);
    V = FFLAS::fflas_new(F, n+1, n);
    Polynomial minP;


    FFPACK::RandomMatrix (F, n, n, A, lda, G);

	FFPACK::NonZeroRandomMatrix(F, 1, n, V, n, G); 

	FFPACK::MatVecMinPoly(F, minP, n, A, lda, V, ldv); //3rd input argument is the matrix order

	/*Check that minP is monic*/

	size_t deg = minP.size() - 1;
	if(!(minP[deg]==F.one))
		return false;

	/*Check that minP(A).V is zero*/


	Element *E, *E_tmp;
    E = FFLAS::fflas_new(F, n+1, n);
	E_tmp = FFLAS::fflas_new(F, n+1, n);
    
    /* Horner's method for polynomial evaluation */
    
	FFLAS::finit(F, n, V, n, E, n); //E <- V

    for(long i = deg; i > 0; --i)
    {
		FFLAS::faxpy(F, F.one, n, minP[i-1], V, n, E, n); //E <- minP[i-1] * V + E
		FFLAS::fassign(F, 1, n, E_tmp, n, E, n); //E_tmp <- E
		FFLAS::fgemv(F, FFLAS::FflasNoTrans, n, n, F.one, A, n, E_tmp, n, F.zero, E, n);//E <- E_tmp * A
    }
	
	FFLAS::fflas_delete(E_tmp);
	FFLAS::fflas_delete(A);
	FFLAS::fflas_delete(V);

    if (!FFLAS::fiszero(F, n, E, n))
	{
		cout<<"NONZERO"<<endl;
		for(long i = 0; i<n; ++i)
			cout<<E[i]<<" ";
		cout<<endl;
		FFLAS::fflas_delete(E);
		return false;
	}

	FFLAS::fflas_delete(E);
	return true;
}

template <class Field>
bool run_with_field (Givaro::Integer q, size_t b, size_t n, size_t iters, uint64_t seed)
{
	bool ok = true;
	int nbiter = iters;

	while (ok && nbiter)
	{
		Field* F = chooseField<Field>(q, b); // F, characteristic q of b bits
		typename Field::RandIter G(*F, 0, seed); //random generator over F
		
		if(F == nullptr)
			return true; //if F is null, nothing to test, just pass

		cout<<"Checking with "; F->write(cout)<<endl;

		ok = ok && check_minpoly(*F, n, G);

		if(!ok)
			cout<<"FAILED"<<endl;
		else
			cout<<"PASS"<<endl;

		delete F;
		nbiter--;
	}


	return ok;
}


int main(int argc, char** argv)
{
    /* Test parameters */
	Givaro::Integer q = -1;
	size_t b = 0;
    size_t n = 128;
	size_t iters = 1;
	bool loop = false;
    uint64_t seed = time(NULL);

	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",	TYPE_INTEGER, &q },
		{ 'b', "-b B", "Set the bitsize of the field characteristic.", TYPE_INT, &b },
		{ 'n', "-n N", "Set the order of the matrix.", TYPE_INT, &b },
		{ 'i', "-i, R", "set the number of repetitions.", TYPE_INT, &iters },
		{ 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
		{ 's', "-s seed", "set seed for the random generator.", TYPE_INT, &seed },
			END_OF_ARGUMENTS
		};

	FFLAS::parseArguments(argc, argv, as);

	bool ok = true;

	do
	{
		ok &= run_with_field<Modular<double>>(q,b,n,iters,seed);
		//more tests?
	} while(ok && loop);

	return !ok ;
}
