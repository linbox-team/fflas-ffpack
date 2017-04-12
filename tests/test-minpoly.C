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
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "test-utils.h"
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include <givaro/modular-integer.h>
#include <givaro/givpoly1factor.h>
#include <givaro/givpoly1.h>

using namespace std;
using namespace FFPACK;
using Givaro::Modular;
using Givaro::ModularBalanced;

typedef Givaro::ModularBalanced<double> Field;
typedef vector<Field::Element> Polynomial;

/* Computes P(A)*V and stores it in E using Horner's scheme for
 * polynomial evaluation. */
template<typename Field, typename Element>
void horner_matrix_vector(const Field &F, size_t n, Element *A, size_t lda, Element *V, 
						  Element *E, Polynomial P)
{

	//TODO: Quite a few copies. Could be improved.
	size_t deg = P.size() - 1;
	Element *E_tmp;
	E_tmp = FFLAS::fflas_new(F, 1, n);

	FFLAS::fassign(F, n, V, 1, E, 1);
	FFLAS::fscalin(F, n, P[deg], E, 1); // E <- P[deg+1] * V

    for(long i = deg; i > 0; --i)
    {
		FFLAS::fassign(F, n, E, 1, E_tmp, 1); //E_tmp <- E
		FFLAS::fgemv(F, FFLAS::FflasNoTrans, n, n, F.one, A, lda, E_tmp, 1, F.zero, E, 1);//E <- A * E_tmp
		FFLAS::faxpy(F, 1, n, P[i-1], V, 1, E, 1); //E <- minP[i-1] * V + E
    }	
	FFLAS::fflas_delete(E_tmp);
}


template<typename Field, class RandIter>
bool check_minpoly(const Field &F, size_t n, RandIter& G)
{
	cout<<"Entering check_minpoly"<<endl;
	typedef typename Field::Element Element;
	size_t lda, ldv;
	Element *A, *V, *Vcst;

	//Default
	lda = n;
	ldv = n; 

	/*Create variables used for testing (matrices, vectors and polynomials) */

    A = FFLAS::fflas_new(F, n, n);
    V = FFLAS::fflas_new(F, n+1, n);
	Vcst = FFLAS::fflas_new(F, 1, n);
    Polynomial minP;


    FFPACK::RandomMatrix (F, n, n, A, lda, G);

	FFPACK::NonZeroRandomMatrix(F, 1, n, V, n, G); 
	FFLAS::fassign(F, n, V, 1, Vcst, 1); //MatVecMinPoly modifies V, we store it in Vcst beforehand

	FFPACK::MatVecMinPoly(F, minP, n, A, lda, V, ldv);
	FFLAS::fflas_delete(V);

	/*Check that minP is monic*/

	size_t deg = minP.size() - 1;
	if(!(minP[deg]==F.one))
		return false;

	/*Check that minP(A).V is zero*/


	Element *E;
    E = FFLAS::fflas_new(F, 1, n);
    
	horner_matrix_vector(F, n, A, lda, Vcst, E, minP);
    
    if (!FFLAS::fiszero(F, n, E, 1))
	{
		cout<<"NONZEROERROR"<<endl;
		FFLAS::fflas_delete(E);
		return false;
	}

	FFLAS::fflas_delete(E);


	/* Check minimality of minP */

	// Krylov matrix computation
	Element *K, *tmp;
	size_t ldk = n;
	K = FFLAS::fflas_new(F, deg+1, ldk);
	tmp = FFLAS::fflas_new(F, 1, n);
	FFLAS::fassign(F, n, K, 1, Vcst, 1);
	Element *Kptr = K;

	for(size_t i = 0; i < deg; ++i, Kptr += ldk)
	{
		FFLAS::fgemv(F, FFLAS::FflasNoTrans, n, n, F.one, A, lda, Kptr, 1, F.zero, tmp, 1);
		FFLAS::fassign(F, n, tmp, 1, Kptr+ldk, 1);
	}
	FFLAS::fflas_delete(tmp);

	// minP factorization
	typedef Givaro::Poly1FactorDom<Field, Givaro::Dense> PolyDom; //defines a polynomial domain for Givaro
	typedef typename PolyDom::Element FieldPoly; //defines an element over this polynomial domain (casting purposes)
	vector<FieldPoly> factors;
	vector<size_t> powers;

	PolyDom PD(F);
	FieldPoly FP = FieldPoly(minP.begin(), minP.end());
	PD.factor(factors, powers, FP);	

	// Factorized minP checks
	//divide minP by each factor, and evaluate it. None shall pass eval==0.
	//call horner_matrix_vector to evaluate.
	


	FFLAS::fflas_delete(A);
	FFLAS::fflas_delete(Vcst);
	FFLAS::fflas_delete(K);
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
