/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2016 FFLAS-FFPACK
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
//          Test for the computations of the LDLT factorization
//--------------------------------------------------------------------------

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"

#include <iostream>
#include <iomanip>
#include <givaro/modular.h>

#include "test-utils.h"
#include "fflas-ffpack/utils/Matio.h"
using namespace FFPACK;
using namespace FFLAS;

template<class Field>
bool run_with_field(Givaro::Integer q, uint64_t b, size_t n, size_t iters, uint64_t seed){
	bool ok = true ;
	int nbit=(int)iters;
	
	while (ok &&  nbit){
		// choose Field 
		Field* F= chooseField<Field>(q,b);
		if (F==nullptr)
			return true;

		typename Field::RandIter G (*F,b,seed);
		std::ostringstream oss;
		F->write(oss);
		
		std::cout.fill('.');
		std::cout<<"Checking ";
		std::cout.width(40);
		std::cout<<oss.str();
		std::cout<<" ... ";

		size_t lda = n+13;
		typename Field::Element_ptr A=FFLAS::fflas_new (*F, n,lda);
		typename Field::Element_ptr B=FFLAS::fflas_new (*F, n,lda);
		RandomSymmetricMatrix (*F, n, true, A, lda, G);

		FFLAS::fassign (*F, n, n, A, lda, B, lda); 

		{ // Testing is B ==  L D^-1 L^T
			FFPACK::fsytrf (*F, FflasLower, n, A, lda);

				// copying L on L^T
			for (size_t i=0; i<n; i++)
				fassign(*F, n-i-1, A+i*(lda+1)+lda, lda, A+i*(lda+1)+1, 1);

				// L^T <- D^-1 L^T 
			for (size_t i=0; i<n; i++){
				F->invin(A[i*(lda+1)]);
				fscalin (*F, n-i-1, A[i*(lda+1)], A+i*(lda+1)+1, 1);
			}
				// A <- L x L^T 
			ftrtrm(*F, FflasRight, FflasNonUnit, n, A, lda);

			ok &= fequal(*F, n, n, A, lda, B, lda);
		}

		{ // Testing is B ==  U^T D^-1 U
			fassign (*F, n, n, B, lda, A, lda);
			FFPACK::fsytrf (*F, FflasUpper, n, A, lda);

				// copying U on U^T
			for (size_t i=0; i<n; i++)
				fassign(*F, n-i-1, A+i*(lda+1)+1, 1, A+i*(lda+1)+lda, lda);

				// U <- D^-1 U 
			for (size_t i=0; i<n; i++){
				F->invin(A[i*(lda+1)]);
				fscalin (*F, n-i-1, A[i*(lda+1)], A+i*(lda+1)+1, 1);
			}
				// A <- U^T x U
			ftrtrm(*F, FflasRight, FflasNonUnit, n, A, lda);

			ok &= fequal(*F, n, n, A, lda, B, lda);
		}
		
		FFLAS::fflas_delete(A);
		FFLAS::fflas_delete(B);
		delete F;

		nbit--;
		if (!ok)
				//std::cout << "\033[1;31mFAILED\033[0m "<<std::endl;		
			std::cout << "FAILED "<<std::endl;
		else
				//std::cout << "\033[1;32mPASSED\033[0m "<<std::endl;
			std::cout << "PASSED "<<std::endl;
	}
	return ok;
	}

int main(int argc, char** argv){
	std::cerr<<std::setprecision(20);

	Givaro::Integer q = -1;
	size_t b = 0;
	size_t n = 280;
	size_t iters = 6 ;
	bool loop=false;
	uint64_t seed=time(NULL);
	Argument as[] = {
		{ 'q', "-q Q", "Set the field cardinality.",         TYPE_INTEGER , &q },
		{ 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
		{ 'n', "-n N", "Set the number of cols in the matrix.", TYPE_INT , &n },
		{ 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
		{ 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
		{ 's', "-s seed", "Set seed for the random generator", TYPE_INT, &seed },
		    // { 'f', "-f file", "Set input file", TYPE_STR, &file },
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(argc,argv,as);

	srand(seed);

	bool ok=true;
	do{
		ok&=run_with_field<Givaro::Modular<float> >           (q,b,n,iters,seed);
		ok&=run_with_field<Givaro::Modular<double> >          (q,b,n,iters,seed);
		ok&=run_with_field<Givaro::ModularBalanced<float> >   (q,b,n,iters,seed);
		ok&=run_with_field<Givaro::ModularBalanced<double> >   (q,b,n,iters,seed);
		ok&=run_with_field<Givaro::Modular<int32_t> >   (q,b,n,iters,seed);
		ok&=run_with_field<Givaro::ModularBalanced<int32_t> >   (q,b,n,iters,seed);
		ok&=run_with_field<Givaro::Modular<int64_t> >   (q,b,n,iters,seed);
		ok&=run_with_field<Givaro::ModularBalanced<int64_t> >   (q,b,n,iters,seed);
			//ok&=run_with_field<Givaro::Modular<Givaro::Integer> >(q,(b?b:128),n/4+1,iters,seed);
	} while (loop && ok);

	return !ok;
}
