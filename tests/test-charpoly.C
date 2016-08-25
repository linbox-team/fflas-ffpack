/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

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
//                        Test for charpoly
//
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------


#include <iomanip>
#include <iostream>
#include "givaro/modular.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "fflas-ffpack/ffpack/ffpack.h"

#include "fflas-ffpack/utils/args-parser.h"


using namespace std;

typedef Givaro::ModularBalanced<double> Field;

typedef vector<Field::Element> Polynomial;

using namespace FFPACK;

template<class Field>
bool launch_test(const Field & F, size_t n, typename Field::Element * A, size_t lda,
				 size_t nbit, FFPACK::FFPACK_CHARPOLY_TAG CT)
{
	Polynomial charp;
	typename Field::Element trace;
	F.init(trace, F.zero);
	for (size_t i = 0; i < n; i++)
		F.subin (trace, A [i*(lda+1)]);

	FFPACK::CharPoly (F, charp, n, A, lda, CT);

		// Checking det(A) == charp[0]
	typename Field::Element det = FFPACK::Det(F, n, n, A, lda);
	if (!F.areEqual(det,charp[0])){
		std::cerr<<"det = "<<det<<" P["<<0<<"] = "<<charp[n]<<std::endl;
		return false;
	}

		// Checking trace(A) == charp[n-1]
	if (!F.areEqual(trace, charp[n-1])){
		std::cerr<<"trace = "<<trace<<" P["<<n-1<<"] = "<<charp[n-1]<<std::endl;
		return false;
	}

	return true ;

}

int main(int argc, char** argv)
{

	size_t       p = 131071; // characteristic
	size_t       nbit = 10; // repetitions
	size_t        n = 500;
	std::string  file = "" ; // file where
	int variant = 6;

	Argument as[] = {
		{ 'p', "-p P", "Set the field characteristic.", TYPE_INT , &p },
		{ 'n', "-n N", "Set the size of the matrix.", TYPE_INT , &n },
		{ 'i', "-i I", "Set number of repetitions.", TYPE_INT , &nbit },
		{ 'f', "-f file", "Set input file", TYPE_STR, &file },
		{ 'a', "-a algorithm", "Set the algorithm variant", TYPE_INT, &variant },
		END_OF_ARGUMENTS
	};
	FFLAS::parseArguments(argc,argv,as);
	size_t      lda = n;

	FFPACK::FFPACK_CHARPOLY_TAG CT;
	switch (variant){
	    case 0: CT = FfpackLUK; break;
	    case 1: CT = FfpackKG; break;
	    case 2: CT = FfpackDanilevski; break;
	    case 3: CT = FfpackKGFast; break;
	    case 4: CT = FfpackKGFastG; break;
	    case 5: CT = FfpackHybrid; break;
	    case 6: CT = FfpackArithProg; break;
	    default: CT = FfpackLUK; break;
	}
	Field F(p);
	Field::Element * A=NULL;
	bool passed = true;
	if (!file.empty()) {
		const char * filestring = file.c_str();
		A = read_field<Field>(F,const_cast<char*>(filestring),&n,&n);
		passed &= launch_test<Field>(F, n, A, lda, nbit, CT);
		FFLAS::fflas_delete( A);
	}

		/* Random matrix test */
	A = FFLAS::fflas_new(F,n,n);
	for(size_t i = 0;i<nbit;++i){
		A = FFPACK::RandomMatrix(F,A,n,n,n);
		passed &= launch_test<Field>(F, n, A, lda, nbit, CT);
	}
	FFLAS::fflas_delete( A);


	return !passed;
}
