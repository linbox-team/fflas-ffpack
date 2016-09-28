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
	std::ostringstream oss;
	switch (CT){
	    case FfpackLUK: oss<<"LUKrylov variant"; break;
	    case FfpackKG: oss<<"Keller-Gehrig variant"; break;
	    case FfpackDanilevski: oss<<"Danilevskii variant"; break;
	    case FfpackKGFast:  oss<<"KGFast variant"; break;
	    case FfpackKGFastG: oss<<"KGFastG variant"; break;
	    case FfpackHybrid: oss<<"Hybrid variant"; break;
	    case FfpackArithProg: oss<<"ArithProg variant"; break;
	    default: oss<<"LUKrylov variant"; break;
	}
	F.write(oss<<" over ");
	std::cout.fill('.');
	std::cout<<"Checking ";
	std::cout.width(65);
	std::cout<<oss.str();
	std::cout<<"...";
	Polynomial charp;

	typename Field::Element_ptr B = FFLAS::fflas_new(F, n,n);
	FFLAS::fassign(F, n, n, A, lda, B, n);

	FFPACK::CharPoly (F, charp, n, A, lda, CT);

	FFLAS::fassign(F, n, n, B, n, A, lda);

		// Checking det(A) == charp[0]
	typename Field::Element det = FFPACK::Det(F, n, n, A, lda);
	if (n&1) // p0 == (-1)^n det
		F.negin(det);
	if (!F.areEqual(det,charp[0])){
			//write_field(F, std::cerr<<"B = "<<std::endl,B, n,n,lda);
		std::cerr<<"FAILED: det = "<<det<<" P["<<0<<"] = "<<charp[0]<<std::endl;
		FFLAS::fflas_delete (B);
		return false;
	}

		// Checking trace(A) == charp[n-1]
	typename Field::Element trace;
	F.init(trace, F.zero);
	for (size_t i = 0; i < n; i++)
		F.subin (trace, B [i*(n+1)]);
	if (!F.areEqual(trace, charp[n-1])){
		std::cerr<<"FAILED: trace = "<<trace<<" P["<<n-1<<"] = "<<charp[n-1]<<std::endl;
		FFLAS::fflas_delete (B);
		return false;
	}
	FFLAS::fflas_delete (B);
	std::cout<<"PASSED"<<std::endl;
	return true ;
}

int main(int argc, char** argv)
{
	size_t       q = 131071; // characteristic
	size_t       nbit = 10; // repetitions
	size_t       n = 500;
	std::string  file = "" ; // file where
	int variant = 0; // default value 0: test all variants
	size_t seed =  time(NULL);

	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic.", TYPE_INT , &q },
		{ 'n', "-n N", "Set the size of the matrix.", TYPE_INT , &n },
		{ 'i', "-i I", "Set number of repetitions.", TYPE_INT , &nbit },
		{ 'f', "-f file", "Set input file", TYPE_STR, &file },
		{ 'a', "-a algorithm", "Set the algorithm variant", TYPE_INT, &variant },
		{ 's', "-s seed", "Set seed for the random generator", TYPE_INT, &seed },
		END_OF_ARGUMENTS
	};
	FFLAS::parseArguments(argc,argv,as);
	size_t lda = n;

	FFPACK::FFPACK_CHARPOLY_TAG CT;
	switch (variant){
	    case 1: CT = FfpackLUK; break;
	    case 2: CT = FfpackKG; break;
	    case 3: CT = FfpackDanilevski; break;
	    case 4: CT = FfpackKGFast; break;
	    case 5: CT = FfpackKGFastG; break;
	    case 6: CT = FfpackHybrid; break;
	    case 7: CT = FfpackArithProg; break;
	    default: CT = FfpackLUK; break;
	}
	
	Field F(q);
	Field::Element * A=NULL;
	bool passed = true;
	Field::RandIter G(F,0,seed);

	for(size_t i = 0;i<nbit;++i){
		if (!file.empty()) {
			const char * filestring = file.c_str();
			A = read_field<Field>(F,const_cast<char*>(filestring),&n,&n);
			passed &= launch_test<Field>(F, n, A, lda, nbit, CT);
			FFLAS::fflas_delete( A);
		} else {
				/* Random matrix test */
			A = FFLAS::fflas_new(F,n,n);
			FFLAS::frand (F,G,n,n,A,n);
			if (variant)
				passed &= launch_test<Field>(F, n, A, lda, nbit, CT);
			else{
				passed &= launch_test<Field>(F, n, A, lda, nbit, FfpackLUK);
				//passed &= launch_test<Field>(F, n, A, lda, nbit, FfpackKG); // fails (variant only implemented for benchmarking comparison
				passed &= launch_test<Field>(F, n, A, lda, nbit, FfpackDanilevski);
				//passed &= launch_test<Field>(F, n, A, lda, nbit, FfpackKGFast); // generic: does not work with any matrix
				//passed &= launch_test<Field>(F, n, A, lda, nbit, FfpackKGFastG); // generic: does not work with any matrix
				//passed &= launch_test<Field>(F, n, A, lda, nbit, FfpackHybrid); // fails with small characteristic
				passed &= launch_test<Field>(F, n, A, lda, nbit, FfpackArithProg);
			}
			FFLAS::fflas_delete( A);
		}
		// if ((i+1)*100 % nbit == 0)
		// 	std::cerr<<double(i+1)/nbit*100<<" % "<<std::endl;
	}
	return !passed;
}
