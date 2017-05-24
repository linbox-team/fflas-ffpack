/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/* Copyright (c) FFLAS-FFPACK
* Written by Clement Pernet <clement.pernet@imag.fr>
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
*/

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "fflas-ffpack/utils/args-parser.h"


using namespace std;

template<class Field>
void run_with_field(int q, size_t bits, size_t n, size_t iter, std::string file, int variant){
	Field F(q);
	typedef typename Field::Element Element;
	FFPACK::FFPACK_CHARPOLY_TAG CT;
	switch (variant){
		case 0: CT = FFPACK::FfpackLUK; break;
		case 1: CT = FFPACK::FfpackKG; break;
		case 2: CT = FFPACK::FfpackDanilevski; break;
		case 3: CT = FFPACK::FfpackKGFast; break;
		case 4: CT = FFPACK::FfpackKGFastG; break;
		case 5: CT = FFPACK::FfpackHybrid; break;
		case 6: CT = FFPACK::FfpackArithProg; break;
		default: CT = FFPACK::FfpackLUK; break;
	}
	FFLAS::Timer chrono;
	Element *A;
	double time_charp=0;
	for (size_t i=0;i<iter;++i){
		if (!file.empty()){
			A = read_field (F, file.c_str(), &n, &n);
		}
		else{
			A = FFLAS::fflas_new (F, n, n);
			typename Field::RandIter G (F, bits);
			FFPACK::RandomMatrix (F, n, n, A, n, G);
		}
		typename Givaro::Poly1Dom<Field>::Element cpol(n+1);
		typename Givaro::Poly1Dom<Field> R(F);
		chrono.clear();
		chrono.start();
		FFPACK::CharPoly (R, cpol, n, A, n, CT);
		chrono.stop();

		time_charp+=chrono.usertime();

		FFLAS::fflas_delete( A);
	}
	// -----------
	// Standard output for benchmark - Alexis Breust 2014/11/14
	std::cerr << "n: "<<n<<" bitsize: "<<bits<<" Time: " << time_charp / double(iter)
			  << " Gflops: " << "Irrelevant";
}

int main(int argc, char** argv) {
  
	size_t iter = 1;
	int    q    = 131071;
	size_t bits = 10;
	size_t    n    = 1000;
	std::string file = "";
	int variant =6;

	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for the ring ZZ).",  TYPE_INT , &q },
		{ 'b', "-b B", "Set the bitsize of the random elements.",         TYPE_INT , &bits},
		{ 'n', "-n N", "Set the dimension of the matrix.",               TYPE_INT , &n },
		{ 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &iter },
		{ 'f', "-f FILE", "Set the input file (empty for random).",  TYPE_STR , &file },
		{ 'a', "-a algorithm", "Set the algorithmic variant", TYPE_INT, &variant },

		END_OF_ARGUMENTS
	};

  FFLAS::parseArguments(argc,argv,as);

  if (q > 0){
	  bits = Givaro::Integer(q).bitsize();
	  run_with_field<Givaro::ModularBalanced<double> >(q, bits, n , iter, file, variant);
  } else
	  run_with_field<Givaro::ZRing<Givaro::Integer> > (q, bits, n , iter, file, variant);

  FFLAS::writeCommandString(std::cerr, as) << std::endl;
  return 0;
}

