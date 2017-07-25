/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/* Copyright (c) FFLAS-FFPACK
 * Written by Philippe LEDENT <philippe.ledent@etu.univ-grenoble-alpes.fr>
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
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/paladin/parallel.h"
#include "fflas-ffpack/paladin/fflas_pfinit.h"
#include "fflas-ffpack/utils/args-parser.h"


using namespace std;
using namespace FFLAS;
using namespace FFPACK;
int main(int argc, char** argv) {

	size_t iter      = 3;
	size_t N      = 5000;
 	size_t BS      = 5000;
    size_t q = 101;
	std::string file = "";
  
	Argument as[] = {
		{ 'n', "-n N", "Set the dimension of the matrix C.",             TYPE_INT , &N },
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INT , &q },
		{ 'i', "-i R", "Set number of repetitions.",                            TYPE_INT , &iter },
		{ 's', "-i S", "Size of the integers.",                            TYPE_INT , &BS },
		END_OF_ARGUMENTS
	};
  
	FFLAS::parseArguments(argc,argv,as);
  
    typedef Givaro::ZRing<Givaro::Integer> Field;
    Field F;
    Givaro::GivRandom generator;
    Givaro::IntegerDom IPD;

// 	typedef Givaro::ModularBalanced<double> Field;
// 	Field F(q);
	
	Field::Element_ptr A, B;
	Field::Element d; F.init(d);
  
	Givaro::OMPTimer chrono, time; time.clear();
  
	for (size_t i=0;i<iter;++i){
		A = fflas_new(F, N);
		B = fflas_new(F, N);
#pragma omp parallel for 
        for (size_t j=0; j<N; ++j) {
            IPD.random(generator,A[j],BS);
            IPD.random(generator,B[j],BS);
        }
// 		Field::RandIter G(F);
// 		RandomMatrix (F, 1, N, A, lda, G);
// 		RandomMatrix (F, 1, N, B, ldb, G);

//         std::cerr << '['; for(size_t i=0; i<N; ++i) F.write(std::cerr,A[i]) << ' '; std::cerr << ']' << std::endl;
//         std::cerr << '['; for(size_t i=0; i<N; ++i) F.write(std::cerr,B[i]) << ' '; std::cerr << ']' << std::endl;
        
        F.assign(d, F.zero);

		chrono.clear();
		chrono.start();
//         PAR_BLOCK {
            F.assign(d, pfdot(F, N, A, 1U, B, 1U) );
//         }
		chrono.stop();

        std::cout << chrono 
                  << " Gfops: " << ((double(N)/1000.)/1000.)/(1000.*chrono.realtime())
                  << std::endl;
	
		F.subin(d, fdot(F, N, A, 1U, B, 1U));
		time+=chrono;
		FFLAS::fflas_delete(A);
		FFLAS::fflas_delete(B);
	}
	F.write(std::cerr, d) << std::endl;
	// -----------
	// Standard output for benchmark - Alexis Breust 2014/11/14
	std::cout << "Time * " << iter << ": " << time 
			  << " Gfops: " << ((double(N)/1000.)/1000.)/(1000.*time.realtime())* double(iter); 
    FFLAS::writeCommandString(std::cout, as) << std::endl;
	return 0;
}
