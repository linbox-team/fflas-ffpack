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
#include "fflas-ffpack/utils/args-parser.h"


using namespace std;
using namespace FFLAS;
using namespace FFPACK;
int main(int argc, char** argv) {

	size_t iter      = 3;
	int    q         = 131071;
	size_t rows      = 5000;
	size_t cols      = 5000;
	double a         = 1.0;
	// size_t threshold = 64;
	std::string file = "";
  
	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INT , &q },
		{ 'm', "-m M", "Set the row dimension of the matrix C.",                TYPE_INT , &rows },
		{ 'n', "-n N", "Set the column dimension of the matrix C.",             TYPE_INT , &cols },
		{ 'a', "-a A", "Set the value of the coefficient alpha for alpha * B.", TYPE_DOUBLE , &a},
		{ 'i', "-i R", "Set number of repetitions.",                            TYPE_INT , &iter },
		END_OF_ARGUMENTS
	};
  
	FFLAS::parseArguments(argc,argv,as);
  
	typedef Givaro::ModularBalanced<double> Field;
	typedef Field::Element Element;
	
	Field F(q);
	Field::Element_ptr A, B;
	Field::Element_ptr C;
	Element alpha;
	F.init(alpha, a);
  
  
	FFLAS::Timer chrono;
	double time=0.0;
  
	for (size_t i=0;i<=iter;++i){
		A = fflas_new(F, rows, cols);
		size_t lda = cols;
		B = fflas_new(F, rows, cols);
		size_t ldb = cols;
		C = fflas_new(F, rows, cols);
		size_t ldc = cols;
		Field::RandIter G(F);
		RandomMatrix (F, rows, cols, A, lda, G);
		RandomMatrix (F, rows, cols, B, ldb, G);
		chrono.clear();
		if (i) chrono.start();
		fadd(F, rows, cols, A, lda, alpha, B, ldb, C, ldc);
		if (i) chrono.stop();
	
		time+=chrono.usertime();
		FFLAS::fflas_delete(A);
		FFLAS::fflas_delete(B);
		FFLAS::fflas_delete(C);
	}

	// -----------
	// Standard output for benchmark - Alexis Breust 2014/11/14
	std::cout << "Time: " << time / double(iter)
			  << " Gfops: " << ((double(rows)/1000)*(double(cols)/1000))/(1000*time)* double(iter); //(n^2/1000^3)/time * iter
	FFLAS::writeCommandString(std::cout, as) << std::endl;
	return 0;
}
