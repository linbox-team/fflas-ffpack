/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
//#include "goto-def.h"

/*
 * Copyright (c) FFLAS-FFPACK
 * Written by Clement Pernet <clement.pernet@imag.fr>
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
 */

#include <iostream>

#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/field/modular-positive.h"
#include "fflas-ffpack/utils/timer.h"
#include "Matio.h"

using namespace std;
using namespace FFPACK;

int main(int argc, char** argv)
{

	// parameter: p, n, iteration, file1, file2

	double    p    = atof(argv[1]);
	int n    = atoi(argv[2]);
	size_t w = atoi (argv[3]);
	size_t iter = atoi(argv[4]);

//	typedef Givaro::Modular<float> Field;
//	typedef Givaro::Modular<double> Field;
//	typedef ModularBalanced<double> Field;
      typedef ModularBalanced<float> Field;
	typedef Field::Element Element;

	Field F((Field::Element)p);
	Element one,zero;
	F.init(one, 1.0);
	F.init(zero,0.0);

 FFLAS::Timer chrono;
	double time=0.0;
	// double time2=0.0;
	// int singular;

	Element * A, * B, * C;

	for (size_t i=0;i<iter;++i){

		Field::RandIter G(F);
		A = FFLAS::fflas_new<Element>(n*n);
		for (size_t i=0; i<(size_t)n*n; ++i)
			G.random (*(A+i));

		B = FFLAS::fflas_new<Element>(n*n);
		for (size_t i=0; i<(size_t)n*n; ++i)
			G.random(*(B+i));

		C = FFLAS::fflas_new<Element>(n*n);

		chrono.clear();
		chrono.start();
		FFLAS::MMHelper<Field, FFLAS::MMHelperAlgo::Winograd> WH (F,w);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n,n,n, one,
			      A, n, B, n, zero, C,n, WH);
		chrono.stop();
		time+=chrono.realtime();

		FFLAS::fflas_delete( A);
		FFLAS::fflas_delete( B);
		FFLAS::fflas_delete( C);
	}

	std::cerr<<"n: "<<n <<" p: "<<p<<" w: "<<w<<std::endl
	<<" time:  "<<time/(double)iter<<" s"<<std::endl
	<<" speed: "<<2.0*n/1000.0*n/1000.0/time*n/1000.0*double(iter)<<" Gffops"
	<<std::endl;

	return 0;
}

