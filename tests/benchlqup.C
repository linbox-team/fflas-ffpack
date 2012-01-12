/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
//
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

#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "timer.h"
#include "Matio.h"

using namespace std;
using namespace FFPACK;

int main(int argc, char** argv) {

	// parameter: p, n, iteration, file

	float    p    = (float)atof(argv[1]);
	int n    = atoi(argv[2]);
	size_t iter = atoi(argv[3]);


	typedef ModularBalanced<double> Field;
	//  typedef ModularBalanced<float> Field;
	typedef Field::Element Element;

	Field F(p);

	Timer chrono;
	double time=0.0;
	// int singular;

	Element *A;

	for (size_t i=0;i<iter;++i){

		A = new Element[n*n];
		Field::RandIter G(F);
		for (size_t i=0; i< (size_t)n*n; ++i)
			G.random(*(A+i));

		size_t * P = new size_t[n];
		size_t * Q = new size_t[n];

		chrono.clear();
		chrono.start();
		FFPACK::LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, n, n, A, n,
				  P, Q);
		chrono.stop();

		time+=chrono.realtime();
		delete[] P;
		delete[] Q;
		delete[] A;

	}

	cerr<<"n: "<<n<<" p: "<<p<<std::endl
	<<" time:  "<<time/(double)iter<<std::endl
	<<" speed: "<<2/3.0*n/1000.0*n/1000.0*n/1000.0/time*double(iter)<<" Gffops"<<std::endl;


	return 0;
}

