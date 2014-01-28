/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
//#include "goto-def.h"

/* Copyright (c) 2012 FFLAS-FFPACK
 * Written by J.G. Dumas <jgdumas@imag.fr>
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

//#define LinBoxSrcOnly
#include <iostream>
#include <fstream>
    //#define _LINBOX_LINBOX_CONFIG_H
// #define __FFLASFFPACK_CONFIGURATION
// #include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/field/modular-positive.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/timer.h"

#define CUBE(x) ((x)*(x)*(x))

template<class Field>
void launch_wino(const Field  &F,
		 const size_t &n,
		 const size_t &NB,
		 const size_t &winomax,
		 const size_t &seed)
{

	typedef typename Field::Element Element ;
	typename Field::RandIter G(F);
	F.write(std::cout<< "Field " ) << std::endl;

	double basetime(0.0), time(0.0);

	Element *A, *C;
	A = new Element[n*n];
	C = new Element[n*n];
	for (size_t i=0; i<n*n;++i)
		G.random(A[i]);


	Timer chrono;
	for(size_t i=0; i<NB; ++i) {
		chrono.start();
		FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
			     n, n, n, F.one, A, n, A, n, F.zero, C, n);
		chrono.stop();
		basetime+= chrono.usertime();
	}
	std::cout << std::endl
	<< "fgemm " << n << "x" << n << ": "
	<< basetime/(double)NB << " s, "
	<< (2.0/basetime*(double)NB*CUBE((double)n/100.0)) << " Mffops"
	<< std::endl;

	for(size_t w=0; w<winomax; ++w) {

		time = 0. ;
		chrono.clear();
		for(size_t i=0; i<NB; ++i) {
			chrono.start();
			FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				     n, n, n, F.one, A, n, A, n, F.zero, C, n, w);
			chrono.stop();
			time+= chrono.usertime();
		}
		std::cout << w << "Wino " << n << "x" << n << ": "
		<< time/(double)NB << " s, "
		<< (2.0/time*(double)NB*CUBE((double)n/100.0)) << " Mffops"
		<< std::endl;
	}

	std::cout << std::endl;
	std::cout << std::endl;

	delete[] A;
	delete[] C;
}

    //using namespace LinBox;
int main (int argc, char ** argv) {
        const size_t p = argc>1 ? atoi(argv[1]) : 65521;
        const size_t n = argc>2 ? atoi(argv[2]) : 1000;
        const size_t NB = argc>3 ? atoi(argv[3]) : 1;
        const size_t winomax = argc>4 ? atoi(argv[4]) : 8;
        const size_t seed = argc>5 ? atoi(argv[5]) : BaseTimer::seed() ;
        srand((unsigned int)seed);

	using namespace FFPACK;
	Modular<double> F1(p);
	Modular<float>  F2(p);
	Modular<int>    F3(p);
	ModularBalanced<double> F4(p);
	ModularBalanced<float>  F5(p);
	ModularBalanced<int>    F6((int)p);
	//! @bug no randiter in UnparametricField !!
	// UnparametricField<double> F7;
	// UnparametricField<float>  F8;
	// UnparametricField<int>    F9;


	launch_wino(F1,n,NB,winomax,seed);
	launch_wino(F2,n,NB,winomax,seed);
	launch_wino(F3,n,NB,winomax,seed);
	launch_wino(F4,n,NB,winomax,seed);
	launch_wino(F5,n,NB,winomax,seed);
	launch_wino(F6,n,NB,winomax,seed);
	// launch_wino(F7,n,NB,winomax,seed);
	// launch_wino(F8,n,NB,winomax,seed);
	// launch_wino(F9,n,NB,winomax,seed);

        return 0;
    }

