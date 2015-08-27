/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2012 FFLAS-FFPACK group.
 *
 * Extirpé form a m4 macro by Brice Boyer (briceboyer) <boyer.brice@gmail.com>.
 *
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
 *
 */


//#define LinBoxSrcOnly
#define DOUBLE_TO_FLOAT_CROSSOVER 0

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <fstream>
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"

#ifndef FLTTYPE
#define FLTTYPE Givaro::Modular<double>
#endif

template<class Field>
bool balanced(const Field & )
{
	return false;
}

template <class T>
bool balanced(const Givaro::ModularBalanced<T>&)
{
	return true;
}

#ifdef __GIVARO_USE_OPENMP
typedef Givaro::OMPTimer TTimer;
#else
typedef Givaro::Timer TTimer;
#endif

#define MFLOPS (2.0*iter/chrono.realtime()*(double)n/100.0*(double)n/100.0*(double)n/100.0)
#define GFLOPS (2.0*iter/chrono.realtime()*(double)n/1000.0*(double)n/1000.0*(double)n/1000.0)

#ifdef __FFLASFFPACK_HAVE_CXX11
#include <ctime>
#endif

//using namespace LinBox;
int main () {
	using namespace std;

	typedef FLTTYPE Field ;
	Field F(17);
	typedef Field::Element Element ;
	size_t n=768, nmax=5000, prec=512, nbest=0, count=0;
    TTimer chrono;
	bool bound=false;
	Field::RandIter G(F); 

	Element *A,*B,*C;
	A = FFLAS::fflas_new<Element>(nmax*nmax);
	B = FFLAS::fflas_new<Element>(nmax*nmax);
	C = FFLAS::fflas_new<Element>(nmax*nmax);
	for (size_t i=0; i<nmax*nmax;++i)
		G.random(A[i]);

	for (size_t i=0; i<nmax*nmax;++i)
		G.random(B[i]);

	for (size_t i=0; i<nmax*nmax;++i)
		G.random(C[i]);
	

	std::ofstream outlog;
	outlog.open("optim.log", std::ofstream::out | std::ofstream::app);
#ifdef __FFLASFFPACK_HAVE_CXX11
    std::time_t result = std::time(NULL);
    outlog << std::endl <<
        "---------------------------------------------------------------------"
           << std::endl << std::asctime(std::localtime(&result));
#endif
	outlog << std::endl
		<< "Threshold for finite field Strassen-Winograd matrix multiplication" ;
	F.write(outlog << "(using ") << ')' << std::endl;
	do {
	double basetime, time;
		FFLAS::MMHelper<Field, FFLAS::MMHelperAlgo::Winograd> ClassicH(F,0, FFLAS::ParSeqHelper::Sequential());
		FFLAS::MMHelper<Field, FFLAS::MMHelperAlgo::Winograd> WinogradH(F,1, FFLAS::ParSeqHelper::Sequential());

		int iter=3;
		    //warm up computation
		FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				n, n, n, F.mOne, A, n, B, n, F.one, C, n, ClassicH);
		chrono.start();
		for (int i=0;i<iter;i++)
			FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n, F.mOne, A, n, B, n, F.one, C, n, ClassicH);
		chrono.stop();
		std::cout << std::endl
			<< "fgemm " << n << "x" << n << ": "
			<< chrono.realtime()/iter << " s, "
			<< GFLOPS << " Gffops"
			<< std::endl;
		outlog << std::endl
			<< "fgemm " << n << "x" << n << ": "
			<< chrono.realtime()/iter << " s, "
			<< GFLOPS << " Gffops"
			<< std::endl;
		basetime= chrono.realtime();
		    //warm up
		FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
			     n, n, n, F.mOne, A, n, B, n, F.one, C, n, WinogradH);
		chrono.clear();
		chrono.start();
		for (int i=0; i<iter; i++)
			FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				     n, n, n, F.mOne, A, n, B,n, F.one, C, n, WinogradH);
		chrono.stop();
		std::cout << "1Wino " << n << "x" << n << ": "
			<< chrono.realtime()/iter << " s, "
			<< GFLOPS  << " Gffops"
			<< std::endl;
		outlog << "1Wino " << n << "x" << n << ": "
			<< chrono.realtime()/iter << " s, "
			<< GFLOPS  << " Gffops"
			<< std::endl;
		time= chrono.realtime();

		if (basetime > time ){
			count++;
			if (count > 1){
				nbest=n;
				bound=true;
				prec=prec>>1;
				n-=prec;
			}
		}
		else{
			count=0;
			if (bound)
				prec=prec>>1;
			n+=prec;
		}
	} while ((prec > 64 ) && (n < nmax));

	std::ofstream out("WinoThreshold");
	if (nbest != 0 ) {
	if (typeid(Element).name() == typeid(double).name()) {
		if ( balanced(F) ) {
			out << "#ifndef __FFLASFFPACK_WINOTHRESHOLD_BAL"  << endl;
			out << "#define __FFLASFFPACK_WINOTHRESHOLD_BAL" << ' ' <<  nbest << endl;
		}
		else {
			out << "#ifndef __FFLASFFPACK_WINOTHRESHOLD"  << endl;
			out << "#define __FFLASFFPACK_WINOTHRESHOLD" << ' ' <<  nbest << endl;
		}
		out << "#endif"                               << endl  << endl;
	}

	if (typeid(Element).name() == typeid(float).name()) {
		if ( balanced(F) ) {
			out << "#ifndef __FFLASFFPACK_WINOTHRESHOLD_BAL_FLT"  << endl;
			out << "#define __FFLASFFPACK_WINOTHRESHOLD_BAL_FLT" << ' ' << nbest << endl;
		}
		else {
			out << "#ifndef __FFLASFFPACK_WINOTHRESHOLD_FLT"  << endl;
			out << "#define __FFLASFFPACK_WINOTHRESHOLD_FLT" << ' ' << nbest << endl;
		}
		out << "#endif"                               << endl  << endl;
	}
	}
	out.close();

	outlog << "defined __FFLASFFPACK_WINOTHRESHOLD to " << nbest << "" << std::endl;
	outlog.close();

	FFLAS::fflas_delete( A);
	FFLAS::fflas_delete( C);

	return 0;
}
