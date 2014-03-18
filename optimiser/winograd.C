/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2012 FFLAS-FFPACK group.
 *
 * Extirpé form a m4 macro by BB <bboyer@imag.fr>.
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
#include <iostream>
#include <fstream>
#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/config.h"
#include "fflas-ffpack/fflas-ffpack-config.h"
// #include "fflas-ffpack/fflas-ffpack-optimise.h"
#include "fflas-ffpack/field/modular-positive.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/timer.h"

#ifndef FLTTYPE
#define FLTTYPE Modular<double>
#endif

template<class Field>
bool balanced(const Field & )
{
	return false;
}

template <class T>
bool balanced(const FFPACK::ModularBalanced<T>&)
{
	return true;
}

//using namespace LinBox;
int main () {
	using namespace std;

	typedef FFPACK:: FLTTYPE Field ;
	Field F(17);
	typedef Field::Element Element ;
	size_t n=1000, nmax=5000, prec=512, nbest=0, count=0;
	Timer chrono;
	bool bound=false;

	Element *A, *C;
	A = new Element[nmax*nmax];
	C = new Element[nmax*nmax];
	for (size_t i=0; i<nmax*nmax;++i){
		A[i]=2.;
	}

	std::ofstream outlog;
	outlog.open("optim.log", std::ofstream::out | std::ofstream::app);
	outlog << std::endl
		<< "Threshold for finite field Strassen-Winograd matrix multiplication" ;
	F.write(outlog << "(using ") << ')' << std::endl;
	do {
	double basetime, time;
		chrono.start();
		FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				n, n, n, F.one, A, n, A, n, F.zero, C, n, 0);
		chrono.stop();
		std::cout << std::endl
			<< "fgemm " << n << "x" << n << ": "
			<< chrono.usertime() << " s, "
			<< (2.0/chrono.usertime()*n/100.0*n/100.0*n/100.0) << " Mffops"
			<< std::endl;
		outlog << std::endl
			<< "fgemm " << n << "x" << n << ": "
			<< chrono.usertime() << " s, "
			<< (2.0/chrono.usertime()*n/100.0*n/100.0*n/100.0) << " Mffops"
			<< std::endl;
		basetime= chrono.usertime();
		chrono.clear();
		chrono.start();
		FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				n, n, n, 1., A, n, A, n, 0., C, n, 1);
		chrono.stop();
		std::cout << "1Wino " << n << "x" << n << ": "
			<< chrono.usertime() << " s, "
			<< (2.0/chrono.usertime()*n/100.0*n/100.0*n/100.0) << " Mffops"
			<< std::endl;
		outlog << "1Wino " << n << "x" << n << ": "
			<< chrono.usertime() << " s, "
			<< (2.0/chrono.usertime()*n/100.0*n/100.0*n/100.0) << " Mffops"
			<< std::endl;
		time= chrono.usertime();

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

	delete[] A;
	delete[] C;

	return 0;
}
