/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

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
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include <iostream>
#include <fstream>
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"

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

#define GFOPS(n,t) (2.0/t*(double)n/1000.0*(double)n/1000.0*(double)n/1000.0)

#include <ctime>

int main () {
	using namespace std;

	typedef FIELD Field;
	Field F(17);
	typedef Field::Element Element ;
	size_t n=512, nmax=4000, prec=512, nbest=0, count=0;
	TTimer chrono;
	bool bound=false;

	Element * A = FFLAS::fflas_new (F,nmax,nmax);
	Element * B = FFLAS::fflas_new (F,nmax,nmax);
	Element * C = FFLAS::fflas_new (F,nmax,nmax);
	FFPACK::RandomMatrix (F, A, nmax,nmax,nmax);
	FFPACK::RandomMatrix (F, B, nmax,nmax,nmax);
	FFPACK::RandomMatrix (F, C, nmax,nmax,nmax);

	time_t result = std::time(NULL);
	cout << std::endl 
		  << "---------------------------------------------------------------------"
		  << std::endl << std::asctime(std::localtime(&result))
		  << std::endl
		  << "Threshold for finite field Strassen-Winograd matrix multiplication" ;
	F.write(cout << " (using ") << ')' << endl << endl;

	cout << "fgemm:  n                   Classic                        Winograd 1 level" << std::endl;
	cout << "                    seconds            Gfops          seconds            Gfops" << std::endl;
	FFLAS::MMHelper<Field, FFLAS::MMHelperAlgo::Winograd> ClassicH(F,0, FFLAS::ParSeqHelper::Sequential());
	FFLAS::MMHelper<Field, FFLAS::MMHelperAlgo::Winograd> WinogradH(F,1, FFLAS::ParSeqHelper::Sequential());
	    //warm up computation
	FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n, F.mOne, A, n, B, n, F.one, C, n, ClassicH);
	do {
		double classicTime, winogradTime;

		int iter=3;
		chrono.start();
		for (int i=0;i<iter;i++)
			FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n, F.mOne, A, n, B, n, F.one, C, n, ClassicH);
		chrono.stop();

		classicTime = chrono.realtime()/iter;

		chrono.clear(); chrono.start();
		for (int i=0; i<iter; i++)
			FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n, F.mOne, A, n, B,n, F.one, C, n, WinogradH);
		chrono.stop();
		
		winogradTime = chrono.realtime()/iter;

		cout << "      ";
		cout.width(4);
		cout << n;
		cout << "  ";
		cout.width(15);
		cout << classicTime;
		cout << "  ";
		cout.width(15);
		cout << GFOPS(n, classicTime) << "  ";
		cout.width(15);
		cout << winogradTime;
		cout << "  ";
		cout.width(15);
		cout << GFOPS(n, winogradTime) << endl;

		if (classicTime > winogradTime ){
			count++;
			if (count > 1){
				nbest = n;
				bound = true;
				prec = prec >> 1;
				n -= prec;
			}
		}
		else{
			count=0;
			if (bound)
				prec=prec>>1;
			n+=prec;
		}
	} while ((prec > 32 ) && (n < nmax));

	cout<<endl;
	if (nbest != 0 ) {
		if (typeid(Element).name() == typeid(double).name()) {
			if ( balanced(F) ) {
				cerr << "#ifndef __FFLASFFPACK_WINOTHRESHOLD_BAL"  << endl;
				cerr << "#define __FFLASFFPACK_WINOTHRESHOLD_BAL" << ' ' <<  nbest << endl;
				cout << "defined __FFLASFFPACK_WINOTHRESHOLD_BAL to " << nbest << "" << std::endl;
		}
			else {
				cerr << "#ifndef __FFLASFFPACK_WINOTHRESHOLD"  << endl;
				cerr << "#define __FFLASFFPACK_WINOTHRESHOLD" << ' ' <<  nbest << endl;
				cout << "defined __FFLASFFPACK_WINOTHRESHOLD to " << nbest << "" << std::endl;
	
			}
			std::cerr << "#endif" << endl  << endl;
		}
		
		if (typeid(Element).name() == typeid(float).name()) {
			if ( balanced(F) ) {
				cerr << "#ifndef __FFLASFFPACK_WINOTHRESHOLD_BAL_FLT"  << endl;
				cerr << "#define __FFLASFFPACK_WINOTHRESHOLD_BAL_FLT" << ' ' << nbest << endl;
				cout << "defined __FFLASFFPACK_WINOTHRESHOLD_BAL_FLT to " << nbest << "" << std::endl;

			}
			else {
				cerr << "#ifndef __FFLASFFPACK_WINOTHRESHOLD_FLT"  << endl;
				cerr << "#define __FFLASFFPACK_WINOTHRESHOLD_FLT" << ' ' << nbest << endl;
				cout << "defined __FFLASFFPACK_WINOTHRESHOLD_FLT to " << nbest << "" << std::endl;
			}
			cerr << "#endif" << endl << endl;
		}
	}

	FFLAS::fflas_delete(A);
	FFLAS::fflas_delete(B);
	FFLAS::fflas_delete(C);

	return 0;
}
