/*
 * Copyright (C) 2014 FFLAS-FFPACK
 * Written by :
 *        Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

// #define SIMD_INT

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <typeinfo>
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/utils/args-parser.h"

#include "fflas-ffpack/utils/test-utils.h"
#include "assert.h"

#include <random>
#include <chrono>
using namespace FFLAS;
template<class Field>
bool test_freduce (const Field & F, size_t m, size_t k, size_t n, bool timing, uint64_t seed)
{
	typedef typename Field::Element T ;
	size_t repet = 3 ;

	T * A = FFLAS::fflas_new<T>(m*n);
	T * B = FFLAS::fflas_new<T>(m*n);

	Givaro::ModularBalanced<T> E(101);
	typename Givaro::ModularBalanced<T>::RandIter G(E,seed);

	if (timing)	std::cout << ">>>" << std::endl ;
	if (timing)	std::cout << "=== inc == 1 ===" << std::endl ;

	FFLAS::Timer chrono, tim, tom ;
	tim.clear(); tom.clear();
	if (timing)		F.write(std::cout << "Field ") << std::endl;
	for (size_t b = 0 ; b < repet ; ++b) {
		FFPACK::RandomMatrix(E, m, k, A, n, G);
		FFLAS::fassign(E,m,k,A,n,B,n);

		chrono.clear();chrono.start();
		for (size_t i = 0 ; i < m ; ++i) {
			for (size_t j = 0 ; j < k ; ++j) {
				F.init(A[i*n+j], A[i*n+j]);
			}
		}
		chrono.stop();
		tim += chrono ;

		chrono.clear();chrono.start();
		FFLAS::freduce (F,m,k,B,n);
		chrono.stop();
		tom += chrono ;

#if 1
		for (size_t i =0 ; i < m ; ++i)
			for (size_t j =0 ; j < k ; ++j)
				if (! F.areEqual(B[i*n+j],A[i*n+j])) {
					F.write(std::cout) << std::endl <<  i << ',' << j  << " : ";
					F.write(std::cout, B[i*n+j]) << "!= (ref)";
					F.write(std::cout, A[i*n+j]) << std::endl;
					return false ;
				}
#endif
	}

	if (timing)	std::cout  << " freduce (___): " << tim.usertime()/(double)repet << 's' << std::endl;
	if (timing)	std::cout  << " freduce (AVX): " << tom.usertime()/(double)repet << 's'<<  std::endl << std::endl;

	if (timing)	std::cout << "=== inc != 1 ===" << std::endl ;

	tim.clear() ; tom.clear();
	if (timing)	F.write(std::cout << "Modular ") << std::endl;
	for (size_t b = 0 ; b < repet ; ++b) {
		FFPACK::RandomMatrix(E, m, n, A, n, G);
		FFLAS::fassign(E,m,n,A,n,B,n);
		size_t incX = 2 ;

		chrono.clear();chrono.start();
		for (size_t i = 1 ; i < m*n ; i += incX) {
			F.init(A[i], A[i]);
		}
		chrono.stop();
		tim += chrono ;

		size_t cnt = (size_t)floor((double)(m*n)/(double)incX) ;

		chrono.clear();chrono.start();
		FFLAS::freduce (F,cnt,B+1,incX);
		chrono.stop();
		tom += chrono ;

#if 1
		for (size_t i =1 ; i < m*n ; i+=incX)
			if (! F.areEqual(B[i],A[i])) {
				F.write(std::cout) << std::endl <<  i << " : ";
				F.write(std::cout, B[i]) << "!= (ref)";
				F.write(std::cout, A[i]) << std::endl;
				return false;
			}
#endif

	}

	if (timing)	std::cout <<  " freduce (___): " << tim.usertime()/(double)repet << 's' << std::endl;
	if (timing)	std::cout <<  " freduce (AVX): " << tom.usertime()/(double)repet << 's'<<  std::endl << std::endl;

	if (timing)	std::cout << "<<<" << std::endl;

	FFLAS::fflas_delete( A );
	FFLAS::fflas_delete( B);

	return true;
}

template<class Field>
bool run_with_field (Givaro::Integer q, size_t b, size_t m, size_t k, size_t n, size_t iters, bool timing, uint64_t seed){
	size_t nbit = iters;
	bool ok = true;
	while (ok && nbit){
		Field* F= FFPACK::chooseField<Field>(q,b,seed);
		if (F==NULL) return true;
		std::ostringstream oss;
		F->write(oss);

		std::cout.fill('.');
		std::cout<<"Checking ";
		std::cout.width(45);
		std::cout<<oss.str();
		std::cout<<"... ";

		ok = ok && test_freduce (*F,m,k,n,timing, seed++);
		if (!ok)
			std::cout << "FAILED "<<std::endl;
		else
			std::cout << "PASSED "<<std::endl;
		delete F;
		nbit--;
	}
	if (!ok)
		std::cout << "with seed = "<< seed << std::endl;
	return ok;
}

int main(int ac, char **av) {
	Givaro::Integer q = -1;
	size_t b = 0;
	size_t m = 297 ;
	size_t n = 301 ;
	size_t k = 299 ;
	size_t iters=3;
	uint64_t seed = getSeed();
	bool loop=false;
	bool timing = false ;

	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic.",  TYPE_INTEGER , &q },
		{ 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
		{ 'n', "-n N", "Set the number of cols in C.",   TYPE_INT , &n },
		{ 'm', "-m N", "Set the number of rows in C.",   TYPE_INT , &m },
		{ 'k', "-k N", "Set the number of rows in B.",   TYPE_INT , &k },
		{ 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
		{ 's', "-s N", "Set the seed                 .", TYPE_UINT64 , &seed },
		{ 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
		{ 't', "-timing", "Output timings"            , TYPE_NONE, &timing},
		END_OF_ARGUMENTS
	};


	FFLAS::parseArguments(ac,av,as);

	if (n < k) {
		std::cout << "Usage : m k n ; matrix of size m x k, lda is n" << std::endl;
		return -1 ;
	}


	bool pass  = true ;
	do{
		pass = pass && run_with_field<Givaro::Modular<float> >(q,b,m,k,n,iters,timing,seed);
		pass = pass && run_with_field<Givaro::Modular<double> >(q,b,m,k,n,iters,timing,seed);
		pass = pass && run_with_field<Givaro::ModularBalanced<float> >(q,b,m,k,n,iters,timing,seed);
		pass = pass && run_with_field<Givaro::ModularBalanced<double> >(q,b,m,k,n,iters,timing,seed);
		pass = pass && run_with_field<Givaro::Modular<int32_t> >(q,b,m,k,n,iters,timing,seed);
		pass = pass && run_with_field<Givaro::ModularBalanced<int32_t> >(q,b,m,k,n,iters,timing,seed);
		pass = pass && run_with_field<Givaro::Modular<int64_t> >(q,b,m,k,n,iters,timing,seed);
		pass = pass && run_with_field<Givaro::ModularBalanced<int64_t> >(q,b,m,k,n,iters,timing,seed);
		pass = pass && run_with_field<Givaro::ZRing<float> >(q,b,m,k,n,iters,timing,seed);
		pass = pass && run_with_field<Givaro::ZRing<double> >(q,b,m,k,n,iters,timing,seed);
		pass = pass && run_with_field<Givaro::ZRing<int32_t> >(q,b,m,k,n,iters,timing,seed);
		pass = pass && run_with_field<Givaro::ZRing<int64_t> >(q,b,m,k,n,iters,timing,seed);
	} while (loop && pass);
	return !pass ;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
