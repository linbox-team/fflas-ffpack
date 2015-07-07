/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

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

#include "Matio.h"
#include "test-utils.h"
#include "assert.h"

template<class Field>
bool test_freduce (const Field & F, size_t m, size_t k, size_t n, bool timing)
{
	typedef typename Field::Element T ;
	size_t repet = 3 ;

	T * A = FFLAS::fflas_new<T>(m*n);
	T * B = FFLAS::fflas_new<T>(m*n);

	Givaro::ModularBalanced<T> E(101);

	if (timing)	std::cout << ">>>" << std::endl ;
	if (timing)	std::cout << "=== inc == 1 ===" << std::endl ;

	FFLAS::Timer chrono, tim, tom ;
	tim.clear(); tom.clear();
	if (timing)		F.write(std::cout << "Field ") << std::endl;
	for (size_t b = 0 ; b < repet ; ++b) {
		FFPACK::RandomMatrix(E,A,m,k,n);
		// RandomMatrix(E,B,m,k,n);
		FFLAS::fassign(E,m,k,A,n,B,n);

		chrono.clear();chrono.start();
		for (size_t i = 0 ; i < m ; ++i)
			for (size_t j = 0 ; j < k ; ++j)
				F.init(A[i*n+j],A[i*n+j]);
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
		FFPACK::RandomMatrix(E,A,m,n,n);
		FFLAS::fassign(E,m,n,A,n,B,n);
		size_t incX = 2 ;

		chrono.clear();chrono.start();
		for (size_t i = 1 ; i < m*n ; i += incX) {
			F.init(A[i],A[i]);
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
				return false ;
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

int main(int ac, char **av) {
	static size_t m = 297 ;
	static size_t n = 301 ;
	static size_t k = 299 ;
	static uint64_t p = 7;
	int seed = (int) time(NULL);
	static bool timing = false ;

	static Argument as[] = {
		{ 'p', "-p P", "Set the field characteristic.",  TYPE_INT , &p },
		{ 'n', "-n N", "Set the number of cols in C.",   TYPE_INT , &n },
		{ 'm', "-m N", "Set the number of rows in C.",   TYPE_INT , &m },
		{ 'k', "-k N", "Set the number of rows in B.",   TYPE_INT , &k },
		{ 's', "-s N", "Set the seed                 .", TYPE_INT , &seed },
		{ 't', "-timing", "Output timings"            , TYPE_NONE, &timing},
		END_OF_ARGUMENTS
	};


	FFLAS::parseArguments(ac,av,as);

	if (n < k) {
		std::cout << "Usage : m k n ; matrix of size m x k, lda is n" << std::endl;
		return -1 ;
	}

	srand(seed);
	srand48(seed);

	bool pass  = true ;
	{ /*  freduce */
		{
			Givaro::Modular<float> F(p) ;
			pass &= test_freduce (F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<float> F(p) ;
			pass &= test_freduce (F,m,k,n,timing);
		}
		{
			Givaro::Modular<double> F(p) ;
			pass &= test_freduce (F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<double> F(p) ;
			pass &= test_freduce (F,m,k,n,timing);
		}
		{
			Givaro::Modular<int32_t> F((int32_t)p) ;
			pass &= test_freduce (F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<int32_t> F((int32_t)p) ;
			pass &= test_freduce (F,m,k,n,timing);
		}
		{
			Givaro::Modular<int64_t> F(p) ;
			pass &= test_freduce (F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<int64_t> F(p) ;
			pass &= test_freduce (F,m,k,n,timing);
		}
#if 1
		{
			Givaro::ZRing<float> F ;
			pass &= test_freduce (F,m,k,n,timing);
		}
		{
			Givaro::ZRing<double> F ;
			pass &= test_freduce (F,m,k,n,timing);
		}
		{
			Givaro::ZRing<int32_t> F;
			pass &= test_freduce (F,m,k,n,timing);
		}
		{
			Givaro::ZRing<int64_t> F ;
			pass &= test_freduce (F,m,k,n,timing);
		}
#endif
	}

	return (pass?0:1) ;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


