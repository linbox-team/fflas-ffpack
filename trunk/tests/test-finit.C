/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2014 FFLAS-FFPACK
 * Written by :
 *        BB <brice.boyer@lip6.fr>
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


#include "fflas-ffpack/utils/timer.h"
#include "Matio.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/field/modular-int64.h"
#include "fflas-ffpack/field/modular-balanced-int64.h"
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "test-utils.h"
#include "assert.h"
#include "fflas-ffpack/utils/args-parser.h"
#include <typeinfo>

// using namespace FFPACK;
template<class Field>
bool test_finit(const Field & F, size_t m, size_t k, size_t n, bool timing)
{
	typedef typename Field::Element T ;

	T * A = FFLAS::fflas_new<T>(m*n);
	T * B = FFLAS::fflas_new<T>(m*n);

	FFPACK::ModularBalanced<T> E(101);

	if (timing)	std::cout << ">>>" << std::endl ;
	if (timing)	std::cout << "=== inc == 1 ===" << std::endl ;

	for (size_t b = 0 ; b < 3 ; ++b) {
		RandomMatrix(E,A,m,k,n);
		// RandomMatrix(E,B,m,k,n);
		FFLAS::fcopy(E,m,k,A,n,B,n);

	 FFLAS::Timer tim;
	 FFLAS::Timer tom;
	if (timing)		F.write(std::cout << "Field ") << std::endl;
		tim.clear();tim.start();
		for (size_t i = 0 ; i < m ; ++i)
			for (size_t j = 0 ; j < k ; ++j)
				F.init(A[i*n+j],A[i*n+j]);
		tim.stop();
		if (timing)	std::cout << "finit (___): " << tim.usertime() << 's' << std::endl;

		tim.clear();tim.start();
		FFLAS::finit(F,m,k,B,n);
		tim.stop();
		if (timing)	std::cout << "finit (AVX): " << tim.usertime() << 's'<<  std::endl << std::endl;

#if 1
		for (size_t i =0 ; i < m ; ++i)
			for (size_t j =0 ; j < k ; ++j)
				if (! F.areEqual(B[i*n+j],A[i*n+j])) {
					if (timing)	std::cout  <<  i << ',' << j << " : " <<  B[i*n+j] << "!= (ref)" << A[i*n+j] << std::endl;
					return false ;
				}
#endif
	}

	if (timing)	std::cout << "=== inc != 1 ===" << std::endl ;


	for (size_t b = 0 ; b < 3 ; ++b) {
		RandomMatrix(E,A,m,n,n);
		FFLAS::fcopy(E,m,n,A,n,B,n);
		size_t incX = 2 ;

	 FFLAS::Timer tim;
	 FFLAS::Timer tom;
		if (timing)	F.write(std::cout << "Modular ") << std::endl;
		tim.clear();tim.start();
		for (size_t i = 1 ; i < m*n ; i += incX) {
				F.init(A[i],A[i]);
		}

		tim.stop();
		size_t cnt = (size_t)floor((double)(m*n)/(double)incX) ;
		if (timing)	std::cout << "finit (___): " << tim.usertime() << 's' << std::endl;
		tim.clear();tim.start();
		FFLAS::finit(F,cnt,B+1,incX);
		tim.stop();
	if (timing)		std::cout << "finit (AVX): " << tim.usertime() << 's'<<  std::endl << std::endl;

#if 1
		for (size_t i =1 ; i < m*n ; i+=incX)
				if (! F.areEqual(B[i],A[i])) {
				if (timing)		std::cout <<  i << " : " << B[i] << "!= (ref)" << A[i] << std::endl;
					return false ;
				}
#endif

	}

	if (timing)	std::cout << "<<<" << std::endl;
	delete[] A ;
	delete[] B;

	return true;
}

int main(int ac, char **av) {
	static size_t m = 297 ;
	static size_t n = 301 ;
	static size_t k = 299 ;
	static unsigned long p = 7;
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

	// std::cout << seed << std::endl;

	bool pass  = true ;
	{ /*  finit */
	{
		FFPACK:: Modular<float> F(p) ;
		pass &= test_finit(F,m,k,n,timing);
	}
	{
		FFPACK:: ModularBalanced<float> F(p) ;
		pass &= test_finit(F,m,k,n,timing);
	}
	{
		FFPACK:: Modular<double> F(p) ;
		pass &= test_finit(F,m,k,n,timing);
	}
	{
		FFPACK:: ModularBalanced<double> F(p) ;
		pass &= test_finit(F,m,k,n,timing);
	}
	{
		FFPACK:: Modular<int32_t> F(p) ;
		pass &= test_finit(F,m,k,n,timing);
	}
	{
		FFPACK:: ModularBalanced<int32_t> F((int)p) ;
		pass &= test_finit(F,m,k,n,timing);
	}
	{
		FFPACK:: Modular<int64_t> F(p) ;
		pass &= test_finit(F,m,k,n,timing);
	}
	{
		FFPACK:: ModularBalanced<int64_t> F(p) ;
		pass &= test_finit(F,m,k,n,timing);
	}
#if 1
	{
		FFPACK:: UnparametricField<float> F ;
		pass &= test_finit(F,m,k,n,timing);
	}
	{
		FFPACK:: UnparametricField<double> F ;
		pass &= test_finit(F,m,k,n,timing);
	}
	{
		FFPACK:: UnparametricField<int32_t> F;
		pass &= test_finit(F,m,k,n,timing);
	}
	{
		FFPACK:: UnparametricField<int64_t> F ;
		pass &= test_finit(F,m,k,n,timing);
	}
#endif
	}

	return (pass?0:1) ;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


