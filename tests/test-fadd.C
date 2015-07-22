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
#include <givaro/modular-balanced.h>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/args-parser.h"

#include "Matio.h"
#include "test-utils.h"
#include "assert.h"

template<class Field>
bool test_fadd(const Field & F, size_t m, size_t k, size_t n, bool timing)
{
	typedef typename Field::Element T ;

	T * A = FFLAS::fflas_new<T>(m*n);
	T * B = FFLAS::fflas_new<T>(m*n);
	T * C = FFLAS::fflas_new<T>(m*n);
	T * D = FFLAS::fflas_new<T>(m*n);

	if (timing)	std::cout << ">>>" << std::endl ;

	size_t iter = 3 ;
 FFLAS::Timer tim, tom, tam ;
	tim.clear() ; tom.clear() ;
		if (timing)	F.write(std::cout << "Field ") << std::endl;
	for (size_t b = 0 ; b < iter ; ++b) {
		FFPACK::RandomMatrix(F,A,m,k,n);
		FFPACK::RandomMatrix(F,B,m,k,n);
		FFPACK::RandomMatrix(F,C,m,k,n);
		FFLAS::fassign(F,m,k,C,n,D,n);

		tam.clear();tam.start();
		for (size_t i = 0 ; i < m ; ++i)
			for (size_t j = 0 ; j < k ; ++j)
				F.add(D[i*n+j],A[i*n+j],B[i*n+j]);
		tam.stop();
		tim += tam ;

		tam.clear();tam.start();
		FFLAS::fadd(F,m,k,A,n,B,n,C,n);
		tam.stop();
		tom += tam ;

#if 1
		for (size_t i =0 ; i < m ; ++i)
			for (size_t j =0 ; j < k ; ++j)
				if (! F.areEqual(C[i*n+j],D[i*n+j])) {
					if (timing)		std::cout  <<  i << ',' << j << " : " <<  C[i*n+j] << "!= (ref)" << D[i*n+j] << std::endl;
					return false ;
				}
#endif
	}
		if (timing)	std::cout << "fadd (___): " << tim.usertime()/(double)iter << 's' << std::endl;
		if (timing)	std::cout << "fadd (AVX): " << tom.usertime()/(double)iter << 's'<<  std::endl;

		if (timing)	std::cout << "<<<" << std::endl;
	FFLAS::fflas_delete( A );
	FFLAS::fflas_delete( B);
	FFLAS::fflas_delete( C );
	FFLAS::fflas_delete( D );

	return true;
}

template<class Field>
bool test_faddin(const Field & F, size_t m, size_t k, size_t n, bool timing)
{
	typedef typename Field::Element T ;

	T * A = FFLAS::fflas_new<T>(m*n);
	T * C = FFLAS::fflas_new<T>(m*n);
	T * D = FFLAS::fflas_new<T>(m*n);

		if (timing)	std::cout << ">>>" << std::endl ;
		if (timing)	F.write(std::cout << "Field ") << std::endl;
	size_t iter = 3 ;
 FFLAS::Timer tim, tom, tam ;
	tim.clear() ; tom.clear() ;

	for (size_t b = 0 ; b < iter ; ++b) {
		FFPACK::RandomMatrix(F,A,m,k,n);
		FFPACK::RandomMatrix(F,C,m,k,n);
		FFLAS::fassign(F,m,k,C,n,D,n);

		tam.clear();tam.start();
		for (size_t i = 0 ; i < m ; ++i)
			for (size_t j = 0 ; j < k ; ++j)
				F.addin(D[i*n+j],A[i*n+j]);
		tam.stop();
		tim += tam ;

		tam.clear();tam.start();
		FFLAS::faddin(F,m,k,A,n,C,n);
		tam.stop();
		tom += tam ;

#if 1
		for (size_t i =0 ; i < m ; ++i)
			for (size_t j =0 ; j < k ; ++j)
				if (! F.areEqual(C[i*n+j],D[i*n+j])) {
						if (timing)	std::cout  <<  i << ',' << j << " : " <<  C[i*n+j] << "!= (ref)" << D[i*n+j] << std::endl;
					return false ;
				}
#endif
	}
		if (timing)	std::cout << "faddin (___): " << tim.usertime()/(double)iter << 's' << std::endl;
		if (timing)	std::cout << "faddin (AVX): " << tom.usertime()/(double)iter << 's'<<  std::endl;


		if (timing)	std::cout << "<<<" << std::endl;
	FFLAS::fflas_delete( A );
	FFLAS::fflas_delete( C );
	FFLAS::fflas_delete( D );

	return true;
}

template<class Field>
bool test_fsub(const Field & F, size_t m, size_t k, size_t n, bool timing)
{
	typedef typename Field::Element T ;

	T * A = FFLAS::fflas_new<T>(m*n);
	T * B = FFLAS::fflas_new<T>(m*n);
	T * C = FFLAS::fflas_new<T>(m*n);
	T * D = FFLAS::fflas_new<T>(m*n);

		if (timing)	std::cout << ">>>" << std::endl ;

	size_t iter = 3 ;
 FFLAS::Timer tim, tom, tam ;
	tim.clear() ; tom.clear() ;
		if (timing)	F.write(std::cout << "Field ") << std::endl;
	for (size_t b = 0 ; b < iter ; ++b) {
		FFPACK::RandomMatrix(F,A,m,k,n);
		FFPACK::RandomMatrix(F,B,m,k,n);
		FFPACK::RandomMatrix(F,C,m,k,n);
		FFLAS::fassign(F,m,k,C,n,D,n);

		tam.clear();tam.start();
		for (size_t i = 0 ; i < m ; ++i)
			for (size_t j = 0 ; j < k ; ++j)
				F.sub(D[i*n+j],A[i*n+j],B[i*n+j]);
		tam.stop();
		tim += tam ;

		tam.clear();tam.start();
		FFLAS::fsub(F,m,k,A,n,B,n,C,n);
		tam.stop();
		tom += tam ;

#if 1
		for (size_t i =0 ; i < m ; ++i)
			for (size_t j =0 ; j < k ; ++j)
				if (! F.areEqual(C[i*n+j],D[i*n+j])) {
						if (timing)	std::cout  <<  i << ',' << j << " : " <<  C[i*n+j] << "!= (ref)" << D[i*n+j] << std::endl;
					return false ;
				}
#endif
	}
		if (timing)	std::cout << "fsub (___): " << tim.usertime()/(double)iter << 's' << std::endl;
		if (timing)	std::cout << "fsub (AVX): " << tom.usertime()/(double)iter << 's'<<  std::endl;

		if (timing)	std::cout << "<<<" << std::endl;
	FFLAS::fflas_delete( A );
	FFLAS::fflas_delete( B);
	FFLAS::fflas_delete( C );
	FFLAS::fflas_delete( D );

	return true;
}

template<class Field>
bool test_fsubin(const Field & F, size_t m, size_t k, size_t n, bool timing)
{
	typedef typename Field::Element T ;

	T * A = FFLAS::fflas_new<T>(m*n);
	T * C = FFLAS::fflas_new<T>(m*n);
	T * D = FFLAS::fflas_new<T>(m*n);

	if (timing)  std::cout << ">>>" << std::endl ;
	if (timing)  F.write(std::cout << "Field ") << std::endl;
	size_t iter = 3 ;
 FFLAS::Timer tim, tom, tam ;
	tim.clear() ; tom.clear() ;

	for (size_t b = 0 ; b < iter ; ++b) {
		FFPACK::RandomMatrix(F,A,m,k,n);
		FFPACK::RandomMatrix(F,C,m,k,n);
		FFLAS::fassign(F,m,k,C,n,D,n);

		tam.clear();tam.start();
		for (size_t i = 0 ; i < m ; ++i)
			for (size_t j = 0 ; j < k ; ++j)
				F.subin(D[i*n+j],A[i*n+j]);
		tam.stop();
		tim += tam ;

		tam.clear();tam.start();
		FFLAS::fsubin(F,m,k,A,n,C,n);
		tam.stop();
		tom += tam ;

#if 1
		for (size_t i =0 ; i < m ; ++i)
			for (size_t j =0 ; j < k ; ++j)
				if (! F.areEqual(C[i*n+j],D[i*n+j])) {
					if (timing)  std::cout  <<  i << ',' << j << " : " <<  C[i*n+j] << "!= (ref)" << D[i*n+j] << std::endl;
					return false ;
				}
#endif
	}
if (timing)	std::cout << "fsubin (___): " << tim.usertime()/(double)iter << 's' << std::endl;
	if (timing) std::cout << "fsubin (AVX): " << tom.usertime()/(double)iter << 's'<<  std::endl;


	if (timing) std::cout << "<<<" << std::endl;
	FFLAS::fflas_delete( A );
	FFLAS::fflas_delete( C );
	FFLAS::fflas_delete( D );

	return true;
}


int main(int ac, char **av) {
	static size_t m = 300 ;
	static size_t n = 301 ;
	static size_t k = 300 ;
	static uint64_t  p = 7;
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
	{ /*  fadd  */
		{
			Givaro::Modular<float> F(p) ;
			pass &= test_fadd(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<float> F(p) ;
			pass &= test_fadd(F,m,k,n,timing);
		}
		{
			Givaro::Modular<double> F(p) ;
			pass &= test_fadd(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<double> F(p) ;
			pass &= test_fadd(F,m,k,n,timing);
		}
		{
			Givaro::Modular<int32_t> F( (int32_t)p ) ;
			pass &= test_fadd(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<int32_t> F((int32_t)p) ;
			pass &= test_fadd(F,m,k,n,timing);
		}
		{
			Givaro::Modular<int64_t> F(p) ;
			pass &= test_fadd(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<int64_t> F(p) ;
			pass &= test_fadd(F,m,k,n,timing);
		}
#if 1
		{
			Givaro::ZRing<float> F ;
			pass &= test_fadd(F,m,k,n,timing);
		}
		{
			Givaro::ZRing<double> F ;
			pass &= test_fadd(F,m,k,n,timing);
		}
		{
			Givaro::ZRing<int32_t> F;
			pass &= test_fadd(F,m,k,n,timing);
		}
		{
			Givaro::ZRing<int64_t> F ;
			pass &= test_fadd(F,m,k,n,timing);
		}
#endif
	}
	{ /*  faddin  */
		{
			Givaro::Modular<float> F(p) ;
			pass &= test_faddin(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<float> F(p) ;
			pass &= test_faddin(F,m,k,n,timing);
		}
		{
			Givaro::Modular<double> F(p) ;
			pass &= test_faddin(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<double> F(p) ;
			pass &= test_faddin(F,m,k,n,timing);
		}
		{
			Givaro::Modular<int32_t> F((int32_t)p) ;
			pass &= test_faddin(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<int32_t> F((int32_t)p) ;
			pass &= test_faddin(F,m,k,n,timing);
		}
		{
			Givaro::Modular<int64_t> F(p) ;
			pass &= test_faddin(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<int64_t> F(p) ;
			pass &= test_faddin(F,m,k,n,timing);
		}
#if 1
		{
			Givaro::ZRing<float> F ;
			pass &= test_faddin(F,m,k,n,timing);
		}
		{
			Givaro::ZRing<double> F ;
			pass &= test_faddin(F,m,k,n,timing);
		}
		{
			Givaro::ZRing<int32_t> F;
			pass &= test_faddin(F,m,k,n,timing);
		}
		{
			Givaro::ZRing<int64_t> F ;
			pass &= test_faddin(F,m,k,n,timing);
		}
#endif
	}
	{ /*  fsub */
		{
			Givaro::Modular<float> F(p) ;
			pass &= test_fsub(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<float> F(p) ;
			pass &= test_fsub(F,m,k,n,timing);
		}
		{
			Givaro::Modular<double> F(p) ;
			pass &= test_fsub(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<double> F(p) ;
			pass &= test_fsub(F,m,k,n,timing);
		}
		{
			Givaro::Modular<int32_t> F((int32_t)p) ;
			pass &= test_fsub(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<int32_t> F((int32_t)p) ;
			pass &= test_fsub(F,m,k,n,timing);
		}
		{
			Givaro::Modular<int64_t> F(p) ;
			pass &= test_fsub(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<int64_t> F(p) ;
			pass &= test_fsub(F,m,k,n,timing);
		}
#if 1
		{
			Givaro::ZRing<float> F ;
			pass &= test_fsub(F,m,k,n,timing);
		}
		{
			Givaro::ZRing<double> F ;
			pass &= test_fsub(F,m,k,n,timing);
		}
		{
			Givaro::ZRing<int32_t> F;
			pass &= test_fsub(F,m,k,n,timing);
		}
		{
			Givaro::ZRing<int64_t> F ;
			pass &= test_fsub(F,m,k,n,timing);
		}
#endif
	}
	{ /*  fsubin */
		{
			Givaro::Modular<float> F(p) ;
			pass &= test_fsubin(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<float> F(p) ;
			pass &= test_fsubin(F,m,k,n,timing);
		}
		{
			Givaro::Modular<double> F(p) ;
			pass &= test_fsubin(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<double> F(p) ;
			pass &= test_fsubin(F,m,k,n,timing);
		}
		{
			Givaro::Modular<int32_t> F((int32_t)p) ;
			pass &= test_fsubin(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<int32_t> F((int32_t)p) ;
			pass &= test_fsubin(F,m,k,n,timing);
		}
		{
			Givaro::Modular<int64_t> F(p) ;
			pass &= test_fsubin(F,m,k,n,timing);
		}
		{
			Givaro::ModularBalanced<int64_t> F(p) ;
			pass &= test_fsubin(F,m,k,n,timing);
		}
#if 1
		{
			Givaro::ZRing<float> F ;
			pass &= test_fsubin(F,m,k,n,timing);
		}
		{
			Givaro::ZRing<double> F ;
			pass &= test_fsubin(F,m,k,n,timing);
		}
		{
			Givaro::ZRing<int32_t> F;
			pass &= test_fsubin(F,m,k,n,timing);
		}
		{
			Givaro::ZRing<int64_t> F ;
			pass &= test_fsubin(F,m,k,n,timing);
		}
#endif
	}

	return (pass?0:1) ;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


