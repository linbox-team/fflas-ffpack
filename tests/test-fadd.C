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

#include "fflas-ffpack/utils/test-utils.h"
#include "assert.h"

using namespace FFLAS;
template<class Field>
bool test_fadd(const Field & F, size_t m, size_t k, size_t n, bool timing, uint64_t seed)
{
	typedef typename Field::Element T ;

	T * A = FFLAS::fflas_new<T>(m*n);
	T * B = FFLAS::fflas_new<T>(m*n);
	T * C = FFLAS::fflas_new<T>(m*n);
	T * D = FFLAS::fflas_new<T>(m*n);

	typename Field::RandIter G(F, seed);
	if (timing)	std::cout << ">>>" << std::endl ;

	size_t iter = 3 ;
	FFLAS::Timer tim, tom, tam ;
	tim.clear() ; tom.clear() ;
	if (timing)	F.write(std::cout << "Field ") << std::endl;
	for (size_t b = 0 ; b < iter ; ++b) {
		FFPACK::RandomMatrix(F,m,k,A,n,G);
		FFPACK::RandomMatrix(F,m,k,B,n,G);
		FFPACK::RandomMatrix(F,m,k,C,n,G);
		FFLAS::fassign(F,m,k,C,n,D,n);

		tam.clear();tam.start();
		for (size_t i = 0 ; i < m ; ++i)
			for (size_t j = 0 ; j < k ; ++j)
				F.add(D[i*n+j],A[i*n+j],B[i*n+j]);
		tam.stop();
		tim += tam ;

		tam.clear();tam.start();
		//WriteMatrix(std::cerr<<"A:=",F,m,k,A,n)<<';'<<std::endl;
		//WriteMatrix(std::cerr<<"B:=",F,m,k,B,n)<<';'<<std::endl;
		FFLAS::fadd(F,m,k,A,n,B,n,C,n);
		//WriteMatrix(std::cerr<<"C:=",F,m,k,C,n)<<';'<<std::endl;
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
bool test_faddin(const Field & F, size_t m, size_t k, size_t n, bool timing, uint64_t seed)
{
	typedef typename Field::Element T ;

	typename Field::RandIter G(F, seed);

	T * A = FFLAS::fflas_new<T>(m*n);
	T * C = FFLAS::fflas_new<T>(m*n);
	T * D = FFLAS::fflas_new<T>(m*n);

	if (timing)	std::cout << ">>>" << std::endl ;
	if (timing)	F.write(std::cout << "Field ") << std::endl;
	size_t iter = 3 ;
	FFLAS::Timer tim, tom, tam ;
	tim.clear() ; tom.clear() ;

	for (size_t b = 0 ; b < iter ; ++b) {
		FFPACK::RandomMatrix(F,m,k,A,n,G);
		FFPACK::RandomMatrix(F,m,k,C,n,G);
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
bool test_fsub(const Field & F, size_t m, size_t k, size_t n, bool timing, uint64_t seed)
{
	typedef typename Field::Element T ;

	typename Field::RandIter G(F, seed);

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
		FFPACK::RandomMatrix(F,m,k,A,n,G);
		FFPACK::RandomMatrix(F,m,k,B,n,G);
		FFPACK::RandomMatrix(F,m,k,C,n,G);
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
bool test_fsubin(const Field & F, size_t m, size_t k, size_t n, bool timing, uint64_t seed)
{
	typedef typename Field::Element T ;

	typename Field::RandIter G(F, seed);

	T * A = FFLAS::fflas_new<T>(m*n);
	T * C = FFLAS::fflas_new<T>(m*n);
	T * D = FFLAS::fflas_new<T>(m*n);

	if (timing)  std::cout << ">>>" << std::endl ;
	if (timing)  F.write(std::cout << "Field ") << std::endl;
	size_t iter = 3 ;
	FFLAS::Timer tim, tom, tam ;
	tim.clear() ; tom.clear() ;

	for (size_t b = 0 ; b < iter ; ++b) {
		FFPACK::RandomMatrix(F,m,k,A,n,G);
		FFPACK::RandomMatrix(F,m,k,C,n,G);
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
	uint64_t seed = getSeed();
	static bool timing = false ;

	static Argument as[] = {
		{ 'p', "-p P", "Set the field characteristic.",  TYPE_INT , &p },
		{ 'n', "-n N", "Set the number of cols in C.",   TYPE_INT , &n },
		{ 'm', "-m N", "Set the number of rows in C.",   TYPE_INT , &m },
		{ 'k', "-k N", "Set the number of rows in B.",   TYPE_INT , &k },
		{ 's', "-s N", "Set the seed                 .", TYPE_UINT64, &seed },
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
			pass = pass && test_fadd(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<float> F(p) ;
			pass = pass && test_fadd(F,m,k,n,timing,seed);
		}
		{
			Givaro::Modular<double> F(p) ;
			pass = pass && test_fadd(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<double> F(p) ;
			pass = pass && test_fadd(F,m,k,n,timing,seed);
		}
		{
			Givaro::Modular<int32_t> F( (int32_t)p ) ;
			pass = pass && test_fadd(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<int32_t> F((int32_t)p) ;
			pass = pass && test_fadd(F,m,k,n,timing,seed);
		}
		{
			Givaro::Modular<int64_t> F(p) ;
			pass = pass && test_fadd(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<int64_t> F(p) ;
			pass = pass && test_fadd(F,m,k,n,timing,seed);
		}
#if 1
		{
			Givaro::ZRing<float> F ;
			pass = pass && test_fadd(F,m,k,n,timing,seed);
		}
		{
			Givaro::ZRing<double> F ;
			pass = pass && test_fadd(F,m,k,n,timing,seed);
		}
		{
			Givaro::ZRing<int32_t> F;
			pass = pass && test_fadd(F,m,k,n,timing,seed);
		}
		{
			Givaro::ZRing<int64_t> F ;
			pass = pass && test_fadd(F,m,k,n,timing,seed);
		}
#endif
	}
	{ /*  faddin  */
		{
			Givaro::Modular<float> F(p) ;
			pass = pass && test_faddin(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<float> F(p) ;
			pass = pass && test_faddin(F,m,k,n,timing,seed);
		}
		{
			Givaro::Modular<double> F(p) ;
			pass = pass && test_faddin(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<double> F(p) ;
			pass = pass && test_faddin(F,m,k,n,timing,seed);
		}
		{
			Givaro::Modular<int32_t> F((int32_t)p) ;
			pass = pass && test_faddin(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<int32_t> F((int32_t)p) ;
			pass = pass && test_faddin(F,m,k,n,timing,seed);
		}
		{
			Givaro::Modular<int64_t> F(p) ;
			pass = pass && test_faddin(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<int64_t> F(p) ;
			pass = pass && test_faddin(F,m,k,n,timing,seed);
		}
#if 1
		{
			Givaro::ZRing<float> F ;
			pass = pass && test_faddin(F,m,k,n,timing,seed);
		}
		{
			Givaro::ZRing<double> F ;
			pass = pass && test_faddin(F,m,k,n,timing,seed);
		}
		{
			Givaro::ZRing<int32_t> F;
			pass = pass && test_faddin(F,m,k,n,timing,seed);
		}
		{
			Givaro::ZRing<int64_t> F ;
			pass = pass && test_faddin(F,m,k,n,timing,seed);
		}
#endif
	}
	{ /*  fsub */
		{
			Givaro::Modular<float> F(p) ;
			pass = pass && test_fsub(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<float> F(p) ;
			pass = pass && test_fsub(F,m,k,n,timing,seed);
		}
		{
			Givaro::Modular<double> F(p) ;
			pass = pass && test_fsub(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<double> F(p) ;
			pass = pass && test_fsub(F,m,k,n,timing,seed);
		}
		{
			Givaro::Modular<int32_t> F((int32_t)p) ;
			pass = pass && test_fsub(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<int32_t> F((int32_t)p) ;
			pass = pass && test_fsub(F,m,k,n,timing,seed);
		}
		{
			Givaro::Modular<int64_t> F(p) ;
			pass = pass && test_fsub(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<int64_t> F(p) ;
			pass = pass && test_fsub(F,m,k,n,timing,seed);
		}
#if 1
		{
			Givaro::ZRing<float> F ;
			pass = pass && test_fsub(F,m,k,n,timing,seed);
		}
		{
			Givaro::ZRing<double> F ;
			pass = pass && test_fsub(F,m,k,n,timing,seed);
		}
		{
			Givaro::ZRing<int32_t> F;
			pass = pass && test_fsub(F,m,k,n,timing,seed);
		}
		{
			Givaro::ZRing<int64_t> F ;
			pass = pass && test_fsub(F,m,k,n,timing,seed);
		}
#endif
	}
	{ /*  fsubin */
		{
			Givaro::Modular<float> F(p) ;
			pass = pass && test_fsubin(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<float> F(p) ;
			pass = pass && test_fsubin(F,m,k,n,timing,seed);
		}
		{
			Givaro::Modular<double> F(p) ;
			pass = pass && test_fsubin(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<double> F(p) ;
			pass = pass && test_fsubin(F,m,k,n,timing,seed);
		}
		{
			Givaro::Modular<int32_t> F((int32_t)p) ;
			pass = pass && test_fsubin(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<int32_t> F((int32_t)p) ;
			pass = pass && test_fsubin(F,m,k,n,timing,seed);
		}
		{
			Givaro::Modular<int64_t> F(p) ;
			pass = pass && test_fsubin(F,m,k,n,timing,seed);
		}
		{
			Givaro::ModularBalanced<int64_t> F(p) ;
			pass = pass && test_fsubin(F,m,k,n,timing,seed);
		}
#if 1
		{
			Givaro::ZRing<float> F ;
			pass = pass && test_fsubin(F,m,k,n,timing,seed);
		}
		{
			Givaro::ZRing<double> F ;
			pass = pass && test_fsubin(F,m,k,n,timing,seed);
		}
		{
			Givaro::ZRing<int32_t> F;
			pass = pass && test_fsubin(F,m,k,n,timing,seed);
		}
		{
			Givaro::ZRing<int64_t> F ;
			pass = pass && test_fsubin(F,m,k,n,timing,seed);
		}
#endif
	}

	return (pass?0:1) ;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
