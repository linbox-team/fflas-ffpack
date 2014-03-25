#include "utils/timer.h"
#include "Matio.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "test-utils.h"
#include "assert.h"
#include "fflas-ffpack/utils/args-parser.h"
#include <typeinfo>

// using namespace FFPACK;
template<class T>
bool test_finit(size_t m, size_t k, size_t n, size_t p)
{

	T * A = new T[m*n];
	T * B = new T[m*n];
	FFPACK::Modular<T> F(p);
	FFPACK::ModularBalanced<T> G((int32_t)p);

	// FFPACK::Modular<T> E(21*p);
	FFPACK::ModularBalanced<T> E(21*(int32_t)p);
	// for (size_t i = 0 ; i < m ; ++i)
	// for (size_t j = 0 ; j < k ; ++j)
	// A[i*k+j] = F.mOne ;
	// for (size_t i = 0 ; i < k ; ++i)
	// for (size_t j = 0 ; j < n ; ++j)
	// B[i*n+j] = F.mOne ;

	std::cout << ">>>" << std::endl ;

	for (size_t b = 0 ; b < 3 ; ++b) {
		RandomMatrix(E,A,m,k,n);
		// RandomMatrix(E,B,m,k,n);
		FFLAS::fcopy(E,m,k,B,n,A,n);

		Timer tim;
		Timer tom;
		F.write(std::cout << "Modular ") << std::endl;
		tim.clear();tim.start();
		for (size_t i = 0 ; i < m ; ++i)
			for (size_t j = 0 ; j < k ; ++j)
				F.init(A[i*n+j],A[i*n+j]);
		tim.stop();
		std::cout << "finit (___): " << tim.usertime() << 's' << std::endl;

		tim.clear();tim.start();
		FFLAS::finit(F,m,k,B,n);
		tim.stop();
		std::cout << "finit (AVX): " << tim.usertime() << 's'<<  std::endl << std::endl;

#if 1
		for (size_t i =0 ; i < m ; ++i)
			for (size_t j =0 ; j < k ; ++j)
				if (! F.areEqual(B[i*n+j],A[i*n+j])) {
					std::cout  <<  i << ',' << j << " : " <<  B[i*n+j] << "!= (ref)" << A[i*n+j] << std::endl;
					return false ;
				}
#endif
	}

	std::cout << "===" << std::endl ;

	for (size_t b = 0 ; b < 3 ; ++b) {
		RandomMatrix(E,A,m,k,n);
		FFLAS::fcopy(E,m,k,B,n,A,n);

		Timer tim;
		Timer tom;
		G.write(std::cout << "Modular ") << std::endl;
		tim.clear();tim.start();
		for (size_t i = 0 ; i < m ; ++i)
			for (size_t j = 0 ; j < k ; ++j)
				G.init(A[i*n+j],A[i*n+j]);

		tim.stop();
		std::cout << "finit (___): " << tim.usertime() << 's' << std::endl;
		tim.clear();tim.start();
		FFLAS::finit(G,m,k,B,n);
		tim.stop();
		std::cout << "finit (AVX): " << tim.usertime() << 's'<<  std::endl << std::endl;

#if 1
		for (size_t i =0 ; i < m ; ++i)
			for (size_t j =0 ; j < k ; ++j)
				if (! G.areEqual(B[i*n+j],A[i*n+j])) {
					std::cout <<  i << ',' << j << " : " << B[i*n+j] << "!= (ref)" << A[i*n+j] << std::endl;
					return false ;
				}
#endif

	}

	std::cout << "===" << std::endl ;

	for (size_t b = 0 ; b < 3 ; ++b) {
		RandomMatrix(E,A,m,n,n);
		FFLAS::fcopy(E,m,n,B,n,A,n);
		size_t incX = 2 ;

		Timer tim;
		Timer tom;
		G.write(std::cout << "Modular ") << std::endl;
		tim.clear();tim.start();
		for (size_t i = 1 ; i < m*n ; i += incX) {
				G.init(A[i],A[i]);
		}

		tim.stop();
		size_t cnt = (size_t)floor((double)(m*n)/(double)incX) ;
		std::cout << "finit (___): " << tim.usertime() << 's' << std::endl;
		tim.clear();tim.start();
		FFLAS::finit(G,cnt,B+1,incX);
		tim.stop();
		std::cout << "finit (AVX): " << tim.usertime() << 's'<<  std::endl << std::endl;

#if 1
		for (size_t i =1 ; i < m*n ; i+=incX)
				if (! G.areEqual(B[i],A[i])) {
					std::cout <<  i << " : " << B[i] << "!= (ref)" << A[i] << std::endl;
					return false ;
				}
#endif

	}

	std::cout << "<<<" << std::endl;
	delete[] A ;
	delete[] B;

	return true;
}

int main(int ac, char **av) {
	static size_t m = 300 ;
	static size_t n = 301 ;
	static size_t k = 300 ;
	static size_t p = 7;
	int seed = (int) time(NULL);

	static Argument as[] = {
		{ 'p', "-p P", "Set the field characteristic.",  TYPE_INT , &p },
		{ 'n', "-n N", "Set the number of cols in C.",   TYPE_INT , &n },
		{ 'm', "-m N", "Set the number of rows in C.",   TYPE_INT , &m },
		{ 'k', "-k N", "Set the number of rows in B.",   TYPE_INT , &k },
		{ 's', "-s N", "Set the seed                 .", TYPE_INT , &seed },
		END_OF_ARGUMENTS
	};
	FFLAS::parseArguments(ac,av,as);

	if (n < k) {
		std::cout << "Usage : m k n ; matrix of size m x k, lda is n" << std::endl;
		return -1 ;
	}

	srand(seed);
	srand48(seed);

	std::cout << seed << std::endl;

	bool pass  = true ;
	// pass &= test_finit<int32_t>(m,k,n,p);
	// pass &= test_finit<int64_t>(m,k,n,p);
	pass &= test_finit<float>  (m,k,n,p);
	pass &= test_finit<double> (m,k,n,p);

	return (pass?0:1) ;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


