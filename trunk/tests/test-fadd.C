#include "fflas-ffpack/utils/timer.h"
#include "Matio.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "test-utils.h"
#include "assert.h"
#include "fflas-ffpack/utils/args-parser.h"
#include <typeinfo>

// using namespace FFPACK;
template<class Field>
bool test_fadd(const Field & F, size_t m, size_t k, size_t n)
{
	typedef typename Field::Element T ;

	T * A = new T[m*n];
	T * B = new T[m*n];
	T * C = new T[m*n];
	T * D = new T[m*n];

	std::cout << ">>>" << std::endl ;

	int iter = 3 ;
	Timer tim, tom, tam ;
	tim.clear() ; tom.clear() ;
	F.write(std::cout << "Field ") << std::endl;
	for (size_t b = 0 ; b < iter ; ++b) {
		RandomMatrix(F,A,m,k,n);
		RandomMatrix(F,B,m,k,n);
		RandomMatrix(F,C,m,k,n);
		FFLAS::fcopy(F,m,k,C,n,D,n);

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
					std::cout  <<  i << ',' << j << " : " <<  C[i*n+j] << "!= (ref)" << D[i*n+j] << std::endl;
					return false ;
				}
#endif
	}
	std::cout << "fadd (___): " << tim.usertime()/iter << 's' << std::endl;
	std::cout << "fadd (AVX): " << tom.usertime()/iter << 's'<<  std::endl;

	std::cout << "<<<" << std::endl;
	delete[] A ;
	delete[] B;
	delete[] C ;
	delete[] D ;

	return true;
}

template<class Field>
bool test_faddin(const Field & F, size_t m, size_t k, size_t n)
{
	typedef typename Field::Element T ;

	T * A = new T[m*n];
	T * C = new T[m*n];
	T * D = new T[m*n];

	std::cout << ">>>" << std::endl ;
	F.write(std::cout << "Field ") << std::endl;
	int iter = 3 ;
	Timer tim, tom, tam ;
	tim.clear() ; tom.clear() ;

	for (size_t b = 0 ; b < iter ; ++b) {
		RandomMatrix(F,A,m,k,n);
		RandomMatrix(F,C,m,k,n);
		FFLAS::fcopy(F,m,k,C,n,D,n);

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
					std::cout  <<  i << ',' << j << " : " <<  C[i*n+j] << "!= (ref)" << D[i*n+j] << std::endl;
					return false ;
				}
#endif
	}
	std::cout << "faddin (___): " << tim.usertime()/iter << 's' << std::endl;
	std::cout << "faddin (AVX): " << tom.usertime()/iter << 's'<<  std::endl;


	std::cout << "<<<" << std::endl;
	delete[] A ;
	delete[] C ;
	delete[] D ;

	return true;
}

template<class Field>
bool test_fsub(const Field & F, size_t m, size_t k, size_t n)
{
	typedef typename Field::Element T ;

	T * A = new T[m*n];
	T * B = new T[m*n];
	T * C = new T[m*n];
	T * D = new T[m*n];

	std::cout << ">>>" << std::endl ;

	int iter = 3 ;
	Timer tim, tom, tam ;
	tim.clear() ; tom.clear() ;
	F.write(std::cout << "Field ") << std::endl;
	for (size_t b = 0 ; b < iter ; ++b) {
		RandomMatrix(F,A,m,k,n);
		RandomMatrix(F,B,m,k,n);
		RandomMatrix(F,C,m,k,n);
		FFLAS::fcopy(F,m,k,C,n,D,n);

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
					std::cout  <<  i << ',' << j << " : " <<  C[i*n+j] << "!= (ref)" << D[i*n+j] << std::endl;
					return false ;
				}
#endif
	}
	std::cout << "fsub (___): " << tim.usertime()/iter << 's' << std::endl;
	std::cout << "fsub (AVX): " << tom.usertime()/iter << 's'<<  std::endl;

	std::cout << "<<<" << std::endl;
	delete[] A ;
	delete[] B;
	delete[] C ;
	delete[] D ;

	return true;
}

template<class Field>
bool test_fsubin(const Field & F, size_t m, size_t k, size_t n)
{
	typedef typename Field::Element T ;

	T * A = new T[m*n];
	T * C = new T[m*n];
	T * D = new T[m*n];

	std::cout << ">>>" << std::endl ;
	F.write(std::cout << "Field ") << std::endl;
	int iter = 3 ;
	Timer tim, tom, tam ;
	tim.clear() ; tom.clear() ;

	for (size_t b = 0 ; b < iter ; ++b) {
		RandomMatrix(F,A,m,k,n);
		RandomMatrix(F,C,m,k,n);
		FFLAS::fcopy(F,m,k,C,n,D,n);

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
					std::cout  <<  i << ',' << j << " : " <<  C[i*n+j] << "!= (ref)" << D[i*n+j] << std::endl;
					return false ;
				}
#endif
	}
	std::cout << "fsubin (___): " << tim.usertime()/iter << 's' << std::endl;
	std::cout << "fsubin (AVX): " << tom.usertime()/iter << 's'<<  std::endl;


	std::cout << "<<<" << std::endl;
	delete[] A ;
	delete[] C ;
	delete[] D ;

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
	{ /*  fadd  */
		{
			FFPACK:: Modular<float> F(p) ;
			pass &= test_fadd(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<float> F(p) ;
			pass &= test_fadd(F,m,k,n);
		}
		{
			FFPACK:: Modular<double> F(p) ;
			pass &= test_fadd(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<double> F(p) ;
			pass &= test_fadd(F,m,k,n);
		}
		{
			FFPACK:: Modular<int32_t> F(p) ;
			pass &= test_fadd(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<int32_t> F((int)p) ;
			pass &= test_fadd(F,m,k,n);
		}
		{
			FFPACK:: Modular<int64_t> F(p) ;
			pass &= test_fadd(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<int64_t> F(p) ;
			pass &= test_fadd(F,m,k,n);
		}
#if 1
		{
			FFPACK:: UnparametricField<float> F ;
			pass &= test_fadd(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<double> F ;
			pass &= test_fadd(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<int32_t> F;
			pass &= test_fadd(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<int64_t> F ;
			pass &= test_fadd(F,m,k,n);
		}
#endif
	}
	{ /*  faddin  */
		{
			FFPACK:: Modular<float> F(p) ;
			pass &= test_faddin(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<float> F(p) ;
			pass &= test_faddin(F,m,k,n);
		}
		{
			FFPACK:: Modular<double> F(p) ;
			pass &= test_faddin(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<double> F(p) ;
			pass &= test_faddin(F,m,k,n);
		}
		{
			FFPACK:: Modular<int32_t> F(p) ;
			pass &= test_faddin(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<int32_t> F((int)p) ;
			pass &= test_faddin(F,m,k,n);
		}
		{
			FFPACK:: Modular<int64_t> F(p) ;
			pass &= test_faddin(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<int64_t> F(p) ;
			pass &= test_faddin(F,m,k,n);
		}
#if 1
		{
			FFPACK:: UnparametricField<float> F ;
			pass &= test_faddin(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<double> F ;
			pass &= test_faddin(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<int32_t> F;
			pass &= test_faddin(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<int64_t> F ;
			pass &= test_faddin(F,m,k,n);
		}
#endif
	}
	{ /*  fsub */
		{
			FFPACK:: Modular<float> F(p) ;
			pass &= test_fsub(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<float> F(p) ;
			pass &= test_fsub(F,m,k,n);
		}
		{
			FFPACK:: Modular<double> F(p) ;
			pass &= test_fsub(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<double> F(p) ;
			pass &= test_fsub(F,m,k,n);
		}
		{
			FFPACK:: Modular<int32_t> F(p) ;
			pass &= test_fsub(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<int32_t> F((int)p) ;
			pass &= test_fsub(F,m,k,n);
		}
		{
			FFPACK:: Modular<int64_t> F(p) ;
			pass &= test_fsub(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<int64_t> F(p) ;
			pass &= test_fsub(F,m,k,n);
		}
#if 1
		{
			FFPACK:: UnparametricField<float> F ;
			pass &= test_fsub(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<double> F ;
			pass &= test_fsub(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<int32_t> F;
			pass &= test_fsub(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<int64_t> F ;
			pass &= test_fsub(F,m,k,n);
		}
#endif
	}
	{ /*  fsubin */
		{
			FFPACK:: Modular<float> F(p) ;
			pass &= test_fsubin(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<float> F(p) ;
			pass &= test_fsubin(F,m,k,n);
		}
		{
			FFPACK:: Modular<double> F(p) ;
			pass &= test_fsubin(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<double> F(p) ;
			pass &= test_fsubin(F,m,k,n);
		}
		{
			FFPACK:: Modular<int32_t> F(p) ;
			pass &= test_fsubin(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<int32_t> F((int)p) ;
			pass &= test_fsubin(F,m,k,n);
		}
		{
			FFPACK:: Modular<int64_t> F(p) ;
			pass &= test_fsubin(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<int64_t> F(p) ;
			pass &= test_fsubin(F,m,k,n);
		}
#if 1
		{
			FFPACK:: UnparametricField<float> F ;
			pass &= test_fsubin(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<double> F ;
			pass &= test_fsubin(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<int32_t> F;
			pass &= test_fsubin(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<int64_t> F ;
			pass &= test_fsubin(F,m,k,n);
		}
#endif
	}

	return (pass?0:1) ;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


