#include "fflas-ffpack/utils/timer.h"
#include "Matio.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "test-utils.h"
#include "assert.h"
#include "fflas-ffpack/utils/args-parser.h"
#include <typeinfo>

// using namespace FFPACK;
using FFPACK::ModularBalanced ;

template<class Field>
bool test_fscal(const Field & F, const typename Field::Element & alpha, size_t m, size_t k, size_t n)
{
	typedef typename Field::Element T ;

	T * A = new T[m*n];
	T * C = new T[m*n];
	T * D = new T[m*n];

	std::cout << ">>>" << std::endl ;

	int iter = 3 ;
	Timer tim, tom, tam ;
	tim.clear() ; tom.clear() ;
	F.write(std::cout << "Field ") << std::endl;
	for (size_t b = 0 ; b < iter ; ++b) {
		RandomMatrix(F,A,m,k,n);
		RandomMatrix(F,C,m,k,n);
		FFLAS::fcopy(F,m,k,C,n,D,n);

		tam.clear();tam.start();
		for (size_t i = 0 ; i < m ; ++i)
			for (size_t j = 0 ; j < k ; ++j)
				F.mul(D[i*n+j],A[i*n+j],alpha);
		tam.stop();
		tim += tam ;

		tam.clear();tam.start();
		FFLAS::fscal(F,m,k,alpha,A,n,C,n);
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
	std::cout << "fscal(___): " << tim.usertime()/iter << 's' << std::endl;
	std::cout << "fscal (AVX): " << tom.usertime()/iter << 's'<<  std::endl;

	std::cout << "<<<" << std::endl;
	delete[] A ;
	delete[] C ;
	delete[] D ;

	return true;
}

template<class Field>
bool test_fscal(const Field & F,  size_t m, size_t k, size_t n)
{
	ModularBalanced<typename Field::Element> G(1234);
	bool pass = true ;
	typename Field::Element  alpha;
	F.init(alpha,F.one);
	pass &= test_fscal(F,alpha,m,k,n);
	F.init(alpha,F.mOne);
	pass &= test_fscal(F,alpha,m,k,n);
	F.init(alpha,F.zero);
	pass &= test_fscal(F,alpha,m,k,n);
	typename ModularBalanced<typename Field::Element>::RandIter RValue( G );
	F.init(alpha,RValue.random(alpha));
	pass &= test_fscal(F,alpha,m,k,n);
	F.init(alpha,RValue.random(alpha));
	pass &= test_fscal(F,alpha,m,k,n);

	return pass ;
}

template<class Field>
bool test_fscalin(const Field & F, const typename Field::Element & alpha, size_t m, size_t k, size_t n)
{
	typedef typename Field::Element T ;

	T * C = new T[m*n];
	T * D = new T[m*n];

	std::cout << ">>>" << std::endl ;

	int iter = 3 ;
	Timer tim, tom, tam ;
	tim.clear() ; tom.clear() ;
	F.write(std::cout << "Field ") << std::endl;
	for (size_t b = 0 ; b < iter ; ++b) {
		RandomMatrix(F,C,m,k,n);
		FFLAS::fcopy(F,m,k,C,n,D,n);

		tam.clear();tam.start();
		for (size_t i = 0 ; i < m ; ++i)
			for (size_t j = 0 ; j < k ; ++j)
				F.mulin(D[i*n+j],alpha);
		tam.stop();
		tim += tam ;

		tam.clear();tam.start();
		FFLAS::fscalin(F,m,k,alpha,C,n);
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
	std::cout << "fscalin(___): " << tim.usertime()/iter << 's' << std::endl;
	std::cout << "fscalin (AVX): " << tom.usertime()/iter << 's'<<  std::endl;

	std::cout << "<<<" << std::endl;
	delete[] C ;
	delete[] D ;

	return true;
}


template<class Field>
bool test_fscalin(const Field & F,  size_t m, size_t k, size_t n)
{
	ModularBalanced<typename Field::Element> G(1234);
	bool pass = true ;
	typename Field::Element  alpha;
	F.init(alpha,F.one);
	pass &= test_fscalin(F,alpha,m,k,n);
	F.init(alpha,F.mOne);
	pass &= test_fscalin(F,alpha,m,k,n);
	F.init(alpha,F.zero);
	pass &= test_fscalin(F,alpha,m,k,n);
	typename ModularBalanced<typename Field::Element>::RandIter RValue( G );
	F.init(alpha,RValue.random(alpha));
	pass &= test_fscalin(F,alpha,m,k,n);
	F.init(alpha,RValue.random(alpha));
	pass &= test_fscalin(F,alpha,m,k,n);

	return pass ;
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
	{ /*  fscal  */
		{
			FFPACK:: Modular<float> F(p) ;
			pass &= test_fscal(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<float> F(p) ;
			pass &= test_fscal(F,m,k,n);
		}
		{
			FFPACK:: Modular<double> F(p) ;
			pass &= test_fscal(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<double> F(p) ;
			pass &= test_fscal(F,m,k,n);
		}
		{
			FFPACK:: Modular<int32_t> F(p) ;
			pass &= test_fscal(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<int32_t> F((int)p) ;
			pass &= test_fscal(F,m,k,n);
		}
		{
			FFPACK:: Modular<int64_t> F(p) ;
			pass &= test_fscal(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<int64_t> F(p) ;
			pass &= test_fscal(F,m,k,n);
		}
#if 1
		{
			FFPACK:: UnparametricField<float> F ;
			pass &= test_fscal(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<double> F ;
			pass &= test_fscal(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<int32_t> F;
			pass &= test_fscal(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<int64_t> F ;
			pass &= test_fscal(F,m,k,n);
		}
#endif
	}
	{ /*  fscalin  */
		{
			FFPACK:: Modular<float> F(p) ;
			pass &= test_fscalin(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<float> F(p) ;
			pass &= test_fscalin(F,m,k,n);
		}
		{
			FFPACK:: Modular<double> F(p) ;
			pass &= test_fscalin(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<double> F(p) ;
			pass &= test_fscalin(F,m,k,n);
		}
		{
			FFPACK:: Modular<int32_t> F(p) ;
			pass &= test_fscalin(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<int32_t> F((int)p) ;
			pass &= test_fscalin(F,m,k,n);
		}
		{
			FFPACK:: Modular<int64_t> F(p) ;
			pass &= test_fscalin(F,m,k,n);
		}
		{
			FFPACK:: ModularBalanced<int64_t> F(p) ;
			pass &= test_fscalin(F,m,k,n);
		}
#if 1
		{
			FFPACK:: UnparametricField<float> F ;
			pass &= test_fscalin(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<double> F ;
			pass &= test_fscalin(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<int32_t> F;
			pass &= test_fscalin(F,m,k,n);
		}
		{
			FFPACK:: UnparametricField<int64_t> F ;
			pass &= test_fscalin(F,m,k,n);
		}
#endif
	}

	return (pass?0:1) ;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


