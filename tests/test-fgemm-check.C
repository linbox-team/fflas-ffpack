#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"

int main(int argc, char** argv) {
	srand (time(NULL));
	typedef Givaro::Modular<double> Field;
	Givaro::Integer q = 131071;
	Field F(q);

	Field::RandIter Rand(F);
	FFLAS::FFLAS_TRANSPOSE ta,tb;

	for (size_t i=0; i<100; ++i) {

		size_t m = rand() % 10000 + 1;
		size_t n = rand() % 10000 + 1;
		size_t k = rand() % 10000 + 1;

		typename Field::Element alpha,beta;
		F.init(alpha, rand()%10000+1);
		F.init(beta,  rand()%10000+1);
		
		ta = rand()%2 ? FFLAS::FflasNoTrans : FFLAS::FflasTrans,
		tb = rand()%2 ? FFLAS::FflasNoTrans : FFLAS::FflasTrans;

		size_t lda = ta == FFLAS::FflasNoTrans ? k : m,
			   ldb = tb == FFLAS::FflasNoTrans ? n : k,
			   ldc = n;

		Field::Element_ptr A = FFLAS::fflas_new(F,m,k);
		Field::Element_ptr B = FFLAS::fflas_new(F,k,n);
		Field::Element_ptr C = FFLAS::fflas_new(F,m,n);

		for( size_t i = 0; i < m*k; ++i )
			Rand.random( *(A+i) );
		for( size_t i = 0; i < k*n; ++i )
			Rand.random( *(B+i) );
		for( size_t i = 0; i < m*n; ++i )
			Rand.random( *(C+i) );

		Checker_fgemm<Field> checker(F,m,n,k,beta,C,ldc);
		FFLAS::fgemm(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
		checker.check(ta,tb,alpha,A,lda,B,ldb,C);

		FFLAS::fflas_delete(A,B,C);
	}



	return 0;
}