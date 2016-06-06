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

	typename Field::Element alpha,tmp;
	Field::RandIter Rand(F);
	Field::NonZeroRandIter NZRand(Rand);

	for (size_t i=0; i<100; ++i) {

		size_t m = rand() % 10000 + 1;
		size_t n = rand() % 10000 + 1;
		F.init(alpha,rand() % 10000 + 1);
		FFLAS::FFLAS_SIDE side = rand()%2?FFLAS::FflasLeft:FFLAS::FflasRight;
		FFLAS::FFLAS_UPLO uplo = rand()%2?FFLAS::FflasLower:FFLAS::FflasUpper;
		FFLAS::FFLAS_TRANSPOSE trans = rand()%2?FFLAS::FflasNoTrans:FFLAS::FflasTrans;
		FFLAS::FFLAS_DIAG diag = rand()%2?FFLAS::FflasNonUnit:FFLAS::FflasUnit;
		size_t k = (side==FFLAS::FflasLeft?m:n);
		//std::cout << alpha << "  " << side << "  " << uplo << "  " << trans << "  " << diag << "  \n";

		Field::Element_ptr X = FFLAS::fflas_new(F,m,n);
		Field::Element_ptr A = FFLAS::fflas_new(F,k,k);

		for( size_t i = 0; i < m*n; ++i )
			Rand.random( *(X+i) );
		//write_field(F,std::cerr<<"X:=",X,m,n,n,true) <<std::endl;

		for (size_t i=0;i<k;++i){
			for (size_t j=0;j<i;++j)
				A[i*k+j]= (uplo == FFLAS::FflasLower)? Rand.random(tmp) : F.zero;
			A[i*k+i]= (diag == FFLAS::FflasNonUnit)? NZRand.random(tmp) : F.one;
			for (size_t j=i+1;j<k;++j)
				A[i*k+j]= (uplo == FFLAS::FflasUpper)? Rand.random(tmp) : F.zero;
		}
		//write_field(F,std::cerr<<"A:=",A,k,k,k,true) <<std::endl;

		Checker_ftrsm<Field> checker(F, m, n, alpha, X, n);
		FFLAS::ftrsm(F, side, uplo, trans, diag, m, n, alpha, A, k, X, n);
		checker.check(side, uplo, trans, diag, A, k, X);

		FFLAS::fflas_delete(X,A);
	}



	return 0;
}