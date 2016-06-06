#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"

int main(int argc, char** argv) {
	srand (time(NULL));
	typedef Givaro::Modular<double> Field;
	Givaro::Integer q = rand()%10000;//131071;
	Field F(q);

	Field::RandIter Rand(F);

	for (size_t i=0; i<10; ++i) {

		size_t m = rand() % 10000 + 1;
		int nullity;

		Field::Element_ptr A = FFLAS::fflas_new(F,m,m);

		for( size_t i = 0; i < m*m; ++i )
			Rand.random( *(A+i) );

		Checker_invert<Field> checker(F,m,A,m);
		FFPACK::Invert(F,m,A,m,nullity);
		checker.check(A,nullity);

		FFLAS::fflas_delete(A);
	}



	return 0;
}