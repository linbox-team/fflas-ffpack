#define ENABLE_CHECKER_charpoly 1
#define ENABLE_CHECKER_PLUQ 1

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"


template <class Field, class Polynomial>
void printPolynomial (const Field &F, Polynomial &v)
{
	for (int i = v.size() - 1; i >= 0; i--) {
		F.write (std::cout, v[i]);
		if (i > 0)
			std::cout << " x^" << i << " + ";
	}
	std::cout << std::endl;
}

int main(int argc, char** argv) {
	srand (time(NULL));
	typedef Givaro::Modular<double> Field;
	Givaro::Integer q = 131071;
	Field F(q);
	typedef std::vector<Field::Element> Polynomial;

	Field::RandIter Rand(F);
	Field::Element_ptr A = FFLAS::fflas_new(F,1000,1000);

	size_t p = 0;
	for (size_t i=0; i<1000; ++i) {

		size_t n = rand() % 1000 + 1;

		Polynomial g(n);

		for( size_t i = 0; i < n*n; ++i )
			Rand.random( *(A+i) );

		//write_field(F,std::cerr<<"A=",A,n,n,n,true) <<std::endl;
		Checker_charpoly<Field,Polynomial> checker(F,n,A);
		FFPACK::CharPoly(F,g,n,A,n,FFPACK::FfpackLUK);
		//printPolynomial(F,g);
		p += checker.check(g);
	}

	std::cout << p << "/1000 successful tests\n";

	return 0;
}