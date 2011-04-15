/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
//--------------------------------------------------------------------------
//                        Test for rank
//
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------

#include <iomanip>
#include <iostream>
#include "fflas-ffpack/field/modular-balanced.h"
#include "timer.h"
#include "Matio.h"
#include "fflas-ffpack/ffpack/ffpack.h"



using namespace std;

typedef Modular<double> Field;

int main(int argc, char** argv){

	int n,m;
	cerr<<setprecision(10);
	if (argc !=  3)	{
		cerr<<"Usage : test-fullranksubmatrix <p> <A> <<i>"
		    <<endl;
		exit(-1);
	}
	Field F (atoi(argv[1]));
	Field::Element * A;
	Field::Element * X;

	A = read_field(F,argv[2],&m ,&n);
	write_field (F, cerr<<"A = "<<endl, A, m, n, n);

	Timer tim,t; t.clear();tim.clear();
	size_t R;

	FFPACK::ColRankProfileSubmatrix (F, m, n, A, n, X, R);

	write_field (F, cerr<<"X = "<<endl, X, R, R, R);

	size_t r2 = FFPACK::Rank(F, R,R, X, R);
	if (r2 != R)
		std::cerr<<"Fail : Rank (X) != Rank (A)"<<std::endl;



	delete[]X;
	delete[]A;
}
