/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
//--------------------------------------------------------------------------
//          Test for the krylov-elimination
//--------------------------------------------------------------------------
// usage: test-krylov-elim p A, to compute the rank profile of the (n+m)xn matrix B
// formed by the n identity vectors and the mxn matrix A over Z/pZ
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
#define DEBUG 0

#include <iostream>
#include <list>
#include <vector>
#include "Matio.h"
#include "timer.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/ffpack/ffpack.h"

#include "fflas-ffpack/utils/args-parser.h"
using namespace std;


typedef NAMESPACE Modular<double> Field;

template<class T>
std::ostream& printvect(std::ostream& o, vector<T>& vect){
	for(size_t i=0; i < vect.size(); ++i)
		o << vect[i] << " " ;
	return o << std::endl;
}

int main(int argc, char** argv)
{


	static Argument as[] = {
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,as);

	// int m,n;

	Field F(65521);
	int N = 17;
	double * A = new double[N*N];
	double * tmp = new double[N*N];
	size_t * deg = new size_t[N];

	for (size_t i=0; i<(size_t)N*N; ++i)
		A[i] = 0;
	for (int i=0; i<3; ++i)
		A[i+i*N] = 1;

	for (int i=3; i<6; ++i)
		A[i+1+i*N] = 1;
	for (int i=6; i<9; ++i)
		A[i+2+i*N] = 1;

	A[12+9*N] = 1;
	A[13+10*N] = 1;
	A[14+12*N] = 1;
	A[15+15*N] = 1;
	A[16+16*N] = 1;
	deg[0]  = 4; deg[1] = 4; deg[2] = 4;deg[3] = 2; deg[4] = 1; deg[5] =2;
	for (size_t i=0; i<size_t(N); ++i)
		A[11+i*N] = A[7+i*N] = A[3+i*N] = i % 10;

	double * B = new double[N*N] ;
	FFLAS::fcopy(F,N*N,B,1,A,1);

	// write_field(F, cerr, A, N, N, N);

	FFPACK::CompressRowsQK (F, N, A+9*N, N, tmp, N, deg+3, 4, 3 );

	// write_field(F, cerr, A, N, N, N);

	FFPACK::DeCompressRowsQK (F, N, N-9, A+9*N, N, tmp, N, deg+3, 4, 3 );

	// write_field(F, cerr, A, N, N, N);

	for (size_t i = 0 ; i < N * N ; ++i)
		if (A[i] != B[i])
			return 1 ;
	return 0 ;

}
