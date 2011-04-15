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
#include "Matio.h"
#include "timer.h"
using namespace std;
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/ffpack/ffpack.h"


typedef Modular<double> Field;

template<class T>
std::ostream& printvect(std::ostream& o, T* vect, size_t dim)
{
	for(size_t i=0; i<dim; ++i)
		o << vect[i] << " " ;
	return o << std::endl;
}

int main(int argc, char** argv){

	size_t m,n;


	if (argc!=3){
		cerr<<"usage : test-lqup <p> <A>"<<endl
		    <<"         to compute the rank profile of the (n+m)xn matrix B formed by the n identity vectors and the mxn matrix A over Z/pZ"
		    <<endl;
		exit(-1);
	}
	Field F(atoi(argv[1]));
	Field::Element one;
	F.init(one, 1UL);
	Field::Element * A = read_field<Field> (F,argv[2],(int*)&m,(int*)&n);

	Field::Element * B = new Field::Element[(m+n)*n];
	for (size_t i=0; i<(n+m)*n;++i) *(B+i)=0;

	size_t deg = (n-1)/m+1;
	size_t curr_row = 0;
	size_t it_idx = 0;
	size_t bk_idx = 0;
	for (size_t i=0; i<m; ++i){
		for (size_t j=0; j<deg; ++j){
			if (curr_row < n+m -1){
				F.assign( *(B + curr_row*n + n-1 - it_idx), one);
				curr_row++;
				it_idx++;
			}
		}
		for (size_t j=0; j<n; ++j)
			*(B + curr_row*n + j) = *(A + bk_idx*n + j);
		bk_idx++;
		curr_row++;
	}
	write_field (F, cout<<"A = "<<endl, A, m,n,n);
	write_field (F, cout<<"B = "<<endl, B, m+n,n,n);

	size_t *rp = new size_t[n];

	FFPACK::SpecRankProfile(F, m, n, A, n, deg,rp);

	size_t * P = new size_t[n];
	size_t * Q = new size_t[n+m];
	FFPACK::LUdivine(F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans,(int)m+n, n, B, n, P, Q, FFPACK::FfpackLQUP);

	printvect (cout<<"RankProfile (A) = "<<endl, rp, n)<<endl;

	printvect (cout<<"RankProfile (B) = "<<endl, Q, n)<<endl;




	return 0;
}
