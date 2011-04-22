/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
//--------------------------------------------------------------------------
//          Test for the krylov-elimination
//--------------------------------------------------------------------------
// usage: test-krylov-elim p A, to compute the rank profile of the (n+m)xn matrix B
// formed by the n identity vectors and the mxn matrix A over Z/pZ
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//#define DEBUG 0

#include <iostream>
#include <iomanip>
#include <list>
#include <vector>
#include "Matio.h"
#include "timer.h"
using namespace std;
#include "fflas-ffpack/field/modular-positive.h"
#include "fflas-ffpack/ffpack/ffpack.h"


typedef Modular<double> Field;

template<class T>
std::ostream& printvect(std::ostream& o, vector<T>& vect){
	for(size_t i=0; i < vect.size()-1; ++i)
		o << vect[i] << " " ;
	return o << vect[vect.size()-1] << std::endl;
	}

int main(int argc, char** argv){

	int m,n;
	cout<<setprecision(20);

	if (argc!=4){
		cerr<<"usage : test-frobenius <p> <A> <c>"<<endl
	 	    <<"         to compute the frobenius normal form of the matrix A over Z/pZ, with conditonning parameter c"
		    <<endl;
		exit(-1);
	}
	Field F(  long(atoi(argv[1])));
	Field::Element one;
	F.init(one, 1UL);
	Field::Element * A = read_field<Field> (F,argv[2],&m,&n);
	size_t c = atoi(argv[3]);

	std::list<vector<Field::Element> > frobForm;
	Timer tim;
	tim.clear();
	tim.start();
	FFPACK::CharpolyArithProg (F, frobForm, n, A, n, c);
	tim.stop();
	std::list<vector<Field::Element> >::iterator it = frobForm.begin();
	while(it != frobForm.end()){
		printvect (cout, *(it++));
	}
	cerr<<c<<" "<<tim.usertime()<<" "<<4.55*n*n/1000000.0*n/tim.usertime()<<endl;
	delete[] A;
	return 0;
}
