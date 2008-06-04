/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
//--------------------------------------------------------------------------
//                        Test for charpoly
//                  
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------

#include <iomanip>
#include <iostream>
#include "fflas-ffpack/modular-balanced.h"
#include "timer.h"
#include "Matio.h"
#include "fflas-ffpack/ffpack.h"




using namespace std;

typedef ModularBalanced<double> Field;
typedef vector<Field::Element> Polynomial;


template <class Field, class Polynomial>
void printPolynomial (const Field &F, const Polynomial &v) 
{
	for (int i = v.size () - 1; i >= 0; i--) {
		F.write (cout, v[i]);
		if (i > 0)
			cout << " x^" << i << " + ";
	}
	cout << endl;
}

int main(int argc, char** argv){

	int n;
	cerr<<setprecision(10);
	if (argc != 4)	{
		cerr<<"Usage : test-charpoly <p> <A> <i>"
		    <<endl
		    <<"        to compute the characteristic polynomial of A mod p (i computations)"
		    <<endl;
		exit(-1);
	}
	int nbit=atoi(argv[3]); // number of times the product is performed
	Field F((long unsigned int)atoi(argv[1]));
	Field::Element * A;
	A = read_field(F,argv[2],&n,&n);
		
	Timer tim,t; t.clear();tim.clear(); 
	list<Polynomial> P_list;
	for(int i = 0;i<nbit;++i){
		P_list.clear();
		t.clear();
		t.start();
		
		FFPACK::CharPoly (F, P_list, n, A, n);
		t.stop();
		tim+=t;
		delete[] A;
		if (i+1<nbit){
		  A = read_field(F,argv[2],&n,&n);
		}
	}
		
	double mflops = (2+(2.0/3.0))*(n*n/1000000.0)*nbit*n/tim.usertime();
	list<Polynomial>::iterator P_it = P_list.begin();
	for (;P_it!=P_list.end();++P_it)
		printPolynomial ( F, *P_it);
	
	cerr<<"n = "<<n<<" #inv. fact = "<<P_list.size()<<" Charpoly (A) over Z/ "<<atoi(argv[1])<<"Z : t= "
		     << tim.usertime()/nbit 
		     << " s, Mffops = "<<mflops
		     << endl;
	
	cout<<n<<" "<<P_list.size()<<" "<<mflops<<" "<<tim.usertime()/nbit<<endl;
}
