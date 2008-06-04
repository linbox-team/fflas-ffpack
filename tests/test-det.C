/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
//--------------------------------------------------------------------------
//                        Test for det
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

int main(int argc, char** argv){

	int n;
	int nbit=atoi(argv[3]); // number of times the product is performed
	cerr<<setprecision(10);
	if (argc != 4)	{
		cerr<<"Usage : test-det <p> <A> <<i>"
		    <<endl
		    <<"         to compute the determinant of A mod p (i computations)"
		    <<endl;
		exit(-1);
	}
	Field F(atof(argv[1]));
	Field::Element * A;
	A = read_field(F,argv[2],&n,&n);
		
	Timer tim,t; t.clear();tim.clear(); 
	Field::Element d=0;
	for(int i = 0;i<nbit;++i){
		t.clear();
		t.start();
		d = FFPACK::Det (F, n, n, A, n);
		t.stop();
		tim+=t;
		if (i+1<nbit){
			delete[] A;
			A = read_field(F,argv[2],&n,&n);
		}
	}
		
	double mflops = 2.0/3.0*(n*n/1000000.0)*nbit*n/tim.usertime();
	F.write (cerr<<"n = "<<n<<" Det (A) = ",d)
		     << " mod "<<atoi(argv[1])<<" : t= "
		     << tim.usertime()/nbit 
		     << " s, Mffops = "<<mflops
		     << endl;
	
	cout<<n<<" "<<mflops<<" "<<tim.usertime()/nbit<<endl;
}
