/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) Fflas-Ffpack
 * Written by Clément Pernet
 * This file is Free Software and part of Fflas-Ffpack.
 * See COPYING for license information.
 */

//--------------------------------------------------------------------------
//                        Test for charpoly
//
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------

#include <iomanip>
#include <iostream>
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/field/modular-positive.h"
// #include "fflas-ffpack/field/modular-int.h"
#include "timer.h"
#include "utils/Matio.h"
#include "fflas-ffpack/ffpack/ffpack.h"

#include "fflas-ffpack/utils/args-parser.h"


using namespace std;

//typedef ModularBalanced<double> Field;
//typedef ModularBalanced<float> Field;
//typedef Modular<double> Field;
//typedef Modular<float> Field;
typedef FFPACK:: Modular<int> Field;
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

template<class Field>
bool launch_test(const Field & F, typename Field::Element * A, int n,
		 size_t p, size_t nbit)
{
	Timer tim,t; t.clear();tim.clear();
	list<Polynomial> P_list;
	for(size_t i = 0;i<nbit;++i){
		P_list.clear();
		t.clear();
		t.start();

		FFPACK::CharPoly (F, P_list, n, A, n);
		t.stop();
		tim+=t;
		/*  test */
		// apply P_list.A.V and check 0 for random V
	}

#ifdef _FF_TIMING
	double mflops = (2+(2.0/3.0))*(n*n/1000000.0)*nbit*n/tim.usertime();
	list<Polynomial>::iterator P_it = P_list.begin();
	for (;P_it!=P_list.end();++P_it)
		printPolynomial ( F, *P_it);

	cerr<<"n = "<<n<<" #inv. fact = "<<P_list.size()<<" Charpoly (A) over Z/ "<<atoi(argv[1])<<"Z : t= "
	<< tim.usertime()/nbit
	<< " s, Mffops = "<<mflops
	<< endl;

	cout<<n<<" "<<P_list.size()<<" "<<mflops<<" "<<tim.usertime()/nbit<<endl;
#endif
	return true ;

}

int main(int argc, char** argv)
{

	cerr<<setprecision(10);

	static size_t       p = 13; // characteristique
	static size_t       nbit = 2; // repetitions
	static int          n = 100;
	static std::string  file = "" ; // file where


	static Argument as[] = {
		{ 'p', "-p P", "Set the field characteristic.", TYPE_INT , &p },
		{ 'n', "-n N", "Set the size of the matrix.", TYPE_INT , &p },
		{ 'r', "-r R", "Set number of repetitions.", TYPE_INT , &nbit },
		{ 'f', "-f file", "Set input file", TYPE_STR, &file },
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,as);

	Field F((long unsigned int)p);
	Field::Element * A;
	if (!file.empty()) {
		const char * filestring = file.c_str();
		A = read_field<Field>(F,const_cast<char*>(filestring),&n,&n);
		bool passed = launch_test<Field>(F,A,n,p,nbit);
		delete[] A;
		return !passed ;
	}
	else {
		std::cerr << std::endl << "##################################"<< std::endl;
		std::cerr << std::endl << "  **** not implemented yet ! ***" << std::endl;
		std::cerr << std::endl << "##################################"<< std::endl;
		// create A random
		// create A diagonal
		// create A nilpotent
		// create A non invertible
		return false ;
	}

}

