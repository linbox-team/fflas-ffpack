/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Clément Pernet
 * This file is Free Software and part of FFLAS-FFPACK.
 *
 * ========LICENCE========
 * This file is part of the library FFLAS-FFPACK.
 *
 * FFLAS-FFPACK is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
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
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/ffpack/ffpack.h"

#include "fflas-ffpack/utils/args-parser.h"


using namespace std;

//typedef ModularBalanced<double> Field;
typedef Givaro::ModularBalanced<double> Field;
//typedef Givaro::Modular<double> Field;
//typedef Givaro::Modular<float> Field;
//typedef Givaro::Modular<int> Field;
//typedef Givaro::Modular<double> Field;

typedef vector<Field::Element> Polynomial;

using namespace FFPACK;
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
		 size_t p, size_t nbit, FFPACK::FFPACK_CHARPOLY_TAG CT)
{
 FFLAS::Timer tim,t; t.clear();tim.clear();
	list<Polynomial> P_list;
	for(size_t i = 0;i<nbit;++i){
		P_list.clear();
		t.clear();
		t.start();

		FFPACK::CharPoly (F, P_list, n, A, n, CT);
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

	F.write(cerr<<"n = "<<n<<" #inv. fact = "<<P_list.size()<<" Charpoly (A) over ") << " : t= "
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

	static size_t       p = 13; // characteristic
	static size_t       nbit = 2; // repetitions
	static int          n = 100;
	static std::string  file = "" ; // file where
	static int variant =0;

	static Argument as[] = {
		{ 'p', "-p P", "Set the field characteristic.", TYPE_INT , &p },
		{ 'n', "-n N", "Set the size of the matrix.", TYPE_INT , &p },
		{ 'r', "-r R", "Set number of repetitions.", TYPE_INT , &nbit },
		{ 'f', "-f file", "Set input file", TYPE_STR, &file },
		{ 'a', "-a algorithm", "Set the algorithm variant", TYPE_INT, &variant },
		END_OF_ARGUMENTS
	};
	FFLAS::parseArguments(argc,argv,as);

	FFPACK::FFPACK_CHARPOLY_TAG CT;
	switch (variant){
	    case 0: CT = FfpackLUK; break;
	    case 1: CT = FfpackKG; break;
	    case 2: CT = FfpackDanilevski; break;
	    case 3: CT = FfpackKGFast; break;
	    case 4: CT = FfpackKGFastG; break;
	    case 5: CT = FfpackHybrid; break;
	    case 6: CT = FfpackArithProg; break;
	    default: CT = FfpackLUK; break;
	}
	Field F((long unsigned int)p);
	Field::Element * A;
	if (!file.empty()) {
		const char * filestring = file.c_str();
		A = read_field<Field>(F,const_cast<char*>(filestring),&n,&n);
		bool passed = launch_test<Field>(F,A,n,p,nbit,CT);
		FFLAS::fflas_delete( A);
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

