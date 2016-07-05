/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/* Copyright (c) FFLAS-FFPACK
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
*/

#include <iostream>
#include <givaro/modular.h>
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/Matio.h"

using namespace std;

int main(int argc, char** argv) {
	
	if (argc > 2){
		std::cerr<<"Usage: 2x2-pluq <p>"<<std::endl;
		return -1;
	}

	int p = (argc>1?atoi(argv[1]):5);
		// Creating the finite field Z/qZ
	Givaro::Modular<double> F(p);

	size_t m(2),n(2);
	double A[4] {1,2,3,4};
	write_field(F,std::cout<<"A = "<<std::endl,A,m,n,n);

    size_t * P = FFLAS::fflas_new<size_t>(m);
    size_t * Q = FFLAS::fflas_new<size_t>(n);

    FFPACK::PLUQ (F, FFLAS::FflasNonUnit, m, n, A, n, P, Q);

	write_perm(std::cout<<"P = "<<std::endl,P,m);
	write_field(F,std::cout<<"LU = "<<std::endl,A,m,n,n);
	write_perm(std::cout<<"Q = "<<std::endl,Q,n);

    FFLAS::fflas_delete( P);
    FFLAS::fflas_delete( Q);
		
	return 0;
}

