/* Copyright (c) FFLAS-FFPACK
* Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
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


#include <fflas-ffpack/fflas-ffpack-config.h>
#include <givaro/modular-balanced.h>
#include <fflas-ffpack/fflas/fflas.h>
#include <fflas-ffpack/utils/timer.h>
#include <fflas-ffpack/utils/Matio.h>
#include <fflas-ffpack/utils/args-parser.h>

#include <iostream>

using namespace FFLAS;

int main(int argc, char** argv) {

	typedef Givaro::Modular<float> Ring;
	Ring F(11);

	Ring::Element A[4]{1,2,3,4}, B[4]{5,6,7,8}, * C;

        size_t m(2),k(2),n(2);

	C = fflas_new(F,m,n);
	
        // A is mxk with leading dimension k
	write_field(F, std::cout << "A:=", A, m, k, k, true) << std::endl;
        // B is kxn with leading dimension n
	write_field(F, std::cout << "B:=", B, k, n, n, true) << std::endl;

	fgemm (F, FflasNoTrans, FflasNoTrans, m, n, k, F.one, A, m, B, n, F.zero, C, n);

        // C is mxn with leading dimension n
	write_field(F, std::cout << "C:=", C, m, n, n, true) << " modulo 11" << std::endl;
	
	fflas_delete( C);

  return 0;
}

