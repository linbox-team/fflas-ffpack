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
#include <vector>
#include <givaro/modular.h>
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/fflas_io.h"

using namespace std;

int main(int argc, char** argv) {

    if (argc > 2){
        std::cerr<<"Usage: 2x2-pluq <p>"<<std::endl;
        return -1;
    }

    int64_t p = (argc>1?atoi(argv[1]):5);
    // Creating the finite field Z/pZ
    Givaro::Modular<double> F(p);

    size_t m(2),n(2);
    double A[4] {1,2,3,4};
    FFLAS::WriteMatrix (std::cout<<"A = "<<std::endl,F,m,n,A,n);

    size_t * P = FFLAS::fflas_new<size_t>(m);
    size_t * Q = FFLAS::fflas_new<size_t>(n);

    FFPACK::PLUQ (F, FFLAS::FflasNonUnit, m, n, A, n, P, Q);

    FFLAS::WritePermutation(std::cout<<"P = "<<std::endl,P,m);
    FFLAS::WriteMatrix (std::cout<<"LU = "<<std::endl,F,m,n,A,n)<< " modulo " << p << std::endl;
    FFLAS::WritePermutation(std::cout<<"Q = "<<std::endl,Q,n);

    FFLAS::fflas_delete( P);
    FFLAS::fflas_delete( Q);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
