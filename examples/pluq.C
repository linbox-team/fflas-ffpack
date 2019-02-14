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
#include "fflas-ffpack/utils/fflas_io.h"

using namespace std;

int main(int argc, char** argv) {

    if (argc != 3){
        std::cerr<<"Usage: pluq <p> <matrix>"<<std::endl;
        return -1;
    }

    int p = atoi(argv[1]);
    std::string file = argv[2];
    size_t m,n;

    // Creating the finite field Z/qZ
    Givaro::Modular<double> F(p);

    // Reading the matrix from a file
    double *A;
    FFLAS::ReadMatrix (file.c_str(),F,m,n,A);

    size_t * P = FFLAS::fflas_new<size_t>(m);
    size_t * Q = FFLAS::fflas_new<size_t>(n);

    size_t r = FFPACK::PLUQ (F, FFLAS::FflasNonUnit, m, n, A, n, P, Q);


    FFLAS::WritePermutation (std::cout<<"P = "<<std::endl,P,m);

    // Displays L and U separately
    double * L = FFLAS::fflas_new<double>(m * r);
    double * U = FFLAS::fflas_new<double>(r * n);
    FFPACK::getTriangular(F, FFLAS::FflasLower, FFLAS::FflasUnit, m, n, r, A, n, L, r, true);
    FFPACK::getTriangular(F, FFLAS::FflasUpper, FFLAS::FflasNonUnit, m, n, r, A, n, U, n, true);
    FFLAS::WriteMatrix(std::cout<<"L = "<<std::endl, F, m, r, L, r);
    std::cout<<"modulo "<<p<<std::endl;
    FFLAS::WriteMatrix(std::cout<<"U = "<<std::endl, F, r, n, U, n);
    std::cout<<"modulo "<<p<<std::endl;

    FFLAS::WritePermutation (std::cout<<"Q = "<<std::endl,Q,n);

    FFLAS::fflas_delete( P);
    FFLAS::fflas_delete( Q);
    FFLAS::fflas_delete( A);
    FFLAS::fflas_delete( L);
    FFLAS::fflas_delete( U);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
