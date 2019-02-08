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

#include "fflas-ffpack/fflas-ffpack.h"

/**
 * This example computes the matrix multiplication
 * over a defined finite field.
 *
 * Outputs the product of the matrix given as input.
 */
int main(int argc, char** argv)
{
    if (argc != 4) {
        std::cerr << "Usage: matmul <p> <matrixA> <matrixB>" << std::endl;
        return -1;
    }

    int p = atoi(argv[1]);
    std::string fileA = argv[2];
    std::string fileB = argv[3];

    // Creating the finite field Z/pZ (assuming p>2, use Modular<float> for p=2)
    Givaro::ModularBalanced<double> F(p);

    // Reading the matrices
    double* A, *B, *C;
    size_t mA, nA, mB, nB;
    FFLAS::ReadMatrix(fileA.c_str(), F, mA, nA, A);
    FFLAS::ReadMatrix(fileB.c_str(), F, mB, nB, B);
    C = FFLAS::fflas_new(F, mA, nB);

    if (nA != mB) {
        std::cerr << "Those matrices can't be multiplied together" << std::endl;
        return -1;
    }

    FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mA, nB, nA,
                 F.one, A, mA, B, mB, F.zero, C, mA);

    FFLAS::WriteMatrix(std::cout << "Result" << std::endl, F, mA, nB, C, mA) << std::endl;

    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);
    FFLAS::fflas_delete(C);

    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
