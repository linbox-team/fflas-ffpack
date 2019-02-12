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
 * This example solve the quare system defined by the input
 * over a defined finite field.
 *
 */
int main(int argc, char** argv)
{
    if (argc != 4) {
        std::cerr << "Usage: solve <p> <matrixA> <matrixB>" << std::endl;
        return -1;
    }

    int p = atoi(argv[1]);
    std::string fileA = argv[2];
    std::string fileB = argv[3];

    // Creating the finite field Z/pZ (assuming p>2, use Modular<float> for p=2)
    Givaro::ModularBalanced<double> F(p);

    // Reading the matrices
    double* A, *B, *X;
    size_t mA, nA, mB, nB;
    int info = 0;
    FFLAS::ReadMatrix(fileA.c_str(), F, mA, nA, A);
    FFLAS::ReadMatrix(fileB.c_str(), F, mB, nB, B);
    X = FFLAS::fflas_new(F, nA, nB);

    if (mA != nA || mA != mB) {
        std::cerr << "At least one input matrix has not the right dimension";
        std::cerr  << std::endl;
        return -1;
    }

    size_t rank = FFPACK::fgesv(F, FFLAS::FflasLeft, mA, nA, nB, A, nA, X, nB,
                                B, nB, &info);

    if (info){
        std::cout << "System is inconsistent" << std::endl;
        return -1;
    }

    std::cout << "rank: " << rank << std::endl;

    FFLAS::WriteMatrix(std::cout << "debug B" << std::endl, F, mB, nB, B, nB) << std::endl;
    FFLAS::WriteMatrix(std::cout << "solution" << std::endl, F, nA, nB, X, nB) << std::endl;

    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);
    FFLAS::fflas_delete(X);

    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
