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

#include <givaro/modular.h>
#include <iostream>

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/fflas_io.h"

/**
 * This example computes the determinant of a matrix
 * over a defined finite field.
 *
 * Outputs the determinant.
 */
int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "Usage: det <p> <matrix>" << std::endl;
        return -1;
    }

    int p = atoi(argv[1]);
    std::string file = argv[2];

    // Creating the finite field Z/qZ
    Givaro::Modular<double> F(p);

    // Reading the matrix from a file
    double* A;
    size_t m, n;
    FFLAS::ReadMatrix(file.c_str(), F, m, n, A);

    double d;
    FFPACK::Det(F, d, n, A, n);
    std::cout << d << std::endl;

    FFLAS::fflas_delete(A);

    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
