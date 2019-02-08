/* Copyright (c) 2017 FFLAS-FFPACK
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
#include <fflas-ffpack/utils/fflas_io.h>
#include <fflas-ffpack/utils/args-parser.h>

#include <iostream>
#include <array>

using namespace FFLAS;

int main(int argc, char** argv) {

    typedef Givaro::Modular<float> Ring;
    Ring F(11);

    Ring::Element L[4]{1,0,2,3};
    Ring::Element B[2]{4,5};

    size_t m(2);

    FFLAS::WriteMatrix (std::cout << "L:=", F, m, m, L, m) << std::endl;
    FFLAS::WriteMatrix (std::cout << "B:=", F, m, 1, B, 1) << std::endl;

    // In place system solve
    ftrsv (F, FflasLower,FflasNoTrans,FflasNonUnit, m, L, m, B, 1);

    FFLAS::WriteMatrix (std::cout << "X:=", F, m, 1, B, 1) << std::endl;
    std::cerr << "0 = L.X - B mod " << F.characteristic() << ';' << std::endl;

    Ring::Element U[4]{3,2,0,5}, C[2]{4,7};


    FFLAS::WriteMatrix (std::cout << "U:=", F, m, m, U, m) << std::endl;
    FFLAS::WriteMatrix (std::cout << "C:=", F, m, 1, C, 1) << std::endl;

    // In place system solve
    ftrsv (F, FflasUpper, FflasNoTrans ,FflasNonUnit, m, U, m, B, 1);

    FFLAS::WriteMatrix (std::cout << "X:=", F, m, 1, C, 1) << std::endl;
    std::cerr << "0 = U.X - C mod " << F.characteristic() << ';' << std::endl;

    Ring::Element D[2]{4,5};

    FFLAS::WriteMatrix (std::cout << "L:=", F, m, m, L, m) << std::endl;
    FFLAS::WriteMatrix (std::cout << "D:=", F, m, 1, D, 1) << std::endl;

    // In place system solve
    ftrsv (F, FflasLower,FflasTrans,FflasNonUnit, m, L, m, D, 1);

    FFLAS::WriteMatrix (std::cout << "X:=", F, m, 1, D, 1) << std::endl;
    std::cerr << "0 = Transpose(X).L - Transpose(D) mod " << F.characteristic() << ';' << std::endl;

    Ring::Element E[2]{4,7};


    FFLAS::WriteMatrix (std::cout << "U:=", F, m, m, U, m) << std::endl;
    FFLAS::WriteMatrix (std::cout << "E:=", F, m, 1, E, 1) << std::endl;

    // In place system solve
    ftrsv (F, FflasUpper, FflasTrans ,FflasNonUnit, m, U, m, E, 1);

    FFLAS::WriteMatrix (std::cout << "X:=", F, m, 1, E, 1) << std::endl;
    std::cerr << "0 = Transpose(X).U - Transpose(E) mod " << F.characteristic() << ';' << std::endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
