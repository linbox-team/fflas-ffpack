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
#include <fflas-ffpack/utils/fflas_io.h>
#include <fflas-ffpack/utils/args-parser.h>

#include <iostream>

using namespace FFLAS;

int main(int argc, char** argv) {

    Givaro::Modular<float> F(11);

    constexpr size_t m(2),k(3),n(1);
    float A[m*k]{1,2,3,4,5,6}, B[k*n]{7,8,9}, * C;

    C = fflas_new(F,m,n);

    // A is mxk with leading dimension k
    FFLAS::WriteMatrix (std::cout << "A:="<<std::endl, F, m, k, A, k) << std::endl;
    // B is kxn with leading dimension n
    FFLAS::WriteMatrix (std::cout << "B:="<<std::endl, F, k, n, B, n) << std::endl;

    // C <-- 1. A.B + 0. C
    fgemm (F, FflasNoTrans, FflasNoTrans, m, n, k, F.one, A, k, B, n, F.zero, C, n);

    // C is mxn with leading dimension n
    FFLAS::WriteMatrix (std::cout << "C:="<<std::endl, F, m, n, C, n) << " modulo 11" << std::endl;

    fflas_delete( C);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
