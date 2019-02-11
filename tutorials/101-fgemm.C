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

    typedef Givaro::ModularBalanced<float> Ring;
    Ring F(101);

    Ring::Element * A, * B, * C;

    A = fflas_new(F,2,3);
    B = fflas_new(F,3,2);
    C = fflas_new(F,2,2);

    F.init(*(A+0));
    F.assign(*(A+0),F.one);
    F.init(*(A+1),2);
    F.init(*(A+2),3);
    F.init(*(A+3),5);
    F.init(*(A+4),7);
    F.init(*(A+5),11);

    Ring::Element t,u,v;
    F.init(t, 2); F.init(u, 4); F.init(v);

    F.assign(*(B+0),F.zero);		// B[0] <- 0
    F.assign(*(B+1),t);			// B[1] <- 2
    F.assign(*(B+2),u);			// B[2] <- 4
    F.add(v,t,u); F.assign(*(B+3),v);	// B[3] <- 2+4
    F.mul(*(B+4),t,u);			// B[4] <- 2*4
    F.add(*(B+5),u,v);			// B[5] <- 4+6

    FFLAS::WriteMatrix (std::cout << "A:=", F, 2, 3, A, 3) << std::endl;
    FFLAS::WriteMatrix (std::cout << "B:=", F, 3, 2, B, 2) << std::endl;

    fgemm (F, FflasNoTrans, FflasNoTrans, 2,2,3, F.one, A, 3, B, 2, F.zero, C, 2 );

    FFLAS::WriteMatrix (std::cout << "C:=", F, 2, 2, C,2) << std::endl;

    fflas_delete( A);
    fflas_delete( B);
    fflas_delete( C);



    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
