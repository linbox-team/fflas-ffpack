/*
 * Copyright (C) 2015 the FFLAS-FFPACK group
 * Written by Ashley Lesdalons <Ashley.Lesdalons@e.ujf-grenoble.fr>
 *
 * This file is Free Software and part of FFLAS-FFPACK.
 *
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
 *
 */

//--------------------------------------------------------------------------
//          Test for Checker_PLUQ
//--------------------------------------------------------------------------

#define ENABLE_ALL_CHECKINGS 1


#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"

using namespace FFLAS;

int main(int argc, char** argv) {
    size_t iter = 3 ;
    Givaro::Integer q = 131071;
    size_t MAXM = 1000;
    size_t MAXN = 1000;
    size_t m=0,n=0;
    uint64_t seed = getSeed();
    bool random_dim = false;

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).", TYPE_INTEGER , &q },
        { 'm', "-m M", "Set the row dimension of A.", TYPE_INT , &m },
        { 'n', "-n N", "Set the col dimension of A.", TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.", TYPE_INT , &iter },
        { 's', "-s N", "Set the seed.", TYPE_UINT64 , &seed },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);
    if (m == 0 || n == 0) random_dim = true;

    srandom (seed);

    typedef Givaro::Modular<double> Field;
    Field F(q);

    Field::RandIter Rand(F,seed);

    size_t pass = 0;	// number of tests that have successfully passed

    FFLAS_DIAG Diag = FflasNonUnit;
    for(size_t it=0; it<iter; ++it) {
        if (random_dim) {
            m = random() % MAXM + 1;
            n = random() % MAXN + 1;
        }

        Field::Element_ptr A = fflas_new(F,m,n);
        size_t *P = fflas_new<size_t>(m);
        size_t *Q = fflas_new<size_t>(n);

        // generate a random matrix A
        PAR_BLOCK { pfrand(F,Rand, m,n,A,m/MAX_THREADS); }

        try {
            FFPACK::PLUQ(F, Diag, m, n, A, n, P, Q);
            std::cerr << m << 'x' << n << ' ' << Diag << " pluq verification PASSED\n";
            pass++;
        } catch(FailurePLUQCheck &e) {
            std::cerr << m << 'x' << n << ' ' << Diag << " pluq verification FAILED!\n";
        }

        fflas_delete(A,P,Q);
        Diag = (Diag == FflasNonUnit) ? FflasUnit : FflasNonUnit;
    }

    std::cerr << pass << "/" << iter << " tests SUCCESSFUL.\n";

    return (iter-pass);
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
