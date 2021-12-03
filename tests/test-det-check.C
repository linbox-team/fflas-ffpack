/*
 * Copyright (C) 2015 the FFLAS-FFPACK group
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
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

//#define ENABLE_ALL_CHECKINGS 1
#define ENABLE_CHECKER_Det 1
#define TIME_CHECKER_Det 1


#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/checkers/checkers_ffpack.h"
#include "fflas-ffpack/checkers/checkers_ffpack.inl"

using namespace FFLAS;

int main(int argc, char** argv) {
    size_t iter = 3 ;
    Givaro::Integer q = 131071;
    size_t MAXN = 1000;
    size_t n=0;
    uint64_t seed = getSeed();
    bool random_dim = false, random_rpm=false;

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).", TYPE_INTEGER , &q },
        { 'm', "-m M", "Set the dimension of A.", TYPE_INT , &n },
        { 'n', "-n N", "Set the dimension of A.", TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.", TYPE_INT , &iter },
        { 's', "-s N", "Set the seed.", TYPE_UINT64 , &seed },
        { 'r', "-r Y/N", "Set random RPM or not.", TYPE_BOOL, &random_rpm },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);
    if (n == 0) random_dim = true;

    typedef Givaro::Modular<double> Field;
    Field F(q);

    Field::RandIter Rand(F,seed);
    srandom(seed);

    size_t pass = 0;	// number of tests that have successfully passed

    for(size_t it=0; it<iter; ++it) {
#ifdef TIME_CHECKER_Det
        FFLAS::Timer init;init.start();
#endif
        if (random_dim) {
            n = random() % MAXN + 1;
        }

        Field::Element_ptr A = FFLAS::fflas_new(F,n,n);
        size_t *P = FFLAS::fflas_new<size_t>(n);
        size_t *Q = FFLAS::fflas_new<size_t>(n);

        // generate a random matrix A
        if (random_rpm)
            FFPACK::RandomMatrixWithRankandRandomRPM(F,n,n,n,A,n,Rand);
        else
            FFLAS::frand(F,Rand, n,n,A,n);
        init.stop();
        std::cerr << "init: " << init << std::endl;

        Field::Element det; F.init(det);
        try {
            FFPACK::ForceCheck_Det<Field> checker (Rand,n,A,n);
            Givaro::Timer chrono; chrono.start();
            FFPACK::Det(F,det,n,A,n,P,Q);
            chrono.stop();
            checker.check(det,A,n,P,Q);
            F.write(std::cerr << n << 'x' << n << ' ' << '(', det) << ')' << " Det verification PASSED\n" ;
#ifdef TIME_CHECKER_Det
            std::cerr << "Det COMPT: " << chrono << std::endl;
#endif
            pass++;
        } catch(FailureDetCheck &e) {
            F.write(std::cerr << n << 'x' << n << ' ' << '(', det) << ')' << " Det verification FAILED!\n";
        }

        FFLAS::fflas_delete(A,P,Q);
    }

    std::cerr << pass << "/" << iter << " tests SUCCESSFUL.\n";

    return (iter-pass);
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

