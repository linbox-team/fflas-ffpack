/* Copyright (c) FFLAS-FFPACK
 * Written by Clément Pernet <clement.pernet@imag.fr> and Philippe Ledent
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

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/args-parser.h"

// declare that the call to openblas_set_numthread will be made here, hence don't do it
// everywhere in the call stack
#define __FFLASFFPACK_OPENBLAS_NT_ALREADY_SET 1

using namespace std;

int main(int argc, char** argv) {

 
#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t iter = 3;
    int    q    = 131071;
    size_t    n    = 2000;
    size_t    t    = 0;
    std::string file = "";

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
        { 'n', "-n N", "Set the dimension of the matrix.",               TYPE_INT , &n },
        { 't', "-t T", "Set the base case threshold.",                   TYPE_INT , &t },
        { 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &iter },
        { 'f', "-f FILE", "Set the input file (empty for random).",  TYPE_STR , &file },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    typedef Givaro::ModularBalanced<double> Field;
    typedef Field::Element Element;

    Field F(q);
    Field::Element * A, *B;

    FFLAS::Timer chrono;
    double time=0.0;

    A = FFLAS::fflas_new<Element>(n*n);
    B = FFLAS::fflas_new<Element>(n*n);
    FFPACK::RandomTriangularMatrix (F, n, n, FFLAS::FflasUpper,FFLAS::FflasNonUnit,true,B,n);

    for (size_t i=0;i<=iter;++i){
        if (!file.empty()){
            FFLAS::ReadMatrix (file.c_str(),F,n,n,A);
        }
        else {
            FFLAS::fassign(F,n,n,B,n,A,n);
        }
        chrono.clear();
        if (i) chrono.start();
        if (t)
            FFPACK::ftrtri(F, FFLAS::FflasUpper, FFLAS::FflasNonUnit, n, A, n, t);
        else
            FFPACK::ftrtri(F, FFLAS::FflasUpper, FFLAS::FflasNonUnit, n, A, n);
        if (i) chrono.stop();

        time+=chrono.realtime();
    }
    FFLAS::fflas_delete (A);
    FFLAS::fflas_delete (B);

    // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14
#define CUBE(x) ((x)*(x)*(x))
    std::cout << "Time: " << time / double(iter)
    << " Gfops: " <<  CUBE(double(n)/1000.)/(3*time) * double(iter);
    FFLAS::writeCommandString(std::cout, as) << std::endl;


    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
