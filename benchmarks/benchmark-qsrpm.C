/* Copyright (c) FFLAS-FFPACK
 * Written by Clement Pernet <clement.pernet@univ-grenoble-alpes.fr>
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

#define __FFLASFFPACK_OPENBLAS_NT_ALREADY_SET 1

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>
#include <givaro/givpoly1.h>

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/args-parser.h"


using namespace std;
using namespace FFLAS;
using namespace FFPACK;

void run (size_t n, size_t t, size_t r, size_t iter){

    FFLAS::Timer chrono;

    double* time = new double[iter];
    for (size_t i=0;i<iter;++i){

        size_t * rows = new size_t[r];
        size_t * cols = new size_t[r];

        chrono.clear();
        chrono.start();
        RandomLTQSRankProfileMatrix (n, r, t, rows, cols);
        chrono.stop();

        time[i] = chrono.usertime();

        FFLAS::fflas_delete(rows,cols);
    }
    std::sort(time, time+iter);
    double mediantime = time[iter/2];
    delete[] time;
    std::cout << "Time: " << mediantime ;

}

int main(int argc, char** argv) {

#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t iter = 10;
    size_t    n    = 2000;
    size_t    t    = 200;
    size_t    r    = 1000;
    uint64_t seed = FFLAS::getSeed();
    Argument as[] = {
        { 'n', "-n N", "Set the order of the square matrix A.",               TYPE_INT , &n },
        { 't', "-t T", "Set the quasiseparability order of A.",  TYPE_INT , &t },
        { 'r', "-r R", "Set the rank of each upper/lower triangular part of A.",  TYPE_INT , &r },
        { 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &iter },
        { 's', "-s S", "Sets seed.", TYPE_INT , &seed },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    srand(seed);

    run (n, t, r, iter);
    
    FFLAS::writeCommandString(std::cout, as) << std::endl;
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
