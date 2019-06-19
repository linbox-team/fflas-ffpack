/* Copyright (c) FFLAS-FFPACK
 * Written by Clément Pernet <clement.pernet@imag.fr>
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

// declare that the call to openblas_set_numthread will be made here, hence don't do it
// everywhere in the call stack
#define __FFLASFFPACK_OPENBLAS_NT_ALREADY_SET 1

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/args-parser.h"


using namespace std;
using namespace FFLAS;
using namespace FFPACK;
int main(int argc, char** argv) {

 #ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t iter = 3;
    int    q    = 131071;
    size_t    n    = 1000;
    size_t    k    = 1000;
    size_t threshold = 64;
    bool up =true;
    std::string file = "";

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
        { 'n', "-n N", "Set the dimension of the matrix C.",               TYPE_INT , &n },
        { 'k', "-k K", "Set the other dimension of the matrix A.",               TYPE_INT , &k },
        { 'u', "-u yes/no", "Updates an upper/lower triangular matrix.",  TYPE_BOOL , &up },
        { 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &iter },
        { 't', "-t T", "Set the threshold to the base case.",                     TYPE_INT , &threshold },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    typedef Givaro::ModularBalanced<double> Field;
    typedef Field::Element Element;

    Field F(q);
    Field::Element * A, *C;

    FFLAS::Timer chrono;
    double time=0.0;

    FFLAS_UPLO uplo = up?FflasUpper:FflasLower;
    for (size_t i=0;i<=iter;++i){
        A = fflas_new<Element>(n*k);
        size_t lda=k;
        C = fflas_new<Element>(n*n);
        size_t ldc=n;
        Field::RandIter G(F);
        RandomMatrix (F, n, k, A, k, G);
        RandomTriangularMatrix (F, n,n,uplo, FflasNonUnit, true, C, ldc, G);
        Field::Element_ptr D = FFLAS::fflas_new(F,k,1);
        Givaro::GeneralRingNonZeroRandIter<Field,Field::RandIter> nzG (G);
        for (size_t i=0; i<k; i++)
            nzG.random(D[i]);
        chrono.clear();
        if (i) chrono.start();
        fsyrk (F, uplo, FflasTrans, n, k, F.mOne, A, lda, D, 1, F.one, C, ldc, threshold);
        if (i) chrono.stop();

        time+=chrono.usertime();
        FFLAS::fflas_delete( A);
        FFLAS::fflas_delete( C);
        FFLAS::fflas_delete( D);
    }

    // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14
#define CUBE(x) ((x)*(x)*(x))
    std::cout << "Time: " << time / double(iter)
    << " Gfops: " << CUBE(double(n)/1000.)/ time * double(iter);
    FFLAS::writeCommandString(std::cout, as) << std::endl;
    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
