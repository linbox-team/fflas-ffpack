/* Copyright (c) FFLAS-FFPACK
 * Written by Hippolyte Signargout <hippolyte.signargout@ens-lyon.fr>
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

// Template from benchmark-quasisep.C

#define __FFLASFFPACK_OPENBLAS_NT_ALREADY_SET 1

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/args-parser.h"


using namespace std;
using namespace FFLAS;
using namespace FFPACK;

template<class Field>
void run_with_field(int q, size_t n, size_t m, size_t s, size_t r, size_t iter, uint64_t seed){

    Field F(q);
    typedef typename Field::Element_ptr Element_ptr;

    FFLAS::Timer chrono;
    Element_ptr A, B, TS;

    double time_gen = 0, time_cbxts =0;
    for (size_t i=0;i<iter;++i){

        A = FFLAS::fflas_new (F, n, n);
        B = FFLAS::fflas_new (F, n, n);
        size_t lda=n;
        TS = FFLAS::fflas_new (F, n, m);
        size_t ldts = m;
        typename Field::RandIter G (F, seed);
        RandomLTQSMatrixWithRankandQSorder (F,n,r,s,A,lda,G);
        RandomLTQSMatrixWithRankandQSorder (F,n,r,s,B,lda,G);
        size_t * p = FFLAS::fflas_new<size_t> (n);
        for (size_t i = 0; i < ceil(n/2.); i++)
            {
                p[i] = n - i - 1;
            }
        applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), A, n, p);
        applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), B, n, p);
        faddin (F, n, n, B, n, A, n);
        RandomMatrix(F, n, m, TS, ldts, G);

        size_t rs = n%s;           // Size of the partial block
        size_t ls = (rs)? rs: s;   // Size of the last block

        Element_ptr D = fflas_new (F, n, s);
        Element_ptr P = fflas_new (F, n - s, s);
        Element_ptr Q = fflas_new (F, n - ls, s);
        Element_ptr R = fflas_new (F, ((n > (s + ls))? (n - s - ls): 0), s);
        Element_ptr U = fflas_new (F, n - ls, s);
        Element_ptr V = fflas_new (F, n - ls, s);
        Element_ptr W = fflas_new (F, ((n > (s + ls))? (n - s - ls): 0), s);

        chrono.clear();
        chrono.start();
        DenseToSSS (F, n, s, P, s, Q, s, R, s, U, s, V, s, W, s, D, s, A, n);
        chrono.stop();

        time_gen+=chrono.usertime();


        Element_ptr Res = fflas_new(F, n, m); // Inadequate name
 
        chrono.clear();
        chrono.start();
        productSSSxTS(F, n, m, s, F.one, P, s, Q, s, R, s, U, s, V, s, W, s,
                      D, s, TS, m, F.zero, Res, m);
        chrono.stop();

        time_cbxts += chrono.usertime();
        FFLAS::fflas_delete(A, D, P, Q, R, U, V, W, B, p); // Could be done once for all iters
        FFLAS::fflas_delete(TS, Res);
    }
    // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14
    std::cout << "Time: " << (time_gen + time_cbxts) / double(iter)
              << " Gfops: Irrelevant Specific times: " << time_gen / double(iter)
              <<" (for construction), " << time_cbxts / double(iter)<<" (for CB x TS)" ;
}

int main(int argc, char** argv) {

#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t iter = 10;
    int    q    = 131071;
    size_t    n    = 2167;
    size_t    m    = 455;
    size_t    t    = 236;
    size_t    r    = 1100;
    uint64_t seed = FFLAS::getSeed();

    Argument as[] = {
                     { 'q', "-q Q", "Set the field characteristic (-1 for the ring ZZ).",     TYPE_INT , &q },
                     { 'n', "-n N", "Set the order of the square matrix A.",                  TYPE_INT , &n },
                     { 'm', "-m M", "Set the column dimension of n x m RHS matrix B.",        TYPE_INT , &m },
                     { 't', "-t T", "Set the quasiseparability order of A.",                  TYPE_INT , &t },
                     { 'r', "-r R", "Set the rank of each upper/lower triangular part of A.", TYPE_INT , &r },
                     { 'i', "-i R", "Set number of repetitions.",                             TYPE_INT , &iter },
                     { 's', "-s S", "Sets seed.",                                             TYPE_INT , &seed },
                     END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    run_with_field<Givaro::ModularBalanced<double> >(q, n, m, t, r, iter, seed);

    std::cout << "( ";
    FFLAS::writeCommandString(std::cout, as) << ")" << std::endl;
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
