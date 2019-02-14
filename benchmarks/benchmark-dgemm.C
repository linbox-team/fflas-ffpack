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

//#include "goto-def.h"

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/args-parser.h"

#ifdef __FFLASFFPACK_USE_OPENMP
typedef FFLAS::OMPTimer TTimer;
#else
typedef FFLAS::Timer TTimer;
#endif

#ifndef __SGEMM__
typedef double Floats;
#define CBLAS_GEMM cblas_dgemm
#else
typedef float Floats;
#define CBLAS_GEMM cblas_sgemm
#endif


using namespace std;

int main(int argc, char** argv) {

    size_t iter = 1;
    int    q    = 1009;
    size_t n    = 2000;
    std::string file1 = "";
    std::string file2 = "";

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
        { 'n', "-n N", "Set the dimension of the matrix.",               TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &iter },
        { 'f', "-f FILE", "Set the first input file (empty for random).",   TYPE_STR , &file1 },
        { 'g', "-g FILE", "Set the second input file (empty for random).",  TYPE_STR , &file2 },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);


    typedef Givaro::ModularBalanced<Floats> Field;
    typedef Field::Element Element;

    Field F(q);

    TTimer chrono;
    double time=0.0;// time2=0.0;

    Element * A, * B, * C;

    if (iter>1) {
        if (!file1.empty()){
            FFLAS::ReadMatrix (file1.c_str(),F,n,n,A);
        }
        else{
            Field::RandIter G(F);
            A = FFLAS::fflas_new<Element>(n*n);
#pragma omp parallel for
            for (size_t i=0; i<n; ++i)
                for (size_t j=0; j<n; ++j)
                    G.random(*(A+i*n+j));
        }

        if (!file2.empty()){
            FFLAS::ReadMatrix (file2.c_str(),F,n,n,B);
        }
        else{
            Field::RandIter G(F);
            B = FFLAS::fflas_new<Element>(n*n);
#pragma omp parallel for
            for (size_t i=0; i<n; ++i)
                for (size_t j=0; j<n; ++j)
                    G.random(*(B+i*n+j));
        }

        C = FFLAS::fflas_new<Element>(n*n);

        CBLAS_GEMM (CblasRowMajor, CblasNoTrans, CblasNoTrans, n,n,n, F.one,
                    A, n, B, n, F.zero, C,n);

        FFLAS::fflas_delete( A);
        FFLAS::fflas_delete( B);
        FFLAS::fflas_delete( C);
    }

    for (size_t it=0;it<iter;++it){

        if (!file1.empty()){
            FFLAS::ReadMatrix (file1.c_str(),F,n,n,A);
        }
        else{
            Field::RandIter G(F);
            A = FFLAS::fflas_new<Element>(n*n);
#pragma omp parallel for
            for (size_t i=0; i<n; ++i)
                for (size_t j=0; j<n; ++j)
                    G.random(*(A+i*n+j));
        }

        if (!file2.empty()){
            FFLAS::ReadMatrix (file2.c_str(),F,n,n,B);
        }
        else{
            Field::RandIter G(F);
            B = FFLAS::fflas_new<Element>(n*n);
#pragma omp parallel for
            for (size_t i=0; i<n; ++i)
                for (size_t j=0; j<n; ++j)
                    G.random(*(B+i*n+j));
        }

        C = FFLAS::fflas_new<Element>(n*n);

        chrono.clear();
        chrono.start();
        CBLAS_GEMM (CblasRowMajor, CblasNoTrans, CblasNoTrans, n,n,n, F.one,
                    A, n, B, n, F.zero, C,n);
        chrono.stop();
        time+=chrono.usertime();

        FFLAS::fflas_delete( A);
        FFLAS::fflas_delete( B);
        FFLAS::fflas_delete( C);
    }

    // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14
    std::cout << "Time: " << time / double(iter)
    << " Gfops: " << (2.*double(n)/1000.*double(n)/1000.*double(n)/1000.0) / time * double(iter);
    FFLAS::writeCommandString(std::cout, as) << std::endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
