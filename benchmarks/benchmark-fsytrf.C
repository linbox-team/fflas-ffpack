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

//#define __FFPACK_FSYTRF_BC_RL
#undef __FFPACK_FSYTRF_BC_RL
#define __FFPACK_FSYTRF_BC_CROUT

#define __FFLASFFPACK_OPENBLAS_NT_ALREADY_SET 1

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/args-parser.h"


using namespace std;

int main(int argc, char** argv) {

#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t iter = 3;
    int    q    = 131071;
    size_t    n    = 1000;
    size_t threshold = 64;
    size_t rank = 500;
    bool up =true;
    bool rpm =true;
    bool grp =true;
    bool par =false;
    int t=MAX_THREADS;
    std::string file = "";

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
        { 'n', "-n N", "Set the dimension of the matrix.",               TYPE_INT , &n },
        { 'u', "-u yes/no", "Computes a UTDU (true) or LDLT decomposition (false).",  TYPE_BOOL , &up },
        { 'm', "-m yes/no", "Use the rank profile matrix revealing algorithm.", TYPE_BOOL , &rpm },
        { 'r', "-r R", "Set the rank (for the RPM version).", TYPE_INT , &rank },
        { 'g', "-g yes/no", "Matrix with generic rank profile (yes) or random rank profile matrix (no).", TYPE_BOOL , &grp },
        { 'i', "-i I", "Set number of repetitions.",                     TYPE_INT , &iter },
        { 'c', "-c C", "Set the cross-over point to the base case.",            TYPE_INT , &threshold },
        { 'f', "-f FILE", "Set the input file (empty for random).",  TYPE_STR , &file },
        { 'p', "-p P", "run the parallel fsytrf (only supported when RPM=N)", TYPE_BOOL , &par },
        { 't', "-t T", "number of virtual threads to drive the partition.", TYPE_INT , &t },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    rank=std::min(n,rank);

    typedef Givaro::ModularBalanced<double> Field;
    typedef Field::Element Element;

    Field F(q);
    Field::Element * A;
    FFLAS::Timer chrono;

    FFLAS::FFLAS_UPLO uplo = up?FFLAS::FflasUpper:FFLAS::FflasLower;
    double *time=new double[iter];
    for (size_t i=0; i<iter; i++) time[i]=0;
    for (size_t i=0;i<=iter;++i){
        if (!file.empty()){
            FFLAS::ReadMatrix (file.c_str(),F,n,n,A);
        }
        else {
            A = FFLAS::fflas_new<Element>(n*n);
            Field::RandIter G(F);
            if (grp){
                size_t * cols = FFLAS::fflas_new<size_t>(n);
                size_t * rows = FFLAS::fflas_new<size_t>(n);
                for (size_t i=0; i<n; ++i)
                    cols[i] = rows[i] = i;
                FFPACK::RandomSymmetricMatrixWithRankandRPM (F, n, rank, A, n, rows, cols, G);
                FFLAS::fflas_delete(cols);
                FFLAS::fflas_delete(rows);
            } else
                FFPACK::RandomSymmetricMatrixWithRankandRandomRPM (F, n, rank, A, n, G);
        }
        size_t*P=FFLAS::fflas_new<size_t>(n);
        if (rpm){
            chrono.clear();
            if (i) chrono.start();
            FFPACK::fsytrf_RPM (F, uplo, n, A, n, P, threshold);
            if (i) chrono.stop();
        }else{
            if (!par){
                chrono.clear();
                if (i) chrono.start();
                FFPACK::fsytrf (F, uplo, n, A, n, threshold);
                if (i) chrono.stop();
            }else{
                chrono.clear();
                if (i) chrono.start();
                PAR_BLOCK{
                    FFPACK::fsytrf (F, uplo, n, A, n, SPLITTER(t),threshold);
                }
                if (i) chrono.stop();
            }
        }
        FFLAS::fflas_delete(P);
        if (i) time[i-1]=chrono.realtime();
        FFLAS::fflas_delete( A);
    }
    std::sort(time, time+iter);
    double mediantime = time[iter/2];
    delete[] time;

    // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14
#define CUBE(x) ((x)*(x)*(x))
    double gfops = 1.0/3.0*CUBE(double(rank)/1000.0) +n/1000.0*n/1000.0*double(rank)/1000.0  - double(rank)/1000.0*double(rank)/1000.0*n/1000;

    std::cout << "Time: " << mediantime << " Gfops: " << gfops / mediantime;
    FFLAS::writeCommandString(std::cout, as) << std::endl;
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
