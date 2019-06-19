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
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/args-parser.h"

using namespace std;

int main(int argc, char** argv) {
 
#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t iter = 3;
    int    q    = 131071;
    size_t    m    = 2000 ;
    size_t    n    = 1000;
    std::string file1 = "";
    std::string file2 = "";
    bool v=false;

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
        { 'm', "-m M", "Set the row dimension of the RHS matrix.", TYPE_INT , &m },
        { 'n', "-n N", "Set the col dimension of the RHS matrix.", TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.",               TYPE_INT , &iter },
        { 'f', "-f FILE", "Set the first input file (empty for random).", TYPE_STR, &file1 },
        { 'g', "-g FILE", "Set the second input file (empty for random).", TYPE_STR, &file2 },
        { 'v', "-v V", "Whether to run the verification of the solution.", TYPE_BOOL, &v },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    typedef Givaro::ModularBalanced<double> Field;

    Field F(q);
    Field::Element_ptr  A, Ac, B, Bc;

    FFLAS::Timer chrono;
    double time=0.0;
    Field::RandIter G(F);

    if (!file1.empty()){
        FFLAS::ReadMatrix (file1.c_str(),F,m,m,A);
    }
    else{
        A = FFLAS::fflas_new (F,m,m,Alignment::CACHE_PAGESIZE);
        FFPACK::RandomMatrixWithRank (F, m,m, m, A, m);
    }

    if (!file2.empty()){
        FFLAS::ReadMatrix (file2.c_str(),F,m,n,B);
    }
    else{
        B = FFLAS::fflas_new(F,m,n,Alignment::CACHE_PAGESIZE);
        FFPACK::RandomMatrix(F, m, n, B, n, G);
    }

    Bc = FFLAS::fflas_new(F,m,n);
    Ac = FFLAS::fflas_new (F,m,m,Alignment::CACHE_PAGESIZE);

    for (size_t i=0;i<=iter;++i){
        FFLAS::fassign (F,m,m,A,m,Ac,m);
        FFLAS::fassign (F,m,n,B,n,Bc,n);
        int info;
        chrono.clear();
        if (i) chrono.start();
        FFPACK::fgesv (F, FFLAS::FflasLeft, m,n, Ac, m, Bc, n, &info);
        if (i) {chrono.stop(); time+=chrono.realtime();}
    }

    // -----------
    // Standard output for benchmark
    std::cout << "Time: " << time / double(iter)
    << " Gfops: " << (double(m*m)/1000000.*(2*double(m)/3000.+2*double(n)/1000.0)) / time * double(iter);
    FFLAS::writeCommandString(std::cout, as) << std::endl;

    if (v){
        FFLAS::fgemm(F,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,m,n,m,F.one,A,m,Bc,n,F.mOne,B,n);
        if (!FFLAS::fiszero(F,m,n,B,n)){
            std::cerr<<"Verification FAILED"<<std::endl;
            // FFLAS::WriteMatrix(std::cerr<<"A= "<<std::endl,F,m,m,A,m);
            // FFLAS::WriteMatrix(std::cerr<<"X= "<<std::endl,F,m,n,Bc,n);
            // FFLAS::WriteMatrix(std::cerr<<"B= "<<std::endl,F,m,n,B,n);
            FFLAS::fflas_delete( A);
            FFLAS::fflas_delete( B);
            return -1;
        }
    }
    FFLAS::fflas_delete( A);
    FFLAS::fflas_delete( B);


    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
