/* Copyright (c) FFLAS-FFPACK 2017
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

int main(int argc, char** argv) {

 #ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t iter = 10;
    int    q    = 131071;
    size_t    n    = 4000;
    bool u = true;  // Upper triangular
    bool d = false; // NonUnit diagonal
    bool t = false; // NoTrans
    bool v = true; // verification

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
        { 'n', "-n N", "Set the col dimension of the RHS matrix.", TYPE_INT , &n },
        { 'u', "-u U", "Upper/Lower triangular.", TYPE_BOOL , &u },
        { 'd', "-d D", "Unit/NonUnit diagonal.", TYPE_BOOL , &d },
        { 't', "-t T", "Trans/NoTrans.", TYPE_BOOL , &t },
        { 'i', "-i R", "Set number of repetitions.",               TYPE_INT , &iter },
        { 'v', "-v V", "Whether to check for correctness (probabilistic).", TYPE_BOOL , &v },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);

    FFLAS_UPLO UpLo = u? FflasUpper : FflasLower;
    FFLAS_TRANSPOSE Trans = t? FflasTrans : FflasNoTrans;
    FFLAS_DIAG Diag = d? FflasUnit : FflasNonUnit;

    typedef Givaro::ModularBalanced<double> Field;
    typedef Field::Element Element;

    Field F(q);
    Element * A;
    Element * b, *c=NULL;

    FFLAS::Timer chrono;
    double time=0.0;
    Field::RandIter G(F);
    Field::Element proj;
    F.init(proj);

    A = fflas_new (F,n,n,Alignment::CACHE_PAGESIZE);

    FFPACK::RandomTriangularMatrix (F, n, n, UpLo, Diag, true, A, n, G);

    b = fflas_new(F,n,1,Alignment::CACHE_PAGESIZE);
    for (size_t i=0;i<=iter;++i){
        chrono.clear();
        frand (F,G,n,1,b,1);
        if (v){
            c = fflas_new(F,n,1,Alignment::CACHE_PAGESIZE);
            frand (F,G,n,1,c,1);
            // proj <- c^T . b
            proj = fdot (F, n, c, 1, b, 1);
            // cerr<<"proj = "<<proj<<endl;
        }
        if (v){
            // c <- c U
            ftrmm (F, FflasRight, UpLo, Trans, Diag,
                   1,n, F.one, A, n, c, n);

        }
        chrono.start();
        // b <- U^-1 b
        ftrsv (F, UpLo, Trans, Diag, n, A, n, b, 1);
        chrono.stop();
        if (i) { time+=chrono.usertime();}

        if (v){
            // check b.c == proj
            Element verif = fdot (F, n, c, 1, b, 1);
            if (!F.areEqual (proj,verif)){
                cerr<<"FAIL: proj = "<<proj<<" verif = "<<verif<<endl;
                fflas_delete (A);
                fflas_delete (b);
                fflas_delete (c);
                return -1;
            }
        }
    }
    if (v) fflas_delete(c);
    fflas_delete (A);
    fflas_delete (b);

    std::cout << "Time: " << time / double(iter)
    << " Gfops: " << (double(n)/1000.*double(n)/1000000.0) / time * double(iter);
    writeCommandString(std::cout, as) << std::endl;

    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
