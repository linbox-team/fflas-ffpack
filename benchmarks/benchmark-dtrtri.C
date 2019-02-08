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

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>

#define __FFLASFFPACK_HAVE_DTRTRI 1

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/args-parser.h"

#ifdef __FFLASFFPACK_USE_OPENMP
typedef FFLAS::OMPTimer TTimer;
#else
typedef FFLAS::Timer TTimer;
#endif


using namespace std;

int main(int argc, char** argv) {

    size_t iter = 1;
    int    q    = 1009;
    size_t    n    = 2000;
    std::string file = "";

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
        { 'n', "-n N", "Set the dimension of the matrix.",               TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &iter },
        { 'f', "-f FILE", "Set the input file (empty for random).",  TYPE_STR , &file },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    typedef Givaro::Modular<double> Field;
    typedef Field::Element Element;

    Field F(q);
    Element * A;

    TTimer chrono;
    double time=0.0;

    Field::RandIter G(F);
    for (size_t i=0;i<iter;++i){
        if (!file.empty()){
            FFLAS::ReadMatrix (file.c_str(),F,n,n,A);
        } else {
            A = FFLAS::fflas_new<Element>(n*n);
            for (size_t j=0; j<(size_t) n*n; ++j)
                G.random(*(A+j));
            for (size_t k=0;k<(size_t)n;++k)
                while (F.isZero( G.random(*(A+k*(n+1)))));
        }

        chrono.clear();
        chrono.start();
        clapack_dtrtri(CblasRowMajor,CblasUpper, CblasNonUnit,n,A,n);
        chrono.stop();

        time+=chrono.usertime();
        FFLAS::fflas_delete( A);

    }

    // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14
    std::cout << "Time: " << time / double(iter)
    << " Gfops: " << (2.*double(n)/1000.*double(n)/1000.*double(n)/1000.0) / time * double(iter) / 3.;
    FFLAS::writeCommandString(std::cout, as) << std::endl;

    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
