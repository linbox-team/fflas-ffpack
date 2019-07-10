/* Copyright (c) FFLAS-FFPACK
 * Written by Clement Pernet <clement.pernet@imag.fr>
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
//#define __FFLASFFPACK_ARITHPROG_PROFILING

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
using namespace FFPACK;

template<class Field>
void run_with_field(int q, uint64_t bits, size_t n, size_t d, size_t iter, std::string file, int variant, uint64_t seed){
    Field F(q);
    typedef typename Field::Element Element;
    FFPACK::FFPACK_CHARPOLY_TAG CT;
    switch (variant){
    case 0: CT = FfpackAuto; break;
    case 1: CT = FfpackDanilevski; break;
    case 2: CT = FfpackLUK; break;
    case 3: CT = FfpackArithProgKrylovPrecond; break;
    case 4: CT = FfpackArithProg; break;
    case 5: CT = FfpackKG; break;
    case 6: CT = FfpackKGFast; break;
    case 7: CT = FfpackHybrid; break;
    case 8: CT = FfpackKGFastG; break;
    default: CT = FfpackAuto; break;
    }
    FFLAS::Timer chrono;
    Element *A;
    double time_charp=0;
    for (size_t i=0;i<iter;++i){
        if (!file.empty()){
            FFLAS::ReadMatrix (file, F, n, n, A);
        }
        else{
            A = FFLAS::fflas_new (F, n, n);
            Givaro::Integer samplesize(1); samplesize <<= bits;
            typename Field::RandIter G (F, seed, samplesize);
            FFPACK::RandomMatrix (F, n, n, A, n, G);
        }
        typename Givaro::Poly1Dom<Field>::Element cpol(n+1);
        typename Givaro::Poly1Dom<Field> R(F);
        chrono.clear();
        chrono.start();
        FFPACK::CharPoly (R, cpol, n, A, n, CT, d);
        chrono.stop();

        time_charp+=chrono.usertime();

        FFLAS::fflas_delete( A);
    }
    // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14
    std::cout << "Time: " << time_charp / double(iter)
    << " Gfops: " << (2.*double(n)/1000.*double(n)/1000.*double(n)/1000.0) / time_charp * double(iter);

}

int main(int argc, char** argv) {

#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t iter = 1;
    int    q    = 131071;
    uint64_t bits = 10;
    size_t    n    = 1000;
    size_t    d    = __FFLASFFPACK_ARITHPROG_THRESHOLD;
    std::string file = "";
    uint64_t seed = FFLAS::getSeed();
    int variant = 0;

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for the ring ZZ).",  TYPE_INT , &q },
        { 'b', "-b B", "Set the bitsize of the random elements.",         TYPE_INT , &bits},
        { 'n', "-n N", "Set the dimension of the matrix.",               TYPE_INT , &n },
        { 'd', "-d D", "Set the degree of the preconditionner (for ArithProg variant).",  TYPE_INT , &d },
        { 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &iter },
        { 'f', "-f FILE", "Set the input file (empty for random).",  TYPE_STR , &file },
        { 'a', "-a algorithm", "Set the algorithmic variant", TYPE_INT, &variant },
        { 's', "-s S", "Sets seed.", TYPE_INT , &seed },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    if (q > 0){
        bits = Givaro::Integer(q).bitsize();
        run_with_field<Givaro::ModularBalanced<double> >(q, bits, n , d, iter, file, variant,seed);
    } else
        run_with_field<Givaro::ZRing<Givaro::Integer> > (q, bits, n , d, iter, file, variant,seed);

    FFLAS::writeCommandString(std::cout, as) << std::endl;
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
