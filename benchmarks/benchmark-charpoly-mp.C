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
#define  __FFLASFFPACK_FORCE_SEQ

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/args-parser.h"


using namespace std;

int main(int argc, char** argv) {

    size_t iter = 1;
    size_t    n    = 100;
    std::string file = "";
    static int variant =0;
    uint64_t b = 150;
    uint64_t seed = FFLAS::getSeed();

    Argument as[] = {
        { 'b', "-b B", "Set the bitsize of the random characteristic.",  TYPE_INT , &b },
        { 'n', "-n N", "Set the dimension of the matrix.",               TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &iter },
        { 'f', "-f FILE", "Set the input file (empty for random).",  TYPE_STR , &file },
        { 'a', "-a algorithm", "Set the algorithmic variant", TYPE_INT, &variant },
        { 's', "-s S", "Sets seed.", TYPE_INT , &seed },

        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);
    typedef Givaro::ZRing<Givaro::Integer> Field;
    FFPACK::FFPACK_CHARPOLY_TAG CT;
    switch (variant){
    case 0: CT = FFPACK::FfpackAuto; break;
    case 1: CT = FFPACK::FfpackLUK; break;
    case 2: CT = FFPACK::FfpackDanilevski; break;
    case 3: CT = FFPACK::FfpackArithProg; break;
    case 4: CT = FFPACK::FfpackKG; break;
    case 5: CT = FFPACK::FfpackKGFast; break;
    case 6: CT = FFPACK::FfpackHybrid; break;
    case 7: CT = FFPACK::FfpackKGFastG; break;
    default: CT = FFPACK::FfpackAuto; break;
    }
    typedef Field::Element Element;

    Field F;
    FFLAS::Timer chrono;
    double time=0.0;

    Element *A;
    uint64_t bs=1;
    uint64_t size=b;
    for (size_t i=0;i<iter;++i){

        if (!file.empty()){
            FFLAS::ReadMatrix (file, F, n, n, A);
        }
        else{
            A = FFLAS::fflas_new<Element>(n*n);
//            typename Field::Residu_t samplesize(1); samplesize <<= size;
            Field::RandIter G(F,seed);
            G.setBitsize(size);
            for (size_t j=0; j< (size_t)n*n; ++j)
                G.random(*(A+j));
        }

        typedef Givaro::Poly1Dom<Field> PolRing;
        PolRing R(F);
        PolRing::Element cpol;
        chrono.clear();
        chrono.start();
        FFPACK::CharPoly (R, cpol, n, A, n, CT);
        chrono.stop();

        time+=chrono.usertime();

        bs = FFLAS::bitsize (F,n,n,A,n);
        FFLAS::fflas_delete( A);
    }

    // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14
    std::cout <<"Time: " << time / double(iter)
    << " Gfops: " << (2.*double(n)/1000.*double(n)/1000.*double(n)/1000.0) / time * double(iter)
    << " bitsize: " << bs;
    FFLAS::writeCommandString(std::cout, as) << std::endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
