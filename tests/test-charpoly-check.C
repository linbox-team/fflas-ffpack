/*
 * Copyright (C) 2015 the FFLAS-FFPACK group
 * Written by Ashley Lesdalons <Ashley.Lesdalons@e.ujf-grenoble.fr>
 * Jean-Guillaume Dumas <Jean-Guillaume.Dumas@univ-grenoble-alpes.fr>
 *
 * This file is Free Software and part of FFLAS-FFPACK.
 *
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
 *
 */

//--------------------------------------------------------------------------
//          Test for Checker_charpoly
//--------------------------------------------------------------------------

#define ENABLE_CHECKER_charpoly 1
#define TIME_CHECKER_CHARPOLY 1


#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/utils/fflas_io.h"


template <class Field, class Polynomial>
void printPolynomial (const Field &F, Polynomial &v)
{
    for (int i = v.size() - 1; i >= 0; i--) {
        F.write (std::cout, v[i]);
        if (i > 0)
            std::cout << " x^" << i << " + ";
    }
    std::cout << std::endl;
}
using namespace FFLAS;
int main(int argc, char** argv) {
    typedef Givaro::ModularBalanced<double> Field;
    Givaro::Integer q = 131071;
    size_t iter = 3;
    size_t MAXN = 100;
    size_t n = 0, N = 0;
    int variant = 6;
    int seed = getSeed();

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).", TYPE_INTEGER , &q },
        { 'i', "-i R", "Set number of repetitions.", TYPE_INT , &iter },
        { 'n', "-n N", "Set the size of the matrix.", TYPE_INT , &n },
        { 's', "-s N", "Set the seed                 .", TYPE_UINT64 , &seed },
        { 'a', "-a algorithm", "Set the algorithmic variant", TYPE_INT, &variant },
        END_OF_ARGUMENTS
    };
    FFLAS::parseArguments(argc,argv,as);
    FFPACK::FFPACK_CHARPOLY_TAG CPalg;
    switch (variant){
    case 0: CPalg = FFPACK::FfpackLUK; break;
    case 1: CPalg = FFPACK::FfpackKG; break;
    case 2: CPalg = FFPACK::FfpackDanilevski; break;
    case 3: CPalg = FFPACK::FfpackKGFast; break;
    case 4: CPalg = FFPACK::FfpackKGFastG; break;
    case 5: CPalg = FFPACK::FfpackHybrid; break;
    case 6: CPalg = FFPACK::FfpackArithProg; break;
    default: CPalg = FFPACK::FfpackLUK; break;
    }

    Field F(q);
    srand (seed);
    Field::RandIter Rand(F,seed);
    typedef Givaro::Poly1Dom<Field> PolRing;
    typedef PolRing::Element Polynomial;
    PolRing R(F);
    size_t pass = 0;
    for (size_t i=0; i<iter; ++i) {

        N = n?n: rand() % MAXN + 1;
        // 		std::cout << "n= " << n << "\n";
        Field::Element_ptr A = FFLAS::fflas_new(F,N,N);

        Polynomial g(n);

        PAR_BLOCK { FFLAS::pfrand(F,Rand,N,N,A,N/MAX_THREADS); }
        try {
            //             FFLAS::WriteMatrix (std::cerr<<"A=",F,N,N,A,N,FflasMaple) <<std::endl;
            //             FFPACK::Checker_charpoly<Field,Polynomial> checker(F,n,A);
            Givaro::Timer charpolytime; charpolytime.start();
            FFPACK::CharPoly (R,g,N,A,N,CPalg);
            charpolytime.stop();
            std::cerr << "CHARPol checked full: " << charpolytime << std::endl;
            //             printPolynomial(F,g);
            //             checker.check(g);
            std::cout << N << 'x' << N << " charpoly verification successful\n";
            pass++;
        } catch(FailureCharpolyCheck &e) {
            std::cout << N << 'x' << N << " charpoly verification failed!\n";
        }
        FFLAS::fflas_delete( A);

    }

    std::cout << pass << "/" << iter << " tests were successful: ";
    FFLAS::writeCommandString(std::cout, as) << std::endl;

    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
