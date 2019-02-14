/* Copyright (c) 2012 FFLAS-FFPACK
 * Written by J.G. Dumas <jgdumas@imag.fr>
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
#include <fstream>
#include <givaro/modular.h>

#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/args-parser.h"

#define CUBE(x) ((x)*(x)*(x))

template<class Field>
void launch_wino(const Field  &F,
                 const size_t &n,
                 const size_t &NB,
                 const size_t &wino,
                 const bool   &asmax,
                 const size_t &seed,
                 const bool   compare)
{

    typedef typename Field::Element Element ;
    typename Field::RandIter G(F);

    if (compare)
        F.write(std::cout << "Field ") << std::endl;

    double basetime(0.0), time(0.0);

    Element *A, *C;
    A = FFLAS::fflas_new<Element>(n*n);
    C = FFLAS::fflas_new<Element>(n*n);
    for (size_t i=0; i<n*n;++i)
        G.random(A[i]);

    // ----- Compare with fgemm
    FFLAS::Timer chrono;
    if (compare) {
        for(size_t i=0; i<NB; ++i) {
            chrono.start();
            FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                         n,n,n, F.one,
                         A, n, A, n, F.zero, C,n);
            chrono.stop();
            basetime+= chrono.usertime();
        }

        std::cout << "Time: " << basetime / double(NB)
        << " Gfops: " << 2. * CUBE(double(n)/1000.0) / basetime * double(NB)
        << " [fgemm result]" << std::endl;
    }

    // ----- Winograd
    for(size_t w = (asmax)? 0 : wino; w <= wino; ++w) {
        FFLAS::MMHelper<Field, FFLAS::MMHelperAlgo::Winograd> WH (F,(int)w);
        time = 0. ;
        chrono.clear();
        for(size_t i=0; i<NB; ++i) {
            chrono.start();
            FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                         n, n, n, F.one, A, n, A, n, F.zero, C, n, WH);
            chrono.stop();
            time+= chrono.usertime();
        }

        // -----------
        // Standard output for benchmark - Alexis Breust 2014/11/14
        std::cout << "Time: " << time / double(NB)
        << " Gfops: " << 2. * CUBE(double(n)/1000.0) / time * double(NB);

        if (compare || asmax)
            std::cout << " [wino" << w << " result]" << std::endl;
    }

    if (compare)
        std::cout << std::endl;

    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(C);
}

int main (int argc, char ** argv) {

    size_t iter = 1;
    int q  = 1009;
    int n       = 1000;
    int w       = 7;
    size_t seed    = 0;
    bool compare = false;
    bool balanced = false;
    std::string type = "double";
    bool levelasmax = false;

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
        { 'n', "-n N", "Set the dimension of the matrix.",               TYPE_INT , &n },
        { 'w', "-w N", "Set the winograd level.",                        TYPE_INT , &w },
        { 'l', "-l {YN}", "Use -w info a max (Yes or No).",              TYPE_BOOL , &levelasmax },
        { 's', "-s S", "Set the seed for randomness (0 for random).",   TYPE_INT , &seed },
        { 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &iter },
        { 'c', "-c {YN}", "Compare mode, overrides -b and -t options (Yes or No).", TYPE_BOOL , &compare },
        { 'b', "-b {YN}", "Use balanced modular (Yes or No).",           TYPE_BOOL , &balanced },
        { 't', "-t TYPE", "Set the field type (double/float/int).",      TYPE_STR , &type },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    if (!seed)
        seed = FFLAS::BaseTimer::seed();
    srand((uint32_t)seed);

    if (compare) {
        Givaro::Modular<double> F1(q);
        Givaro::Modular<float>  F2(q);
        Givaro::Modular<int>    F3(q);
        Givaro::ModularBalanced<double> F4(q);
        Givaro::ModularBalanced<float>  F5(q);
        Givaro::ModularBalanced<int>    F6(q);
        // ZZ<double> F7;
        // ZZ<float>  F8;
        // ZZ<int>    F9;

        launch_wino(F1,n,iter,w,levelasmax,seed,true);
        launch_wino(F2,n,iter,w,levelasmax,seed,true);
        launch_wino(F3,n,iter,w,levelasmax,seed,true);
        launch_wino(F4,n,iter,w,levelasmax,seed,true);
        launch_wino(F5,n,iter,w,levelasmax,seed,true);
        launch_wino(F6,n,iter,w,levelasmax,seed,true);
        // launch_wino(F7,n,iter,winomax,seed);
        // launch_wino(F8,n,iter,winomax,seed);
        // launch_wino(F9,n,iter,winomax,seed);
    }
    else {
        if (balanced) {
            if (type == "double")     launch_wino(Givaro::ModularBalanced<double>(q),n,iter,w,levelasmax,seed,false);
            else if (type == "float") launch_wino(Givaro::ModularBalanced<float>(q),n,iter,w,levelasmax,seed,false);
            else if (type == "int")   launch_wino(Givaro::ModularBalanced<int>(q),n,iter,w,levelasmax,seed,false);
        }
        else {
            if (type == "double")     launch_wino(Givaro::Modular<double>(q),n,iter,w,levelasmax,seed,false);
            else if (type == "float") launch_wino(Givaro::Modular<float>(q),n,iter,w,levelasmax,seed,false);
            else if (type == "int")   launch_wino(Givaro::Modular<int>(q),n,iter,w,levelasmax,seed,false);
        }
    }

    if (compare || levelasmax)
        std::cout << "Lauch with:";
    FFLAS::writeCommandString(std::cout, as) << std::endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
