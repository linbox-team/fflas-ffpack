/*
 * Copyright (C) 2020 the FFLAS-FFPACK group
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@univ-grenoble-alpes.fr>
 *
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
 *.
 */


#include "fflas-ffpack/fflas-ffpack-config.h"

#include <iostream>
#include <givaro/modular-balanced.h>

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/args-parser.h"

using namespace std;
using namespace FFLAS;


int main(int argc, char** argv) {

    size_t iter = 3 ;
    Givaro::Integer q = 131071 ;
    size_t k = 2000 ;
    size_t n = 2000 ;
    int t=MAX_THREADS;
    int NBK = -1;
    bool up =true;

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'k', "-k K", "Set the col dimension of A.",      TYPE_INT , &k },
        { 'n', "-n N", "Set the col dimension of B.",      TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.",       TYPE_INT , &iter },
        { 't', "-t T", "number of virtual threads to drive the partition.", TYPE_INT , &t },
        { 'b', "-b B", "number of numa blocks per dimension for the numa placement", TYPE_INT , &NBK },
        { 'u', "-u yes/no", "Updates an upper/lower triangular matrix.",  TYPE_BOOL , &up },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);

    FFLAS_UPLO uplo = up?FflasUpper:FflasLower;

    if (NBK==-1) NBK = t;
    //  typedef Givaro::Modular<Givaro::Integer> Field;
    //  typedef Givaro::Modular<int64_t> Field;
    typedef Givaro::Modular<double> Field;
    //  typedef Givaro::Modular<float> Field;
    //  typedef Givaro::ModularBalanced<float> Field;
    //  typedef Givaro::ModularBalanced<double> Field;
    //  typedef Givaro::ModularBalanced<int64_t> Field;
    //  typedef Givaro::Modular<Givaro::Integer> Field;
    typedef Field::Element Element;

    Field F(q);

    Timer chrono, TimFreivalds;
    double time=0.0;

    Element * A, * B, * C, * D;

    Field::RandIter G(F);
    A = fflas_new(F,n,k,Alignment::CACHE_PAGESIZE);
    //#pragma omp parallel for collapse(2) schedule(runtime)
    PAR_BLOCK { pfrand(F,G, n,k,A,n/size_t(NBK)); }

    B = fflas_new(F,k,n,Alignment::CACHE_PAGESIZE);
    //#pragma omp parallel for collapse(2) schedule(runtime)
    PAR_BLOCK { pfrand(F,G, k,n,B,k/NBK); }

    C = fflas_new(F,n,n,Alignment::CACHE_PAGESIZE);
    //#pragma omp parallel for collapse(2) schedule(runtime)
    PAR_BLOCK { pfzero(F, n,n,C,n/NBK); }

    D = fflas_new(F,n,n,Alignment::CACHE_PAGESIZE);
    //#pragma omp parallel for collapse(2) schedule(runtime)
    PAR_BLOCK { pfzero(F, n,n,D,n/NBK); }

    for (size_t i=0;i<=iter;++i){

        chrono.clear();
        if (i) chrono.start();

        fsyr2k(F, uplo, FflasNoTrans, n,k, F.one, A, k, B, n, F.zero, C,n);

        if (i) chrono.stop();
        time+=chrono.usertime();

    }
    fflas_delete( A);
    fflas_delete( B);
    fflas_delete( C);

    // -----------
    std::cout << "Time: " << time / double(iter)
              << " Gfops: " << 2.*(double(n)/1000.)*(double(n)/1000.)*(double(k)/1000.)/ time * double(iter);
    FFLAS::writeCommandString(std::cout, as) << std::endl;
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
