/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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
 *.
 */
#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "givaro/modular-integer.h"

int main(int argc, char** argv){
    srand((int)time(NULL));
    srand48(time(NULL));

    static size_t iters = 3 ;
    static Givaro::Integer q = -1 ;
    static unsigned long b = 512 ;
    static size_t m = 512 ;
    static size_t n = 512 ;
    static Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the random characteristic.",         TYPE_INT , &b },
        { 'm', "-m M", "Set the dimension m of the matrix.",                    TYPE_INT , &m },
        { 'n', "-n N", "Set the dimension n of the matrix.",                    TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.",                            TYPE_INT , &iters },
        END_OF_ARGUMENTS
    };
    FFLAS::parseArguments(argc,argv,as);

    size_t seed= time(NULL);
    typedef Givaro::Modular<Givaro::Integer> Field;
    FFLAS::Timer chrono;
    double time=0.;
    Givaro::Integer p;
    Givaro::IntPrimeDom IPD;

    for (size_t i=0;i<iters;i++) {

        Givaro::Integer::random_exact_2exp(p, b);
        IPD.prevprimein(p);
        Field F(p);
        size_t lda;
        lda=n;

        typename Field::RandIter Rand(F,seed);
        Field::Element_ptr A;
        A= FFLAS::fflas_new(F,m,lda);
        size_t * P = FFLAS::fflas_new<size_t>(n) ;
        size_t * Q = FFLAS::fflas_new<size_t>(m) ;

        for (size_t ii=0;ii<m*lda;++ii)
            Rand.random(A[ii]);

        Givaro::Integer alpha;
        alpha=1;
        chrono.clear();chrono.start();
        FFPACK::LUdivine (F, FFLAS::FflasUnit, FFLAS::FflasNoTrans, m, n, A, lda, P, Q);
        chrono.stop();
        time+=chrono.usertime();

        FFLAS::fflas_delete(A);
        FFLAS::fflas_delete(P);
        FFLAS::fflas_delete(Q);
    }
    double Gflops=(2./3.*double(m)/1000.*double(m)/1000.*double(n)/1000.0) / chrono.usertime() * double(iters);
    Gflops*=p.bitsize()/16.;
    cout<<"Time: "<<time/iters<<"  Gfops: "<<Gflops<<endl;


    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
