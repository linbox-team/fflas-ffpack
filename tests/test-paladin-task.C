/*
 * Copyright (C) the FFLAS-FFPACK group
 * Written by Ziad Sultan <ziad.sultan@imag.fr>
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

#undef  __FFLASFFPACK_USE_OPENMP
#define __FFLASFFPACK_USE_TBB

#include <string>
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/args-parser.h"



size_t add(const size_t x, const size_t y)
{
    return x+y;

}

size_t seq_fib(size_t n) {
    if (n < 2)
        return n;
    else
        return seq_fib(n-1) + seq_fib(n-2);
}

size_t par_fib(const size_t n, const size_t cutoff)
{

    size_t x=0, y=0, z=0;
    //	if (n < 2)
    //		return n;
    if (n < cutoff) // The bigger the cutoff the bigger is the parallel speed-up
        return seq_fib(n);
    else{
        SYNCH_GROUP(
                    TASK(MODE(READ(n) WRITE(x) CONSTREFERENCE(x)),
                         x = par_fib(n-1, cutoff);
                        );

                    TASK(MODE(READ(n) WRITE(y) CONSTREFERENCE(y)),
                         y = par_fib(n-2, cutoff);
                        );

                    CHECK_DEPENDENCIES;

                    TASK(MODE(READ(x,y) WRITE(z) CONSTREFERENCE(z,x,y)),
                         z=add(x,y);
                        );
                   );//end SYNCH_GROUP
        return z;
    }
}


int main(int argc, char** argv)
{

    size_t n = 20;
    bool p = true;
    size_t iters = 3;
    size_t cutoff = 2;
    //    int64_t q = 131071 ;
    //    int proc = MAX_THREADS;

    Argument as[] = {
        { 'n', "-n N", "Set the nth number of fibonacci to compute",      TYPE_INT , &n },
        { 'c', "-c N", "Set the Cutoff at which the sequential base case is called (the bigger the cuttof is the better is the parallel speed-up)",      TYPE_INT , &cutoff },
        { 'i', "-i N", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'p', "-p Y/N", "run the parallel program using Parallel(Y)/Sequential(N).", TYPE_BOOL , &p },
        END_OF_ARGUMENTS
    };
    FFLAS::parseArguments(argc,argv,as);
    //       { 't', "-t N", "Set number of processors.",            TYPE_INT , &proc },

    size_t f=0;


    // time
    FFLAS::Timer chrono;
    double *time=new double[iters];

    // parallel add using PARFOR1D
    for (size_t it=0;it<=iters;++it){
        chrono.clear();
        if (it) chrono.start();

        if(p){
            PAR_BLOCK{
                f=par_fib(n, cutoff);
            }// end of PAR_BLOCK
        }
        else
            f=seq_fib(n);
        if (it) {chrono.stop(); time[it-1]=chrono.realtime();}
    }

    std::sort(time, time+iters);
    double meantime = time[iters/2];
    delete[] time;

    // sequential add
    chrono.clear();
    chrono.start();
    size_t l=seq_fib(n);
    chrono.stop();
    double timeseq = chrono.realtime();


    // verification of the parallel result
    if (f!=l)
        std::cout<<"FAIL: Par_Fib("<<n<<") = "<<f<<" and Seq_Fib("<<n<<") = "<<l<<std::endl;
    else
        std::cout<<"PASS"<<std::endl;


    std::cout<<" n: "<<n;
    std::cout<<" SeqTime: "<<timeseq;
    std::cout<<" ParTime: " << meantime;

    std::string dataflow;

#ifdef __FFLASFFPACK_USE_DATAFLOW // OMP/KAAPI dataflow option
    dataflow  = " with dataflow synch!";
#else
    dataflow  = " with explicit synch!";
#endif

    std::cout<<dataflow<<std::endl;

    return 0;

}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
