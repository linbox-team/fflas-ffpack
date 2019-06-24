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


//#define __FFLASFFPACK_FORCE_SEQ

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iomanip>
#include <iostream>
#include <givaro/modular.h>
#include <givaro/udl.h>
#include <recint/rint.h>
#include <string>
#include <givaro/givintprime.h>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"

#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"

typedef Givaro::ModularBalanced<double> Field;

using namespace FFLAS;
template<class CutStrat, class StratParam>
bool tmain(int argc, char** argv, std::string printStrat)
{

    std::cerr << "tmain: " << printStrat << std::endl;


    size_t n = 2000;
    bool p = true;
    size_t iters = 3;
    int64_t q = 131071 ;
    bool dataPar = true;
    int proc = MAX_THREADS;
    uint64_t seed=getSeed();
    int strat = 1;

    Argument as[] = {
        { 'n', "-n N", "Set the dimension of the matrix.",      TYPE_INT , &n },
        { 'i', "-i N", "Set number of repetitions.",            TYPE_INT , &iters },
        { 't', "-t N", "Set number of processors.",            TYPE_INT , &proc },
        { 'u', "-u N", "Set the strategy parameter using t: 1 for (t, BLOCK, THREADS), 2 for (t, BLOCK, GRAIN), 3 for (t, BLOCK, FIXED), 4 for (t, ROW, THREADS), 5 for (t, ROW, GRAIN), 6 for (t, ROW, FIXED), 7 for (t, COLUMN, THREADS), 8 for (t, COLUMN, GRAIN), 9 for (t, COLUMN, FIXED), 10 for SINGLE strategy.",            TYPE_INT , &strat },
        { 'p', "-p Y/N", "run the parallel program using Parallel(Y)/Sequential(N).", TYPE_BOOL , &p },
        { 'd', "-d Y/N", "run the parallel program using data parallelism(Y)/task parallelism(N).", TYPE_BOOL , &dataPar },
        { 's', "-s seed", "Set seed for the random generator", TYPE_INT, &seed },
        END_OF_ARGUMENTS
    };
    parseArguments(argc,argv,as);

    size_t m = n; // matrices are square in this test

    Field F(q);
    Field::RandIter G(F,seed);

    // Allocate matrices
    typename Field::Element_ptr A = fflas_new (F, m, n);
    typename Field::Element_ptr B = fflas_new (F, m, n);
    typename Field::Element_ptr C = fflas_new (F, m, n);
    typename Field::Element_ptr Acop = fflas_new (F, m, n);


    auto CUTTER = SPLITTER(proc, CutStrat, StratParam);

    // initialize
    if(dataPar){
        PARFOR1D(i, m, CUTTER,
                 for (size_t j=0; j<(size_t)n; ++j)
                 G.random (*(A+i*n+j));
                );

        PARFOR1D(i, m, CUTTER,
                 for (size_t j=0; j<(size_t)n; ++j)
                 G.random (*(B+i*n+j));
                );

        PARFOR1D(i, m, CUTTER,
                 for (size_t j=0; j<(size_t)n; ++j)
                 G.random (*(C+i*n+j));
                );
    }
    else{ // initialize with tasks using FORBLOCK1D
        PAR_BLOCK{
            SYNCH_GROUP(
                        FORBLOCK1D(itt, m*n, CUTTER,
                                   TASK(MODE(WRITE(A)),
                                        for(size_t i=itt.begin(); i!=itt.end(); ++i)
                                        G.random (*(A+i)););

                                   TASK(MODE(WRITE(B)),
                                        for(size_t i=itt.begin(); i!=itt.end(); ++i)
                                        G.random (*(B+i)););

                                   TASK(MODE(WRITE(C)),
                                        for(size_t i=itt.begin(); i!=itt.end(); ++i)
                                        G.random (*(C+i)););
                                  );// end of FORBLOCK1D
                       );// end of SYNCH_GROUP
        }// end of PAR_BLOCK
    }

    // copy A for verification
    fassign(F,m,n,A,n,Acop,n);

    // time
    Timer chrono;
    double *time=new double[iters];

    // parallel add using PARFOR1D
    for (size_t it=0;it<=iters;++it){
        chrono.clear();
        if (it) chrono.start();

        if(dataPar){

            PARFOR1D(i, m, CUTTER,
                     for (size_t j=0; j<(size_t)n; ++j)
                     A[i*n+j] = B[i*n+j] + C[i*n+j];
                    );
        }
        else{
            PAR_BLOCK{
                FORBLOCK1D(itt, m*n, CUTTER,
                           TASK(MODE(READ(B,C) WRITE(A)),
                                for(size_t i=itt.begin(); i!=itt.end(); ++i)
                                A[i] = B[i] + C[i];
                               );
                          );
            }
        }


        if (it) {chrono.stop(); time[it-1]=chrono.realtime();}

    }
    std::sort(time, time+iters);
    double meantime = time[iters/2];
    delete[] time;

    // sequential add
    chrono.clear();
    chrono.start();
    for(size_t i=0; i<m*n; ++i)
        Acop[i]=B[i]+C[i];
    chrono.stop();
    double timeseq = chrono.usertime();


    // verification of the parallel result
    bool fail = false;
    for(size_t i=0; i<m; ++i)
        for (size_t j=0; j<n; ++j)
            if (!F.areEqual (*(Acop+i*n+j), *(A+i*n+j))){
                std::cout << " Seq["<<i<<","<<j<<"] = " << (*(Acop+i*n+j))
                << " Par["<<i<<","<<j<<"] = " << (*(A+i*n+j))
                << std::endl;
                fail=true;
            }

    if (fail)
        std::cout<<"FAIL"<<std::endl;
    else
        std::cout<<"PASS"<<std::endl;


    std::cout<<"m: "<<m<<" n: "<<n;
    std::cout<<" SeqTime: "<<timeseq;
    std::cout<<" ParTime: " << meantime;

    //       std::cout<<" Strategy:("<<proc<<", "<<CutStrat<<", "<<StratParam<<")";
    std::cout<<" Strategy:("<<proc<<", "<<printStrat<<")";

    std::string dataflow;


#ifdef __FFLASFFPACK_USE_DATAFLOW // OMP/KAAPI dataflow option
    dataflow  = " with dataflow synch!";
#else
    dataflow  = " with explicit synch!";
#endif

    if(dataPar)
        std::cout<<" Data parallelism is used!"<<std::endl;
    else
        std::cout<<" TASK parallelism is used"<<dataflow<<std::endl;

    fflas_delete(A);
    fflas_delete(Acop);
    fflas_delete(B);
    fflas_delete(C);

    return fail;

}



int main(int argc, char** argv)
{

    size_t n = 2000;
    bool p = true;
    size_t iters = 3;
    bool dataPar = true;
    int proc = MAX_THREADS;

    int strat = 1;

    Argument as[] = {
        { 'n', "-n N", "Set the dimension of the matrix.",      TYPE_INT , &n },
        { 'i', "-i N", "Set number of repetitions.",            TYPE_INT , &iters },
        { 't', "-t N", "Set number of processors.",            TYPE_INT , &proc },
        { 's', "-s N", "Set the strategy parameter using t: 1 for (t, BLOCK, THREADS), 2 for (t, BLOCK, GRAIN), 3 for (t, BLOCK, FIXED), 4 for (t, ROW, THREADS), 5 for (t, ROW, GRAIN), 6 for (t, ROW, FIXED), 7 for (t, COLUMN, THREADS), 8 for (t, COLUMN, GRAIN), 9 for (t, COLUMN, FIXED), 10 for SINGLE strategy.",            TYPE_INT , &strat },
        { 'p', "-p Y/N", "run the parallel program using Parallel(Y)/Sequential(N).", TYPE_BOOL , &p },
        { 'd', "-d Y/N", "run the parallel program using data parallelism(Y)/task parallelism(N).", TYPE_BOOL , &dataPar },
        END_OF_ARGUMENTS
    };
    parseArguments(argc,argv,as);



    bool fail = false;

    switch (strat){
    case 1: fail |= tmain<CuttingStrategy::Block,StrategyParameter::Threads>(argc,argv,std::string("BLOCK, THREADS"));
    case 2: fail |= tmain<CuttingStrategy::Block,StrategyParameter::Grain>(argc,argv,std::string("BLOCK, GRAIN"));
    case 3: fail |= tmain<CuttingStrategy::Block,StrategyParameter::Fixed>(argc,argv,std::string("BLOCK, FIXED"));
    case 4: fail |= tmain<CuttingStrategy::Row,StrategyParameter::Threads>(argc,argv,std::string("ROW, THREADS"));
    case 5: fail |= tmain<CuttingStrategy::Row,StrategyParameter::Grain>(argc,argv,std::string("ROW, GRAIN"));
    case 6: fail |= tmain<CuttingStrategy::Row,StrategyParameter::Fixed>(argc,argv,std::string("ROW, FIXED"));
    case 7: fail |= tmain<CuttingStrategy::Column,StrategyParameter::Threads>(argc,argv,std::string("COLUMN, THREADS"));
    case 8: fail |= tmain<CuttingStrategy::Column,StrategyParameter::Grain>(argc,argv,std::string("COLUMN, GRAIN"));
    case 9: fail |= tmain<CuttingStrategy::Column,StrategyParameter::Fixed>(argc,argv,std::string("COLUMN, FIXED"));
    case 10: fail |= tmain<CuttingStrategy::Single,StrategyParameter::Threads>(argc,argv,std::string("SINGLE, THREADS"));
    }

    return fail;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
