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

// declare that the call to openblas_set_numthread will be made here, hence don't do it
// everywhere in the call stack
#define __FFLASFFPACK_OPENBLAS_NT_ALREADY_SET 1

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/args-parser.h"

using namespace std;

int main(int argc, char** argv) {

#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t iter = 3;
    int    q    = 1009;
    size_t    m    = 2000 ;
    size_t    n    = 2000;
    std::string file1 = "";
    std::string file2 = "";
    int t=MAX_THREADS;
    int NBK = -1;
    int p = 0; // 0 for sequential 1 for pIter-sRec ; 2 for pRec; 3 for hybrid

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
        { 'm', "-m M", "Set the row dimension of the RHS matrix.", TYPE_INT , &m },
        { 'n', "-n N", "Set the col dimension of the RHS matrix.", TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.",               TYPE_INT , &iter },
        { 'f', "-f FILE", "Set the first input file (empty for random).", TYPE_STR, &file1 },
        { 'g', "-g FILE", "Set the second input file (empty for random).", TYPE_STR, &file2 },
        { 't', "-t T", "number of virtual threads to drive the partition.", TYPE_INT, &t },
        { 'b', "-b B", "number of numa blocks per dimension for the numa placement", TYPE_INT, &NBK },
        { 'p', "-p P", "0 for sequential, 1 for Iterative, 2 for Recursive, 3 for Hybrid.", TYPE_INT , &p },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);
    if (NBK==-1) NBK = t;

    typedef Givaro::ModularBalanced<double> Field;
    typedef Field::Element Element;

    Field F(q);
    Element * A;
    Element * B;

    FFLAS::Timer chrono;
    double time=0.0;
    Field::RandIter G(F);

    if (!file1.empty()){
        FFLAS::ReadMatrix (file1.c_str(),F,m,m,A);
    }
    else{
        A = FFLAS::fflas_new (F,m,m,Alignment::CACHE_PAGESIZE);
        PAR_BLOCK{ FFLAS::pfrand(F,G,m,m,A,m/NBK); }

        for (size_t k=0;k<(size_t)m;++k)
            while (F.isZero( G.random(*(A+k*(m+1)))));
    }

    if (!file2.empty()){
        FFLAS::ReadMatrix (file2.c_str(),F,m,n,B);
    }
    else{
        B = FFLAS::fflas_new(F,m,n,Alignment::CACHE_PAGESIZE);
        PAR_BLOCK{ FFLAS::pfrand(F,G,m,n,B,m/NBK); }
    }
    //}

for (size_t i=0;i<=iter;++i){
    chrono.clear();
    if (i) chrono.start();

    if (!p){
        FFLAS::ParSeqHelper::Sequential H;
        FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower,
                      FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
                      m,n, F.one, A, m, B, n, H);
    }
    else{
        FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads> PSH(t);
        PAR_BLOCK{
            switch (p) {
            case 1: {
                        FFLAS::TRSMHelper<FFLAS::StructureHelper::Iterative,
                        FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads> >
                        PH (PSH);
                        FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower,
                                      FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
                                      m,n, F.one, A, m, B, n, PH);
                        break;}
            case 2: {FFLAS::TRSMHelper<FFLAS::StructureHelper::Recursive,
                        FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads> >
                        PH (PSH);
                        FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower,
                                      FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
                                      m,n, F.one, A, m, B, n, PH);
                        break;}
            case 3:
                    FFLAS::TRSMHelper<FFLAS::StructureHelper::Hybrid, FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads> >
                    PH (PSH);
                    FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower,
                                  FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
                                  m,n, F.one, A, m, B, n, PH);
                    break;
            }

        }
    }
    if (i) {chrono.stop(); time+=chrono.realtime();}
}

FFLAS::fflas_delete( A);
FFLAS::fflas_delete( B);

// -----------
// Standard output for benchmark - Alexis Breust 2014/11/14
std::cout << "Time: " << time / double(iter)
<< " Gfops: " << (double(m)/1000.*double(m)/1000.*double(n)/1000.0) / time * double(iter);
FFLAS::writeCommandString(std::cout, as) << std::endl;

return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
