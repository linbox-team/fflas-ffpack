/* Copyright (c) FFLAS-FFPACK
 * Written by Hippolyte Signargout <hippolyte.signargout@ens-lyon.fr>
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

// Template from benchmark-quasisep.C

#define __FFLASFFPACK_OPENBLAS_NT_ALREADY_SET 1

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/ffpack/ffpack_rrrgen.inl"



using namespace std;
using namespace FFLAS;
using namespace FFPACK;

template<class Field>
void run_with_field(int q, size_t n, size_t m, size_t t, size_t r, size_t iter, uint64_t seed){

    Field F(q);
    typedef typename Field::Element_ptr Element_ptr;
    typename Field::RandIter G (F, seed);
    FFLAS::Timer chrono;
    Element_ptr A, B;
    A = FFLAS::fflas_new (F, n, n);
    B = FFLAS::fflas_new (F, n, n);

    RRgen<Field>* RRA;
    RRgen<Field>* RRB;
    RRgen<Field>* RRres;

    RRRgen<Field>* RRRA;
    RRRgen<Field>* RRRB;
    RRRgen<Field>* RRRres;

    Element_ptr A2 = fflas_new (F, n, n);
    Element_ptr B2 = fflas_new (F, n, n);

    size_t lda=n;
    size_t ldb = n;

    double time_rrxrr = 0, time_rraddrr = 0, time_rrraddrr = 0, time_rrrxrr = 0, time_rrrxrrr = 0, time_invert = 0;
    for (size_t i=0;i<iter;++i){
        // generate A a t qsmatrix
        RandomLTQSMatrixWithRankandQSorder (F,n,r,t,A, n,G);
        RandomLTQSMatrixWithRankandQSorder (F,n,r,t, A2, n,G);
        size_t * p = FFLAS::fflas_new<size_t> (ceil(n/2.));
        for (size_t i = 0; i < ceil(n/2.); i++)
        {
            p[i] = n - i - 1;
        }
        applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), A, n, p);
        applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), A2, n, p);
        faddin (F, n, n, A2, n, A, n);
        FFLAS::fflas_delete(A2, p);

        // generate B a t qsmatrix
        RandomLTQSMatrixWithRankandQSorder (F,n,r,t,B, n,G);
        RandomLTQSMatrixWithRankandQSorder (F,n,r,t, B2, n,G);
        p = FFLAS::fflas_new<size_t> (ceil(n/2.));
        for (size_t i = 0; i < ceil(n/2.); i++)
        {
            p[i] = n - i - 1;
        }
        applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), B, n, p);
        applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), B2, n, p);
        faddin (F, n, n, B2, n, B, n);
        FFLAS::fflas_delete(B2, p);
        
        // Computes the RRgen of A and B
        RRA = new RRgen<Field>(F, n, n, A, lda);
        RRB = new RRgen<Field>(F, n, n, B, ldb);
        
        // Computes the RRRgen of A and B
        RRRA = new RRRgen<Field>(F, n, t, A, lda,true,true);
        RRRB = new RRRgen<Field>(F, n, t, B, lda,true,true);
        

        // RRxRR
        chrono.clear();
        chrono.start();
        RRres = RRxRR(F,RRA,RRB);
        chrono.stop();
        delete  RRres;
        time_rrxrr+=chrono.usertime();

        // RRaddRR
        chrono.clear();
        chrono.start();
        RRres = RRaddRR(F,RRA,RRB);
        chrono.stop();
        delete  RRres;
        time_rraddrr+=chrono.usertime();
        
        // RRRaddRR
        chrono.clear();
        chrono.start();
        RRRres = RRRaddRR(F,RRRA,RRB);
        chrono.stop();
        delete  RRRres;
        time_rrraddrr+=chrono.usertime();

        // RRRxRR
        chrono.clear();
        chrono.start();
        RRres = RRRxRR(F,RRRA,RRB);
        chrono.stop();
        delete  RRres;
        time_rrrxrr+=chrono.usertime();

        // RRRxRRR
        chrono.clear();
        chrono.start();
        RRRres = RRRxRRR(F,RRRA,RRRB);
        chrono.stop();
        delete  RRRres;
        time_rrrxrrr+=chrono.usertime();

        // RRR invert
        chrono.clear();
        chrono.start();
        RRRres = RRRinvert(F,RRRA);
        chrono.stop();
        delete  RRRres;
        time_invert+=chrono.usertime();

        
        delete RRA;
        delete RRB;
        delete RRRA;
        delete RRRB;
    }
    FFLAS::fflas_delete(A, B); 
    std::cout << "Total Time: " << (time_rrxrr + time_rraddrr + time_rrraddrr + time_rrrxrr + time_rrrxrrr + time_invert) / double(iter) << std::endl
              << " RRxRR : " << time_rrxrr / double(iter) << std::endl
              << " RR+RR : " << time_rraddrr / double(iter) << std::endl
              << " RRR+RR : " << time_rrraddrr / double(iter) << std::endl
              << " RRRxRR : " << time_rrrxrr / double(iter) << std::endl
              << " RRRxRRR : " << time_rrrxrrr / double(iter) << std::endl
              << " RRR invert : " << time_invert / double(iter) << std::endl;
}

int main(int argc, char** argv) {

#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t iter = 10;
    int    q    = 131071;
    size_t    n    = 2167;
    size_t    m    = 455;
    size_t    t    = 236;
    size_t    r    = 1100;
    uint64_t seed = FFLAS::getSeed();

    Argument as[] = {
                     { 'q', "-q Q", "Set the field characteristic (-1 for the ring ZZ).",     TYPE_INT , &q },
                     { 'n', "-n N", "Set the order of the square matrix A.",                  TYPE_INT , &n },
                     { 'm', "-m M", "Set the column dimension of n x m RHS matrix B.",        TYPE_INT , &m },
                     { 't', "-t T", "Set the quasiseparability order of A.",                  TYPE_INT , &t },
                     { 'r', "-r R", "Set the rank of each upper/lower triangular part of A.", TYPE_INT , &r },
                     { 'i', "-i R", "Set number of repetitions.",                             TYPE_INT , &iter },
                     { 's', "-s S", "Sets seed.",                                             TYPE_INT , &seed },
                     END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    run_with_field<Givaro::ModularBalanced<double> >(q, n, m, t, r, iter, seed);

    std::cout << "( ";
    FFLAS::writeCommandString(std::cout, as) << ")" << std::endl;
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
