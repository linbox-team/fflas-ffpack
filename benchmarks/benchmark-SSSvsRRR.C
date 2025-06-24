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
    
    size_t lda = n;
    size_t ldts = m;
    size_t rs = n%t;           // Size of the partial block
    size_t ls = (rs)? rs: t;   // Size of the last block

    Element_ptr A, B, TS;
    A = FFLAS::fflas_new (F, n, n);
    B = FFLAS::fflas_new (F, n, m);
    TS = FFLAS::fflas_new (F, n, m);
    Element_ptr D = fflas_new (F, n, t);
    Element_ptr P = fflas_new (F, n - t, t);
    Element_ptr Q = fflas_new (F, n - ls, t);
    Element_ptr R = fflas_new (F, ((n > (t + ls))? (n - t - ls): 0), t);
    Element_ptr U = fflas_new (F, n - ls, t);
    Element_ptr V = fflas_new (F, n - ls, t);
    Element_ptr W = fflas_new (F, ((n > (t + ls))? (n - t - ls): 0), t);
    Element_ptr Res = fflas_new(F, n, m); // Inadequate name
    RRRgen<Field>* RRRA;


    Element_ptr A2 = fflas_new (F, n, n);


    double time_RRRxTS = 0, time_SSSxTS = 0;
    size_t * p = FFLAS::fflas_new<size_t> (ceil(n/2.));
    for (size_t i = 0; i < ceil(n/2.); i++)
    {
        p[i] = n - i - 1;
    }
    
    for (size_t i=0; i<iter;++i){

        std::cout << "Generation start"<< std::endl;
        // generates random matrices
        typename Field::RandIter G (F, seed);
        RandomMatrix(F, n, m, TS, ldts, G);
        
        // generate A a t qsmatrix
        RandomLTQSMatrixWithRankandQSorder (F,n,r,t,A, n,G);
        RandomLTQSMatrixWithRankandQSorder (F,n,r,t, A2, n,G);
        applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), A, n, p);
        applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), A2, n, p);
        faddin (F, n, n, A2, n, A, n);
        
        std::cout << "Generation ended"<< std::endl;

        // create RRR
        RRRA = new RRRgen<Field>(F, n, t, A, lda,true,true);
        
        // RRRxTS product
        chrono.clear();
        chrono.start(); 
        RRRxTS(F,n,m,RRRA,TS,m,B,m);
        chrono.stop();
        time_RRRxTS += chrono.usertime();
        delete RRRA;
        
        
        
        
        // create SSS
        DenseToSSS (F, n, t, P, t, Q, t, R, t, U, t, V, t, W, t, D, t, A, n);
        
        // SSSxTS product
        chrono.clear();
        chrono.start();
        productSSSxTS(F, n, m, t, F.one, P, t, Q, t, R, t, U, t, V, t, W, t,
            D, t, TS, m, F.zero, Res, m);
            chrono.stop();
            
        time_SSSxTS += chrono.usertime();
    }
    
    FFLAS::fflas_delete(A, D, P, Q, R, U, V, W, B,TS, Res); 
    FFLAS::fflas_delete(A2, p);
    
    double mean_time_RRRxTS = time_RRRxTS / double(iter);
    double mean_time_SSSxTS = time_SSSxTS / double(iter);
    std::cout << "Finished for n = " << n << " for RRR product took " << mean_time_RRRxTS << " and SSS product took " << mean_time_SSSxTS << std::endl;
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
