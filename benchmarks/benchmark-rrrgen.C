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
#include <fstream>

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
    // size_t ldts = m;


    // Element_ptr A, B, TS;
    Element_ptr A;
    A = FFLAS::fflas_new (F, n, n);
    // B = FFLAS::fflas_new (F, n, n);
    // TS = FFLAS::fflas_new (F, n, m);
    // Element_ptr Res = fflas_new(F, n, m); // Inadequate name
    RRRgen<Field>* RRRA;
    // RRRgen<Field>* RRRB;
    RRRgen<Field>* RRRres;
    RRRgen<Field>* RRRL;
    RRRgen<Field>* RRRU;


    Element_ptr A2 = fflas_new (F, n, n);
    // Element_ptr B2 = fflas_new (F, n, n);
    double time_invert = 0, time_LU = 0;
    // double time_RRRxTS = 0, time_RRRxRRR = 0,time_gen_qs = 0, time_gen_rrr = 0;
    size_t * p = FFLAS::fflas_new<size_t> (ceil(n/2.));
    for (size_t i = 0; i < ceil(n/2.); i++)
    {
        p[i] = n - i - 1;
    }

    Givaro::GeneralRingNonZeroRandIter<Field,typename Field::RandIter> nzG (G);
    for (size_t i=0; i<iter;++i){

        // cout << "Generation start"<< endl;
        // generates random matrices
        typename Field::RandIter G (F, seed);
        // RandomMatrix(F, n, m, TS, ldts, G);
        
        // generate A a t qsmatrix
        
        // chrono.clear();
        // chrono.start(); 


        RandomLTQSMatrixWithRankandQSorder (F,n,r,t,A, n,G);
        RandomLTQSMatrixWithRankandQSorder (F,n,r,t, A2, n,G);
        applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), A, n, p);
        applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), A2, n, p);
        for (size_t i=0; i< n; ++i)
                nzG.random (A2[i*(n+1)]);
        faddin (F, n, n, A2, n, A, n);


        // chrono.stop();
        // time_gen_qs = chrono.usertime();
        
        // cout << "Generation of A ended in "<< time_gen_qs << endl;

        // generate B a t qsmatrix
        // chrono.clear();
        // chrono.start(); 

        // RandomLTQSMatrixWithRankandQSorder (F,n,r,t, B, n,G);
        // RandomLTQSMatrixWithRankandQSorder (F,n,r,t, B2, n,G);
        // applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), B, n, p);
        // applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), B2, n, p);
        // faddin (F, n, n, B2, n, B, n);
        // chrono.stop();
        // time_gen_qs += chrono.usertime();
        
        // cout << "Generation of B ended in "<< time_gen_qs << endl;

        // create RRR
        // chrono.clear();
        // chrono.start(); 
        RRRA = new RRRgen<Field>(F, n, t, A, lda,true,true);
        // chrono.stop();
        // time_gen_rrr += chrono.usertime();

        // create RRR
        // chrono.clear();
        // chrono.start(); 
        // RRRB = new RRRgen<Field>(F, n, t, B, lda,true,true);
        // chrono.stop();
        // time_gen_rrr += chrono.usertime();

        // RRRxTS product
        // chrono.clear();
        // chrono.start(); 
        // RRRxTS(F,n,m,RRRA,TS,m, Res,m);
        // chrono.stop();
        // time_RRRxTS += chrono.usertime();
        
        // RRRxRRR
        // chrono.clear();
        // chrono.start(); 
        // RRRres = RRRxRRR(F,RRRA,RRRB);
        // chrono.stop();
        // time_RRRxRRR += chrono.usertime();
        // delete RRRres;
        
        
        
        RRRL = nullptr;
        RRRU = nullptr;
        chrono.clear();
        try {
            chrono.start(); 
            LUfactRRR(F,RRRA,RRRL,RRRU);
            chrono.stop();
        }
        catch(...){
            delete RRRU;
            delete RRRL;
            delete RRRA;
            FFLAS::fflas_delete(A);
            FFLAS::fflas_delete(A2, p);
            throw std::runtime_error("RRR Error: LU with non invertible matrix ");
        }
        time_LU += chrono.usertime();
        
        chrono.clear();
        chrono.start(); 
        RRRres = RRRinvert(F,RRRA);
        chrono.stop();
        time_invert += chrono.usertime();
        delete RRRres;



        delete RRRA;
        delete RRRU;
        delete RRRL;
        // delete RRRB;
    }
    
    FFLAS::fflas_delete(A);
    // FFLAS::fflas_delete(B);
    // FFLAS::fflas_delete(TS);
    // FFLAS::fflas_delete(Res); 
    // FFLAS::fflas_delete(A2,B2, p);
    FFLAS::fflas_delete(A2, p);
    
    // double mean_time_RRRxTS = time_RRRxTS / double(iter);
    // double mean_time_gen_RRR = time_gen_rrr / double(2*iter);
    // double mean_time_RRRxRRR = time_RRRxRRR / double(iter);
    double mean_time_invert = time_invert / double(iter);
    double mean_time_LU_fact = time_LU / double(iter);
    double time = mean_time_invert + mean_time_LU_fact;
    #define GFOPS_LU(x) (double(10*n*t*t)*log(double(n))/x) //Gfops = 10*n*s^2*log(n/s) / Time
    
    std::cout << "Time: " << time 
    << " Gfops: " << GFOPS_LU(time)  
    << " | details { Time : " 
    << "LU_RRR : " << mean_time_LU_fact << "; " 
    << "Invert_RRR : " << mean_time_invert << "; "  
    << "Gfops : "
    << "Gfops_LU_RRR : " << GFOPS_LU(mean_time_LU_fact) << "; "
    << "Gfops_Invert_RRR : " << GFOPS_LU(mean_time_invert) << "} ";
    return ;
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
    FFLAS::writeCommandString(std::cout, as)  << std::endl;
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
