/* Copyright (c) FFLAS-FFPACK
 * Written by 
 * Clement Pernet <clement.pernet@univ-grenoble-alpes.fr>
 * Hippolyte Signargout <hippolyte.signargout@ens-lyon.fr>
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


using namespace std;
using namespace FFLAS;
using namespace FFPACK;

template<class Field>
void run_with_field(int q, size_t n, size_t m, size_t t, size_t r, size_t iter, uint64_t seed){

    Field F(q);
    typedef typename Field::Element_ptr Element_ptr;

    FFLAS::Timer chrono;
    Element_ptr A, TS;

    double time_gen = 0, time_cbxts =0;
    double time_RPMGen_Tom = 0, time_RPMGen = 0;
    for (size_t i=0;i<iter;++i){

        if ( r >= 1000 && r <= 1400){
            size_t * rows1 = FFLAS::fflas_new<size_t> (r);
            size_t * cols1 = FFLAS::fflas_new<size_t> (r);

            chrono.clear();
            chrono.start();
            RandomLTQSRankProfileMatrix (n, r, t, rows1, cols1);
            chrono.stop();

            time_RPMGen += chrono.usertime();
            FFLAS::fflas_delete(rows1,cols1);
        }


        //std::cout << "rows :";
        //for ( size_t i = 0; i<r; i++){
        //    std::cout << rows1[i] << ",";
        //}
        //std::cout << std::endl;

        //std::cout << "cols :";
        //for ( size_t j = 0; j<r; j++){
        //    std::cout << cols1[j] << ",";
        //}
        //std::cout << std::endl;

        size_t * rows2 = FFLAS::fflas_new<size_t> (r);
        size_t * cols2 = FFLAS::fflas_new<size_t> (r);

        chrono.clear();
        chrono.start();
        RandomLTQSRankProfileMatrix_Tom (n, r, t, rows2, cols2);
        chrono.stop();

        time_RPMGen_Tom += chrono.usertime();

        FFLAS::fflas_delete(rows2,cols2);

        A = FFLAS::fflas_new (F, n, n);
        size_t lda=n;
        TS = FFLAS::fflas_new (F, n, m);
        size_t ldts = m;
        typename Field::RandIter G (F, seed);

        RandomLTQSMatrixWithRankandQSorder (F,n,r,t,A,lda,G);
        RandomMatrix(F, n, m, TS, ldts, G);
        
        size_t * P = fflas_new<size_t> (n);
        size_t * Q = fflas_new<size_t> (n);
        Element_ptr L = fflas_new(F,n,n);
            //      Element_ptr R = fflas_new(F,n,n);
        Element_ptr U = fflas_new(F,n,n);
        
        Element_ptr Xu = fflas_new(F, 2*t, n);
        size_t * Ku = fflas_new<size_t> (r+1);
        size_t * Mu = fflas_new<size_t> (n);
        size_t * Tu = fflas_new<size_t>(r);
        Element_ptr Xl = fflas_new(F, n, 2*t);
        size_t * Kl = fflas_new<size_t> (r+1);
        size_t * Ml = fflas_new<size_t> (n);
        size_t * Tl = fflas_new<size_t>(r);
  
        size_t r2;
        chrono.clear();
        chrono.start();
        r2 =  LTBruhatGen (F, FflasNonUnit, n, A, lda, P, Q);
//        getLTBruhatGen(F, n, r, P, Q, R, n);
        getLTBruhatGen(F, FflasLower, FflasUnit, n, r, P, Q, A, lda, L,n);
        size_t NbBlocksL = CompressToBlockBiDiagonal(F, FflasLower, n, t, r, P, Q, L,n ,Xl,2*t,Kl,Ml,Tl);
        getLTBruhatGen(F, FflasUpper, FflasNonUnit, n, r, P, Q, A, lda, U, n);
        size_t NbBlocksU = CompressToBlockBiDiagonal(F, FflasUpper, n, t, r, P, Q, U,n ,Xu,n,Ku,Mu,Tu);
        chrono.stop();

        if (r2!=r){ std::cerr<<"ERROR: r != r2"<<std::endl; exit(-1);}

        time_gen+=chrono.usertime();
        FFLAS::fflas_delete(A,L,U);

        Element_ptr CBruhat = fflas_new(F, n, m);
 
        chrono.clear();
        chrono.start();
        productBruhatxTS(F, n, t, r, m, P, Q, Xu, n, NbBlocksU, Ku, Tu, Mu,Xl, 2*t, NbBlocksL, Kl, Tl, Ml,TS, ldts, F.zero, CBruhat, m);
        chrono.stop();

        time_cbxts += chrono.usertime();
        
        FFLAS::fflas_delete(TS,P,Q,Xu,Ku,Mu,Tu,Xl,Ml,Tl);
    }
    // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14
    // std::cout << "Time: " << (time_gen + time_cbxts) / double(iter)  << " Gfops: Irrelevant (Generator) Specific times: " << time_gen / double(iter)<<" (for construction)" << time_cbxts / double(iter)<<" (for CB x TS)" ;
    std::cout << "Time_stand : " << time_RPMGen / double(iter) << std::endl;
    std::cout << "Time_Tom : " << time_RPMGen_Tom / double(iter) << std::endl;
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
        { 'q', "-q Q", "Set the field characteristic (-1 for the ring ZZ).",  TYPE_INT , &q },
        { 'n', "-n N", "Set the order of the square matrix A.",               TYPE_INT , &n },
        { 'm', "-m M", "Set the column dimension of n x m RHS matrix B.",               TYPE_INT , &m },
        { 't', "-t T", "Set the quasiseparability order of A.",  TYPE_INT , &t },
        { 'r', "-r R", "Set the rank of each upper/lower triangular part of A.",  TYPE_INT , &r },
        { 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &iter },
        { 's', "-s S", "Sets seed.", TYPE_INT , &seed },
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
