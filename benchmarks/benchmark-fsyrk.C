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
#include "fflas-ffpack/utils/args-parser.h"


using namespace std;
using namespace FFLAS;
using namespace FFPACK;
int main(int argc, char** argv) {

 #ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t iter = 3;
    size_t    q    = 131071;
    size_t    n    = 1000;
    size_t    k    = 1000;
    size_t algo    = 0;
    size_t threshold = 64;
    bool up =true;
    std::string file = "";

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
        { 'n', "-n N", "Set the dimension of the matrix C.",               TYPE_INT , &n },
        { 'k', "-k K", "Set the other dimension of the matrix A.",               TYPE_INT , &k },
        { 'u', "-u yes/no", "Updates an upper/lower triangular matrix.",  TYPE_BOOL , &up },
        { 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &iter },
        { 'a', "-a A", "Set algorithmic variant.",                     TYPE_INT , &algo },    
        { 't', "-t T", "Set the threshold to the base case.",                     TYPE_INT , &threshold },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    typedef Givaro::ModularBalanced<double> Field;
    typedef Field::Element Element;

    Field F(q);
    Field::Element * A, *C;

    FFLAS::Timer chrono;
    double time=0.0;

    FFLAS_UPLO uplo = up?FflasUpper:FflasLower;
    for (size_t i=0;i<=iter;++i){
        A = fflas_new<Element>(n*k);
        size_t lda=k;
        C = fflas_new<Element>(n*n);
        size_t ldc=n;
        Field::RandIter G(F);
        RandomMatrix (F, n, k, A, k, G);
        RandomTriangularMatrix (F, n,n,uplo, FflasNonUnit, true, C, ldc, G);
        Field::Element_ptr D = FFLAS::fflas_new(F,k,1);
        Givaro::GeneralRingNonZeroRandIter<Field,Field::RandIter> nzG (G);
        std::vector<bool> twoBlocks(k);
        for (size_t i=0; i<k; i++)
            nzG.random(D[i]);
        chrono.clear();
        if (i) chrono.start();
        switch (algo){
            case 0: // fsyrk with no diagonal scaling
                fsyrk (F, uplo, FflasNoTrans, n, k, F.mOne, A, lda, F.one, C, ldc);
                break;
            case 1: // fsyrk with diagonal scaling
                fsyrk (F, uplo, FflasNoTrans, n, k, F.mOne, A, lda, D, 1, F.one, C, ldc, threshold);
                break;
            case 2: // fsyrk with diagonal scaling and 2x2 diagonal blocks
                fsyrk (F, uplo, FflasNoTrans, n, k, F.one, A, lda, D, 1, twoBlocks, F.zero, C, ldc, threshold);
                break;
            case 3: // fsyrk with Strassen and no diagonal scaling
            {
                size_t reclevel = 0;
                size_t dim = n;
                while(dim > threshold) {reclevel++; dim>>=1;}
                MMHelper<Field, MMHelperAlgo::Winograd> H(F,reclevel);
                fsyrk (F, uplo, FflasNoTrans, n, k, F.one, A, lda, F.zero, C, ldc, H);
                break;
            }
            case 4: // cblas_dsyrk
                cblas_dsyrk (CblasRowMajor, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) FflasNoTrans, n, k, 1.0, A, lda, 0.0, C, ldc);
                break;
            case 5: // fsyrk with the classic Divide and conquer algorithm
                size_t reclevel = 0;
                size_t dim = n;
                while(dim > threshold) {reclevel++; dim>>=1;}
                MMHelper<Field, MMHelperAlgo::DivideAndConquer> H(F,reclevel);
                fsyrk (F, uplo, FflasNoTrans, n, k, F.one, A, lda, F.zero, C, ldc, H);
                break;
        }
        if (i) chrono.stop();

        time+=chrono.usertime();
        FFLAS::fflas_delete( A);
        FFLAS::fflas_delete( C);
        FFLAS::fflas_delete( D);
    }

    // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14
    std::cout << "Time: " << time / double(iter)
              << " Gfops: " << (double(n)/1000.)*(double(n)/1000.)*(double(k)/1000.)/ time * double(iter);
    FFLAS::writeCommandString(std::cout, as) << std::endl;
    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
