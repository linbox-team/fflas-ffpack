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

// Template from benchmark-sss.C

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
void run_with_field(int q, size_t n, size_t m, size_t s, size_t r, size_t iter, uint64_t seed){
    Field F(q);
    typedef typename Field::Element_ptr Element_ptr;
    FFLAS::Timer chrono;
    size_t lda=n;
    size_t ldts = m;
    size_t rs = n%s;           // Size of the partial block
    size_t ls = (rs)? rs: s;   // Size of the last block
    
    double time_gens = 0, time_cbxtss =0;
    Element_ptr A = FFLAS::fflas_new (F, n, n);
    Element_ptr A2 = FFLAS::fflas_new (F, n, n);
    Element_ptr B = FFLAS::fflas_new (F, n, n);
    Element_ptr B2 = FFLAS::fflas_new (F, n, n);
    Element_ptr TSS = FFLAS::fflas_new (F, n, m);
    Element_ptr D = fflas_new (F, n, s);
    Element_ptr P = fflas_new (F, n - s, s);
    Element_ptr Q = fflas_new (F, n - ls, s);
    Element_ptr R = fflas_new (F, ((n > (s + ls))? (n - s - ls): 0), s);
    Element_ptr U = fflas_new (F, n - ls, s);
    Element_ptr V = fflas_new (F, n - ls, s);
    Element_ptr W = fflas_new (F, ((n > (s + ls))? (n - s - ls): 0), s);
    Element_ptr Res = fflas_new(F, n, m); // Inadequate name
    size_t * p = FFLAS::fflas_new<size_t> (n);
        for (size_t i = 0; i < ceil(n/2.); i++)
            {
                p[i] = n - i - 1;
            }
    
    double time_genb = 0, time_cbxtsb =0;
    Element_ptr H = FFLAS::fflas_new (F, n, 1);
    Element_ptr TSB = FFLAS::fflas_new (F, n, m);
    size_t * pa = fflas_new<size_t> (n);
    size_t * qa = fflas_new<size_t> (n);
    Element_ptr L = fflas_new(F,n,n);
    Element_ptr Ua = fflas_new(F,n,n);
    Element_ptr Xu = fflas_new(F, 2*s, n);
    size_t * Ku = fflas_new<size_t> (r+1);
    size_t * Mu = fflas_new<size_t> (n);
    size_t * Tu = fflas_new<size_t>(r);
    Element_ptr Xl = fflas_new(F, n, 2*s);
    size_t * Kl = fflas_new<size_t> (r+1);
    size_t * Ml = fflas_new<size_t> (n);
    size_t * Tl = fflas_new<size_t>(r);
    size_t * pb = fflas_new<size_t> (n);
    size_t * qb = fflas_new<size_t> (n);
    Element_ptr Lb= fflas_new(F,n,n);
    Element_ptr Ub= fflas_new(F,n,n);
    Element_ptr Xub= fflas_new(F, 2*s, n);
    size_t * Kub= fflas_new<size_t> (r+1);
    size_t * Mub= fflas_new<size_t> (n);
    size_t * Tub= fflas_new<size_t>(r);
    Element_ptr Xlb= fflas_new(F, n, 2*s);
    size_t * Klb= fflas_new<size_t> (r+1);
    size_t * Mlb= fflas_new<size_t> (n);
    size_t * Tlb= fflas_new<size_t>(r);
    size_t r2;
        size_t r3;
    Element_ptr CBruhat = fflas_new(F, n, m);
    
    for (size_t i=0;i<iter;++i){
        std::cout << "." <<std::flush;
        typename Field::RandIter G (F, seed + i);
        RandomMatrix(F, n, 1, H, 1, G);
        RandomLTQSMatrixWithRankandQSorder (F,n,r,s,A,lda,G);
        fassign (F, n, n, A, lda, A2, lda);
        RandomLTQSMatrixWithRankandQSorder (F,n,r,s,B,lda,G);
        fassign (F, n, n, B, lda, B2, lda);
        
        applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), A, n, p);
        applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), B, n, p);
        faddin (F, n, n, B, n, A, n);
        faddin (F, n, H, 1, A, n+1);
        RandomMatrix(F, n, m, TSS, ldts, G);
        fassign(F, n, m, TSS, ldts, TSB, ldts);
                
        chrono.clear();
        chrono.start();
        DenseToSSS (F, n, s, P, s, Q, s, R, s, U, s, V, s, W, s, D, s, A, n);
        chrono.stop();
        time_gens+=chrono.usertime();

        chrono.clear();
        chrono.start();
        r2 =  LTBruhatGen (F, FflasNonUnit, n, A2, lda, pa, qa);
        r3 =  LTBruhatGen (F, FflasNonUnit, n, B2, lda, pb, qb);
        chrono.stop();
        if ((r2!=r)||(r3 != r))
            {
                std::cerr<<"ERROR: r != r2 or r3"<<std::endl;
                exit(-1);
            }
        time_genb+=chrono.usertime();
 
        chrono.clear();
        chrono.start();
        productSSSxTS(F, n, m, s, F.one, P, s, Q, s, R, s, U, s, V, s, W, s,
                      D, s, TSS, m, F.zero, Res, m);
        chrono.stop();
        time_cbxtss += chrono.usertime();
 
        chrono.clear();
        chrono.start();
        getLTBruhatGen(F, FflasLower, FflasUnit, n, r, pa, qa, A2, lda, L,n);
        size_t NbBlocksL = CompressToBlockBiDiagonal(F, FflasLower, n, s, r, pa, qa, L,n ,Xl,2*s,Kl,Ml,Tl);
        getLTBruhatGen(F, FflasUpper, FflasNonUnit, n, r, pa, qa, A2, lda, Ua, n);
        size_t NbBlocksU = CompressToBlockBiDiagonal(F, FflasUpper, n, s, r, pa, qa, Ua,n ,Xu,n,Ku,Mu,Tu);
        getLTBruhatGen(F, FflasLower, FflasUnit, n, r, pb, qb, B2, lda, Lb,n);
        size_t NbBlocksLb = CompressToBlockBiDiagonal(F, FflasLower, n, s, r, pb, qb, Lb,n ,Xlb,2*s,Klb,Mlb,Tlb);
        getLTBruhatGen(F, FflasUpper, FflasNonUnit, n, r, pb, qb, B2, lda, Ub, n);
        size_t NbBlocksUb = CompressToBlockBiDiagonal(F, FflasUpper, n, s, r, pb, qb, Ub,n ,Xub,n,Kub,Mub,Tub);
        productBruhatxTS(F, n, s, r, m, pa, qa, Xu, n, NbBlocksU, Ku, Tu, Mu,Xl, 2*s,
                         NbBlocksL, Kl, Tl, Ml,TSB, ldts, F.zero, CBruhat, m);
        applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, m, 0, ceil(n/2.), CBruhat, m, p);
        for (size_t i = 0 ; i < n ; ++i)
            faxpy(F, 1, m, H[i], TSB + ldts * i, ldts, CBruhat + ldts*i, ldts);
        applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, m, 0, ceil(n/2.), TSB, m, p);
        productBruhatxTS(F, n, s, r, m, pb, qb, Xub, n, NbBlocksUb, Kub, Tub, Mub,Xlb, 2*s,
                         NbBlocksLb, Klb, Tlb, Mlb,TSB, ldts, F.one, CBruhat, m);
        chrono.stop();
        time_cbxtsb += chrono.usertime();
    }
    FFLAS::fflas_delete(A, A2, B, B2, D, P, Q, R, U, V, W, p, H, L, Ua, Lb, Klb, CBruhat, TSS, Res, pa,qa,Xu,Ku);
    FFLAS::fflas_delete(Mu,Tu,Xl,Ml,Tl, Kl, TSB,pb,qb, Kub, Mub,Tub, Mlb,Tlb);
    FFLAS::fflas_delete( Ub);
    FFLAS::fflas_delete(Xub);
    FFLAS::fflas_delete(Xlb);
        // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14
    std::cout << std:: endl << "Time: " << (time_gens + time_genb + time_cbxtss + time_cbxtsb) / double(iter)
              << " Gfops: Irrelevant. Specific times: "
              << " DenseToSSS: " << time_gens / double(iter)
              << "; AppSSS: " << time_cbxtss / double(iter)
              << " DenseToBruhat: " << time_genb / double(iter)
              << "; AppBruhat: " << time_cbxtsb / double(iter) << ". r ="
              << r << " s = " << s << std::endl ;
}

int main(int argc, char** argv) {

#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t iter = 50;
    int    q    = 131071;
    size_t    n    = 3000;
    size_t    m    = 500;
    size_t    t    = 300;
    size_t    r    = 1500;
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
    for (size_t rank = r;  rank < 2100; rank += 500)
        {
            for (size_t order = t;  order < rank / 2; order += 50)
                {
                    run_with_field<Givaro::ModularBalanced<double> >(q, n, m, order, rank, iter, seed);
                }
        }
    std::cout << "( ";
    FFLAS::writeCommandString(std::cout, as) << ")" << std::endl;
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
