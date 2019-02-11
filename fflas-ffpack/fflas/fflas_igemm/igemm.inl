/*
 * Copyright (C) 2013,2014  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 * the code is inspired and adapted from the Eigen library
 * modified by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

#ifndef __FFLASFFPACK_fflas_igemm_igemm_INL
#define __FFLASFFPACK_fflas_igemm_igemm_INL

#include "fflas-ffpack/utils/fflas_memory.h"

namespace FFLAS { namespace Protected {


    // Assume matrices A,B,C are stored in column major order
    template<enum FFLAS_TRANSPOSE tA, enum FFLAS_TRANSPOSE tB>
    void igemm_colmajor(size_t rows, size_t cols, size_t depth,
                        const int64_t alpha,
                        const int64_t* A, size_t lda, const int64_t* B, size_t ldb,
                        int64_t* C, size_t ldc)
    {
        FFLASFFPACK_check(alpha != 0);
        switch(alpha) {
        case 1:
            igemm_colmajor<tA,tB,number_kind::one>(rows,cols,depth,
                                                   alpha,A,lda,B,ldb,
                                                   C,ldc);
            break;
        case -1:
            igemm_colmajor<tA,tB,number_kind::mone>(rows,cols,depth,
                                                    alpha,A,lda,B,ldb,
                                                    C,ldc);
            break;
        default:
            igemm_colmajor<tA,tB,number_kind::other>(rows,cols,depth,
                                                     alpha,A,lda,B,ldb,
                                                     C,ldc);


        }
    }

    template<enum FFLAS_TRANSPOSE tA, enum FFLAS_TRANSPOSE tB, enum number_kind alpha_kind>
    void igemm_colmajor(size_t rows, size_t cols, size_t depth,
                        const int64_t alpha,
                        const int64_t* A, size_t lda, const int64_t* B, size_t ldb,
                        int64_t* C, size_t ldc)
    {
        using simd = Simd<int64_t> ;
        size_t mc,kc,nc;
        mc=rows;
        nc=cols;
        kc=depth;
        FFLAS::details::BlockingFactor(mc,nc,kc);
        size_t sizeA = mc*kc;
        size_t sizeB = kc*cols;

        // these data must be simd::alignment byte aligned
        int64_t *blockA, *blockB;


        blockA = fflas_new<int64_t>(sizeA, (Alignment)simd::alignment);
        blockB = fflas_new<int64_t>(sizeB, (Alignment)simd::alignment);

        // For each horizontal panel of B, and corresponding vertical panel of A
        for(size_t k2=0; k2<depth; k2+=kc){

            const size_t actual_kc = std::min(k2+kc,depth)-k2;
            FFLASFFPACK_check(kc <= depth);

            // pack horizontal panel of B into sequential memory (L2 cache)
            if (tB == FflasNoTrans)
                FFLAS::details::pack_rhs<_nr,false>(blockB, B+k2, ldb, actual_kc, cols);
            else
                FFLAS::details::pack_lhs<_nr,true>(blockB, B+k2*ldb, ldb, cols, actual_kc);

            // For each mc x kc block of the lhs's vertical panel...
            for(size_t i2=0; i2<rows; i2+=mc){

                const size_t actual_mc = std::min(i2+mc,rows)-i2;


                FFLASFFPACK_check(mc <= rows);
                // pack a chunk of the vertical panel of A into a sequential memory (L1 cache)
                if (tA == FflasNoTrans)
                    FFLAS::details::pack_lhs<_mr,false>(blockA, A+i2+k2*lda, lda, actual_mc, actual_kc);
                else
                    FFLAS::details::pack_rhs<_mr,true>(blockA, A+i2*lda+k2, lda, actual_kc, actual_mc);

                // call block*panel kernel
                FFLAS::details::igebp<alpha_kind>(actual_mc, cols, actual_kc
                                                  , alpha
                                                  , blockA, actual_kc, blockB, actual_kc
                                                  , C+i2, ldc);
            }
        }

        fflas_delete(blockA);
        fflas_delete(blockB);
    }

    void igemm( const enum FFLAS_TRANSPOSE TransA, const enum FFLAS_TRANSPOSE TransB,
                size_t rows, size_t cols, size_t depth
                , const int64_t alpha
                , const int64_t* A, size_t lda, const int64_t* B, size_t ldb
                , const int64_t beta
                , int64_t* C, size_t ldc
              )
    {
        if (!rows || !cols) {
            return ;
        }

        //! @todo use primitive (no Field()) and  specialise for int64.
        // CP: fscalin assumes C in row major mode and we are here in col major mode
        // hence let's transpose the arguments.
        fscalin(Givaro::ZRing<int64_t>(),cols,rows, beta,C,ldc);
        if (!depth || alpha == 0) {
            return ;
        }
        if (TransA == FflasNoTrans) {
            if (TransB == FflasNoTrans) {
                igemm_colmajor<FflasNoTrans,FflasNoTrans>(rows, cols, depth, alpha, A, lda, B, ldb, C, ldc);
            }
            else {
                igemm_colmajor<FflasNoTrans,FflasTrans>(rows, cols, depth, alpha, A, lda, B, ldb, C, ldc);
            }
        }
        else {
            if (TransB == FflasNoTrans) {
                igemm_colmajor<FflasTrans,FflasNoTrans>(rows, cols, depth, alpha, A, lda, B, ldb, C, ldc);
            }
            else {
                igemm_colmajor<FflasTrans,FflasTrans>(rows, cols, depth, alpha, A, lda, B, ldb, C, ldc);
            }
        }
    }

} // Protected
} // FFLAS


// igemm

namespace FFLAS {
    inline void igemm_(const enum FFLAS_ORDER Order, const enum FFLAS_TRANSPOSE TransA, const enum FFLAS_TRANSPOSE TransB,
                       const size_t M, const size_t N, const size_t K,
                       const int64_t alpha,
                       const int64_t *A, const size_t lda,
                       const int64_t *B, const size_t ldb,
                       const int64_t beta,
                       int64_t *C, const size_t ldc)
    {

        if (Order == FflasColMajor)
            Protected::igemm(TransA,TransB,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);
        else
            Protected::igemm(TransB,TransA,N,M,K,alpha,B,ldb,A,lda,beta,C,ldc);
    }


} // FFLAS

#endif // __FFLASFFPACK_fflas_igemm_igemm_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
