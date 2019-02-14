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

#ifndef __FFLASFFPACK_fflas_igemm_igemm_H
#define __FFLASFFPACK_fflas_igemm_igemm_H

namespace FFLAS {

    enum number_kind {
        zero  =0,
        one   =1,
        mone  =-1,
        other =2
    }  ;

} // FFLAS

#include "igemm_kernels.h"
#include "igemm_tools.h"
#include "fflas-ffpack/utils/fflas_memory.h"

namespace FFLAS { namespace Protected {



    template<enum FFLAS_TRANSPOSE tA, enum FFLAS_TRANSPOSE tB>
    void igemm_colmajor( size_t rows, size_t cols, size_t depth
                         , const int64_t alpha
                         , const int64_t* A, size_t lda, const int64_t* B, size_t ldb
                         , int64_t* C, size_t ldc
                       ) ;

    template<enum FFLAS_TRANSPOSE tA, enum FFLAS_TRANSPOSE tB, enum number_kind alpha_kind>
    void igemm_colmajor( size_t rows, size_t cols, size_t depth
                         , const int64_t alpha
                         , const int64_t* A, size_t lda, const int64_t* B, size_t ldb
                         , int64_t* C, size_t ldc
                       ) ;



    inline void igemm(const enum FFLAS_TRANSPOSE TransA, const enum FFLAS_TRANSPOSE TransB
                      , size_t rows, size_t cols, size_t depth
                      , const int64_t alpha
                      , const int64_t* A, size_t lda, const int64_t* B, size_t ldb
                      , const int64_t beta
                      , int64_t* C, size_t ldc
                     ) ;

} // Protected
} // FFLAS

namespace FFLAS { /*  igemm */

    inline void igemm_(const enum FFLAS_ORDER Order, const enum FFLAS_TRANSPOSE TransA, const enum FFLAS_TRANSPOSE TransB
                       , const size_t M, const size_t N, const size_t K
                       , const int64_t alpha
                       , const int64_t *A, const size_t lda, const int64_t *B, const size_t ldb
                       , const int64_t beta
                       , int64_t *C, const size_t ldc);


} // FFLAS

#include "igemm.inl"
#endif // __FFLASFFPACK_fflas_igemm_igemm_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
