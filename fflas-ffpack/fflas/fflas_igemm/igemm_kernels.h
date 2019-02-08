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

#ifndef __FFLASFFPACK_fflas_igemm_igemm_kernels_H
#define __FFLASFFPACK_fflas_igemm_igemm_kernels_H

namespace FFLAS { namespace details {


    /* ************* */
    /*  GEBP KERNELS */
    /* ************* */

    template<enum number_kind K>
    inline void igebb44(size_t i, size_t j, size_t depth, size_t pdeth
                        , const int64_t alpha
                        , const int64_t *blA, const int64_t* blB
                        , int64_t* C, size_t ldc
                       );

    template<enum number_kind K>
    inline void igebb24(size_t i, size_t j, size_t depth, size_t pdeth
                        , const int64_t alpha
                        , const int64_t *blA, const int64_t* blB
                        , int64_t* C, size_t ldc
                       );

    template<enum number_kind K>
    inline void igebb14(size_t i, size_t j, size_t depth, size_t pdeth
                        , const int64_t alpha
                        , const int64_t *blA, const int64_t* blB
                        , int64_t* C, size_t ldc
                       );

    template<enum number_kind K>
    inline void igebb41(size_t i, size_t j, size_t depth, size_t pdeth
                        , const int64_t alpha
                        , const int64_t *blA, const int64_t* blB
                        , int64_t* C, size_t ldc
                       );

    template<enum number_kind K>
    inline void igebb21(size_t i, size_t j, size_t depth, size_t pdeth
                        , const int64_t alpha
                        , const int64_t *blA, const int64_t* blB
                        , int64_t* C, size_t ldc
                       );

    template<enum number_kind K>
    inline void igebb11(size_t i, size_t j, size_t depth, size_t pdeth
                        , const int64_t alpha
                        , const int64_t *blA, const int64_t* blB
                        , int64_t* C, size_t ldc
                       );


    /*************************
     *  MAIN GEBP OPERATION  *
     ************************/

    template<enum number_kind K>
    void igebp( size_t rows, size_t cols, size_t depth
                , const int64_t alpha
                , const int64_t* blockA, size_t lda,
                const int64_t* blockB, size_t ldb,
                int64_t* C, size_t ldc);

} // details
} // FFLAS

#include "igemm_kernels.inl" // could be .C

#endif // __FFLASFFPACK_fflas_igemm_igemm_kernels_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
