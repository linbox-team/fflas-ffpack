/*
 * Copyright (C) 2015 FFLAS-FFPACK
 *
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *
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
 * Lesser General Public License for more detAils.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/** @file fflas_lvl3.C
 * @author  Brice Boyer
 * @brief C functions calls for level 3 FFLAS in fflas-c.h
 * @see fflas/fflas_level3.inl
 */

#include "fflas-ffpack/interfaces/libs/fflas_c.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "givaro//modular-balanced.h"
#include "givaro//modular.h"

using Givaro::Modular ;
using Givaro::ModularBalanced ;
using namespace FFLAS ;


#ifdef __cplusplus
extern "C" {
#endif

    void
    ftrsm_3_modular_double (const double p, const enum FFLAS_C_SIDE Side,
                            const enum FFLAS_C_UPLO Uplo,
                            const enum FFLAS_C_TRANSPOSE tA,
                            const enum FFLAS_C_DIAG Diag,
                            const size_t m, const size_t n,
                            const double alpha,
                            const double * A,
                            const size_t ldA,
                            double * B, const size_t ldB
                            , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            ftrsm(F,(enum FFLAS_SIDE)Side,(enum FFLAS_UPLO)Uplo,(FFLAS_TRANSPOSE)tA,(enum FFLAS_DIAG)Diag,m,n,alpha,A,ldA,B,ldB);
        } else {
            ModularBalanced<double> F(p);
            ftrsm(F,(enum FFLAS_SIDE)Side,(enum FFLAS_UPLO)Uplo,(FFLAS_TRANSPOSE)tA,(enum FFLAS_DIAG)Diag,m,n,alpha,A,ldA,B,ldB);
        }
    }


    void
    ftrmm_3_modular_double (const double p, const enum FFLAS_C_SIDE Side,
                            const enum FFLAS_C_UPLO Uplo,
                            const enum FFLAS_C_TRANSPOSE tA,
                            const enum FFLAS_C_DIAG Diag,
                            const size_t m, const size_t n,
                            const double alpha,
                            double * A, const size_t ldA,
                            double * B, const size_t ldB
                            , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            ftrmm(F,(enum FFLAS_SIDE)Side,(enum FFLAS_UPLO)Uplo,(FFLAS_TRANSPOSE)tA,(enum FFLAS_DIAG)Diag,m,n,alpha,A,ldA,B,ldB);
        } else {
            ModularBalanced<double> F(p);
            ftrmm(F,(enum FFLAS_SIDE)Side,(enum FFLAS_UPLO)Uplo,(FFLAS_TRANSPOSE)tA,(enum FFLAS_DIAG)Diag,m,n,alpha,A,ldA,B,ldB);
        }
    }

    double *
    fgemm_3_modular_double( const double p,
                            const enum FFLAS_C_TRANSPOSE tA,
                            const enum FFLAS_C_TRANSPOSE tB,
                            const size_t m,
                            const size_t n,
                            const size_t k,
                            const double alpha,
                            const double * A, const size_t ldA,
                            const double * B, const size_t ldB,
                            const double betA,
                            double * C, const size_t ldC,
                            bool positive )

    {
        if (positive) {
            Modular<double> F(p);
            return fgemm(F,(FFLAS_TRANSPOSE)tA,(FFLAS_TRANSPOSE)tB,m,n,k,alpha,A,ldA,B,ldB,betA,C,ldC);
        } else {
            ModularBalanced<double> F(p);
            return fgemm(F,(FFLAS_TRANSPOSE)tA,(FFLAS_TRANSPOSE)tB,m,n,k,alpha,A,ldA,B,ldB,betA,C,ldC);
        }
        return nullptr;
    }


    double *
    fsquare_3_modular_double (const double p,
                              const enum FFLAS_C_TRANSPOSE tA,
                              const size_t n,
                              const double alpha,
                              const double * A,
                              const size_t ldA,
                              const double betA,
                              double * C,
                              const size_t ldC
                              , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            return fsquare(F,(FFLAS_TRANSPOSE)tA,n,alpha,A,ldA,betA,C,ldC);
        } else {
            ModularBalanced<double> F(p);
            return fsquare(F,(FFLAS_TRANSPOSE)tA,n,alpha,A,ldA,betA,C,ldC);
        }
        return nullptr;
    }



#ifdef __cplusplus
}
#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
