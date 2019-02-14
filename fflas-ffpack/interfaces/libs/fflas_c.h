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

/** @file fflas-c.h
 * @author  Brice Boyer
 * @brief C functions calls for FFLAS
 * @see fflas/fflas.h
 */

#ifndef __FFLASFFPACK_interfaces_libs_fflas_c_H
#define __FFLASFFPACK_interfaces_libs_fflas_c_H
//#include "fflas-ffpack/fflas-ffpack-config.h"

#ifndef FFLAS_COMPILED
#define FFLAS_COMPILED
#endif

#include <stdbool.h>
#include <stdlib.h>
#include <inttypes.h>


#ifdef __cplusplus
extern "C" {
#endif

    /// Storage by row or col ?
    enum FFLAS_C_ORDER{
        FflasRowMajor=101, /**< row major */
        FflasColMajor=102  /**< col major */
    };
    // public:
    /// Is matrix transposed ?
    enum FFLAS_C_TRANSPOSE {
        FflasNoTrans = 111, /**< Matrix is not transposed */
        FflasTrans   = 112  /**< Matrix is transposed */
    };
    /// Is triangular matrix's shape upper ?
    enum FFLAS_C_UPLO {
        FflasUpper = 121, /**< Triangular matrix is Upper triangular (if \f$i>j\f$ then \f$T_{i,j} = 0\f$)*/
        FflasLower = 122  /**< Triangular matrix is Lower triangular (if \f$i<j\f$ then \f$T_{i,j} = 0\f$)*/
    };

    /// Is the triangular matrix implicitly unit diagonal ?
    enum FFLAS_C_DIAG {
        FflasNonUnit = 131, /**< Triangular matrix has an explicit arbitrary diagonal */
        FflasUnit    = 132 /**< Triangular matrix has an implicit unit diagonal (\f$T_{i,i} = 1\f$)*/ /**< */
    };

    /// On what side ?
    enum FFLAS_C_SIDE {
        FflasLeft  = 141,/**< Operator applied on the left */
        FflasRight = 142 /**< Operator applied on the rigth*/
    };

    /** \p FFLAS_C_BASE  determines the type of the element representation for Matrix Mult kernel. (deprecated, should not be used) */
    enum FFLAS_C_BASE {
        FflasDouble  = 151,  /**<  to use the double precision BLAS */
        FflasFloat   = 152,  /**<  to use the single precison BLAS */
        FflasGeneric = 153   /**< for any other domain, that can not be converted to floating point integers */
    };

    /* ******** *
     * LEVEL1   *
     * ******** */

    /*  Modular<double>          */
    /*  ModularBalanced<double>  */

    void
    freducein_1_modular_double (const double p, const size_t n,
                                double * X, const size_t incX
                                , bool positive );

    void
    freduce_1_modular_double (const double  F, const size_t n,
                              const double * Y, const size_t incY,
                              double * X, const size_t incX
                              , bool positive );


    void
    fnegin_1_modular_double (const double  F, const size_t n,
                             double * X, const size_t incX
                             , bool positive );


    void
    fneg_1_modular_double (const double p, const size_t n,
                           const double * Y, const size_t incY,
                           double * X, const size_t incX
                           , bool positive );

    void
    fzero_1_modular_double (const double p, const size_t n,
                            double * X, const size_t incX
                            , bool positive );


    bool
    fiszero_1_modular_double (const double p, const size_t n,
                              const double * X, const size_t incX
                              , bool positive );

    bool
    fequal_1_modular_double (const double p, const size_t n,
                             const double * X, const size_t incX,
                             const double * Y, const size_t incY
                             , bool positive );


    void
    fassign_1_modular_double (const double p, const size_t n,
                              const double * Y, const size_t incY ,
                              double * X, const size_t incX
                              , bool positive );


    void
    fscalin_1_modular_double (const double p, const size_t n, const double alpha,
                              double * X, const size_t incX
                              , bool positive );


    void
    fscal_1_modular_double (const double p, const size_t n
                            , const double alpha
                            , const double * X, const size_t incX
                            , double * Y, const size_t incY
                            , bool positive );


    void
    faxpy_1_modular_double (const double p, const size_t n,
                            const double alpha,
                            const double * X, const size_t incX,
                            double * Y, const size_t incY
                            , bool positive );

#if 0
    void
    faxpby_1_modular_double (const double p, const size_t n,
                             const double alpha,
                             const double * X, const size_t incX,
                             const double betA,
                             double * Y, const size_t incY
                             , bool positive );
#endif



    double
    fdot_1_modular_double (const double p, const size_t n,
                           const double * X, const size_t incX,
                           const double * Y, const size_t incY
                           , bool positive );


    void
    fswap_1_modular_double (const double p, const size_t n,
                            double * X, const size_t incX,
                            double * Y, const size_t incY
                            , bool positive );


    void
    fadd_1_modular_double (const double p,  const size_t n,
                           const double * A, const size_t incA,
                           const double * B, const size_t incB,
                           double * C, const size_t incC
                           , bool positive );

    void
    fsub_1_modular_double (const double p,  const size_t n,
                           const double * A, const size_t incA,
                           const double * B, const size_t incB,
                           double * C, const size_t incC
                           , bool positive );

    void
    faddin_1_modular_double (const double p,  const size_t n,
                             const double * B, const size_t incB,
                             double * C, const size_t incC
                             , bool positive );

    void
    fsubin_1_modular_double (const double p,  const size_t n,
                             const double * B, const size_t incB,
                             double * C, const size_t incC
                             , bool positive );

    /* ******** *
     * LEVEL1.5 *
     * ******** */

    // fspmv

    /* ******** *
     * LEVEL2   *
     * ******** */


    /*  Modular<double>          */
    /*  ModularBalanced<double>  */


    void
    fassign_2_modular_double (const double p, const size_t m, const size_t n,
                              const double * B, const size_t ldB ,
                              double * A, const size_t ldA
                              , bool positive  );



    void
    fzero_2_modular_double (const double p, const size_t m, const size_t n,
                            double * A, const size_t ldA
                            , bool positive  );


    bool
    fequal_2_modular_double (const double p, const size_t m, const size_t n,
                             const double * A, const size_t ldA,
                             const double * B, const size_t ldB
                             , bool positive  );


    bool
    fiszero_2_modular_double (const double p, const size_t m, const size_t n,
                              const double * A, const size_t ldA
                              , bool positive  );


    void
    fidentity_2_modular_double (const double p, const size_t m, const size_t n,
                                double * A, const size_t ldA,
                                const double d
                                , bool positive  );



    void
    freducein_2_modular_double (const double p, const size_t m , const size_t n,
                                double * A, const size_t ldA
                                , bool positive  );


    void
    freduce_2_modular_double (const double p, const size_t m , const size_t n,
                              const double * B, const size_t ldB,
                              double * A, const size_t ldA
                              , bool positive  );

    void
    fnegin_2_modular_double (const double p, const size_t m , const size_t n,
                             double * A, const size_t ldA
                             , bool positive  );


    void
    fneg_2_modular_double (const double p, const size_t m , const size_t n,
                           const double * B, const size_t ldB,
                           double * A, const size_t ldA
                           , bool positive  );


    void
    fscalin_2_modular_double (const double p, const size_t m , const size_t n,
                              const double alpha,
                              double * A, const size_t ldA
                              , bool positive  );


    void
    fscal_2_modular_double (const double p, const size_t m , const size_t n,
                            const double alpha,
                            const double * A, const size_t ldA,
                            double * B, const size_t ldB
                            , bool positive  );


    void
    faxpy_2_modular_double (const double p, const size_t m, const size_t n
                            , const double alpha,
                            const double * X, const size_t ldX,
                            double * Y, const size_t ldY
                            , bool positive  );


#if 0
    void
    faxpby_2_modular_double (const double p, const size_t m, const size_t n,
                             const double alpha,
                             const double * X, const size_t ldX,
                             const double betA,
                             double * Y, const size_t ldY
                             , bool positive  );
#endif


    void
    fmove_2_modular_double (const double p, const size_t m, const size_t n,
                            double * A, const size_t ldA,
                            double * B, const size_t ldB
                            , bool positive  );


    void
    fadd_2_modular_double (const double p, const size_t m, const size_t n,
                           const double * A, const size_t ldA,
                           const double * B, const size_t ldB,
                           double * C, const size_t ldC
                           , bool positive  );



    void
    fsub_2_modular_double (const double p, const size_t m, const size_t n,
                           const double * A, const size_t ldA,
                           const double * B, const size_t ldB,
                           double * C, const size_t ldC
                           , bool positive  );


    void
    fsubin_2_modular_double (const double p, const size_t m, const size_t n,
                             const double * B, const size_t ldB,
                             double * C, const size_t ldC
                             , bool positive  );



    void
    faddin_2_modular_double (const double p, const size_t m, const size_t n,
                             const double * B, const size_t ldB,
                             double * C, const size_t ldC
                             , bool positive  );



    double *
    fgemv_2_modular_double (const double p, const enum FFLAS_C_TRANSPOSE TransA,
                            const size_t m, const size_t n,
                            const double alpha,
                            const double * A, const size_t ldA,
                            const double * X, const size_t incX,
                            const  double betA,
                            double * Y, const size_t incY
                            , bool positive  );


    void
    fger_2_modular_double (const double p, const size_t m, const size_t n,
                           const double alpha,
                           const double * x, const size_t incX,
                           const double * y, const size_t incY,
                           double * A, const size_t ldA
                           , bool positive  );


    void
    ftrsv_2_modular_double (const double p, const enum FFLAS_C_UPLO Uplo,
                            const enum FFLAS_C_TRANSPOSE TransA, const enum FFLAS_C_DIAG Diag,
                            const size_t n,const double * A, const size_t ldA,
                            double * X, int incX
                            , bool positive  );

    /* ******** *
     * LEVEL2.5 *
     * ******** */

    // fspmm

    /* ******** *
     * LEVEL3   *
     * ******** */


    void
    ftrsm_3_modular_double (const double p, const enum FFLAS_C_SIDE Side,
                            const enum FFLAS_C_UPLO Uplo,
                            const enum FFLAS_C_TRANSPOSE TransA,
                            const enum FFLAS_C_DIAG Diag,
                            const size_t m, const size_t n,
                            const double alpha,
                            const double * A,
                            const size_t ldA,
                            double * B, const size_t ldB
                            , bool positive  );


    void
    ftrmm_3_modular_double (const double p, const enum FFLAS_C_SIDE Side,
                            const enum FFLAS_C_UPLO Uplo,
                            const enum FFLAS_C_TRANSPOSE TransA,
                            const enum FFLAS_C_DIAG Diag,
                            const size_t m, const size_t n,
                            const double alpha,
                            double * A, const size_t ldA,
                            double * B, const size_t ldB
                            , bool positive  );


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
                            double * C, const size_t ldC
                            , bool positive  );


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
                              , bool positive  );

#ifdef __cplusplus
}
#endif

#endif // __FFLASFFPACK_interfaces_libs_fflas_c_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
