/*
 * Copyright_2_modular_double (C) 2015 FFLAS-FFPACK
 *
 * Written by Brice Boyer_2_modular_double (briceboyer) <boyer.brice@gmail.com>
 *
 *
 * ========LICENCE========
 * This file is part of the library FFLAS-FFPACK.
 *
 * FFLAS-FFPACK is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or_2_modular_double (at your option) any later version.
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

/** @file fflas_lvl2.C
 * @author  Brice Boyer
 * @brief C functions calls for level 2 FFLAS in fflas-c.h
 * @see fflas/fflas_level2.inl
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
    fassign_2_modular_double (const double p, const size_t m, const size_t n,
                              const double * A, const size_t lda ,
                              double * B, const size_t ldb
                              , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            fassign(F,m,n,A,lda,B,ldb);
        } else {
            ModularBalanced<double> F(p);
            fassign(F,m,n,A,lda,B,ldb);
        }
    }



    void
    fzero_2_modular_double (const double p, const size_t m, const size_t n,
                            double * A, const size_t lda
                            , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            fzero(F,m,n,A,lda);
        } else {
            ModularBalanced<double> F(p);
            fzero(F,m,n,A,lda);
        }
    }

    bool
    fequal_2_modular_double (const double p, const size_t m, const size_t n,
                             const double * A, const size_t lda,
                             const double * B, const size_t ldb
                             , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            return fequal(F,m,n,A,lda,B,ldb);
        } else {
            ModularBalanced<double> F(p);
            return fequal(F,m,n,A,lda,B,ldb);
        }
    }

    bool
    fiszero_2_modular_double (const double p, const size_t m, const size_t n,
                              const double * A, const size_t lda
                              , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            return fiszero(F,m,n,A,lda);
        } else {
            ModularBalanced<double> F(p);
            return	fiszero(F,m,n,A,lda);
        }
    }

    void
    fidentity_2_modular_double (const double p, const size_t m, const size_t n,
                                double * A, const size_t lda,
                                const double  d
                                , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            fidentity(F,m,n,A,lda,d);
        } else {
            ModularBalanced<double> F(p);
            fidentity(F,m,n,A,lda,d);
        }
    }


    void
    freducein_2_modular_double (const double p, const size_t m , const size_t n,
                                double * A, const size_t lda
                                , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            freduce(F,m,n,A,lda);
        } else {
            ModularBalanced<double> F(p);
            freduce(F,m,n,A,lda);
        }
    }

    void
    freduce_2_modular_double (const double p, const size_t m , const size_t n,
                              const double * A, const size_t lda,
                              double * B, const size_t ldb
                              , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            freduce(F,m,n,A,lda,B,ldb);
        } else {
            ModularBalanced<double> F(p);
            freduce(F,m,n,A,lda,B,ldb);
        }
    }
    void
    fnegin_2_modular_double (const double p, const size_t m , const size_t n,
                             double * A, const size_t lda
                             , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            fnegin(F,m,n,A,lda);
        } else {
            ModularBalanced<double> F(p);
            fnegin(F,m,n,A,lda);
        }
    }

    void
    fneg_2_modular_double (const double p, const size_t m , const size_t n,
                           const double * A, const size_t lda,
                           double * B, const size_t ldb
                           , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            fneg(F,m,n,A,lda,B,ldb);
        } else {
            ModularBalanced<double> F(p);
            fneg(F,m,n,A,lda,B,ldb);
        }
    }

    void
    fscalin_2_modular_double (const double p, const size_t m , const size_t n,
                              const double alpha,
                              double * A, const size_t lda
                              , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            fscalin(F,m,n,alpha,A,lda);
        } else {
            ModularBalanced<double> F(p);
            fscalin(F,m,n,alpha,A,lda);
        }
    }

    void
    fscal_2_modular_double (const double p, const size_t m , const size_t n,
                            const double alpha,
                            const double * A, const size_t lda,
                            double * B, const size_t ldb
                            , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            fscal(F,m,n,alpha,A,lda,B,ldb);
        } else {
            ModularBalanced<double> F(p);
            fscal(F,m,n,alpha,A,lda,B,ldb);
        }
    }

    void
    faxpy_2_modular_double (const double p, const size_t m, const size_t n
                            , const double alpha,
                            const double * A, const size_t lda,
                            double * B, const size_t ldb
                            , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            faxpy(F,m,n,alpha,A,lda,B,ldb);
        } else {
            ModularBalanced<double> F(p);
            faxpy(F,m,n,alpha,A,lda,B,ldb);
        }
    }

#if 0
    void
    faxpby_2_modular_double (const double p, const size_t m, const size_t n,
                             const double alpha,
                             const double * A, const size_t lda,
                             const double beta,
                             double * B, const size_t ldb
                             , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            faxpby(F,m,n,alpha,A,lda,beta,B,ldb);
        } else {
            ModularBalanced<double> F(p);
            faxpby(F,m,n,alpha,A,lda,beta,B,ldb);
        }
    }
#endif

    void
    fmove_2_modular_double (const double p, const size_t m, const size_t n,
                            double * A, const size_t lda,
                            double * B, const size_t ldb
                            , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            fmove(F,m,n,A,lda,B,ldb);
        } else {
            ModularBalanced<double> F(p);
            fmove(F,m,n,A,lda,B,ldb);
        }
    }

    void
    fadd_2_modular_double (const double p, const size_t m, const size_t n,
                           const double * A, const size_t lda,
                           const double * B, const size_t ldb,
                           double * C, const size_t ldc
                           , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            fadd(F,m,n,A,lda,B,ldb,C,ldc);
        } else {
            ModularBalanced<double> F(p);
            fadd(F,m,n,A,lda,B,ldb,C,ldc);
        }
    }


    void
    fsub_2_modular_double (const double p, const size_t m, const size_t n,
                           const double * A, const size_t lda,
                           const double * B, const size_t ldb,
                           double * C, const size_t ldc
                           , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            fsub(F,m,n,A,lda,B,ldb,C,ldc);
        } else {
            ModularBalanced<double> F(p);
            fsub(F,m,n,A,lda,B,ldb,C,ldc);
        }
    }

    void
    fsubin_2_modular_double (const double p, const size_t m, const size_t n,
                             const double * B, const size_t ldb,
                             double * C, const size_t ldc
                             , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            fsubin(F,m,n,B,ldb,C,ldc);
        } else {
            ModularBalanced<double> F(p);
            fsubin(F,m,n,B,ldb,C,ldc);
        }
    }

    void
    faddin_2_modular_double (const double p, const size_t m, const size_t n,
                             const double * B, const size_t ldb,
                             double * C, const size_t ldc
                             , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            faddin(F,m,n,B,ldb,C,ldc);
        } else {
            ModularBalanced<double> F(p);
            faddin(F,m,n,B,ldb,C,ldc);
        }
    }


    double *
    fgemv_2_modular_double (const double p, const enum FFLAS_C_TRANSPOSE TransA,
                            const size_t m, const size_t n,
                            const double alpha,
                            const double * A, const size_t lda,
                            const double * X, const size_t incX,
                            const  double beta,
                            double * Y, const size_t incY
                            , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            return fgemv(F,(enum FFLAS::FFLAS_TRANSPOSE)TransA,m,n,alpha,A,lda,X,incX,beta,Y,incY);
        } else {
            ModularBalanced<double> F(p);
            return fgemv(F,(enum FFLAS::FFLAS_TRANSPOSE)TransA,m,n,alpha,A,lda,X,incX,beta,Y,incY);
        }
        return nullptr;
    }

    void
    fger_2_modular_double (const double p, const size_t m, const size_t n,
                           const double alpha,
                           const double * X, const size_t incX,
                           const double * Y, const size_t incY,
                           double * A, const size_t lda
                           , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            fger(F,m,n,alpha,X,incX,Y,incY,A,lda);
        } else {
            ModularBalanced<double> F(p);
            fger(F,m,n,alpha,X,incX,Y,incY,A,lda);
        }
    }

    void
    ftrsv_2_modular_double (const double p, const enum FFLAS_C_UPLO Uplo,
                            const enum FFLAS_C_TRANSPOSE TransA, const enum FFLAS_C_DIAG Diag,
                            const size_t n,const double * A, const size_t lda,
                            double * X, int incX
                            , bool positive  )
    {
        if (positive) {
            Modular<double> F(p);
            ftrsv(F,(enum FFLAS::FFLAS_UPLO)Uplo,(enum FFLAS::FFLAS_TRANSPOSE)TransA,(enum FFLAS::FFLAS_DIAG)Diag,n,A,lda,X,incX);
        } else {
            ModularBalanced<double> F(p);
            ftrsv(F,(enum FFLAS::FFLAS_UPLO)Uplo,(enum FFLAS::FFLAS_TRANSPOSE)TransA,(enum FFLAS::FFLAS_DIAG)Diag,n,A,lda,X,incX);
        }
    }

#ifdef __cplusplus
}
#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
