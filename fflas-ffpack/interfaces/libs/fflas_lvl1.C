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
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/** @file fflas_lvl1.C
 * @author  Brice Boyer
 * @brief C functions calls for level 1 FFLAS in fflas-c.h
 * @see fflas/fflas_level1.inl
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
    /*
     * level 1
     */

    void
    freducein_1_modular_double (const double p, const size_t n,
                                double * X, const size_t incX
                                , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            freduce(F,n,X,incX);
        } else {
            ModularBalanced<double> F(p);
            freduce(F,n,X,incX);
        }
    }

    void
    freduce_1_modular_double (const double p, const size_t n,
                              const double * Y, const size_t incY,
                              double * X, const size_t incX
                              , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            freduce(F,n,Y,incY,X,incX);
        } else {
            ModularBalanced<double> F(p);
            freduce(F,n,Y,incY,X,incX);
        }
    }


    void
    fnegin_1_modular_double (const double p, const size_t n,
                             double * X, const size_t incX
                             , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            fnegin(F,n,X,incX);
        } else {
            ModularBalanced<double> F(p);
            fnegin(F,n,X,incX);
        }
    }


    void
    fneg_1_modular_double (const double p, const size_t n,
                           const double * Y, const size_t incY,
                           double * X, const size_t incX
                           , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            fneg(F,n,Y,incY,X,incX);
        } else {
            ModularBalanced<double> F(p);
            fneg(F,n,Y,incY,X,incX);
        }
    }

    void
    fzero_1_modular_double (const double p, const size_t n,
                            double * X, const size_t incX
                            , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            fzero(F,n,X,incX);
        } else {
            ModularBalanced<double> F(p);
            fzero(F,n,X,incX);
        }
    }


    bool
    fiszero_1_modular_double (const double p, const size_t n,
                              const double * X, const size_t incX
                              , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            return fiszero(F,n,X,incX);
        } else {
            ModularBalanced<double> F(p);
            return fiszero(F,n,X,incX);
        }
    }

    bool
    fequal_1_modular_double (const double p, const size_t n,
                             const double * X, const size_t incX,
                             const double * Y, const size_t incY
                             , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            return fequal(F,n,Y,incY,X,incX);
        } else {
            ModularBalanced<double> F(p);
            return fequal(F,n,Y,incY,X,incX);
        }
    }


    void
    fassign_1_modular_double (const double p, const size_t n,
                              const double * Y, const size_t incY ,
                              double * X, const size_t incX
                              , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            fassign(F,n,Y,incY,X,incX);
        } else {
            ModularBalanced<double> F(p);
            fassign(F,n,Y,incY,X,incX);
        }
    }


    void
    fscalin_1_modular_double (const double p, const size_t n, const double alpha,
                              double * X, const size_t incX
                              , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            fscalin(F,n,alpha,X,incX);
        } else {
            ModularBalanced<double> F(p);
            fscalin(F,n,alpha,X,incX);
        }
    }



    void
    fscal_1_modular_double (const double p, const size_t n
                            , const double alpha
                            , const double * X, const size_t incX
                            , double * Y, const size_t incY
                            , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            fscal(F,n,alpha,X,incX,Y,incY);
        } else {
            ModularBalanced<double> F(p);
            fscal(F,n,alpha,X,incX,Y,incY);
        }
    }


    void
    faxpy_1_modular_double (const double p, const size_t n,
                            const double alpha,
                            const double * X, const size_t incX,
                            double * Y, const size_t incY
                            , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            faxpy(F,n,alpha,X,incX,Y,incY);
        } else {
            ModularBalanced<double> F(p);
            faxpy(F,n,alpha,X,incX,Y,incY);
        }
    }

#if 0
    void
    faxpby_1_modular_double (const double p, const size_t n,
                             const double alpha,
                             const double * X, const size_t incX,
                             const double beta,
                             double * Y, const size_t incY
                             , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            faxpby(F,n,alpha,X,incX,beta,Y,incY);
        } else {
            ModularBalanced<double> F(p);
            faxpby(F,n,alpha,X,incX,beta,Y,incY);
        }
    }
#endif


    double
    fdot_1_modular_double (const double p, const size_t n,
                           const double * X, const size_t incX,
                           const double * Y, const size_t incY
                           , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            return fdot(F,n,Y,incY,X,incX);
        } else {
            ModularBalanced<double> F(p);
            return fdot(F,n,Y,incY,X,incX);
        }
    }


    void
    fswap_1_modular_double (const double p, const size_t n,
                            double * X, const size_t incX,
                            double * Y, const size_t incY
                            , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            fswap(F,n,Y,incY,X,incX);
        } else {
            ModularBalanced<double> F(p);
            fswap(F,n,Y,incY,X,incX);
        }
    }


    void
    fadd_1_modular_double (const double p,  const size_t n,
                           const double * A, const size_t incA,
                           const double * B, const size_t incB,
                           double * C, const size_t incC
                           , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            fadd(F,n,A,incA,B,incB,C,incC);
        } else {
            ModularBalanced<double> F(p);
            fadd(F,n,A,incA,B,incB,C,incC);
        }
    }

    void
    fsub_1_modular_double (const double p,  const size_t n,
                           const double * A, const size_t incA,
                           const double * B, const size_t incB,
                           double * C, const size_t incC
                           , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            fsub(F,n,A,incA,B,incB,C,incC);
        } else {
            ModularBalanced<double> F(p);
            fsub(F,n,A,incA,B,incB,C,incC);
        }
    }

    void
    faddin_1_modular_double (const double p,  const size_t n,
                             const double * B, const size_t incB,
                             double * C, const size_t incC
                             , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            faddin(F,n,B,incB,C,incC);
        } else {
            ModularBalanced<double> F(p);
            faddin(F,n,B,incB,C,incC);
        }
    }

    void
    fsubin_1_modular_double (const double p,  const size_t n,
                             const double * B, const size_t incB,
                             double * C, const size_t incC
                             , bool positive )
    {
        if (positive) {
            Modular<double> F(p);
            fsubin(F,n,B,incB,C,incC);
        } else {
            ModularBalanced<double> F(p);
            fsubin(F,n,B,incB,C,incC);
        }
    }

#ifdef __cplusplus
}
#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
