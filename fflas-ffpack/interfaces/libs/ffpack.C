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

/** @file ffpack.C
 * @author  Brice Boyer
 * @brief C functions calls for FFPACK in ffpack-c.h
 * @see ffpack/ffpack.h
 */

#include "fflas-ffpack/interfaces/libs/ffpack_c.h"

#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "givaro//modular-balanced.h"
#include "givaro//modular.h"

using Givaro::Modular ;
using Givaro::ModularBalanced ;
using namespace FFLAS ;
using namespace FFPACK;

/*****************/
/* PERMUTATIONS  */
/*****************/


void LAPACKPerm2MathPerm (size_t * MathP, const size_t * LapackP,
                          const size_t N)
{
    FFPACK::LAPACKPerm2MathPerm(MathP,LapackP,N);
}

void MathPerm2LAPACKPerm (size_t * LapackP, const size_t * MathP,
                          const size_t N)
{
    FFPACK::MathPerm2LAPACKPerm(LapackP, MathP, N);
}

void MatrixApplyS_modular_double (const double p, double * A, const size_t lda, const size_t width,
                                  const size_t M2,
                                  const size_t R1, const size_t R2,
                                  const size_t R3, const size_t R4
                                  , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        MatrixApplyS(F,A,lda,width,M2,R1,R2,R3,R4);
    } else {
        ModularBalanced<double> F(p);
        MatrixApplyS(F,A,lda,width,M2,R1,R2,R3,R4);
    }
}

void PermApplyS_double (double * A, const size_t lda, const size_t width,
                        const size_t M2,
                        const size_t R1, const size_t R2,
                        const size_t R3, const size_t R4)
{
    PermApplyS<double>(A,lda,width,M2,R1,R2,R3,R4);
}


void MatrixApplyT_modular_double (const double p, double * A, const size_t lda, const size_t width,
                                  const size_t N2,
                                  const size_t R1, const size_t R2,
                                  const size_t R3, const size_t R4
                                  , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        MatrixApplyT(F,A,lda,width,N2,R1,R2,R3,R4);
    } else {
        ModularBalanced<double> F(p);
        MatrixApplyT(F,A,lda,width,N2,R1,R2,R3,R4);
    }
}


void PermApplyT_double (double * A, const size_t lda, const size_t width,
                        const size_t N2,
                        const size_t R1, const size_t R2,
                        const size_t R3, const size_t R4)
{
    PermApplyT<double>(A,lda,width,N2,R1,R2,R3,R4);
}

void composePermutationsLLM (size_t * MathP,
                             const size_t * P1,
                             const size_t * P2,
                             const size_t R, const size_t N)
{
    FFPACK::composePermutationsLLM(MathP,P1,P2,R,N);
}


void composePermutationsLLL (size_t * P1,
                             const size_t * P2,
                             const size_t R, const size_t N)
{
    FFPACK::composePermutationsLLL(P1,P2,R,N);
}

void composePermutationsMLM (size_t * MathP1,
                             const size_t * P2,
                             const size_t R, const size_t N)
{
    FFPACK::composePermutationsMLM(MathP1,P2,R,N);
}

void cyclic_shift_mathPerm (size_t * P,  const size_t s)
{
    FFPACK::cyclic_shift_mathPerm(P,s);
}

#if 0
template<typename Base_t>
void cyclic_shift_row_col(Base_t * A, size_t m, size_t n, size_t lda);
#endif


void cyclic_shift_row_modular_double(const double p, double * A, size_t m, size_t n, size_t lda
                                     , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        cyclic_shift_row(F,A,m,n,lda);
    } else {
        ModularBalanced<double> F(p);
        cyclic_shift_row(F,A,m,n,lda);
    }
}


void cyclic_shift_col_modular_double(const double p, double * A, size_t m, size_t n, size_t lda
                                     , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        cyclic_shift_col(F,A,m,n,lda);
    } else {
        ModularBalanced<double> F(p);
        cyclic_shift_col(F,A,m,n,lda);
    }
}




void
applyP_modular_double( const double p,
                       const enum FFLAS::FFLAS_SIDE Side,
                       const enum FFLAS::FFLAS_TRANSPOSE Trans,
                       const size_t M, const size_t ibeg, const size_t iend,
                       double * A, const size_t lda, const size_t * P
                       , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        applyP(F,(enum FFLAS::FFLAS_SIDE)Side,(enum FFLAS::FFLAS_TRANSPOSE)Trans,M,ibeg,iend,A,lda,P);
    } else {
        ModularBalanced<double> F(p);
        applyP(F,(enum FFLAS::FFLAS_SIDE)Side,(enum FFLAS::FFLAS_TRANSPOSE)Trans,M,ibeg,iend,A,lda,P);
    }
}

/* fgetrs, fgesv */

void
fgetrsin_modular_double (const double p,
                         const enum FFLAS::FFLAS_SIDE Side,
                         const size_t M, const size_t N, const size_t R,
                         double * A, const size_t lda,
                         const size_t *P, const size_t *Q,
                         double * B, const size_t ldb,
                         int * info
                         , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        fgetrs(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,R,A,lda,P,Q,B,ldb,info);
    } else {
        ModularBalanced<double> F(p);
        fgetrs(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,R,A,lda,P,Q,B,ldb,info);
    }
}


double *
fgetrsv_modular_double (const double p,
                        const enum FFLAS::FFLAS_SIDE Side,
                        const size_t M, const size_t N, const size_t NRHS, const size_t R,
                        double * A, const size_t lda,
                        const size_t *P, const size_t *Q,
                        double * X, const size_t ldx,
                        const double * B, const size_t ldb,
                        int * info
                        , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return fgetrs(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,NRHS,R,A,lda,P,Q,X,ldx,B,ldb,info);
    } else {
        ModularBalanced<double> F(p);
        return fgetrs(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,NRHS,R,A,lda,P,Q,X,ldx,B,ldb,info);
    }
}



size_t
fgesvin_modular_double (const double p,
                        const enum FFLAS::FFLAS_SIDE Side,
                        const size_t M, const size_t N,
                        double * A, const size_t lda,
                        double * B, const size_t ldb,
                        int * info
                        , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return fgesv(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,A,lda,B,ldb,info);
    } else {
        ModularBalanced<double> F(p);
        return fgesv(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,A,lda,B,ldb,info);
    }
}



size_t
fgesv_modular_double (const double p,
                      const enum FFLAS::FFLAS_SIDE Side,
                      const size_t M, const size_t N, const size_t NRHS,
                      double * A, const size_t lda,
                      double * X, const size_t ldx,
                      const double * B, const size_t ldb,
                      int * info
                      , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return fgesv(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,NRHS,A,lda,X,ldx,B,ldb,info);
    } else {
        ModularBalanced<double> F(p);
        return fgesv(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,NRHS,A,lda,X,ldx,B,ldb,info);
    }
}

/* ftrtr */


void
ftrtri_modular_double (const double p, const enum FFLAS::FFLAS_UPLO Uplo, const enum FFLAS::FFLAS_DIAG Diag,
                       const size_t N, double * A, const size_t lda
                       , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        ftrtri(F,(enum FFLAS::FFLAS_UPLO)Uplo,(enum FFLAS::FFLAS_DIAG)Diag,N,A,lda);
    } else {
        ModularBalanced<double> F(p);
        ftrtri(F,(enum FFLAS::FFLAS_UPLO)Uplo,(enum FFLAS::FFLAS_DIAG)Diag,N,A,lda);
    }
}


void trinv_left_modular_double( const double p, const size_t N, const double * L, const size_t ldl,
                                double * X, const size_t ldx
                                , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        trinv_left(F,N,L,ldl,X,ldx);
    } else {
        ModularBalanced<double> F(p);
        trinv_left(F,N,L,ldl,X,ldx);
    }
}

void
ftrtrm_modular_double (const double p, const FFLAS::FFLAS_SIDE side, const enum FFLAS::FFLAS_DIAG Diag,
                       const size_t N, double * A, const size_t lda, bool positive)
{
    if (positive) {
        Modular<double> F(p);
        ftrtrm(F,side,Diag,N,A,lda);
    } else {
        ModularBalanced<double> F(p);
        ftrtrm(F,side,Diag,N,A,lda);
    }
}



/* PLUQ */

size_t
PLUQ_modular_double (const double p, const enum FFLAS::FFLAS_DIAG Diag,
                     const size_t M, const size_t N,
                     double * A, const size_t lda,
                     size_t*P, size_t *Q
                     , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return PLUQ(F,(enum FFLAS::FFLAS_DIAG)Diag,M,N,A,lda,P,Q);
    } else {
        ModularBalanced<double> F(p);
        return PLUQ(F,(enum FFLAS::FFLAS_DIAG)Diag,M,N,A,lda,P,Q);
    }
}

size_t
LUdivine_modular_double (const double p, const enum FFLAS::FFLAS_DIAG Diag,  const enum FFLAS::FFLAS_TRANSPOSE Trans,
                         const size_t M, const size_t N,
                         double * A, const size_t lda,
                         size_t* P, size_t* Qt,
                         const enum FFPACK_C_LU_TAG LuTag,
                         const size_t cutoff
                         , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return LUdivine(F,(enum FFLAS::FFLAS_DIAG)Diag,(enum FFLAS::FFLAS_TRANSPOSE)Trans,M,N,A,lda,P,Qt,(enum FFPACK::FFPACK_LU_TAG)LuTag,cutoff);
    } else {
        ModularBalanced<double> F(p);
        return LUdivine(F,(enum FFLAS::FFLAS_DIAG)Diag,(enum FFLAS::FFLAS_TRANSPOSE)Trans,M,N,A,lda,P,Qt,(enum FFPACK::FFPACK_LU_TAG)LuTag,cutoff);
    }
}


#if 0 /*  UTILE ?? */

size_t
LUdivine_small_modular_double (const double p, const enum FFLAS::FFLAS_DIAG Diag,  const enum FFLAS::FFLAS_TRANSPOSE Trans,
                               const size_t M, const size_t N,
                               double * A, const size_t lda,
                               size_t* P, size_t* Q,
                               const enum FFPACK_C_LU_TAG LuTag);


size_t
LUdivine_gauss_modular_double (const double p, const enum FFLAS::FFLAS_DIAG Diag,
                               const size_t M, const size_t N,
                               double * A, const size_t lda,
                               size_t* P, size_t* Q,
                               const enum FFPACK_C_LU_TAG LuTag);
#endif



/*****************/
/* ECHELON FORMS */
/*****************/


size_t
ColumnEchelonForm_modular_double (const double p, const size_t M, const size_t N,
                                  double * A, const size_t lda,
                                  size_t* P, size_t* Qt, bool transform ,
                                  const enum FFPACK_C_LU_TAG LuTag
                                  , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return ColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        return ColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}


size_t
RowEchelonForm_modular_double (const double p, const size_t M, const size_t N,
                               double * A, const size_t lda,
                               size_t* P, size_t* Qt, const bool transform,
                               const enum FFPACK_C_LU_TAG LuTag
                               , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return RowEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        return RowEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}



size_t
ReducedColumnEchelonForm_modular_double (const double p, const size_t M, const size_t N,
                                         double * A, const size_t lda,
                                         size_t* P, size_t* Qt, const bool transform,
                                         const enum FFPACK_C_LU_TAG LuTag
                                         , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return ReducedColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        return ReducedColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}



size_t
ReducedRowEchelonForm_modular_double (const double p, const size_t M, const size_t N,
                                      double * A, const size_t lda,
                                      size_t* P, size_t* Qt, const bool transform,
                                      const enum FFPACK_C_LU_TAG LuTag
                                      , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return ReducedRowEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        return ReducedRowEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}

size_t
ColumnEchelonForm_modular_float (const float p, const size_t M, const size_t N,
                                 float * A, const size_t lda,
                                 size_t* P, size_t* Qt, bool transform ,
                                 const enum FFPACK_C_LU_TAG LuTag
                                 , bool positive)
{
    if (positive) {
        Modular<float> F(p);
        return ColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<float> F(p);
        return ColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}


size_t
RowEchelonForm_modular_float (const float p, const size_t M, const size_t N,
                              float * A, const size_t lda,
                              size_t* P, size_t* Qt, const bool transform,
                              const enum FFPACK_C_LU_TAG LuTag
                              , bool positive)
{
    if (positive) {
        Modular<float> F(p);
        return RowEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<float> F(p);
        return RowEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}



size_t
ReducedColumnEchelonForm_modular_float (const float p, const size_t M, const size_t N,
                                        float * A, const size_t lda,
                                        size_t* P, size_t* Qt, const bool transform,
                                        const enum FFPACK_C_LU_TAG LuTag
                                        , bool positive)
{
    if (positive) {
        Modular<float> F(p);
        return ReducedColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<float> F(p);
        return ReducedColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}



size_t
ReducedRowEchelonForm_modular_float (const float p, const size_t M, const size_t N,
                                     float * A, const size_t lda,
                                     size_t* P, size_t* Qt, const bool transform,
                                     const enum FFPACK_C_LU_TAG LuTag
                                     , bool positive)
{
    if (positive) {
        Modular<float> F(p);
        return ReducedRowEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<float> F(p);
        return ReducedRowEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}

size_t
ColumnEchelonForm_modular_int32_t (const int32_t p, const size_t M, const size_t N,
                                   int32_t * A, const size_t lda,
                                   size_t* P, size_t* Qt, bool transform ,
                                   const enum FFPACK_C_LU_TAG LuTag
                                   , bool positive)
{
    if (positive) {
        Modular<int32_t> F(p);
        return ColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<int32_t> F(p);
        return ColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}


size_t
RowEchelonForm_modular_int32_t (const int32_t p, const size_t M, const size_t N,
                                int32_t * A, const size_t lda,
                                size_t* P, size_t* Qt, const bool transform,
                                const enum FFPACK_C_LU_TAG LuTag
                                , bool positive)
{
    if (positive) {
        Modular<int32_t> F(p);
        return RowEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<int32_t> F(p);
        return RowEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}



size_t
ReducedColumnEchelonForm_modular_int32_t (const int32_t p, const size_t M, const size_t N,
                                          int32_t * A, const size_t lda,
                                          size_t* P, size_t* Qt, const bool transform,
                                          const enum FFPACK_C_LU_TAG LuTag
                                          , bool positive)
{
    if (positive) {
        Modular<int32_t> F(p);
        return ReducedColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<int32_t> F(p);
        return ReducedColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}



size_t
ReducedRowEchelonForm_modular_int32_t (const int32_t p, const size_t M, const size_t N,
                                       int32_t * A, const size_t lda,
                                       size_t* P, size_t* Qt, const bool transform,
                                       const enum FFPACK_C_LU_TAG LuTag
                                       , bool positive)
{
    if (positive) {
        Modular<int32_t> F(p);
        return ReducedRowEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<int32_t> F(p);
        return ReducedRowEchelonForm(F,M,N,A,lda,P,Qt,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}

size_t
pColumnEchelonForm_modular_double (const double p, const size_t M, const size_t N,
                                  double * A, const size_t lda,
                                  size_t* P, size_t* Qt, bool transform ,
                                  const enum FFPACK_C_LU_TAG LuTag
                                  , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return pColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,0,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        return pColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,0,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}


size_t
pRowEchelonForm_modular_double (const double p, const size_t M, const size_t N,
                               double * A, const size_t lda,
                               size_t* P, size_t* Qt, const bool transform,
                               const enum FFPACK_C_LU_TAG LuTag
                               , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return pRowEchelonForm(F,M,N,A,lda,P,Qt,transform,0,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        return pRowEchelonForm(F,M,N,A,lda,P,Qt,transform,0,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}



size_t
pReducedColumnEchelonForm_modular_double (const double p, const size_t M, const size_t N,
                                         double * A, const size_t lda,
                                         size_t* P, size_t* Qt, const bool transform,
                                         const enum FFPACK_C_LU_TAG LuTag
                                         , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return pReducedColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,0,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        return pReducedColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,0,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}



size_t
pReducedRowEchelonForm_modular_double (const double p, const size_t M, const size_t N,
                                      double * A, const size_t lda,
                                      size_t* P, size_t* Qt, const bool transform,
                                      const enum FFPACK_C_LU_TAG LuTag
                                      , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return pReducedRowEchelonForm(F,M,N,A,lda,P,Qt,transform,0,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        return pReducedRowEchelonForm(F,M,N,A,lda,P,Qt,transform,0,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}

size_t
pColumnEchelonForm_modular_float (const float p, const size_t M, const size_t N,
                                 float * A, const size_t lda,
                                 size_t* P, size_t* Qt, bool transform ,
                                 const enum FFPACK_C_LU_TAG LuTag
                                 , bool positive)
{
    if (positive) {
        Modular<float> F(p);
        return pColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,0,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<float> F(p);
        return pColumnEchelonForm(F,M,N,A,lda,P,Qt,transform,0,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}


size_t
pRowEchelonForm_modular_float (const float p, const size_t M, const size_t N,
                              float * A, const size_t lda,
                              size_t* P, size_t* Qt, const bool transform,
                              const enum FFPACK_C_LU_TAG LuTag
                              , bool positive)
{
    if (positive) {
        Modular<float> F(p);
        return pRowEchelonForm(F,M,N,A,lda,P,Qt,transform,0,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<float> F(p);
        return pRowEchelonForm(F,M,N,A,lda,P,Qt,transform,0,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}



size_t
pReducedColumnEchelonForm_modular_float (const float p, const size_t M, const size_t N,
                                        float * A, const size_t lda,
                                        size_t* P, size_t* Qt, const bool transform,
                                        const enum FFPACK_C_LU_TAG LuTag
                                        , bool positive)
{
    if (positive) {
        Modular<float> F(p);
        return pReducedColumnEchelonForm(F,M,N,A,lda,P,Qt,0,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<float> F(p);
        return pReducedColumnEchelonForm(F,M,N,A,lda,P,Qt,0,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}



size_t
pReducedRowEchelonForm_modular_float (const float p, const size_t M, const size_t N,
                                     float * A, const size_t lda,
                                     size_t* P, size_t* Qt, const bool transform,
                                     const enum FFPACK_C_LU_TAG LuTag
                                     , bool positive)
{
    if (positive) {
        Modular<float> F(p);
        return pReducedRowEchelonForm(F,M,N,A,lda,P,Qt,0,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<float> F(p);
        return pReducedRowEchelonForm(F,M,N,A,lda,P,Qt,0,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}

size_t
pColumnEchelonForm_modular_int32_t (const int32_t p, const size_t M, const size_t N,
                                   int32_t * A, const size_t lda,
                                   size_t* P, size_t* Qt, bool transform ,
                                   const enum FFPACK_C_LU_TAG LuTag
                                   , bool positive)
{
    if (positive) {
        Modular<int32_t> F(p);
        return pColumnEchelonForm(F,M,N,A,lda,P,Qt,0,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<int32_t> F(p);
        return pColumnEchelonForm(F,M,N,A,lda,P,Qt,0,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}


size_t
pRowEchelonForm_modular_int32_t (const int32_t p, const size_t M, const size_t N,
                                int32_t * A, const size_t lda,
                                size_t* P, size_t* Qt, const bool transform,
                                const enum FFPACK_C_LU_TAG LuTag
                                , bool positive)
{
    if (positive) {
        Modular<int32_t> F(p);
        return pRowEchelonForm(F,M,N,A,lda,P,Qt,0,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<int32_t> F(p);
        return pRowEchelonForm(F,M,N,A,lda,P,Qt,0,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}



size_t
pReducedColumnEchelonForm_modular_int32_t (const int32_t p, const size_t M, const size_t N,
                                          int32_t * A, const size_t lda,
                                          size_t* P, size_t* Qt, const bool transform,
                                          const enum FFPACK_C_LU_TAG LuTag
                                          , bool positive)
{
    if (positive) {
        Modular<int32_t> F(p);
        return pReducedColumnEchelonForm(F,M,N,A,lda,P,Qt,0,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<int32_t> F(p);
        return pReducedColumnEchelonForm(F,M,N,A,lda,P,Qt,0,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}



size_t
pReducedRowEchelonForm_modular_int32_t (const int32_t p, const size_t M, const size_t N,
                                       int32_t * A, const size_t lda,
                                       size_t* P, size_t* Qt, const bool transform,
                                       const enum FFPACK_C_LU_TAG LuTag
                                       , bool positive)
{
    if (positive) {
        Modular<int32_t> F(p);
        return pReducedRowEchelonForm(F,M,N,A,lda,P,Qt,0,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<int32_t> F(p);
        return pReducedRowEchelonForm(F,M,N,A,lda,P,Qt,0,transform,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}




/*****************/
/*   INVERSION   */
/*****************/


double *
Invertin_modular_double (const double p, const size_t M,
                         double * A, const size_t lda,
                         int * nullity
                         , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return Invert(F,M,A,lda,*nullity);
    } else {
        ModularBalanced<double> F(p);
        return Invert(F,M,A,lda,*nullity);
    }
}



double *
Invert_modular_double (const double p, const size_t M,
                       const double * A, const size_t lda,
                       double * X, const size_t ldx,
                       int* nullity
                       , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return Invert(F,M,A,lda,X,ldx,*nullity);
    } else {
        ModularBalanced<double> F(p);
        return Invert(F,M,A,lda,X,ldx,*nullity);
    }
}


double *
Invert2_modular_double( const double p, const size_t M,
                        double * A, const size_t lda,
                        double * X, const size_t ldx,
                        int* nullity
                        , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return Invert2(F,M,A,lda,X,ldx,*nullity);
    } else {
        ModularBalanced<double> F(p);
        return Invert2(F,M,A,lda,X,ldx,*nullity);
    }
}


/*****************************/
/* CHARACTERISTIC POLYNOMIAL */
/*****************************/


#if 0 /*  pas pour le moment */
template <class Polynomial, class>
std::list<Polynomial>&
CharPoly( const double p, std::list<Polynomial>& charp, const size_t N,
          double * A, const size_t lda,
          const enum FFPACK_C_CHARPOLY_TAG CharpTag= FfpackArithProg);

template<class Polynomial>
Polynomial & mulpoly_modular_double(const double p, Polynomial &res, const Polynomial & P1, const Polynomial & P2);

template <class Polynomial>
Polynomial&
CharPoly_modular_double( const double p, Polynomial& charp, const size_t N,
                         double * A, const size_t lda,
                         const enum FFPACK_C_CHARPOLY_TAG CharpTag= FfpackArithProg);



template <class Polynomial>
std::list<Polynomial>&
CharpolyArithProg_modular_double (const double p, std::list<Polynomial>& frobeniusForm,
                                  const size_t N, double * A, const size_t lda, const size_t c);
#endif



/**********************/
/* MINIMAL POLYNOMIAL */
/**********************/

#if 0 /*  pas pour le moment */
template <class Polynomial>
Polynomial&
MinPoly_modular_double( const double p, Polynomial& minP, const size_t N,
                        const double * A, const size_t lda,
                        double * X, const size_t ldx, size_t* P,
                        const enum FFPACK_C_MINPOLY_TAG MinTag= FFPACK::FfpackDense,
                        const size_t kg_mc=0, const size_t kg_mb=0, const size_t kg_j=0 );
#endif


/* Krylov Elim */


size_t KrylovElim_modular_double( const double p, const size_t M, const size_t N,
                                  double * A, const size_t lda, size_t*P,
                                  size_t *Q, const size_t deg, size_t *iterates, size_t * inviterates, const size_t maxit,size_t virt
                                  , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return KrylovElim(F,M,N,A,lda,P,Q,deg,iterates,inviterates,maxit, virt);
    } else {
        ModularBalanced<double> F(p);
        return KrylovElim(F,M,N,A,lda,P,Q,deg,iterates,inviterates,maxit, virt);
    }
}



size_t  SpecRankProfile_modular_double (const double p, const size_t M, const size_t N,
                                        double * A, const size_t lda, const size_t deg, size_t *rankProfile
                                        , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return SpecRankProfile(F,M,N,A,lda,deg,rankProfile);
    } else {
        ModularBalanced<double> F(p);
        return SpecRankProfile(F,M,N,A,lda,deg,rankProfile);
    }
}



/********/
/* RANK */
/********/


size_t
Rank_modular_double( const double p, const size_t M, const size_t N,
                     double * A, const size_t lda
                     , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return Rank(F,M,N,A,lda);
    } else {
        ModularBalanced<double> F(p);
        return Rank(F,M,N,A,lda);
    }
}


/********/
/* DET  */
/********/


bool
IsSingular_modular_double( const double p, const size_t M, const size_t N,
                           double * A, const size_t lda
                           , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return IsSingular(F,M,N,A,lda);
    } else {
        ModularBalanced<double> F(p);
        return IsSingular(F,M,N,A,lda);
    }
}



double
Det_modular_double( const double p, const size_t N,
                    double * A, const size_t lda
                    , bool positive)
{
    double d;
    if (positive) {
        Modular<double> F(p);
        return Det(F,d,N,A,lda);
    } else {
        ModularBalanced<double> F(p);
        return Det(F,d,N,A,lda);
    }
}



/*********/
/* SOLVE */
/*********/



double *
Solve_modular_double( const double p, const size_t M,
                      double * A, const size_t lda,
                      double * x, const int incx,
                      const double * b, const int incb
                      , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return Solve(F,M,A,lda,x,incx,b,incb);
    } else {
        ModularBalanced<double> F(p);
        return Solve(F,M,A,lda,x,incx,b,incb);
    }
}




void
solveLB_modular_double( const double p, const enum FFLAS::FFLAS_SIDE Side,
                        const size_t M, const size_t N, const size_t R,
                        double * L, const size_t ldl,
                        const size_t * Q,
                        double * B, const size_t ldb
                        , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        solveLB(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,R,L,ldl,Q,B,ldb);
    } else {
        ModularBalanced<double> F(p);
        solveLB(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,R,L,ldl,Q,B,ldb);
    }
}



void
solveLB2_modular_double( const double p, const enum FFLAS::FFLAS_SIDE Side,
                         const size_t M, const size_t N, const size_t R,
                         double * L, const size_t ldl,
                         const size_t * Q,
                         double * B, const size_t ldb
                         , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        solveLB2(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,R,L,ldl,Q,B,ldb);
    } else {
        ModularBalanced<double> F(p);
        solveLB2(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,R,L,ldl,Q,B,ldb);
    }
}



/*************/
/* NULLSPACE */
/*************/


void RandomNullSpaceVector_modular_double (const double p, const enum FFLAS::FFLAS_SIDE Side,
                                           const size_t M, const size_t N,
                                           double * A, const size_t lda,
                                           double * X, const size_t incX
                                           , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        RandomNullSpaceVector(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,A,lda,X,incX);
    } else {
        ModularBalanced<double> F(p);
        RandomNullSpaceVector(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,A,lda,X,incX);
    }
}



size_t NullSpaceBasis_modular_double (const double p, const enum FFLAS::FFLAS_SIDE Side,
                                      const size_t M, const size_t N,
                                      double * A, const size_t lda,
                                      double ** NS, size_t* ldn,
                                      size_t * NSdim
                                      , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return NullSpaceBasis(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,A,lda,*NS,*ldn,*NSdim);
    } else {
        ModularBalanced<double> F(p);
        return NullSpaceBasis(F,(enum FFLAS::FFLAS_SIDE)Side,M,N,A,lda,*NS,*ldn,*NSdim);
    }
}


/*****************/
/* RANK PROFILES */
/*****************/


size_t RowRankProfile_modular_double (const double p, const size_t M, const size_t N,
                                      double * A, const size_t lda,
                                      size_t ** rkprofile,
                                      const enum FFPACK_C_LU_TAG LuTag
                                      , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return RowRankProfile(F,M,N,A,lda,*rkprofile,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        return RowRankProfile(F,M,N,A,lda,*rkprofile,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}




size_t ColumnRankProfile_modular_double (const double p, const size_t M, const size_t N,
                                         double * A, const size_t lda,
                                         size_t ** rkprofile,
                                         const enum FFPACK_C_LU_TAG LuTag
                                         , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return ColumnRankProfile(F,M,N,A,lda,*rkprofile,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        return ColumnRankProfile(F,M,N,A,lda,*rkprofile,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}

void RankProfileFromLU (const size_t* P, const size_t N, const size_t R,
                        size_t* rkprofile, const enum FFPACK_C_LU_TAG LuTag)
{
    FFPACK::RankProfileFromLU(P,N,R,rkprofile,(enum FFPACK::FFPACK_LU_TAG)LuTag);
}

size_t LeadingSubmatrixRankProfiles (const size_t M, const size_t N, const size_t R,
                                     const size_t LSm, const size_t LSn,
                                     const size_t* P, const size_t* Q,
                                     size_t* RRP, size_t* CRP)
{
    return FFPACK::LeadingSubmatrixRankProfiles(M,N,R,LSm,LSn,P,Q,RRP,CRP);
}



size_t RowRankProfileSubmatrixIndices_modular_double (const double p,
                                                      const size_t M, const size_t N,
                                                      double * A,
                                                      const size_t lda,
                                                      size_t ** rowindices,
                                                      size_t ** colindices,
                                                      size_t * R
                                                      , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return RowRankProfileSubmatrixIndices(F,M,N,A,lda,*rowindices,*colindices,*R);
    } else {
        ModularBalanced<double> F(p);
        return RowRankProfileSubmatrixIndices(F,M,N,A,lda,*rowindices,*colindices,*R);
    }
}



size_t ColRankProfileSubmatrixIndices_modular_double (const double p,
                                                      const size_t M, const size_t N,
                                                      double * A,
                                                      const size_t lda,
                                                      size_t** rowindices,
                                                      size_t** colindices,
                                                      size_t* R
                                                      , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return ColRankProfileSubmatrixIndices(F,M,N,A,lda,*rowindices,*colindices,*R);
    } else {
        ModularBalanced<double> F(p);
        return ColRankProfileSubmatrixIndices(F,M,N,A,lda,*rowindices,*colindices,*R);
    }
}



size_t RowRankProfileSubmatrix_modular_double (const double p,
                                               const size_t M, const size_t N,
                                               double * A,
                                               const size_t lda,
                                               double ** X, size_t* R
                                               , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return RowRankProfileSubmatrix(F,M,N,A,lda,*X,*R);
    } else {
        ModularBalanced<double> F(p);
        return RowRankProfileSubmatrix(F,M,N,A,lda,*X,*R);
    }
}



size_t ColRankProfileSubmatrix_modular_double (const double p, const size_t M, const size_t N,
                                               double * A, const size_t lda,
                                               double ** X, size_t* R
                                               , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        return ColRankProfileSubmatrix(F,M,N,A,lda,*X,*R);
    } else {
        ModularBalanced<double> F(p);
        return ColRankProfileSubmatrix(F,M,N,A,lda,*X,*R);
    }
}


/*********************************************/
/* Accessors to Triangular and Echelon forms */
/*********************************************/


void
getTriangular_modular_double (const double p, const enum FFLAS::FFLAS_UPLO Uplo,
                              const enum FFLAS::FFLAS_DIAG Diag,
                              const size_t M, const size_t N, const size_t R,
                              const double * A, const size_t lda,
                              double * T, const size_t ldt,
                              const bool OnlyNonZeroVectors
                              , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        getTriangular(F,(enum FFLAS::FFLAS_UPLO)Uplo,(enum FFLAS::FFLAS_DIAG)Diag,M,N,R,A,lda,T,ldt,OnlyNonZeroVectors);
    } else {
        ModularBalanced<double> F(p);
        getTriangular(F,(enum FFLAS::FFLAS_UPLO)Uplo,(enum FFLAS::FFLAS_DIAG)Diag,M,N,R,A,lda,T,ldt,OnlyNonZeroVectors);
    }
}



void
getTriangularin_modular_double (const double p, const enum FFLAS::FFLAS_UPLO Uplo,
                                const enum FFLAS::FFLAS_DIAG Diag,
                                const size_t M, const size_t N, const size_t R,
                                double * A, const size_t lda
                                , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        getTriangular(F,(enum FFLAS::FFLAS_UPLO)Uplo,(enum FFLAS::FFLAS_DIAG)Diag,M,N,R,A,lda);
    } else {
        ModularBalanced<double> F(p);
        getTriangular(F,(enum FFLAS::FFLAS_UPLO)Uplo,(enum FFLAS::FFLAS_DIAG)Diag,M,N,R,A,lda);
    }
}



void
getEchelonForm_modular_double (const double p, const enum FFLAS::FFLAS_UPLO Uplo,
                               const enum FFLAS::FFLAS_DIAG Diag,
                               const size_t M, const size_t N, const size_t R, const size_t* P,
                               const double * A, const size_t lda,
                               double * T, const size_t ldt,
                               const bool OnlyNonZeroVectors,
                               const enum FFPACK_C_LU_TAG LuTag
                               , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        getEchelonForm(F,(enum FFLAS::FFLAS_UPLO)Uplo,(enum FFLAS::FFLAS_DIAG)Diag,M,N,R,P,A,lda,T,ldt,OnlyNonZeroVectors,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        getEchelonForm(F,(enum FFLAS::FFLAS_UPLO)Uplo,(enum FFLAS::FFLAS_DIAG)Diag,M,N,R,P,A,lda,T,ldt,OnlyNonZeroVectors,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}



void
getEchelonFormin_modular_double (const double p, const enum FFLAS::FFLAS_UPLO Uplo,
                                 const enum FFLAS::FFLAS_DIAG Diag,
                                 const size_t M, const size_t N, const size_t R, const size_t* P,
                                 double * A, const size_t lda,
                                 const enum FFPACK_C_LU_TAG LuTag
                                 , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        getEchelonForm(F,(enum FFLAS::FFLAS_UPLO)Uplo,(enum FFLAS::FFLAS_DIAG)Diag,M,N,R,P,A,lda,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        getEchelonForm(F,(enum FFLAS::FFLAS_UPLO)Uplo,(enum FFLAS::FFLAS_DIAG)Diag,M,N,R,P,A,lda,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}


void
getEchelonTransform_modular_double (const double p, const enum FFLAS::FFLAS_UPLO Uplo,
                                    const enum FFLAS::FFLAS_DIAG Diag,
                                    const size_t M, const size_t N, const size_t R, const size_t* P, const size_t* Q,
                                    const double * A, const size_t lda,
                                    double * T, const size_t ldt,
                                    const enum FFPACK_C_LU_TAG LuTag
                                    , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        getEchelonTransform(F,(enum FFLAS::FFLAS_UPLO)Uplo,(enum FFLAS::FFLAS_DIAG)Diag,M,N,R,P,Q,A,lda,T,ldt,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        getEchelonTransform(F,(enum FFLAS::FFLAS_UPLO)Uplo,(enum FFLAS::FFLAS_DIAG)Diag,M,N,R,P,Q,A,lda,T,ldt,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}




void
getReducedEchelonForm_modular_double (const double p, const enum FFLAS::FFLAS_UPLO Uplo,
                                      const size_t M, const size_t N, const size_t R, const size_t* P,
                                      const double * A, const size_t lda,
                                      double * T, const size_t ldt,
                                      const bool OnlyNonZeroVectors,
                                      const enum FFPACK_C_LU_TAG LuTag
                                      , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        getReducedEchelonForm(F,(enum FFLAS::FFLAS_UPLO)Uplo,M,N,R,P,A,lda,T,ldt,OnlyNonZeroVectors,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        getReducedEchelonForm(F,(enum FFLAS::FFLAS_UPLO)Uplo,M,N,R,P,A,lda,T,ldt,OnlyNonZeroVectors,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}



void
getReducedEchelonFormin_modular_double (const double p, const enum FFLAS::FFLAS_UPLO Uplo,
                                        const size_t M, const size_t N, const size_t R, const size_t* P,
                                        double * A, const size_t lda,
                                        const enum FFPACK_C_LU_TAG LuTag
                                        , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        getReducedEchelonForm(F,(enum FFLAS::FFLAS_UPLO)Uplo,M,N,R,P,A,lda,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        getReducedEchelonForm(F,(enum FFLAS::FFLAS_UPLO)Uplo,M,N,R,P,A,lda,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}



void
getReducedEchelonTransform_modular_double (const double p, const enum FFLAS::FFLAS_UPLO Uplo,
                                           const size_t M, const size_t N, const size_t R, const size_t* P, const size_t* Q,
                                           const double * A, const size_t lda,
                                           double * T, const size_t ldt,
                                           const enum FFPACK_C_LU_TAG LuTag
                                           , bool positive)
{
    if (positive) {
        Modular<double> F(p);
        getReducedEchelonTransform(F,(enum FFLAS::FFLAS_UPLO)Uplo,M,N,R,P,Q,A,lda,T,ldt,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    } else {
        ModularBalanced<double> F(p);
        getReducedEchelonTransform(F,(enum FFLAS::FFLAS_UPLO)Uplo,M,N,R,P,Q,A,lda,T,ldt,(enum FFPACK::FFPACK_LU_TAG)LuTag);
    }
}


void
PLUQtoEchelonPermutation (const size_t N, const size_t R, const size_t * P, size_t * outPerm)
{
    FFPACK::PLUQtoEchelonPermutation(N,R,P,outPerm);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
