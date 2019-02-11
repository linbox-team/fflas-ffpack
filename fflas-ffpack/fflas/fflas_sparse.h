/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *              Bastien Vialla <bastien.vialla@lirmm.fr>
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 * ========LICENCE========
 *.
 */

/** @file fflas/fflas_sparse.h
*/

#ifndef __FFLASFFPACK_fflas_fflas_sparse_H
#define __FFLASFFPACK_fflas_fflas_sparse_H

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/config.h"
#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/paladin/parallel.h"

#include <recint/recint.h>
#include <givaro/udl.h>

#ifndef index_t
#define index_t uint32_t
#endif

#ifdef __FFLASFFPACK_HAVE_MKL
#ifndef _MKL_H_ // temporary
#error "MKL (mkl.h) not present, while you have MKL enabled"
#endif
#undef index_t
#define index_t MKL_INT

#endif // __FFLASFFPACK_HAVE_MKL

// Bigger multiple of s lesser or equal than x, s must be a power of two
#ifndef ROUND_DOWN
#define ROUND_DOWN(x, s) ((x) & ~((s)-1))
#endif

#ifndef __FFLASFFPACK_CACHE_LINE_SIZE
#define __FFLASFFPACK_CACHE_LINE_SIZE 64
#endif

#if (__GNUC_MAJOR > 4 || (__GNUC_MAJOR == 4 &&__GNUC_MINOR__ >= 7)) || defined(__clang__)
#define assume_aligned(pout, pin, v) decltype(pin) __restrict__ pout = static_cast<decltype(pin)>(__builtin_assume_aligned(pin, v));
#elif defined(__INTEL_COMPILER)
#define assume_aligned(pout, pin, v) \
decltype(pin) __restrict pout = pin; \
__assume_aligned(pout)
#else
#define assume_aligned(pout, pin, v) decltype(pin) pout = pin;
#endif

#define DENSE_THRESHOLD 0.5

#include "fflas-ffpack/fflas/fflas.h"

#include "fflas-ffpack/field/field-traits.h"
#include "fflas-ffpack/fflas/fflas_bounds.inl"
#include "fflas-ffpack/utils/fflas_memory.h"
#include "fflas-ffpack/paladin/parallel.h"

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
#include "fflas-ffpack/fflas/fflas_simd.h"
#endif

#include <type_traits>
#include <vector>
#include <iostream>


namespace MKL_CONFIG {
    static const double dalpha = 1;
    static const float  salpha = 1;
    static const double dbeta = 0;
    static const float  sbeta = 0;
    static const char metaChar[4] = {'G', ' ', ' ', 'C'};
    static const char trans[1] = {'N'};
}

namespace FFLAS {

    enum class SparseMatrix_t {
        CSR,
        CSR_ZO,
        CSC,
        CSC_ZO,
        COO,
        COO_ZO,
        ELL,
        ELL_ZO,
        SELL,
        SELL_ZO,
        ELL_simd,
        ELL_simd_ZO,
        CSR_HYB,
        HYB_ZO
    };

    template <class Field, SparseMatrix_t, class IdxT = index_t, class PtrT = index_t> struct Sparse;

} // FFLAS
#include "fflas-ffpack/fflas/fflas_sparse/sparse_matrix_traits.h"
#include "fflas-ffpack/fflas/fflas_sparse/utils.h"
#include "fflas-ffpack/fflas/fflas_sparse/csr.h"
#include "fflas-ffpack/fflas/fflas_sparse/coo.h"
#include "fflas-ffpack/fflas/fflas_sparse/ell.h"
#include "fflas-ffpack/fflas/fflas_sparse/sell.h"
#include "fflas-ffpack/fflas/fflas_sparse/csr_hyb.h"
#include "fflas-ffpack/fflas/fflas_sparse/ell_simd.h"
#include "fflas-ffpack/fflas/fflas_sparse/hyb_zo.h"
// #include "fflas-ffpack/fflas/fflas_sparse/sparse_matrix.h"

namespace FFLAS {

    /*********************************************************************************************************************
     *
     *    Sparse Details
     *
     *********************************************************************************************************************/

    namespace sparse_details {

        template <class Field>
        inline void init_y(const Field &F, const size_t m, const typename Field::Element b, typename Field::Element_ptr y);

        template <class Field>
        inline void init_y(const Field &F, const size_t m, const size_t n, const typename Field::Element b,
                           typename Field::Element_ptr y, const int ldy);

        /*************************************
          fspmv
         **************************************/

        template <class Field, class SM, class FC, class MZO>
        inline typename std::enable_if<
        !(std::is_same<typename ElementTraits<typename Field::Element>::value, ElementCategories::MachineFloatTag>::value ||
          std::is_same<typename ElementTraits<typename Field::Element>::value,
          ElementCategories::MachineIntTag>::value)>::type
        fspmv_dispatch(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y, FC fc,
                       MZO mzo);

        template <class Field, class SM, class FC, class MZO>
        inline typename std::enable_if<
        std::is_same<typename ElementTraits<typename Field::Element>::value, ElementCategories::MachineFloatTag>::value ||
        std::is_same<typename ElementTraits<typename Field::Element>::value, ElementCategories::MachineIntTag>::value>::type
        fspmv_dispatch(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y, FC fc,
                       MZO mzo);

        // non ZO matrix
        template <class Field, class SM>
        inline void fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                          FieldCategories::GenericTag, NotZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<!isSparseMatrixSimdFormat<Field, SM>::value>::type
        fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
              FieldCategories::UnparametricTag, NotZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<isSparseMatrixSimdFormat<Field, SM>::value>::type
        fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
              FieldCategories::UnparametricTag, NotZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<!isSparseMatrixSimdFormat<Field, SM>::value>::type
        fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
              FieldCategories::ModularTag, NotZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<isSparseMatrixSimdFormat<Field, SM>::value>::type
        fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
              FieldCategories::ModularTag, NotZOSparseMatrix);

        // ZO matrix

        template <class Field, class SM>
        inline void fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                          FieldCategories::GenericTag, ZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<!isSparseMatrixSimdFormat<Field, SM>::value>::type
        fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
              FieldCategories::UnparametricTag, ZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<isSparseMatrixSimdFormat<Field, SM>::value>::type
        fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
              FieldCategories::UnparametricTag, ZOSparseMatrix);

        template <class Field, class SM>
        inline void fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                          FieldCategories::ModularTag, std::true_type);

        /*************************************
          fspmm
         **************************************/

        template <class Field, class SM, class FCat, class MZO>
        inline typename std::enable_if<
        !(std::is_same<typename ElementTraits<typename Field::Element>::value, ElementCategories::MachineFloatTag>::value ||
          std::is_same<typename ElementTraits<typename Field::Element>::value,
          ElementCategories::MachineIntTag>::value)>::type
        fspmm_dispatch(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
                       typename Field::Element_ptr y, int ldy, FCat, MZO);

        template <class Field, class SM, class FCat, class MZO>
        inline typename std::enable_if<
        std::is_same<typename ElementTraits<typename Field::Element>::value, ElementCategories::MachineFloatTag>::value ||
        std::is_same<typename ElementTraits<typename Field::Element>::value, ElementCategories::MachineIntTag>::value>::type
        fspmm_dispatch(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
                       typename Field::Element_ptr y, int ldy, FCat, MZO);

        template <class Field, class SM>
        inline void fspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
                          typename Field::Element_ptr y, int ldy, FieldCategories::GenericTag, NotZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<support_simd<typename Field::Element>::value>::type
        fspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
              typename Field::Element_ptr y, int ldy, FieldCategories::UnparametricTag, NotZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<!support_simd<typename Field::Element>::value>::type
        fspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
              typename Field::Element_ptr y, int ldy, FieldCategories::UnparametricTag, NotZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<support_simd<typename Field::Element>::value>::type
        fspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
              typename Field::Element_ptr y, int ldy, FieldCategories::ModularTag, NotZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<!support_simd<typename Field::Element>::value>::type
        fspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
              typename Field::Element_ptr y, int ldy, FieldCategories::ModularTag, NotZOSparseMatrix);

        // ZO matrix
        template <class Field, class SM>
        inline void fspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
                          typename Field::Element_ptr y, int ldy, FieldCategories::GenericTag, ZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<support_simd<typename Field::Element>::value>::type
        fspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
              typename Field::Element_ptr y, int ldy, FieldCategories::UnparametricTag, ZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<!support_simd<typename Field::Element>::value>::type
        fspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
              typename Field::Element_ptr y, int ldy, FieldCategories::UnparametricTag, ZOSparseMatrix);

        template <class Field, class SM>
        inline void fspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
                          typename Field::Element_ptr y, int ldy, FieldCategories::ModularTag, ZOSparseMatrix);

        /*************************************
          pfspmm
         **************************************/

        template <class Field, class SM, class FCat, class MZO>
        inline typename std::enable_if<
        !(std::is_same<typename ElementTraits<typename Field::Element>::value, ElementCategories::MachineFloatTag>::value ||
          std::is_same<typename ElementTraits<typename Field::Element>::value,
          ElementCategories::MachineIntTag>::value)>::type
        pfspmm_dispatch(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
                        typename Field::Element_ptr y, int ldy, FCat, MZO);

        template <class Field, class SM, class FCat, class MZO>
        inline typename std::enable_if<
        std::is_same<typename ElementTraits<typename Field::Element>::value, ElementCategories::MachineFloatTag>::value ||
        std::is_same<typename ElementTraits<typename Field::Element>::value, ElementCategories::MachineIntTag>::value>::type
        pfspmm_dispatch(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
                        typename Field::Element_ptr y, int ldy, FCat, MZO);

        template <class Field, class SM>
        inline void pfspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
                           typename Field::Element_ptr y, int ldy, FieldCategories::GenericTag, NotZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<support_simd<typename Field::Element>::value>::type
        pfspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
               typename Field::Element_ptr y, int ldy, FieldCategories::UnparametricTag, NotZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<!support_simd<typename Field::Element>::value>::type
        pfspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
               typename Field::Element_ptr y, int ldy, FieldCategories::UnparametricTag, NotZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<support_simd<typename Field::Element>::value>::type
        pfspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
               typename Field::Element_ptr y, int ldy, FieldCategories::ModularTag, NotZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<!support_simd<typename Field::Element>::value>::type
        pfspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
               typename Field::Element_ptr y, int ldy, FieldCategories::ModularTag, NotZOSparseMatrix);

        // ZO matrix
        template <class Field, class SM>
        inline void pfspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
                           typename Field::Element_ptr y, int ldy, FieldCategories::GenericTag, ZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<support_simd<typename Field::Element>::value>::type
        pfspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
               typename Field::Element_ptr y, int ldy, FieldCategories::UnparametricTag, ZOSparseMatrix);

        template <class Field, class SM>
        inline typename std::enable_if<!support_simd<typename Field::Element>::value>::type
        pfspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
               typename Field::Element_ptr y, int ldy, FieldCategories::UnparametricTag, ZOSparseMatrix);

        template <class Field, class SM>
        inline void pfspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
                           typename Field::Element_ptr y, int ldy, FieldCategories::ModularTag, ZOSparseMatrix);

        /*************************************
          pfspmv
         **************************************/
        template <class Field, class SM>
        inline void pfspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                           FieldCategories::GenericTag, std::false_type);

        template <class Field, class SM>
        inline void pfspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                           FieldCategories::UnparametricTag, std::false_type);

        template <class Field, class SM>
        inline void pfspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                           FieldCategories::ModularTag, std::false_type);

        template <class Field, class SM>
        inline void pfspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                           FieldCategories::GenericTag, std::true_type);

        template <class Field, class SM>
        inline void pfspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                           FieldCategories::UnparametricTag, std::true_type);

        template <class Field, class SM>
        inline void pfspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                           FieldCategories::ModularTag, std::true_type);

    } // sparse_details

    /*********************************************************************************************************************
     *
     *    SpMV, SpMM, pSpMV, pSpMM
     *
     *********************************************************************************************************************/

    template <class Field, class SM>
    inline void fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, const typename Field::Element &beta,
                      typename Field::Element_ptr y);

    template <class Field, class SM>
    inline void fspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
                      const typename Field::Element &beta, typename Field::Element_ptr y, int ldy);

#if defined(__FFLASFFPACK_USE_OPENMP)
    template <class Field, class SM>
    inline void pfspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, const typename Field::Element &beta,
                       typename Field::Element_ptr y);

    template <class Field, class SM>
    inline void pfspmm(const Field &F, const SM &A, size_t blockSize, typename Field::ConstElement_ptr x, int ldx,
                       const typename Field::Element &beta, typename Field::Element_ptr y, int ldy);
#endif
}

#include "fflas-ffpack/fflas/fflas_sparse.inl"

#include "fflas-ffpack/fflas/fflas_sparse/read_sparse.h"


namespace FFLAS {
    struct HelperFlag {
        static constexpr uint64_t none = 0_ui64;
        static constexpr uint64_t coo  = 1_ui64;
        static constexpr uint64_t csr  = 1_ui64 << 1;
        static constexpr uint64_t ell  = 1_ui64 << 2;
        static constexpr uint64_t aut  = 1_ui64 << 32;
        static constexpr uint64_t pm1  = 1_ui64 << 33;
    };

    template<class Field>
    struct CsrMat {
        typename FFLAS::Sparse<Field,SparseMatrix_t::CSR,int16_t> * _csr16 = nullptr;
        typename FFLAS::Sparse<Field,SparseMatrix_t::CSR,int32_t> * _csr32 = nullptr ;
        typename FFLAS::Sparse<Field,SparseMatrix_t::CSR,int64_t> * _csr64 = nullptr ;

        typename FFLAS::Sparse<Field,SparseMatrix_t::CSR_ZO,int16_t> * _csr16_zo = nullptr ;
        typename FFLAS::Sparse<Field,SparseMatrix_t::CSR_ZO,int32_t> * _csr32_zo = nullptr ;
        typename FFLAS::Sparse<Field,SparseMatrix_t::CSR_ZO,int64_t> * _csr64_zo = nullptr ;
    };

    template<class Field>
    struct CooMat {
        typename FFLAS::Sparse<Field,SparseMatrix_t::COO,int16_t> * _coo16 = nullptr;
        typename FFLAS::Sparse<Field,SparseMatrix_t::COO,int32_t> * _coo32 = nullptr ;
        typename FFLAS::Sparse<Field,SparseMatrix_t::COO,int64_t> * _coo64 = nullptr ;

        typename FFLAS::Sparse<Field,SparseMatrix_t::COO_ZO,int16_t> * _coo16_zo = nullptr ;
        typename FFLAS::Sparse<Field,SparseMatrix_t::COO_ZO,int32_t> * _coo32_zo = nullptr ;
        typename FFLAS::Sparse<Field,SparseMatrix_t::COO_ZO,int64_t> * _coo64_zo = nullptr ;
    };

    template<class Field>
    struct EllMat {
        typename FFLAS::Sparse<Field,SparseMatrix_t::ELL,int16_t> * _ell16 = nullptr;
        typename FFLAS::Sparse<Field,SparseMatrix_t::ELL,int32_t> * _ell32 = nullptr ;
        typename FFLAS::Sparse<Field,SparseMatrix_t::ELL,int64_t> * _ell64 = nullptr ;

        typename FFLAS::Sparse<Field,SparseMatrix_t::ELL_ZO,int16_t> * _ell16_zo = nullptr ;
        typename FFLAS::Sparse<Field,SparseMatrix_t::ELL_ZO,int32_t> * _ell32_zo = nullptr ;
        typename FFLAS::Sparse<Field,SparseMatrix_t::ELL_ZO,int64_t> * _ell64_zo = nullptr ;
    };


    template<class Field, int flag = HelperFlag::none >
    struct SpMat {
        typename FFLAS::CooMat<Field> * _coo = nullptr ;
        typename FFLAS::CsrMat<Field> * _csr = nullptr ;
        typename FFLAS::EllMat<Field> * _ell = nullptr ;
    };
}

#undef ROUND_DOWN
#undef DENSE_THRESHOLD
#undef assume_aligned

#endif // __FFLASFFPACK_fflas_fflas_sparse_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
