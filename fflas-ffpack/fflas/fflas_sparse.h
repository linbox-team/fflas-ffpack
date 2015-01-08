/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   BB <bbboyer@ncsu.edu>
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

#ifndef index_t
#define index_t uint64_t
#endif

// Bigger multiple of s lesser or equal than x, s must be a power of two
#ifndef ROUND_DOWN
#define ROUND_DOWN(x, s) ((x) & ~((s) - 1))
#endif

#ifndef __FFLASFFPACK_CACHE_LINE_SIZE
#define __FFLASFFPACK_CACHE_LINE_SIZE 64
#endif

#if __GNUC_MINOR__ >= 7
#define assume_aligned(pout, pin, v)                                                                                   \
    const decltype(pin) pout = static_cast<decltype(pin)>(__builtin_assume_aligned(pin, v));
#else
#define assume_aligned(pout, pin, v) const decltype(pin) pout = pin;
#endif

#include "fflas-ffpack/config.h"
#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/field/field-traits.h"
#include <type_traits>
#include <vector>

#ifdef __FFLASFFPACK_HAVE_MKL
#ifndef _MKL_H_ // temporary
#error "MKL (mkl.h) not present, while you have MKL enabled"
#endif
#undef index_t
#define index_t MKL_INT
#endif

namespace FFLAS {

enum class SparseMatrix_t {
    CSR,
    CSR_ZO,
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

template <class Field, SparseMatrix_t> struct Sparse;
} // FFLAS

#include "fflas-ffpack/fflas/fflas_sparse/csr.h"
#include "fflas-ffpack/fflas/fflas_sparse/coo.h"
#include "fflas-ffpack/fflas/fflas_sparse/ell.h"
#include "fflas-ffpack/fflas/fflas_sparse/sell.h"
#include "fflas-ffpack/fflas/fflas_sparse/csr_hyb.h"
#include "fflas-ffpack/fflas/fflas_sparse/ell_simd.h"
#include "fflas-ffpack/fflas/fflas_sparse/hyb_zo.h"

namespace FFLAS {
/****************************
 *
 *  SparseMatrix Traits
 *
 ****************************/

template <class F, class M> struct isSparseMatrix : public std::false_type {};

template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::CSR>> : public std::true_type {};

template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::CSR_ZO>> : public std::true_type {};

template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::COO>> : public std::true_type {};

template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::COO_ZO>> : public std::true_type {};

template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::ELL>> : public std::true_type {};

template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::ELL_ZO>> : public std::true_type {};

template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::SELL>> : public std::true_type {};

template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::SELL_ZO>> : public std::true_type {};

template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::ELL_simd>> : public std::true_type {};

template <class Field>
struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::ELL_simd_ZO>> : public std::true_type {};

template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::CSR_HYB>> : public std::true_type {};

template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::HYB_ZO>> : public std::true_type {};

template <class F, class M> struct isZOSparseMatrix : public std::false_type {};

template <class Field> struct isZOSparseMatrix<Field, Sparse<Field, SparseMatrix_t::CSR_ZO>> : public std::true_type {};

template <class Field> struct isZOSparseMatrix<Field, Sparse<Field, SparseMatrix_t::COO_ZO>> : public std::true_type {};

template <class Field> struct isZOSparseMatrix<Field, Sparse<Field, SparseMatrix_t::ELL_ZO>> : public std::true_type {};

template <class Field>
struct isZOSparseMatrix<Field, Sparse<Field, SparseMatrix_t::SELL_ZO>> : public std::true_type {};

template <class Field>
struct isZOSparseMatrix<Field, Sparse<Field, SparseMatrix_t::ELL_simd_ZO>> : public std::true_type {};

/*********************************************************************************************************************
 *
 *    Sparse Details
 *
 *********************************************************************************************************************/

namespace sparse_details {
template <class Field>
inline void init_y(const Field &F, const size_t m, const typename Field::Element b, typename Field::Element_ptr y,
                   FieldCategories::ModularTag);

template <class Field>
inline void init_y(const Field &F, const size_t m, const typename Field::Element b, typename Field::Element_ptr y,
                   FieldCategories::UnparametricTag);

template <class Field>
inline void init_y(const Field &F, const size_t m, const typename Field::Element b, typename Field::Element_ptr y,
                   FieldCategories::GenericTag);

template <class Field>
inline void init_y(const Field &F, const size_t m, const size_t n, const typename Field::Element b,
                   typename Field::Element_ptr y, const int ldy, FieldCategories::UnparametricTag);

template <class Field>
inline void init_y(const Field &F, const size_t m, const size_t n, const typename Field::Element b,
                   typename Field::Element_ptr y, const int ldy, FieldCategories::GenericTag);

template <class Field>
inline void init_y(const Field &F, const size_t m, const size_t n, const typename Field::Element b,
                   typename Field::Element_ptr y, const int ldy, FieldCategories::ModularTag);
}

/*********************************************************************************************************************
 *
 *    SpMV, SpMM, pSpMV, pSpMM
 *
 *********************************************************************************************************************/

template <class Field, class SM>
inline void fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, const typename Field::Element &beta,
                  typename Field::Element_ptr y);

template <class Field, class SM>
inline void pfspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, const typename Field::Element &beta,
                   typename Field::Element_ptr y);

template <class Field, class SM>
inline void fspmm(const Field &F, const SM &A, typename Field::ConstElement_ptr x, const typename Field::Element &beta,
                  typename Field::Element_ptr y);

template <class Field, class SM>
inline void pfspmm(const Field &F, const SM &A, typename Field::ConstElement_ptr x, const typename Field::Element &beta,
                   typename Field::Element_ptr y);
}

#include "fflas-ffpack/fflas/fflas_sparse.inl"

#endif // __FFLASFFPACK_fflas_fflas_sparse_H
