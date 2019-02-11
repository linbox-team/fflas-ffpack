/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   Bastien Vialla <bastien.vialla@lirmm.fr>
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

/** @file fflas/fflas_fspmv_ELL_R.inl
 * NO DOC
 */

#ifndef __FFLASFFPACK_fflas_sparse_ELL_R_H
#define __FFLASFFPACK_fflas_sparse_ELL_R_H

namespace FFLAS { /*  ELL_R */

    template <class _Field> struct Sparse<_Field, SparseMatrix_t::ELL_R> {
        bool delayed = false;
        uint64_t kmax = 0;
        index_t m = 0;
        index_t n = 0;
        index_t ld = 0;
        uint64_t nnz = 0;
        uint64_t maxrow = 0;
        uint64_t mRow = 0;
        index_t *col = nullptr;
        index_t *row = nullptr;
        typename _Field::Element_ptr dat;
    };

    template <class _Field>
    struct Sparse<_Field, SparseMatrix_t::ELL_R_ZO>
    : public Sparse<_Field, SparseMatrix_t::ELL_R> {
        typename _Field::Element cst = 1;
    };

    template <class Field>
    void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_R> &A,
               typename Field::ConstElement_ptr x,
               const typename Field::Element &beta, typename Field::Element_ptr y);

    template <class Field>
    void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_R_ZO> &A,
               typename Field::ConstElement_ptr x,
               const typename Field::Element &beta, typename Field::Element_ptr y);

    template <class Field>
    void fspmm(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_R> &A,
               const size_t blockSize, const typename Field::Element_ptr &x,
               const int ldx, const typename Field::Element &beta,
               typename Field::Element_ptr &y, const int ldy);

    template <class Field>
    void fspmm(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_R_ZO> &A,
               const size_t blockSize, const typename Field::Element_ptr &x,
               const int ldx, const typename Field::Element &beta,
               typename Field::Element_ptr &y, const int ldy);
} // FFLAS

#include "fflas-ffpack/fflas/fflas_sparse/ell_r_spmv.inl"
// #include "fflas-ffpack/fflas/fflas_sparse/ell_r_spmm.inl"

#endif // __FFLASFFPACK_fflas_sparse_ELL_R_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
