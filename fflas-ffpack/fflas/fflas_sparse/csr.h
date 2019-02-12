/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by  Bastien Vialla <bastien.vialla@lirmm.fr>
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

/** @file fflas/fflas_fspmv_CSR.inl
 * NO DOC
 */

#ifndef __FFLASFFPACK_fflas_sparse_CSR_H
#define __FFLASFFPACK_fflas_sparse_CSR_H

namespace FFLAS { /*  CSR */

    template <class _Field> struct Sparse<_Field, SparseMatrix_t::CSR> {
        using Field = _Field;
        bool delayed = false;
        uint64_t kmax = 0;
        index_t m = 0;
        index_t n = 0;
        uint64_t nnz = 0;
        uint64_t nElements = 0;
        uint64_t maxrow = 0;
        index_t *col = nullptr;
        index_t *st = nullptr;
        index_t *stend = nullptr;
        typename _Field::Element_ptr dat;
    };

    template <class _Field>
    struct Sparse<_Field, SparseMatrix_t::CSR_ZO>
    : public Sparse<_Field, SparseMatrix_t::CSR> {
        using Field = _Field;
        int64_t cst = 1;
        bool delayed = false;
    };

    template <class Field, class IndexT>
    inline void sparse_init(const Field &F, Sparse<Field, SparseMatrix_t::CSR> &A,
                            const IndexT *row, const IndexT *col,
                            typename Field::ConstElement_ptr dat, uint64_t rowdim,
                            uint64_t coldim, uint64_t nnz);

    template <class Field, class IndexT>
    inline void sparse_init(const Field &F,
                            Sparse<Field, SparseMatrix_t::CSR_ZO> &A,
                            const IndexT *row, const IndexT *col,
                            typename Field::ConstElement_ptr dat, uint64_t rowdim,
                            uint64_t coldim, uint64_t nnz);

    template <class Field>
    inline void sparse_delete(const Sparse<Field, SparseMatrix_t::CSR> &A);

    template <class Field>
    inline void sparse_delete(const Sparse<Field, SparseMatrix_t::CSR_ZO> &A);

} // FFLAS

#include "fflas-ffpack/fflas/fflas_sparse/csr/csr_utils.inl"
#include "fflas-ffpack/fflas/fflas_sparse/csr/csr_spmv.inl"
#include "fflas-ffpack/fflas/fflas_sparse/csr/csr_spmm.inl"

#if defined(__FFLASFFPACK_USE_OPENMP) || defined(__FFLASFFPACK_USE_TBB)

#include "fflas-ffpack/fflas/fflas_sparse/csr/csr_pspmv.inl"
#include "fflas-ffpack/fflas/fflas_sparse/csr/csr_pspmm.inl"

#endif

#endif // __FFLASFFPACK_fflas_sparse_CSR_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
