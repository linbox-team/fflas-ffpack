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

/** @file fflas/fflas_fspmv_CSR_HYB.inl
 * NO DOC
 */

#ifndef __FFLASFFPACK_fflas_sparse_CSR_HYB_H
#define __FFLASFFPACK_fflas_sparse_CSR_HYB_H

namespace FFLAS { /*  CSR_HYB */

    template <class _Field> struct Sparse<_Field, SparseMatrix_t::CSR_HYB> {
        using Field = _Field;
        bool delayed = false;
        index_t *col = nullptr;
        index_t *st = nullptr;
        typename _Field::Element_ptr dat;
        uint64_t kmax = 0;
        index_t m = 0;
        index_t n = 0;
        uint64_t nnz = 0;
        uint64_t nElements = 0;
        uint64_t maxrow = 0;
        uint64_t nOnes = 0;
        uint64_t nMOnes = 0;
        uint64_t nOthers = 0;
    };

    template <class Field>
    void sparse_delete(const Sparse<Field, SparseMatrix_t::CSR_HYB> &A);

    template <class Field, class IndexT>
    void sparse_init(const Field &F, Sparse<Field, SparseMatrix_t::CSR_HYB> &A,
                     const IndexT *row, const IndexT *col,
                     typename Field::ConstElement_ptr dat, uint64_t rowdim,
                     uint64_t coldim, uint64_t nnz);
} // FFLAS

#include "fflas-ffpack/fflas/fflas_sparse/csr_hyb/csr_hyb_utils.inl"
#include "fflas-ffpack/fflas/fflas_sparse/csr_hyb/csr_hyb_spmv.inl"
#if defined(__FFLASFFPACK_USE_OPENMP)
#include "fflas-ffpack/fflas/fflas_sparse/csr_hyb/csr_hyb_pspmv.inl"
#endif
#include "fflas-ffpack/fflas/fflas_sparse/csr_hyb/csr_hyb_spmm.inl"
// #include "fflas-ffpack/fflas/fflas_sparse/csr_hyb/csr_hyb_pspmm.inl"

#endif // __FFLASFFPACK_fflas_sparse_CSR_HYB_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
