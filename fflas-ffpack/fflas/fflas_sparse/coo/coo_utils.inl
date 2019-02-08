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

#ifndef __FFLASFFPACK_fflas_sparse_coo_utils_INL
#define __FFLASFFPACK_fflas_sparse_coo_utils_INL

namespace FFLAS {

    template <class Field> inline void sparse_delete(const Sparse<Field, SparseMatrix_t::COO> &A) {
        fflas_delete(A.dat);
        fflas_delete(A.col);
        fflas_delete(A.row);
    }

    template <class Field> inline void sparse_delete(const Sparse<Field, SparseMatrix_t::COO_ZO> &A) {
        fflas_delete(A.col);
        fflas_delete(A.row);
    }

    template <class Field, class IndexT>
    inline void sparse_init(const Field &F, Sparse<Field, SparseMatrix_t::COO> &A, const IndexT *row, const IndexT *col,
                            typename Field::ConstElement_ptr dat, uint64_t rowdim, uint64_t coldim, uint64_t nnz) {
        A.kmax = Protected::DotProdBoundClassic(F, F.one);
        A.m = rowdim;
        A.n = coldim;
        A.nnz = nnz;
        A.nElements = nnz;
        std::vector<uint64_t> rows(rowdim, 0);
        for (uint64_t i = 0; i < A.nnz; ++i)
            rows[row[i]]++;
        A.maxrow = *(std::max_element(rows.begin(), rows.end()));
        if (A.kmax > A.maxrow)
            A.delayed = true;
        A.col = fflas_new<index_t>(nnz, Alignment::CACHE_LINE);
        A.row = fflas_new<index_t>(nnz, Alignment::CACHE_LINE);
        A.dat = fflas_new(F, nnz, Alignment::CACHE_LINE);

        for (uint64_t i = 0; i < nnz; ++i) {
            A.col[i] = (index_t)col[i];
            A.row[i] = (index_t)row[i];
            A.dat[i] = dat[i];
        }
    }

    template <class Field, class IndexT>
    inline void sparse_init(const Field &F, Sparse<Field, SparseMatrix_t::COO_ZO> &A, const IndexT *row, const IndexT *col,
                            typename Field::ConstElement_ptr dat, uint64_t rowdim, uint64_t coldim, uint64_t nnz) {
        A.kmax = Protected::DotProdBoundClassic(F, F.one);
        A.m = rowdim;
        A.n = coldim;
        A.nnz = nnz;
        A.nElements = nnz;
        std::vector<uint64_t> rows(A.m, 0);
        for (uint64_t i = 0; i < A.nnz; ++i)
            rows[row[i]]++;
        A.maxrow = *(std::max_element(rows.begin(), rows.end()));
        if (A.kmax > A.maxrow)
            A.delayed = true;

        A.col = fflas_new<index_t>(nnz, Alignment::CACHE_LINE);
        A.row = fflas_new<index_t>(nnz, Alignment::CACHE_LINE);

        for (uint64_t i = 0; i < nnz; ++i) {
            A.col[i] = (index_t)col[i];
            A.row[i] = (index_t)row[i];
        }
    }
}

#endif // __FFLASFFPACK_fflas_sparse_coo_spmv_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
