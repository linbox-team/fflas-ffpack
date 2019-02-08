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

#ifndef __FFLASFFPACK_fflas_sparse_ELL_utils_INL
#define __FFLASFFPACK_fflas_sparse_ELL_utils_INL

#include <vector>

namespace FFLAS {
    template <class Field> inline void sparse_delete(const Sparse<Field, SparseMatrix_t::ELL> &A) {
        fflas_delete(A.dat);
        fflas_delete(A.col);
    }

    template <class Field> inline void sparse_delete(const Sparse<Field, SparseMatrix_t::ELL_ZO> &A) {
        fflas_delete(A.col);
    }

    template <class Field, class IndexT>
    inline void sparse_init(const Field &F, Sparse<Field, SparseMatrix_t::ELL> &A, const IndexT *row, const IndexT *col,
                            typename Field::ConstElement_ptr dat, uint64_t rowdim, uint64_t coldim, uint64_t nnz) {
        A.kmax = Protected::DotProdBoundClassic(F, F.one);
        A.m = rowdim;
        A.n = coldim;
        A.nnz = nnz;
        // cout << A.m << " " << A.n << endl;
        std::vector<uint64_t> rows(A.m, 0);
        for (uint64_t i = 0; i < A.nnz; ++i)
            rows[row[i]]++;
        A.maxrow = *(std::max_element(rows.begin(), rows.end()));
        // cout << "maxrow : " << A.maxrow << endl;
        A.ld = A.maxrow;
        if (A.kmax > A.maxrow)
            A.delayed = true;
        // cout << A.ld << " " << rowdim << " " << nnz << " " << A.ld*rowdim << endl;
        A.nElements = A.m * A.ld;
        A.col = fflas_new<index_t>(rowdim * A.ld, Alignment::CACHE_LINE);
        A.dat = fflas_new(F, rowdim * A.ld, Alignment::CACHE_LINE);

        for (size_t i = 0; i < rowdim * A.ld; ++i) {
            A.col[i] = 0;
            F.assign(A.dat[i], F.zero);
        }

        size_t currow = row[0], it = 0;
        for (size_t i = 0; i < nnz; ++i) {
            if (row[i] != currow) {
                it = 0;
                currow = row[i];
            }
            A.col[row[i] * A.ld + it] = col[i];
            A.dat[row[i] * A.ld + it] = dat[i];
            ++it;
        }
    }

    template <class Field, class IndexT>
    inline void sparse_init(const Field &F, Sparse<Field, SparseMatrix_t::ELL_ZO> &A, const IndexT *row, const IndexT *col,
                            typename Field::ConstElement_ptr dat, uint64_t rowdim, uint64_t coldim, uint64_t nnz) {
        A.kmax = Protected::DotProdBoundClassic(F, F.one);
        A.m = rowdim;
        A.n = coldim;
        A.nnz = nnz;
        std::vector<uint64_t> rows(A.m, 0);
        for (uint64_t i = 0; i < A.nnz; ++i)
            rows[row[i]]++;
        A.maxrow = *(std::max_element(rows.begin(), rows.end()));
        A.ld = A.maxrow;
        if (A.kmax > A.maxrow)
            A.delayed = true;
        A.nElements = A.m * A.ld;
        A.col = fflas_new<index_t>(rowdim * A.ld, Alignment::CACHE_LINE);

        for (size_t i = 0; i < rowdim * A.ld; ++i) {
            A.col[i] = 0;
        }

        size_t currow = row[0], it = 0;

        for (size_t i = 0; i < nnz; ++i) {
            if (row[i] != currow) {
                it = 0;
                currow = row[i];
            }
            A.col[row[i] * A.ld + it] = col[i];
            ++it;
        }
    }
}

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
