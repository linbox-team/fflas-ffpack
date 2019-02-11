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

#ifndef __FFLASFFPACK_fflas_sparse_HYB_ZO_utils_INL
#define __FFLASFFPACK_fflas_sparse_HYB_ZO_utils_INL

namespace FFLAS {

    // #define HYB_ZO_DEBUG 1

    template <class Field> inline void sparse_delete(const Sparse<Field, SparseMatrix_t::HYB_ZO> &A) {
        if (A.dat != nullptr)
            sparse_delete(*(A.dat));
        if (A.one != nullptr)
            sparse_delete(*(A.one));
        if (A.mone != nullptr)
            sparse_delete(*(A.mone));
    }

    template <class Field, class IndexT>
    inline void sparse_init(const Field &F, Sparse<Field, SparseMatrix_t::HYB_ZO> &A, const IndexT *row, const IndexT *col,
                            typename Field::ConstElement_ptr dat, uint64_t rowdim, uint64_t coldim, uint64_t nnz) {

        A.m = rowdim;
        A.n = coldim;
        A.nnz = nnz;
        A.delayed = true;
        A.nElements = nnz;
        uint64_t nOnes = 0, nMOnes = 0, nOthers = 0;
        for (uint64_t i = 0; i < nnz; ++i) {
            if (F.isOne(dat[i]))
                nOnes++;
            else if (F.isMOne(dat[i]))
                nMOnes++;
            else
                nOthers++;
        }

        typename Field::Element_ptr dat2(0);
        index_t *colOne = nullptr, *colMOne = nullptr, *colOther = nullptr, *rowOne = nullptr, *rowMOne = nullptr,
                *rowOther = nullptr;
        if (nOnes) {
            colOne = fflas_new<index_t>(nOnes, Alignment::CACHE_LINE);
            rowOne = fflas_new<index_t>(nOnes, Alignment::CACHE_LINE);
        }
        if (nMOnes) {
            colMOne = fflas_new<index_t>(nMOnes, Alignment::CACHE_LINE);
            rowMOne = fflas_new<index_t>(nMOnes, Alignment::CACHE_LINE);
        }
        if (nOthers) {
            dat2 = fflas_new(F, nOthers, Alignment::CACHE_LINE);
            colOther = fflas_new<index_t>(nOthers, Alignment::CACHE_LINE);
            rowOther = fflas_new<index_t>(nOthers, Alignment::CACHE_LINE);
        }

        uint64_t itOne = 0, itMOne = 0, itOther = 0;
        for (uint64_t i = 0; i < nnz; ++i) {
            if (F.isOne(dat[i])) {
                colOne[itOne] = col[i];
                rowOne[itOne] = row[i];
                ++itOne;
            } else if (F.isMOne(dat[i])) {
                colMOne[itMOne] = col[i];
                rowMOne[itMOne] = row[i];
                ++itMOne;
            } else {
                dat2[itOther] = dat[i];
                colOther[itOther] = col[i];
                rowOther[itOther] = row[i];
                ++itOther;
            }
        }

        if (nOnes) {
            A.one = new Sparse<Field, SparseMatrix_t::CSR_ZO>();
            sparse_init(F, *(A.one), rowOne, colOne, nullptr, rowdim, coldim, nOnes);
        }
        if (nMOnes) {
            A.mone = new Sparse<Field, SparseMatrix_t::CSR_ZO>();
            sparse_init(F, *(A.mone), rowMOne, colMOne, nullptr, rowdim, coldim, nMOnes);
            A.mone->cst = -1;
        }
        if (nOthers) {
            A.dat = new Sparse<Field, SparseMatrix_t::CSR>();
            sparse_init(F, *(A.dat), rowOther, colOther, dat2, rowdim, coldim, nOthers);
        }

        if (nOnes) {
            fflas_delete(colOne);
            fflas_delete(rowOne);
        }
        if (nMOnes) {
            fflas_delete(colMOne);
            fflas_delete(rowMOne);
        }
        if (nOthers) {
            fflas_delete(colOther);
            fflas_delete(rowOther);
            fflas_delete(dat2);
        }
    }

    template<typename _Field>
    std::ostream& operator<<(std::ostream& os, const Sparse<_Field, SparseMatrix_t::HYB_ZO>& A) {
        return sparse_print(os << "non-ones: ", *(A.dat));
    }


}

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
