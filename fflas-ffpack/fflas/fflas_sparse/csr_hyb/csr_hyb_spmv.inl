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

#ifndef __FFLASFFPACK_fflas_sparse_CSR_HYB_spmv_INL
#define __FFLASFFPACK_fflas_sparse_CSR_HYB_spmv_INL

namespace FFLAS {
    namespace sparse_details_impl {
        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_HYB> &A, typename Field::ConstElement_ptr x_,
                          typename Field::Element_ptr y_, FieldCategories::GenericTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (uint64_t i = 0; i < A.m; ++i) {
                index_t start = st[4 * i], stop = st[4 * i + 1];
                for (uint64_t j = start; j < stop; ++j) {
                    F.subin(y[i], x[col[j]]);
                }
                start = st[4 * i + 1], stop = st[4 * i + 2];
                for (uint64_t j = start; j < stop; ++j) {
                    F.addin(y[i], x[col[j]]);
                }
                start = st[4 * i + 2], stop = st[4 * (i + 1)];
                index_t startDat = st[4 * i + 3];
                for (uint64_t j = start, k = 0; j < stop; ++j, ++k) {
                    F.axpyin(y[i], dat[startDat + k], x[col[j]]);
                }
            }
        }

        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_HYB> &A, typename Field::ConstElement_ptr x_,
                          typename Field::Element_ptr y_, FieldCategories::UnparametricTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (uint64_t i = 0; i < A.m; ++i) {
                index_t start = st[4 * i], stop = st[4 * i + 1];
                index_t diff = stop - start;
                typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
                uint64_t j = 0;
                for (; j < ROUND_DOWN(diff, 4); j += 4) {
                    y1 += x[col[start + j]];
                    y2 += x[col[start + j + 1]];
                    y3 += x[col[start + j + 2]];
                    y4 += x[col[start + j + 3]];
                }
                for (; j < diff; ++j) {
                    y1 += x[col[start + j]];
                }
                y[i] -= y1 + y2 + y3 + y4;
                y1 = 0;
                y2 = 0;
                y3 = 0;
                y4 = 0;
                start = st[4 * i + 1], stop = st[4 * i + 2];
                diff = stop - start;
                j = 0;
                for (; j < ROUND_DOWN(diff, 4); j += 4) {
                    y1 += x[col[start + j]];
                    y2 += x[col[start + j + 1]];
                    y3 += x[col[start + j + 2]];
                    y4 += x[col[start + j + 3]];
                }
                for (; j < diff; ++j) {
                    y1 += x[col[start + j]];
                }
                y[i] += y1 + y2 + y3 + y4;
                y1 = 0;
                y2 = 0;
                y3 = 0;
                y4 = 0;
                start = st[4 * i + 2], stop = st[4 * (i + 1)];
                diff = stop - start;
                index_t startDat = st[4 * i + 3];
                j = 0;
                for (; j < ROUND_DOWN(diff, 4); j += 4) {
                    y1 += dat[startDat + j] * x[col[start + j]];
                    y2 += dat[startDat + j + 1] * x[col[start + j + 1]];
                    y3 += dat[startDat + j + 2] * x[col[start + j + 2]];
                    y4 += dat[startDat + j + 3] * x[col[start + j + 3]];
                }
                for (; j < diff; ++j) {
                    y1 += dat[startDat + j] * x[col[start + j]];
                }
                y[i] += y1 + y2 + y3 + y4;
            }
        }

        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_HYB> &A, typename Field::ConstElement_ptr x_,
                          typename Field::Element_ptr y_, const uint64_t kmax) {
            return;
            // TODO
        }

    } // CSR_HYB_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_CSR_HYB_spmv_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
