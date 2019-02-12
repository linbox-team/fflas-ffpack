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

#ifndef __FFLASFFPACK_fflas_sparse_ELL_spmv_INL
#define __FFLASFFPACK_fflas_sparse_ELL_spmv_INL

namespace FFLAS {
    namespace sparse_details_impl {
        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::ELL> &A, typename Field::ConstElement_ptr x_,
                          typename Field::Element_ptr y_, FieldCategories::GenericTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);

            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t start = 0;
            for (index_t i = 0; i < A.m; ++i, start += A.ld) {
                index_t j = 0;
                typename Field::Element y1, y2, y3, y4;
                F.assign(y1, F.zero);
                F.assign(y2, F.zero);
                F.assign(y3, F.zero);
                F.assign(y4, F.zero);
                for (; j < ROUND_DOWN(A.ld, 4); j += 4) {
                    F.axpyin(y1, dat[start + j], x[col[start + j]]);
                    F.axpyin(y2, dat[start + j + 1], x[col[start + j + 1]]);
                    F.axpyin(y3, dat[start + j + 2], x[col[start + j + 2]]);
                    F.axpyin(y4, dat[start + j + 3], x[col[start + j + 3]]);
                }
                for (; j < A.ld; ++j) {
                    F.axpyin(y1, dat[start + j], x[col[start + j]]);
                }
                F.addin(y[i], y1);
                F.addin(y[i], y2);
                F.addin(y[i], y3);
                F.addin(y[i], y4);
            }
        }

        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::ELL> &A, typename Field::ConstElement_ptr x_,
                          typename Field::Element_ptr y_, FieldCategories::UnparametricTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);

            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t start = 0;
            for (index_t i = 0; i < A.m; ++i, start += A.ld) {
                index_t j = 0;
                typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
                for (; j < ROUND_DOWN(A.ld, 4); j += 4) {
                    y1 += dat[start + j] * x[col[start + j]];
                    y2 += dat[start + j + 1] * x[col[start + j + 1]];
                    y3 += dat[start + j + 2] * x[col[start + j + 2]];
                    y4 += dat[start + j + 3] * x[col[start + j + 3]];
                }
                for (; j < A.ld; ++j) {
                    y1 += dat[start + j] * x[col[start + j]];
                }
                y[i] += y1 + y2 + y3 + y4;
            }
        }

        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::ELL> &A, typename Field::ConstElement_ptr x_,
                          typename Field::Element_ptr y_, const uint64_t kmax) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);

            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t start = 0;
            index_t block = (A.ld) / kmax;
            for (index_t i = 0; i < A.m; ++i, start += A.ld) {
                index_t j_loc = 0, j = 0;
                for (index_t l = 0; l < (index_t)block; ++l) {
                    j_loc += kmax;
                    for (; j < j_loc; ++j) {
                        y[i] += dat[start + j] * x[col[start + j]];
                    }
                    F.reduce(y[i]);
                }
                for (; j < A.ld; ++j) {
                    y[i] += dat[start + j] * x[col[start + j]];
                }
                F.reduce(y[i]);
            }
        }

        // template <class Field, class Func>
        // inline void
        // fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_ZO> &A,
        //       typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
        //       Func &&func, FieldCategories::GenericTag) {
        //     index_t start = 0;
        //     for (index_t i = 0; i < A.m; ++i, start += A.ld) {
        //         index_t j = 0;
        //         typename Field::Element y1, y2, y3, y4;
        //         F.assign(y1, F.zero);
        //         F.assign(y2, F.zero);
        //         F.assign(y3, F.zero);
        //         F.assign(y4, F.zero);
        //         for (; j < ROUND_DOWN(A.ld, 4); j += 4) {
        //             func(y1, x[col[start + j]]);
        //             func(y2, x[col[start + j + 1]]);
        //             func(y3, x[col[start + j + 2]]);
        //             func(y4, x[col[start + j + 3]]);
        //         }
        //         for (; j < A.ld; ++j) {
        //             func(y1, x[col[start + j]]);
        //         }
        //         F.addin(y[i], y1);
        //         F.addin(y[i], y2);
        //         F.addin(y[i], y3);
        //         F.addin(y[i], y4);
        //     }
        // }

        template <class Field>
        inline void fspmv_one(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_ZO> &A,
                              typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                              FieldCategories::GenericTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);

            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t start = 0;
            for (index_t i = 0; i < A.m; ++i, start += A.ld) {
                index_t j = 0;
                typename Field::Element y1, y2, y3, y4;
                F.assign(y1, F.zero);
                F.assign(y2, F.zero);
                F.assign(y3, F.zero);
                F.assign(y4, F.zero);
                for (; j < ROUND_DOWN(A.ld, 4); j += 4) {
                    F.addin(y1, x[col[start + j]]);
                    F.addin(y2, x[col[start + j + 1]]);
                    F.addin(y3, x[col[start + j + 2]]);
                    F.addin(y4, x[col[start + j + 3]]);
                }
                for (; j < A.ld; ++j) {
                    F.addin(y1, x[col[start + j]]);
                }
                F.addin(y[i], y1);
                F.addin(y[i], y2);
                F.addin(y[i], y3);
                F.addin(y[i], y4);
            }
        }

        template <class Field>
        inline void fspmv_mone(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_ZO> &A,
                               typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                               FieldCategories::GenericTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);

            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t start = 0;
            for (index_t i = 0; i < A.m; ++i, start += A.ld) {
                index_t j = 0;
                typename Field::Element y1, y2, y3, y4;
                F.assign(y1, F.zero);
                F.assign(y2, F.zero);
                F.assign(y3, F.zero);
                F.assign(y4, F.zero);
                for (; j < ROUND_DOWN(A.ld, 4); j += 4) {
                    F.addin(y1, x[col[start + j]]);
                    F.addin(y2, x[col[start + j + 1]]);
                    F.addin(y3, x[col[start + j + 2]]);
                    F.addin(y4, x[col[start + j + 3]]);
                }
                for (; j < A.ld; ++j) {
                    F.addin(y1, x[col[start + j]]);
                }
                F.subin(y[i], y1);
                F.subin(y[i], y2);
                F.subin(y[i], y3);
                F.subin(y[i], y4);
            }
        }

        template <class Field>
        inline void fspmv_one(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_ZO> &A,
                              typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                              FieldCategories::UnparametricTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);

            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t start = 0;
            for (index_t i = 0; i < A.m; ++i, start += A.ld) {
                index_t j = 0;
                typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
                for (; j < ROUND_DOWN(A.ld, 4); j += 4) {
                    y1 += x[col[start + j]];
                    y2 += x[col[start + j + 1]];
                    y3 += x[col[start + j + 2]];
                    y4 += x[col[start + j + 3]];
                }
                for (; j < A.ld; ++j) {
                    y1 += x[col[start + j]];
                }
                y[i] += y1 + y2 + y3 + y4;
            }
        }

        template <class Field>
        inline void fspmv_mone(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_ZO> &A,
                               typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                               FieldCategories::UnparametricTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);

            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t start = 0;
            for (index_t i = 0; i < A.m; ++i, start += A.ld) {
                index_t j = 0;
                typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
                for (; j < ROUND_DOWN(A.ld, 4); j += 4) {
                    y1 += x[col[start + j]];
                    y2 += x[col[start + j + 1]];
                    y3 += x[col[start + j + 2]];
                    y4 += x[col[start + j + 3]];
                }
                for (; j < A.ld; ++j) {
                    y1 += x[col[start + j]];
                }
                y[i] -= y1 + y2 + y3 + y4;
            }
        }

    } // ELL_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_ELL_spmv_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
