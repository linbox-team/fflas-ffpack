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

#ifndef __FFLASFFPACK_fflas_sparse_CSR_spmv_INL
#define __FFLASFFPACK_fflas_sparse_CSR_spmv_INL


namespace FFLAS {
    namespace sparse_details_impl {
        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A, typename Field::ConstElement_ptr x_,
                          typename Field::Element_ptr y_, FieldCategories::GenericTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.m; ++i) {
                auto start = st[i], stop = st[i + 1];
                index_t j = 0;
                index_t diff = stop - start;
                typename Field::Element y1, y2, y3, y4;
                F.assign(y1, F.zero);
                F.assign(y2, F.zero);
                F.assign(y3, F.zero);
                F.assign(y4, F.zero);
                for (; j < ROUND_DOWN(diff, 4); j += 4) {
                    F.axpyin(y1, dat[start + j], x[col[start + j]]);
                    F.axpyin(y2, dat[start + j + 1], x[col[start + j + 1]]);
                    F.axpyin(y3, dat[start + j + 2], x[col[start + j + 2]]);
                    F.axpyin(y4, dat[start + j + 3], x[col[start + j + 3]]);
                }
                for (; j < diff; ++j) {
                    F.axpyin(y1, dat[start + j], x[col[start + j]]);
                }
                F.addin(y[i], y1);
                F.addin(y[i], y2);
                F.addin(y[i], y3);
                F.addin(y[i], y4);
            }
        }

#if 0
        template <class Field>
        inline void fspmv_task(const Field &F, const index_t start_, const index_t size_ const Sparse<Field, SparseMatrix_t::CSR> &A, typename Field::ConstElement_ptr x_,
                               typename Field::Element_ptr y_, FieldCategories::GenericTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = start_; i < start_+size_; ++i) {
                auto start = st[i], stop = st[i + 1];
                index_t j = 0;
                index_t diff = stop - start;
                typename Field::Element y1, y2, y3, y4;
                F.assign(y1, F.zero);
                F.assign(y2, F.zero);
                F.assign(y3, F.zero);
                F.assign(y4, F.zero);
                for (; j < ROUND_DOWN(diff, 4); j += 4) {
                    F.axpyin(y1, dat[start + j], x[col[start + j]]);
                    F.axpyin(y2, dat[start + j + 1], x[col[start + j + 1]]);
                    F.axpyin(y3, dat[start + j + 2], x[col[start + j + 2]]);
                    F.axpyin(y4, dat[start + j + 3], x[col[start + j + 3]]);
                }
                for (; j < diff; ++j) {
                    F.axpyin(y1, dat[start + j], x[col[start + j]]);
                }
                F.addin(y[i], y1);
                F.addin(y[i], y2);
                F.addin(y[i], y3);
                F.addin(y[i], y4);
            }
        }
#endif

        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A, typename Field::ConstElement_ptr x_,
                          typename Field::Element_ptr y_, FieldCategories::UnparametricTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);


            for (index_t i = 0; i < A.m; ++i) {
                auto start = st[i], stop = st[i + 1];
                index_t j = 0;
                index_t diff = stop - start;
                typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
                for (; j < ROUND_DOWN(diff, 4); j += 4) {
                    y1 += dat[start + j] * x[col[start + j]];
                    y2 += dat[start + j + 1] * x[col[start + j + 1]];
                    y3 += dat[start + j + 2] * x[col[start + j + 2]];
                    y4 += dat[start + j + 3] * x[col[start + j + 3]];
                }
                for (; j < diff; ++j) {
                    y1 += dat[start + j] * x[col[start + j]];
                }
                y[i] += y1 + y2 + y3 + y4;
            }
        }

#ifdef __FFLASFFPACK_HAVE_MKL
        inline void fspmv_mkl(const Givaro::DoubleDomain &F, const Sparse<Givaro::DoubleDomain, SparseMatrix_t::CSR> &A,
                              Givaro::DoubleDomain::ConstElement_ptr x_,
                              Givaro::DoubleDomain::Element_ptr y_, FieldCategories::UnparametricTag) {
            // assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            // assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            // assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            // assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            // assume_aligned(y, y_, (size_t)Alignment::DEFAULT);

            mkl_dcsrmv(MKL_CONFIG::trans, &A.m , &A.n, &MKL_CONFIG::dalpha, MKL_CONFIG::metaChar,
                       A.dat, A.col, A.st, A.st+1, x_,  &MKL_CONFIG::dbeta, y_ );

            // void mkl_dcsrmv (char *transa, MKL_INT *m, MKL_INT *k, double *alpha, char *matdescra, double *val, MKL_INT *indx, MKL_INT *pntrb, MKL_INT *pntre, double *x, double *beta, double *y);

        }

        inline void fspmv_mkl(const Givaro::FloatDomain &F, const Sparse<Givaro::FloatDomain, SparseMatrix_t::CSR> &A,
                              Givaro::FloatDomain::ConstElement_ptr x_,
                              Givaro::FloatDomain::Element_ptr y_, FieldCategories::UnparametricTag) {
            // assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            // assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            // assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            // assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            // assume_aligned(y, y_, (size_t)Alignment::DEFAULT);

            mkl_scsrmv(MKL_CONFIG::trans, &A.m , &A.n, &MKL_CONFIG::salpha, MKL_CONFIG::metaChar,
                       A.dat, A.col, A.st, A.st+1, x_,  &MKL_CONFIG::sbeta, y_ );

            // void mkl_scsrmv (char *transa, MKL_INT *m, MKL_INT *k, float *alpha, char *matdescra, float *val, MKL_INT *indx, MKL_INT *pntrb, MKL_INT *pntre, float *x, float *beta, float *y);

        }
#endif // __FFLASFFPACK_HAVE_MKL

#if 0
        template <class Field>
        inline void fspmv_task(const Field &F, const index_t start_, const index_t size_, const Sparse<Field, SparseMatrix_t::CSR> &A, typename Field::ConstElement_ptr x_,
                               typename Field::Element_ptr y_, FieldCategories::UnparametricTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);

            for (index_t i = start_; i < start_+size_; ++i) {
                auto start = st[i], stop = st[i + 1];
                index_t j = 0;
                index_t diff = stop - start;
                typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
                for (; j < ROUND_DOWN(diff, 4); j += 4) {
                    y1 += dat[start + j] * x[col[start + j]];
                    y2 += dat[start + j + 1] * x[col[start + j + 1]];
                    y3 += dat[start + j + 2] * x[col[start + j + 2]];
                    y4 += dat[start + j + 3] * x[col[start + j + 3]];
                }
                for (; j < diff; ++j) {
                    y1 += dat[start + j] * x[col[start + j]];
                }
                y[i] += y1 + y2 + y3 + y4;
            }
        }
#endif

        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A, typename Field::ConstElement_ptr x_,
                          typename Field::Element_ptr y_, const int64_t kmax) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);

            for (index_t i = 0; i < A.m; ++i) {
                index_t j = st[i];
                index_t j_loc = j;
                index_t j_end = st[i + 1];
                index_t block = (j_end - j_loc) / kmax;
                for (index_t l = 0; l < (index_t)block; ++l) {
                    j_loc += kmax;
                    for (; j < j_loc; ++j) {
                        y[i] += dat[j] * x[col[j]];
                    }
                    F.reduce(y[i]);
                }
                for (; j < j_end; ++j) {
                    y[i] += dat[j] * x[col[j]];
                }
                F.reduce(y[i]);
            }
        }

        template <class Field>
        inline void fspmv_one(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A,
                              typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                              FieldCategories::GenericTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.m; ++i) {
                auto start = st[i], stop = st[i + 1];
                index_t j = 0;
                index_t diff = stop - start;
                typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
                for (; j < ROUND_DOWN(diff, 4); j += 4) {
                    F.addin(y1, x[col[start + j]]);
                    F.addin(y2, x[col[start + j + 1]]);
                    F.addin(y3, x[col[start + j + 2]]);
                    F.addin(y4, x[col[start + j + 3]]);
                }
                for (; j < diff; ++j) {
                    F.addin(y1, x[col[start + j]]);
                }
                F.addin(y[i], y1 + y2 + y3 + y4);
            }
        }

        template <class Field>
        inline void fspmv_mone(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A,
                               typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                               FieldCategories::GenericTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.m; ++i) {
                auto start = st[i], stop = st[i + 1];
                index_t j = 0;
                index_t diff = stop - start;
                typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
                for (; j < ROUND_DOWN(diff, 4); j += 4) {
                    F.addin(y1, x[col[start + j]]);
                    F.addin(y2, x[col[start + j + 1]]);
                    F.addin(y3, x[col[start + j + 2]]);
                    F.addin(y4, x[col[start + j + 3]]);
                }
                for (; j < diff; ++j) {
                    F.addin(y1, x[col[start + j]]);
                }
                F.subin(y[i], y1 + y2 + y3 + y4);
            }
        }

        template <class Field>
        inline void fspmv_one(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A,
                              typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                              FieldCategories::UnparametricTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.m; ++i) {
                auto start = st[i], stop = st[i + 1];
                index_t j = 0;
                index_t diff = stop - start;
                typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
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
            }
        }

        template <class Field>
        inline void fspmv_mone(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A,
                               typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                               FieldCategories::UnparametricTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.m; ++i) {
                auto start = st[i], stop = st[i + 1];
                index_t j = 0;
                index_t diff = stop - start;
                typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
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
            }
        }

    } // sparse_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_CSR_spmv_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
