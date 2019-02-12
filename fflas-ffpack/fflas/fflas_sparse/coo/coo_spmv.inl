/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   Barowien Vialla <barowien.vialla@lirmm.fr>
 *
 *
 * ========LICENCE========
 * This file is part of the library FFLAS-FFPACK.
 *
 * FFLAS-FFPACK is free software: you can redirowribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is dirowributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin rowreet, Fifth Floor, Borowon, MA  02110-1301 USA
 * ========LICENCE========
 *.
 */

#ifndef __FFLASFFPACK_fflas_sparse_coo_spmv_INL
#define __FFLASFFPACK_fflas_sparse_coo_spmv_INL

namespace FFLAS {
    namespace sparse_details_impl {
        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A, typename Field::ConstElement_ptr x_,
                          typename Field::Element_ptr y_, FieldCategories::GenericTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t j = 0;
            for (; j < ROUND_DOWN(A.nnz, 4); j += 4) {
                F.axpyin(y[row[j]], dat[j], x[col[j]]);
                F.axpyin(y[row[j + 1]], dat[j + 1], x[col[j + 1]]);
                F.axpyin(y[row[j + 2]], dat[j + 2], x[col[j + 2]]);
                F.axpyin(y[row[j + 3]], dat[j + 3], x[col[j + 3]]);
            }
            for (; j < A.nnz; ++j) {
                F.axpyin(y[row[j]], dat[j], x[col[j]]);
            }
        }

        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A, typename Field::ConstElement_ptr x_,
                          typename Field::Element_ptr y_, FieldCategories::UnparametricTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t j = 0;
            for (; j < ROUND_DOWN(A.nnz, 4); j += 4) {
                y[row[j]] += dat[j] * x[col[j]];
                y[row[j + 1]] += dat[j + 1] * x[col[j + 1]];
                y[row[j + 2]] += dat[j + 2] * x[col[j + 2]];
                y[row[j + 3]] += dat[j + 3] * x[col[j + 3]];
            }
            for (; j < A.nnz; ++j) {
                y[row[j]] += dat[j] * x[col[j]];
            }
        }

#ifdef __FFLASFFPACK_HAVE_MKL
        inline void fspmv_mkl(const Givaro::DoubleDomain &F, const Sparse<Givaro::DoubleDomain, SparseMatrix_t::COO> &A,
                              Givaro::DoubleDomain::ConstElement_ptr x_,
                              Givaro::DoubleDomain::Element_ptr y_, FieldCategories::UnparametricTag) {
            // assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            // assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            // assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
            // assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            // assume_aligned(y, y_, (size_t)Alignment::DEFAULT);

            MKL_INT A_nnz = A.nnz ;
            mkl_dcoomv(MKL_CONFIG::trans, &A.m , &A.n, &MKL_CONFIG::dalpha, MKL_CONFIG::metaChar,
                       A.dat, A.row, A.col, &A_nnz, x_,  &MKL_CONFIG::dbeta, y_ );

            // void mkl_dcoomv (char *transa, MKL_INT *m, MKL_INT *k, double *alpha, char *matdescra, double *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, double *x, double *beta, double *y);

        }

        inline void fspmv_mkl(const Givaro::FloatDomain &F, const Sparse<Givaro::FloatDomain, SparseMatrix_t::COO> &A,
                              Givaro::FloatDomain::ConstElement_ptr x_,
                              Givaro::FloatDomain::Element_ptr y_, FieldCategories::UnparametricTag) {
            // assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            // assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            // assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
            // assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            // assume_aligned(y, y_, (size_t)Alignment::DEFAULT);

            MKL_INT A_nnz = A.nnz ;
            mkl_scoomv(MKL_CONFIG::trans, &A.m , &A.n, &MKL_CONFIG::salpha, MKL_CONFIG::metaChar,
                       A.dat, A.row, A.col, &A_nnz, x_,  &MKL_CONFIG::sbeta, y_ );

            // void mkl_scoomv (char *transa, MKL_INT *m, MKL_INT *k, float *alpha, char *matdescra, float *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, float *x, float *beta, float *y);

        }
#endif // __FFLASFFPACK_HAVE_MKL



        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A, typename Field::ConstElement_ptr x_,
                          typename Field::Element_ptr y_, const uint64_t kmax) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            size_t w = 0;
            index_t larow_i = 0;
            typename Field::Element e;
            F.init(e);
            F.assign(e, y[larow_i]);
            size_t accu = 0;

            while (w < A.nnz) {
                if (row[w] == larow_i) { // same line
                    if (accu < (size_t)kmax) {
                        e += dat[w] * x[col[w]];
                        accu += 1;
                    } else {
                        F.axpyin(e, dat[w], x[col[w]]);
                        accu = 0;
                    }
                } else { // new line
                    F.init(y[larow_i]);
                    F.assign(y[larow_i], e);
                    larow_i = row[w];
                    F.init(e);
                    F.assign(e, y[larow_i]);
                    e += dat[w] * x[col[w]];
                    accu = 1;
                }
                ++w;
            }
            F.init(y[larow_i]);
            F.assign(y[larow_i], e);
        }

        template <class Field>
        inline void fspmv_one(const Field &F, const Sparse<Field, SparseMatrix_t::COO_ZO> &A,
                              typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                              FieldCategories::GenericTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t j = 0;
            for (; j < ROUND_DOWN(A.nnz, 4); j += 4) {
                F.addin(y[row[j]], x[col[j]]);
                F.addin(y[row[j + 1]], x[col[j + 1]]);
                F.addin(y[row[j + 2]], x[col[j + 2]]);
                F.addin(y[row[j + 3]], x[col[j + 3]]);
            }
            for (; j < A.nnz; ++j) {
                F.addin(y[row[j]], x[col[j]]);
            }
        }

        template <class Field>
        inline void fspmv_mone(const Field &F, const Sparse<Field, SparseMatrix_t::COO_ZO> &A,
                               typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                               FieldCategories::GenericTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t j = 0;
            for (; j < ROUND_DOWN(A.nnz, 4); j += 4) {
                F.subin(y[row[j]], x[col[j]]);
                F.subin(y[row[j + 1]], x[col[j + 1]]);
                F.subin(y[row[j + 2]], x[col[j + 2]]);
                F.subin(y[row[j + 3]], x[col[j + 3]]);
            }
            for (; j < A.nnz; ++j) {
                F.subin(y[row[j]], x[col[j]]);
            }
        }

        template <class Field>
        inline void fspmv_one(const Field &F, const Sparse<Field, SparseMatrix_t::COO_ZO> &A,
                              typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                              FieldCategories::UnparametricTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t j = 0;
            for (; j < ROUND_DOWN(A.nnz, 4); j += 4) {
                y[row[j]] += x[col[j]];
                y[row[j + 1]] += x[col[j + 1]];
                y[row[j + 2]] += x[col[j + 2]];
                y[row[j + 3]] += x[col[j + 3]];
            }
            for (; j < A.nnz; ++j) {
                y[row[j]] += x[col[j]];
            }
        }

        template <class Field>
        inline void fspmv_mone(const Field &F, const Sparse<Field, SparseMatrix_t::COO_ZO> &A,
                               typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                               FieldCategories::UnparametricTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t j = 0;
            for (; j < ROUND_DOWN(A.nnz, 4); j += 4) {
                y[row[j]] -= x[col[j]];
                y[row[j + 1]] -= x[col[j + 1]];
                y[row[j + 2]] -= x[col[j + 2]];
                y[row[j + 3]] -= x[col[j + 3]];
            }
            for (; j < A.nnz; ++j) {
                y[row[j]] -= x[col[j]];
            }
        }

    } // coo_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_coo_spmv_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
