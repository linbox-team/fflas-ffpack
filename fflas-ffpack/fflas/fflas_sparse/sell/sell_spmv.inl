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

#ifndef __FFLASFFPACK_fflas_sparse_sell_spmv_INL
#define __FFLASFFPACK_fflas_sparse_sell_spmv_INL

// #define SELL_DEBUG 1

namespace FFLAS {
    namespace sparse_details_impl {
        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::SELL> &A, typename Field::ConstElement_ptr x_,
                          typename Field::Element_ptr y_, FieldCategories::GenericTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(chunkSize, A.chunkSize, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.nChunks; ++i) {
                index_t start = st[i];
                index_t size = chunkSize[i];
                index_t j = 0;
                for (; j < size; j++) {
                    for (index_t k = 0; k < A.chunk; ++k) {
                        F.axpyin(y[i * A.chunk + k], dat[start + j * A.chunk + k], x[col[start + j * A.chunk + k]]);
                    }
                }
            }
        }

        // #ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

        template <class Field>
        inline void fspmv_simd(const Field &F, const Sparse<Field, SparseMatrix_t::SELL> &A,
                               typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                               FieldCategories::UnparametricTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(chunkSize, A.chunkSize, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            using simd = Simd<typename Field::Element>;
            using vect_t = typename simd::vect_t;
            for (index_t i = 0; i < A.nChunks; ++i) {
                index_t start = st[i];
                index_t size = chunkSize[i];
                vect_t x1, x2, y1, y2, dat1, dat2;
                y1 = simd::zero();
                y2 = simd::zero();
                index_t j = 0;
                for (; j < ROUND_DOWN(size, 2); j += 2) {
                    dat1 = simd::load(dat + start + j * A.chunk);
                    dat2 = simd::load(dat + start + (j + 1) * A.chunk);
                    x1 = simd::gather(x, col + start + j * A.chunk);
                    x2 = simd::gather(x, col + start + (j + 1) * A.chunk);
                    y1 = simd::fmadd(y1, dat1, x1);
                    y2 = simd::fmadd(y2, dat2, x2);
                }
                if (size % 2 != 0) {
                    dat1 = simd::load(dat + start + j * A.chunk);
                    x1 = simd::gather(x, col + start + j * A.chunk);
                    y1 = simd::fmadd(y1, dat1, x1);
                }
                simd::store(y + i * A.chunk, simd::add(simd::load(y + i * A.chunk), simd::add(y1, y2)));
            }
        }

        // #endif

        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::SELL> &A, typename Field::ConstElement_ptr x_,
                          typename Field::Element_ptr y_, FieldCategories::UnparametricTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(chunkSize, A.chunkSize, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.nChunks; ++i) {
                index_t start = st[i];
                index_t size = chunkSize[i];
                for (index_t j = 0; j < size; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(A.chunk, 4); k += 4) {
                        y[i * A.chunk + k] += dat[start + j * A.chunk + k] * x[col[start + j * A.chunk + k]];
                        y[i * A.chunk + k + 1] += dat[start + j * A.chunk + k + 1] * x[col[start + j * A.chunk + k + 1]];
                        y[i * A.chunk + k + 2] += dat[start + j * A.chunk + k + 2] * x[col[start + j * A.chunk + k + 2]];
                        y[i * A.chunk + k + 3] += dat[start + j * A.chunk + k + 3] * x[col[start + j * A.chunk + k + 3]];
                    }
                    for (; k < size; ++k) {
                        y[i * A.chunk + k] += dat[start + j * A.chunk + k] * x[col[start + j * A.chunk + k]];
                    }
                }
            }
        }

        // #ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
        template <class Field>
        inline void fspmv_simd(const Field &F, const Sparse<Field, SparseMatrix_t::SELL> &A,
                               typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_, const uint64_t kmax) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(chunkSize, A.chunkSize, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t chunk = A.chunk;
            using simd = Simd<typename Field::Element>;
            using vect_t = typename simd::vect_t;

            vect_t X, Y, D, C, Q, TMP, NEGP, INVP, MIN, MAX, P;
            double p = (typename Field::Element)F.characteristic();

            P = simd::set1(p);
            NEGP = simd::set1(-p);
            INVP = simd::set1(1 / p);
            MIN = simd::set1(F.minElement());
            MAX = simd::set1(F.maxElement());

            for (size_t i = 0; i < A.nChunks; ++i) {
                index_t j = 0;
                index_t j_loc = 0;
                Y = simd::load(y + i * chunk);
                index_t size = chunkSize[i];
                index_t start = st[i];
                index_t block = size / kmax;
                for (size_t l = 0; l < block; ++l) {
                    j_loc += kmax;
                    for (; j < j_loc; ++j) {
                        D = simd::load(dat + start + j * A.chunk);
                        X = simd::gather(x, col + start + j * A.chunk);
                        Y = simd::fmadd(Y, D, X);
                    }
                    simd::mod(Y, P, INVP, NEGP, MIN, MAX, Q, TMP);
                }
                for (; j < size; ++j) {
                    D = simd::load(dat + start + j * A.chunk);
                    X = simd::gather(x, col + start + j * A.chunk);
                    Y = simd::fmadd(Y, D, X);
                }
                simd::mod(Y, P, INVP, NEGP, MIN, MAX, Q, TMP);
                simd::store(y + i * A.chunk, Y);
            }
        }
        // #endif
        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::SELL> &A, typename Field::ConstElement_ptr x_,
                          typename Field::Element_ptr y_, const uint64_t kmax) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(chunkSize, A.chunkSize, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t chunk = A.chunk;
            for (size_t i = 0; i < A.nChunks; ++i) {
                index_t j = 0;
                index_t j_loc = 0;
                index_t size = chunkSize[i];
                index_t start = st[i];
                index_t block = size / kmax;
                for (size_t l = 0; l < block; ++l) {
                    j_loc += kmax;
                    for (; j < j_loc; ++j) {
                        size_t k = 0;
                        for (; k < ROUND_DOWN(A.chunk, 4); k += 4) {
                            y[i * A.chunk + k] += dat[start + j * chunk + k] * x[col[start + j * chunk + k]];
                            y[i * A.chunk + k + 1] += dat[start + j * chunk + k + 1] * x[col[start + j * chunk + k + 1]];
                            y[i * A.chunk + k + 2] += dat[start + j * chunk + k + 2] * x[col[start + j * chunk + k + 2]];
                            y[i * A.chunk + k + 3] += dat[start + j * chunk + k + 3] * x[col[start + j * chunk + k + 3]];
                        }
                        for (; k < size; ++k) {
                            y[i * A.chunk + k] += dat[start + j * chunk + k] * x[col[start + j * chunk + k]];
                        }
                    }
                    for (size_t k = 0; k < size; ++k) {
                        F.reduce(y[i * A.chunk + k]);
                    }
                }
                for (; j < size; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(A.chunk, 4); k += 4) {
                        y[i * A.chunk + k] += dat[start + j * chunk + k] * x[col[start + j * chunk + k]];
                        y[i * A.chunk + k + 1] += dat[start + j * chunk + k + 1] * x[col[start + j * chunk + k + 1]];
                        y[i * A.chunk + k + 2] += dat[start + j * chunk + k + 2] * x[col[start + j * chunk + k + 2]];
                        y[i * A.chunk + k + 3] += dat[start + j * chunk + k + 3] * x[col[start + j * chunk + k + 3]];
                    }
                    for (; k < size; ++k) {
                        y[i * A.chunk + k] += dat[start + j * chunk + k] * x[col[start + j * chunk + k]];
                    }
                }
                for (size_t k = 0; k < size; ++k) {
                    F.reduce(y[i * A.chunk + k]);
                }
            }
        }

        template <class Field>
        inline void fspmv_one(const Field &F, const Sparse<Field, SparseMatrix_t::SELL_ZO> &A,
                              typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                              FieldCategories::GenericTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(chunkSize, A.chunkSize, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.nChunks; ++i) {
                index_t start = st[i];
                index_t size = chunkSize[i];
                index_t j = 0;
                for (; j < size; j++) {
                    for (index_t k = 0; k < A.chunk; ++k) {
                        F.addin(y[i * A.chunk + k], x[col[start + j * A.chunk + k]]);
                    }
                }
            }
        }

        template <class Field>
        inline void fspmv_mone(const Field &F, const Sparse<Field, SparseMatrix_t::SELL_ZO> &A,
                               typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                               FieldCategories::GenericTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(chunkSize, A.chunkSize, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.nChunks; ++i) {
                index_t start = st[i];
                index_t size = chunkSize[i];
                index_t j = 0;
                for (; j < size; j++) {
                    for (index_t k = 0; k < A.chunk; ++k) {
                        F.subin(y[i * A.chunk + k], x[col[start + j * A.chunk + k]]);
                    }
                }
            }
        }

        // #ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
        template <class Field>
        inline void fspmv_one_simd(const Field &F, const Sparse<Field, SparseMatrix_t::SELL_ZO> &A,
                                   typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                                   FieldCategories::UnparametricTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(chunkSize, A.chunkSize, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);

            using simd = Simd<typename Field::Element>;
            using vect_t = typename simd::vect_t;
            for (index_t i = 0; i < A.nChunks; ++i) {
                index_t start = st[i];
                index_t size = chunkSize[i];
                vect_t x1, x2, y1, y2;
                y1 = simd::zero();
                y2 = simd::zero();
                index_t j = 0;
                for (; j < ROUND_DOWN(size, 2); j += 2) {
                    x1 = simd::gather(x, col + start + j * A.chunk);
                    x2 = simd::gather(x, col + start + (j + 1) * A.chunk);
                    y1 = simd::add(y1, x1);
                    y2 = simd::add(y2, x2);
                }
                if (size % 2 != 0) {
                    x1 = simd::gather(x, col + start + j * A.chunk);
                    y1 = simd::add(y1, x1);
                }
                simd::store(y + i * A.chunk, simd::add(simd::load(y + i * A.chunk), simd::add(y1, y2)));
            }
        }

        template <class Field>
        inline void fspmv_mone_simd(const Field &F, const Sparse<Field, SparseMatrix_t::SELL_ZO> &A,
                                    typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                                    FieldCategories::UnparametricTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(chunkSize, A.chunkSize, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);

            using simd = Simd<typename Field::Element>;
            using vect_t = typename simd::vect_t;
            for (index_t i = 0; i < A.nChunks; ++i) {
                index_t start = st[i];
                index_t size = chunkSize[i];
                vect_t x1, x2, y1, y2;
                y1 = simd::zero();
                y2 = simd::zero();
                index_t j = 0;
                for (; j < ROUND_DOWN(size, 2); j += 2) {
                    x1 = simd::gather(x, col + start + j * A.chunk);
                    x2 = simd::gather(x, col + start + (j + 1) * A.chunk);
                    y1 = simd::add(y1, x1);
                    y2 = simd::add(y2, x2);
                }
                if (size % 2 != 0) {
                    x1 = simd::gather(x, col + start + j * A.chunk);
                    y1 = simd::add(y1, x1);
                }
                simd::store(y + i * A.chunk, simd::sub(simd::load(y + i * A.chunk), simd::add(y1, y2)));
            }
        }

        // #endif

        template <class Field>
        inline void fspmv_one(const Field &F, const Sparse<Field, SparseMatrix_t::SELL_ZO> &A,
                              typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                              FieldCategories::UnparametricTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(chunkSize, A.chunkSize, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            auto chunk = A.chunk;
            for (index_t i = 0; i < A.nChunks; ++i) {
                index_t start = st[i];
                index_t size = chunkSize[i];
                for (index_t j = 0; j < size; j++) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(A.chunk, 4); k += 4) {
                        y[i * A.chunk + k] += x[col[start + j * chunk + k]];
                        y[i * A.chunk + k + 1] += x[col[start + j * chunk + k + 1]];
                        y[i * A.chunk + k + 2] += x[col[start + j * chunk + k + 2]];
                        y[i * A.chunk + k + 3] += x[col[start + j * chunk + k + 3]];
                    }
                    for (; k < size; ++k) {
                        y[i * A.chunk + k] += x[col[start + j * chunk + k]];
                    }
                }
            }
        }

        template <class Field>
        inline void fspmv_mone(const Field &F, const Sparse<Field, SparseMatrix_t::SELL_ZO> &A,
                               typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                               FieldCategories::UnparametricTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(chunkSize, A.chunkSize, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            auto chunk = A.chunk;
            for (index_t i = 0; i < A.nChunks; ++i) {
                index_t start = st[i];
                index_t size = chunkSize[i];
                for (index_t j = 0; j < size; j++) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(A.chunk, 4); k += 4) {
                        y[i * A.chunk + k] -= x[col[start + j * chunk + k]];
                        y[i * A.chunk + k + 1] -= x[col[start + j * chunk + k + 1]];
                        y[i * A.chunk + k + 2] -= x[col[start + j * chunk + k + 2]];
                        y[i * A.chunk + k + 3] -= x[col[start + j * chunk + k + 3]];
                    }
                    for (; k < size; ++k) {
                        y[i * A.chunk + k] -= x[col[start + j * chunk + k]];
                    }
                }
            }
        }

    } // SELL_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_SELL_spmv_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
