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

#ifndef __FFLASFFPACK_fflas_sparse_CSR_pspmv_INL
#define __FFLASFFPACK_fflas_sparse_CSR_pspmv_INL

#ifdef __FFLASFFPACK_USE_TBB
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#endif

#include <thread>

namespace FFLAS {
    namespace sparse_details_impl {
        template <class Field>
        inline void pfspmv(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A, typename Field::ConstElement_ptr x_,
                           typename Field::Element_ptr y_, FieldCategories::GenericTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
#if defined(__FFLASFFPACK_USE_TBB)
            int step = __FFLASFFPACK_CACHE_LINE_SIZE / sizeof(typename Field::Element);
            tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m, step),
                              [&F, &A, x, y, dat, col, st](const tbb::blocked_range<index_t> &r) {
                              for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
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
                              });
#else
            // The minimum size has to be a multiple of cache_line/sizeof(Element) to avoid
            // cache coherency problem (ex: 8 for double, 16 for float)
#pragma omp parallel for
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
#endif
        }

        template<class Field>
        inline void pfspmv_task(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A, typename Field::ConstElement_ptr x_,
                                typename Field::Element_ptr y_, const index_t iStart, const index_t iStop, FieldCategories::UnparametricTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for(index_t i = iStart ; i < iStop ; ++i){
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

        template <class Field>
        inline void pfspmv(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A, typename Field::ConstElement_ptr x_,
                           typename Field::Element_ptr y_, FieldCategories::UnparametricTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
#if defined(__FFLASFFPACK_USE_TBB)
            int step = __FFLASFFPACK_CACHE_LINE_SIZE / sizeof(typename Field::Element);
            tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m, step),
                              [&F, &A, x, y, dat, col, st](const tbb::blocked_range<index_t> &r) {
                              for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
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
                              });
#else
            // #pragma omp parallel for schedule(static, 8)
            //     for (index_t i = 0; i < A.m; ++i) {
            //         auto start = st[i], stop = st[i + 1];
            //         index_t j = 0;
            //         index_t diff = stop - start;
            //         typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
            //         for (; j < ROUND_DOWN(diff, 4); j += 4) {
            //             y1 += dat[start + j] * x[col[start + j]];
            //             y2 += dat[start + j + 1] * x[col[start + j + 1]];
            //             y3 += dat[start + j + 2] * x[col[start + j + 2]];
            //             y4 += dat[start + j + 3] * x[col[start + j + 3]];
            //         }
            //         for (; j < diff; ++j) {
            //             y1 += dat[start + j] * x[col[start + j]];
            //         }
            //         y[i] += y1 + y2 + y3 + y4;
            //     }
            std::vector<std::thread> pool(6);

#endif
        }

        template <class Field>
        inline void pfspmv(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A, typename Field::ConstElement_ptr x_,
                           typename Field::Element_ptr y_, const int64_t kmax) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
#if defined(__FFLASFFPACK_USE_TBB)
            int step = __FFLASFFPACK_CACHE_LINE_SIZE / sizeof(typename Field::Element);
            tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m, step),
                              [&F, &A, x, y, kmax, dat, col, st](const tbb::blocked_range<index_t> &r) {
                              for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
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
                              });
#else
#pragma omp parallel for
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
#endif
        }

        template <class Field>
        inline void pfspmv_one(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A,
                               typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                               FieldCategories::GenericTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            size_t am = A.m;
            SYNCH_GROUP(
                        FORBLOCK1D(it, am, SPLITTER(NUM_THREADS),
                                   TASK(MODE(CONSTREFERENCE(F) READ(col, st, x) READWRITE(y)),
                                        for (index_t i = it.begin(); i < it.end(); ++i) {
                                        auto start = st[i];
                                        auto stop = st[i + 1];
                                        index_t j = 0;
                                        index_t diff = stop - start;
                                        typename Field::Element y1;
                                        typename Field::Element y2;
                                        typename Field::Element y3;
                                        typename Field::Element y4;
                                        F.assign(y1, F.zero);
                                        F.assign(y2, F.zero);
                                        F.assign(y3, F.zero);
                                        F.assign(y4, F.zero);
                                        for (; j < ROUND_DOWN(diff, 4); j += 4) {
                                        F.addin(y1, x[col[start + j]]);
                                        F.addin(y2, x[col[start + j + 1]]);
                                        F.addin(y3, x[col[start + j + 2]]);
                                        F.addin(y4, x[col[start + j + 3]]);
                                        }
                                        for (; j < diff; ++j) {
                                            F.addin(y1, x[col[start + j]]);
                                        }
                                        F.addin(y[i], y1);
                                        F.addin(y[i], y2);
                                        F.addin(y[i], y3);
                                        F.addin(y[i], y4);
                                        }
            );
            );
            );
        }

        template <class Field>
        inline void pfspmv_mone(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A,
                                typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                                FieldCategories::GenericTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            size_t am = A.m;
            SYNCH_GROUP(
                        FORBLOCK1D(it, am, SPLITTER(NUM_THREADS),
                                   TASK(MODE(CONSTREFERENCE(F) READ(col, st, x) READWRITE(y)),
                                        for (index_t i = it.begin(); i < it.end(); ++i) {
                                        auto start = st[i];
                                        auto stop = st[i + 1];
                                        index_t j = 0;
                                        index_t diff = stop - start;
                                        typename Field::Element y1;
                                        typename Field::Element y2;
                                        typename Field::Element y3;
                                        typename Field::Element y4;
                                        F.assign(y1, F.zero);
                                        F.assign(y2, F.zero);
                                        F.assign(y3, F.zero);
                                        F.assign(y4, F.zero);
                                        for (; j < ROUND_DOWN(diff, 4); j += 4) {
                                        F.addin(y1, x[col[start + j]]);
                                        F.addin(y2, x[col[start + j + 1]]);
                                        F.addin(y3, x[col[start + j + 2]]);
                                        F.addin(y4, x[col[start + j + 3]]);
                                        }
                                        for (; j < diff; ++j) {
                                            F.addin(y1, x[col[start + j]]);
                                        }
                                        F.subin(y[i], y1);
                                        F.subin(y[i], y2);
                                        F.subin(y[i], y3);
                                        F.subin(y[i], y4);
                                        }
            );
            );
            );
        }

        template <class Field>
        inline void pfspmv_one(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A,
                               typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                               FieldCategories::UnparametricTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
#if defined(__FFLASFFPACK_USE_TBB)
            int step = __FFLASFFPACK_CACHE_LINE_SIZE / sizeof(typename Field::Element);
            tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m, step),
                              [&F, &A, x, y, col, st](const tbb::blocked_range<index_t> &r) {
                              for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
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
                              });
#else
#pragma omp parallel for
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
#endif
        }

        template <class Field>
        inline void pfspmv_mone(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A,
                                typename Field::ConstElement_ptr x_, typename Field::Element_ptr y_,
                                FieldCategories::UnparametricTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
#if defined(__FFLASFFPACK_USE_TBB)
            int step = __FFLASFFPACK_CACHE_LINE_SIZE / sizeof(typename Field::Element);
            tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m, step),
                              [&F, &A, x, y, col, st](const tbb::blocked_range<index_t> &r) {
                              for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
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
                              });
#else
#pragma omp parallel for
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
#endif
        }

    } // CSR_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_CSR_pspmv_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
