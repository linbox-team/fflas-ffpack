/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

#ifndef __FFLASFFPACK_fflas_sparse_CSR_pspmm_INL
#define __FFLASFFPACK_fflas_sparse_CSR_pspmm_INL

#ifdef __FFLASFFPACK_USE_TBB
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#endif

namespace FFLAS {
namespace sparse_details_impl {
template <class Field>
inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A,
                   int blockSize, typename Field::ConstElement_ptr x,
                   typename Field::Element_ptr y, FieldCategories::GenericTag) {
#ifdef __FFLASFFPACK_USE_TBB
    tbb::parallel_for(
        tbb::blocked_range<index_t>(0, A.m),
        [&F, &A, &x, &y, blockSize](const tbb::blocked_range<index_t> &r) {
            for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
                auto start = A.st[i], stop = A.st[i + 1];
                for (index_t j = start; j < stop; ++j) {
                    int k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        F.axpyin(y[i * blockSize + k], A.dat[j],
                                 x[A.col[j] * blockSize + k]);
                        F.axpyin(y[i * blockSize + k + 1], A.dat[j],
                                 x[A.col[j] * blockSize + k + 1]);
                        F.axpyin(y[i * blockSize + k + 2], A.dat[j],
                                 x[A.col[j] * blockSize + k + 2]);
                        F.axpyin(y[i * blockSize + k + 3], A.dat[j],
                                 x[A.col[j] * blockSize + k + 3]);
                    }
                    for (; k < blockSize; ++k)
                        F.axpyin(y[i * blockSize + k], A.dat[j],
                                 x[A.col[j] * blockSize + k]);
                }
            }
        });
#else
#pragma omp parallel for
    for (index_t i = 0; i < A.m; ++i) {
        auto start = A.st[i], stop = A.st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                F.axpyin(y[i * blockSize + k], A.dat[j],
                         x[A.col[j] * blockSize + k]);
                F.axpyin(y[i * blockSize + k + 1], A.dat[j],
                         x[A.col[j] * blockSize + k + 1]);
                F.axpyin(y[i * blockSize + k + 2], A.dat[j],
                         x[A.col[j] * blockSize + k + 2]);
                F.axpyin(y[i * blockSize + k + 3], A.dat[j],
                         x[A.col[j] * blockSize + k + 3]);
            }
            for (; k < blockSize; ++k)
                F.axpyin(y[i * blockSize + k], A.dat[j],
                         x[A.col[j] * blockSize + k]);
        }
    }
#endif
}

template <class Field>
inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A,
                   int blockSize, typename Field::ConstElement_ptr x, int ldx,
                   typename Field::Element_ptr y, int ldy,
                   FieldCategories::GenericTag) {
#ifdef __FFLASFFPACK_USE_TBB
    tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m),
                      [&F, &A, &x, &y, ldx, ldy, blockSize](
                          const tbb::blocked_range<index_t> &r) {
        for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
            auto start = A.st[i], stop = A.st[i + 1];
            for (index_t j = start; j < stop; ++j) {
                typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
                int k = 0;
                for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                    F.axpyin(y[i * ldy + k], A.dat[j], x[A.col[j] * ldx + k]);
                    F.axpyin(y[i * ldy + k + 1], A.dat[j],
                             x[A.col[j] * ldx + k + 1]);
                    F.axpyin(y[i * ldy + k + 2], A.dat[j],
                             x[A.col[j] * ldx + k + 2]);
                    F.axpyin(y[i * ldy + k + 3], A.dat[j],
                             x[A.col[j] * ldx + k + 3]);
                }
                for (; k < blockSize; ++k)
                    F.axpyin(y[i * ldy + k], A.dat[j], x[A.col[j] * ldx + k]);
            }
        }
    });
#else
#pragma omp parallel for
    for (index_t i = 0; i < A.m; ++i) {
        auto start = A.st[i], stop = A.st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                F.axpyin(y[i * ldy + k], A.dat[j], x[A.col[j] * ldx + k]);
                F.axpyin(y[i * ldy + k + 1], A.dat[j],
                         x[A.col[j] * ldx + k + 1]);
                F.axpyin(y[i * ldy + k + 2], A.dat[j],
                         x[A.col[j] * ldx + k + 2]);
                F.axpyin(y[i * ldy + k + 3], A.dat[j],
                         x[A.col[j] * ldx + k + 3]);
            }
            for (; k < blockSize; ++k)
                F.axpyin(y[i * ldy + k], A.dat[j], x[A.col[j] * ldx + k]);
        }
    }
#endif
}

template <class Field>
inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A,
                   int blockSize, typename Field::ConstElement_ptr x,
                   typename Field::Element_ptr y,
                   FieldCategories::UnparametricTag) {
#ifdef __FFLASFFPACK_USE_TBB
    tbb::parallel_for(
        tbb::blocked_range<index_t>(0, A.m),
        [&F, &A, &x, &y, blockSize](const tbb::blocked_range<index_t> &r) {
            for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
                auto start = A.st[i], stop = A.st[i + 1];
                for (index_t j = start; j < stop; ++j) {
                    int k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        y[i * blockSize + k] +=
                            A.dat[j] * x[A.col[j] * blockSize + k];
                        y[i * blockSize + k + 1] +=
                            A.dat[j] * x[A.col[j] * blockSize + k + 1];
                        y[i * blockSize + k + 2] +=
                            A.dat[j] * x[A.col[j] * blockSize + k + 2];
                        y[i * blockSize + k + 3] +=
                            A.dat[j] * x[A.col[j] * blockSize + k + 3];
                    }
                    for (; k < blockSize; ++k)
                        y[i * blockSize + k] +=
                            A.dat[j] * x[A.col[j] * blockSize + k];
                }
            }
        });
#else
#pragma omp parallel for
    for (index_t i = 0; i < A.m; ++i) {
        auto start = A.st[i], stop = A.st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                y[i * blockSize + k] += A.dat[j] * x[A.col[j] * blockSize + k];
                y[i * blockSize + k + 1] +=
                    A.dat[j] * x[A.col[j] * blockSize + k + 1];
                y[i * blockSize + k + 2] +=
                    A.dat[j] * x[A.col[j] * blockSize + k + 2];
                y[i * blockSize + k + 3] +=
                    A.dat[j] * x[A.col[j] * blockSize + k + 3];
            }
            for (; k < blockSize; ++k)
                y[i * blockSize + k] += A.dat[j] * x[A.col[j] * blockSize + k];
        }
    }
#endif
}

template <class Field>
inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A,
                   int blockSize, typename Field::ConstElement_ptr x, int ldx,
                   typename Field::Element_ptr y, int ldy,
                   FieldCategories::UnparametricTag) {
#ifdef __FFLASFFPACK_USE_TBB
    tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m),
                      [&F, &A, &x, &y, ldx, ldy, blockSize](
                          const tbb::blocked_range<index_t> &r) {
        for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
            auto start = A.st[i], stop = A.st[i + 1];
            for (index_t j = start; j < stop; ++j) {
                typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
                int k = 0;
                for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                    y[i * ldy + k] += A.dat[j] * x[A.col[j] * ldx + k];
                    y[i * ldy + k + 1] += A.dat[j] * x[A.col[j] * ldx + k + 1];
                    y[i * ldy + k + 2] += A.dat[j] * x[A.col[j] * ldx + k + 2];
                    y[i * ldy + k + 3] += A.dat[j] * x[A.col[j] * ldx + k + 3];
                }
                for (; k < blockSize; ++k)
                    y[i * ldy + k] += A.dat[j] * x[A.col[j] * ldx + k];
            }
        }
    });
#else
#pragma omp parallel for
    for (index_t i = 0; i < A.m; ++i) {
        auto start = A.st[i], stop = A.st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                y[i * ldy + k] += A.dat[j] * x[A.col[j] * ldx + k];
                y[i * ldy + k + 1] += A.dat[j] * x[A.col[j] * ldx + k + 1];
                y[i * ldy + k + 2] += A.dat[j] * x[A.col[j] * ldx + k + 2];
                y[i * ldy + k + 3] += A.dat[j] * x[A.col[j] * ldx + k + 3];
            }
            for (; k < blockSize; ++k)
                y[i * ldy + k] += A.dat[j] * x[A.col[j] * ldx + k];
        }
    }
#endif
}

#ifdef __FFLASFFPACK_USE_SIMD
template <class Field, class LFunc, class SFunc>
inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A,
                   int blockSize, typename Field::ConstElement_ptr x,
                   typename Field::Element_ptr y, LFunc &&lfunc, SFunc &&sfunc,
                   FieldCategories::UnparametricTag) {
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_t;

#ifdef __FFLASFFPACK_USE_TBB
    tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m),
                      [&F, &A, &x, &y, lfunc, sfunc, blockSize](
                          const tbb::blocked_range<index_t> &r) {
        for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
            auto start = A.st[i], stop = A.st[i + 1];
            for (index_t j = start; j < stop; ++j) {
                vect_t y1, x1, y2, x2, dat;
                y1 = simd::zero();
                y2 = simd::zero();
                int k = 0;
                dat = simd::set1(A.dat[j]);
                for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size);
                     k += 2 * simd::vect_size) {
                    x1 = lfunc(x + A.col[j] * blockSize + k);
                    x2 = lfunc(x + A.col[j] * blockSize + k + simd::vect_size);
                    y1 = simd::fmadd(y1, x1, dat);
                    y2 = simd::fmadd(y2, x2, dat);
                    sfunc(y + i * blockSize + k, y1);
                    sfunc(y + i * blockSize + k + simd::vect_size, y2);
                }
                for (; k < ROUND_DOWN(blockSize, simd::vect_size);
                     k += simd::vect_size) {
                    x1 = lfunc(x + A.col[j] * blockSize + k);
                    y1 = simd::fmadd(y1, x1, dat);
                    sfunc(y + i * blockSize + k, y1);
                }
                for (; k < blockSize; ++k) {
                    y[i * blockSize + k] +=
                        A.dat[j] * x[A.col[j] * blockSize + k];
                }
            }
        }
    });
#else
#pragma omp parallel for
    for (index_t i = 0; i < A.m; ++i) {
        auto start = A.st[i], stop = A.st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            vect_t y1, x1, y2, x2, dat;
            y1 = simd::zero();
            y2 = simd::zero();
            int k = 0;
            dat = simd::set1(A.dat[j]);
            for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size);
                 k += 2 * simd::vect_size) {
                x1 = lfunc(x + A.col[j] * blockSize + k);
                x2 = lfunc(x + A.col[j] * blockSize + k + simd::vect_size);
                y1 = simd::fmadd(y1, x1, dat);
                y2 = simd::fmadd(y2, x2, dat);
                sfunc(y + i * blockSize + k, y1);
                sfunc(y + i * blockSize + k + simd::vect_size, y2);
            }
            for (; k < ROUND_DOWN(blockSize, simd::vect_size);
                 k += simd::vect_size) {
                x1 = lfunc(x + A.col[j] * blockSize + k);
                y1 = simd::fmadd(y1, x1, dat);
                sfunc(y + i * blockSize + k, y1);
            }
            for (; k < blockSize; ++k) {
                y[i * blockSize + k] += A.dat[j] * x[A.col[j] * blockSize + k];
            }
        }
    }
#endif
}

template <class Field, class LFunc, class SFunc>
inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A,
                   int blockSize, typename Field::ConstElement_ptr x, int ldx,
                   typename Field::Element_ptr y, int ldy, LFunc &&lfunc,
                   SFunc &&sfunc, FieldCategories::UnparametricTag) {
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_t;

#ifdef __FFLASFFPACK_USE_TBB
    tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m),
                      [&F, &A, &x, &y, ldx, ldy, lfunc, sfunc, blockSize](
                          const tbb::blocked_range<index_t> &r) {
        for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
            auto start = A.st[i], stop = A.st[i + 1];
            for (index_t j = start; j < stop; ++j) {
                vect_t y1, x1, y2, x2, dat;
                y1 = simd::zero();
                y2 = simd::zero();
                int k = 0;
                dat = simd::set1(A.dat[j]);
                for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size);
                     k += 2 * simd::vect_size) {
                    x1 = lfunc(x + A.col[j] * ldx + k);
                    x2 = lfunc(x + A.col[j] * ldx + k + simd::vect_size);
                    y1 = simd::fmadd(y1, x1, dat);
                    y2 = simd::fmadd(y2, x2, dat);
                    sfunc(y + i * ldy + k, y1);
                    sfunc(y + i * ldy + k + simd::vect_size, y2);
                }
                for (; k < ROUND_DOWN(blockSize, simd::vect_size);
                     k += simd::vect_size) {
                    x1 = lfunc(x + A.col[j] * ldx + k);
                    y1 = simd::fmadd(y1, x1, dat);
                    sfunc(y + i * ldy + k, y1);
                }
                for (; k < blockSize; ++k) {
                    y[i * ldy + k] += A.dat[j] * x[A.col[j] * ldx + k];
                }
            }
        }
    });
#else
#pragma omp parallel for
    for (index_t i = 0; i < A.m; ++i) {
        auto start = A.st[i], stop = A.st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            vect_t y1, x1, y2, x2, dat;
            y1 = simd::zero();
            y2 = simd::zero();
            int k = 0;
            dat = simd::set1(A.dat[j]);
            for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size);
                 k += 2 * simd::vect_size) {
                x1 = lfunc(x + A.col[j] * ldx + k);
                x2 = lfunc(x + A.col[j] * ldx + k + simd::vect_size);
                y1 = simd::fmadd(y1, x1, dat);
                y2 = simd::fmadd(y2, x2, dat);
                sfunc(y + i * ldy + k, y1);
                sfunc(y + i * ldy + k + simd::vect_size, y2);
            }
            for (; k < ROUND_DOWN(blockSize, simd::vect_size);
                 k += simd::vect_size) {
                x1 = lfunc(x + A.col[j] * ldx + k);
                y1 = simd::fmadd(y1, x1, dat);
                sfunc(y + i * ldy + k, y1);
            }
            for (; k < blockSize; ++k) {
                y[i * ldy + k] += A.dat[j] * x[A.col[j] * ldx + k];
            }
        }
    }
#endif
}
#endif

template <class Field>
inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A,
                   int blockSize, typename Field::ConstElement_ptr x,
                   typename Field::Element_ptr y, const int64_t kmax) {
#ifdef __FFLASFFPACK_USE_TBB
    tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m),
                      [&F, &A, &x, &y, kmax, blockSize](
                          const tbb::blocked_range<index_t> &r) {
        for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
            index_t j = A.st[i];
            index_t j_loc = j;
            index_t j_end = A.st[i + 1];
            index_t block = (j_end - j_loc) / kmax;
            for (index_t l = 0; l < (index_t)block; ++l) {
                j_loc += kmax;
                for (; j < j_loc; ++j) {
                    for (int k = 0; k < blockSize; ++k) {
                        y[i * blockSize + k] +=
                            A.dat[j] * x[A.col[j] * blockSize + k];
                    }
                }
                // TODO : replace with freduce
                for (int k = 0; k < blockSize; ++k) {
                    F.reduce(y[i * blockSize + k]);
                }
            }
            for (; j < j_end; ++j) {
                for (int k = 0; k < blockSize; ++k) {
                    y[i * blockSize + k] +=
                        A.dat[j] * x[A.col[j] * blockSize + k];
                }
            }
            for (int k = 0; k < blockSize; ++k) {
                F.reduce(y[i * blockSize + k]);
            }
        }
    });
#else
#pragma omp parallel for
    for (index_t i = 0; i < A.m; ++i) {
        index_t j = A.st[i];
        index_t j_loc = j;
        index_t j_end = A.st[i + 1];
        index_t block = (j_end - j_loc) / kmax;
        for (index_t l = 0; l < (index_t)block; ++l) {
            j_loc += kmax;
            for (; j < j_loc; ++j) {
                for (int k = 0; k < blockSize; ++k) {
                    y[i * blockSize + k] +=
                        A.dat[j] * x[A.col[j] * blockSize + k];
                }
            }
            // TODO : replace with freduce
            for (int k = 0; k < blockSize; ++k) {
                F.reduce(y[i * blockSize + k]);
            }
        }
        for (; j < j_end; ++j) {
            for (int k = 0; k < blockSize; ++k) {
                y[i * blockSize + k] += A.dat[j] * x[A.col[j] * blockSize + k];
            }
        }
        for (int k = 0; k < blockSize; ++k) {
            F.reduce(y[i * blockSize + k]);
        }
    }
#endif
}

template <class Field>
inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A,
                   int blockSize, typename Field::ConstElement_ptr x, int ldx,
                   typename Field::Element_ptr y, int ldy, const int64_t kmax) {
#ifdef __FFLASFFPACK_USE_TBB
    tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m),
                      [&F, &A, &x, &y, kmax, ldx, ldy, blockSize](
                          const tbb::blocked_range<index_t> &r) {
        for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
            index_t j = A.st[i];
            index_t j_loc = j;
            index_t j_end = A.st[i + 1];
            index_t block = (j_end - j_loc) / kmax;
            for (index_t l = 0; l < (index_t)block; ++l) {
                j_loc += kmax;
                for (; j < j_loc; ++j) {
                    for (int k = 0; k < blockSize; ++k) {
                        y[i * ldy + k] += A.dat[j] * x[A.col[j] * ldx + k];
                    }
                }
                // TODO : replace with freduce
                for (int k = 0; k < blockSize; ++k) {
                    F.reduce(y[i * ldy + k]);
                }
            }
            for (; j < j_end; ++j) {
                for (int k = 0; k < blockSize; ++k) {
                    y[i * ldy + k] += A.dat[j] * x[A.col[j] * ldx + k];
                }
            }
            for (int k = 0; k < blockSize; ++k) {
                F.reduce(y[i * ldy + k]);
            }
        }
    });
#else
#pragma omp parallel for
    for (index_t i = 0; i < A.m; ++i) {
        index_t j = A.st[i];
        index_t j_loc = j;
        index_t j_end = A.st[i + 1];
        index_t block = (j_end - j_loc) / kmax;
        for (index_t l = 0; l < (index_t)block; ++l) {
            j_loc += kmax;
            for (; j < j_loc; ++j) {
                for (int k = 0; k < blockSize; ++k) {
                    y[i * ldy + k] += A.dat[j] * x[A.col[j] * ldx + k];
                }
            }
            // TODO : replace with freduce
            for (int k = 0; k < blockSize; ++k) {
                F.reduce(y[i * ldy + k]);
            }
        }
        for (; j < j_end; ++j) {
            for (int k = 0; k < blockSize; ++k) {
                y[i * ldy + k] += A.dat[j] * x[A.col[j] * ldx + k];
            }
        }
        for (int k = 0; k < blockSize; ++k) {
            F.reduce(y[i * ldy + k]);
        }
    }
#endif
}

template <class Field, class Func>
inline void pfspmm_zo(const Field &F,
                      const Sparse<Field, SparseMatrix_t::CSR_ZO> &A,
                      int blockSize, typename Field::ConstElement_ptr x,
                      typename Field::Element_ptr y, Func &&func) {
#ifdef __FFLASFFPACK_USE_TBB
    tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m),
                      [&F, &A, &x, &y, blockSize, func](
                          const tbb::blocked_range<index_t> &r) {
        for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
            auto start = A.st[i], stop = A.st[i + 1];
            for (index_t j = start; j < stop; ++j) {
                int k = 0;
                for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                    func(y[i * blockSize + k], x[A.col[j] * blockSize + k]);
                    func(y[i * blockSize + k + 1],
                         x[A.col[j] * blockSize + k + 1]);
                    func(y[i * blockSize + k + 2],
                         x[A.col[j] * blockSize + k + 2]);
                    func(y[i * blockSize + k + 3],
                         x[A.col[j] * blockSize + k + 3]);
                }
                for (; k < blockSize; ++k)
                    func(y[i * blockSize + k], x[A.col[j] * blockSize + k]);
            }
        }
    });
#else
#pragma omp parallel for
    for (index_t i = 0; i < A.m; ++i) {
        auto start = A.st[i], stop = A.st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                func(y[i * blockSize + k], x[A.col[j] * blockSize + k]);
                func(y[i * blockSize + k + 1], x[A.col[j] * blockSize + k + 1]);
                func(y[i * blockSize + k + 2], x[A.col[j] * blockSize + k + 2]);
                func(y[i * blockSize + k + 3], x[A.col[j] * blockSize + k + 3]);
            }
            for (; k < blockSize; ++k)
                func(y[i * blockSize + k], x[A.col[j] * blockSize + k]);
        }
    }
#endif
}

template <class Field, class Func>
inline void
pfspmm_zo(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A,
          int blockSize, typename Field::ConstElement_ptr x, int ldx,
          typename Field::Element_ptr y, int ldy, Func &&func) {
#ifdef __FFLASFFPACK_USE_TBB
    tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m),
                      [&F, &A, &x, &y, blockSize, ldx, ldy, func](
                          const tbb::blocked_range<index_t> &r) {
        for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
            auto start = A.st[i], stop = A.st[i + 1];
            for (index_t j = start; j < stop; ++j) {
                int k = 0;
                for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                    func(y[i * ldy + k], x[A.col[j] * ldx + k]);
                    func(y[i * ldy + k + 1], x[A.col[j] * ldx + k + 1]);
                    func(y[i * ldy + k + 2], x[A.col[j] * ldx + k + 2]);
                    func(y[i * ldy + k + 3], x[A.col[j] * ldx + k + 3]);
                }
                for (; k < blockSize; ++k)
                    func(y[i * ldy + k], x[A.col[j] * ldx + k]);
            }
        }
    });
#else
#pragma omp parallel for
    for (index_t i = 0; i < A.m; ++i) {
        auto start = A.st[i], stop = A.st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                func(y[i * ldy + k], x[A.col[j] * ldx + k]);
                func(y[i * ldy + k + 1], x[A.col[j] * ldx + k + 1]);
                func(y[i * ldy + k + 2], x[A.col[j] * ldx + k + 2]);
                func(y[i * ldy + k + 3], x[A.col[j] * ldx + k + 3]);
            }
            for (; k < blockSize; ++k)
                func(y[i * ldy + k], x[A.col[j] * ldx + k]);
        }
    }
#endif
}

#ifdef __FFLASFFPACK_USE_SIMD

template <class Field, class LFunc, class SFunc, class FuncVect, class FuncScal>
inline void pfspmm_zo(const Field &F,
                      const Sparse<Field, SparseMatrix_t::CSR_ZO> &A,
                      int blockSize, typename Field::ConstElement_ptr x,
                      typename Field::Element_ptr y, FuncVect &&funcv,
                      FuncScal &&funcs, LFunc &&lfunc, SFunc &&sfunc) {
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_t;
#ifdef __FFLASFFPACK_USE_TBB
    tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m),
                      [&F, &A, &x, &y, blockSize, funcv, funcs, lfunc, sfunc](
                          const tbb::blocked_range<index_t> &r) {
        for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
            auto start = A.st[i], stop = A.st[i + 1];
            for (index_t j = start; j < stop; ++j) {
                vect_t y1, x1, y2, x2;
                y1 = simd::zero();
                y2 = simd::zero();
                int k = 0;
                for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size);
                     k += 2 * simd::vect_size) {
                    x1 = lfunc(x + A.col[j] * blockSize + k);
                    x2 = lfunc(x + A.col[j] * blockSize + k + simd::vect_size);
                    sfunc(y + i * blockSize + k, funcv(y1, x1));
                    sfunc(y + i * blockSize + k + simd::vect_size,
                          funcv(y2, x2));
                }
                for (; k < ROUND_DOWN(blockSize, simd::vect_size);
                     k += simd::vect_size) {
                    x1 = lfunc(x + A.col[j] * blockSize + k);
                    sfunc(y + i * blockSize + k, funcv(y1, x1));
                }
                for (; k < blockSize; ++k) {
                    funs(y[i * blockSize + k], x[A.col[j] * blockSize + k]);
                }
            }
        }
    });
#else
#pragma omp parallel for
    for (index_t i = 0; i < A.m; ++i) {
        auto start = A.st[i], stop = A.st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            vect_t y1, x1, y2, x2;
            y1 = simd::zero();
            y2 = simd::zero();
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size);
                 k += 2 * simd::vect_size) {
                x1 = lfunc(x + A.col[j] * blockSize + k);
                x2 = lfunc(x + A.col[j] * blockSize + k + simd::vect_size);
                sfunc(y + i * blockSize + k, funcv(y1, x1));
                sfunc(y + i * blockSize + k + simd::vect_size, funcv(y2, x2));
            }
            for (; k < ROUND_DOWN(blockSize, simd::vect_size);
                 k += simd::vect_size) {
                x1 = lfunc(x + A.col[j] * blockSize + k);
                sfunc(y + i * blockSize + k, funcv(y1, x1));
            }
            for (; k < blockSize; ++k) {
                funs(y[i * blockSize + k], x[A.col[j] * blockSize + k]);
            }
        }
    }
#endif
}

template <class Field, class LFunc, class SFunc, class FuncVect, class FuncScal>
inline void
pfspmm_zo(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A,
          int blockSize, typename Field::ConstElement_ptr x, int ldx,
          typename Field::Element_ptr y, int ldy, FuncVect &&funcv,
          FuncScal &&funs, LFunc &&lfunc, SFunc &&sfunc) {
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_t;
#ifdef __FFLASFFPACK_USE_TBB
    tbb::parallel_for(
        tbb::blocked_range<index_t>(0, A.m),
        [&F, &A, &x, &y, blockSize, ldx, ldy, funcv, funs, lfunc, sfunc](
            const tbb::blocked_range<index_t> &r) {
            for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
                auto start = A.st[i], stop = A.st[i + 1];
                for (index_t j = start; j < stop; ++j) {
                    vect_t y1, x1, y2, x2;
                    y1 = simd::zero();
                    y2 = simd::zero();
                    int k = 0;
                    for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size);
                         k += 2 * simd::vect_size) {
                        x1 = lfunc(x + A.col[j] * ldx + k);
                        x2 = lfunc(x + A.col[j] * ldx + k + simd::vect_size);
                        sfunc(y + i * ldy + k, funcv(y1, x1));
                        sfunc(y + i * ldy + k + simd::vect_size, funcv(y2, x2));
                    }
                    for (; k < ROUND_DOWN(blockSize, simd::vect_size);
                         k += simd::vect_size) {
                        x1 = lfunc(x + A.col[j] * ldx + k);
                        sfunc(y + i * ldy + k, funcv(y1, x1));
                    }
                    for (; k < blockSize; ++k) {
                        funcs(y[i * ldy + k], x[A.col[j] * ldx + k]);
                    }
                }
            }
        });
#else
#pragma omp parallel for
    for (index_t i = 0; i < A.m; ++i) {
        auto start = A.st[i], stop = A.st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            vect_t y1, x1, y2, x2;
            y1 = simd::zero();
            y2 = simd::zero();
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size);
                 k += 2 * simd::vect_size) {
                x1 = lfunc(x + A.col[j] * ldx + k);
                x2 = lfunc(x + A.col[j] * ldx + k + simd::vect_size);
                sfunc(y + i * ldy + k, funcv(y1, x1));
                sfunc(y + i * ldy + k + simd::vect_size, funcv(y2, x2));
            }
            for (; k < ROUND_DOWN(blockSize, simd::vect_size);
                 k += simd::vect_size) {
                x1 = lfunc(x + A.col[j] * ldx + k);
                sfunc(y + i * ldy + k, funcv(y1, x1));
            }
            for (; k < blockSize; ++k) {
                funcs(y[i * ldy + k], x[A.col[j] * ldx + k]);
            }
        }
    }
#endif
}

#endif //__FFLASFFPACK_USE_SIMD

} // CSR_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_CSR_pspmm_INL