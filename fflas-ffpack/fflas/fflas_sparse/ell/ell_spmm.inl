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

#ifndef __FFLASFFPACK_fflas_sparse_ELL_spmm_INL
#define __FFLASFFPACK_fflas_sparse_ELL_spmm_INL

namespace FFLAS {
    namespace sparse_details_impl {

        template <class Field>
        inline void fspmm(const Field &F, const Sparse<Field, SparseMatrix_t::ELL> &A, size_t blockSize,
                          typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                          FieldCategories::GenericTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.m; ++i) {
                for (index_t j = 0; j < A.ld; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        F.axpyin(y[i * ldy + k], dat[i * A.ld + j], x[col[i * A.ld + j] * ldx + k]);
                        F.axpyin(y[i * ldy + k + 1], dat[i * A.ld + j], x[col[i * A.ld + j] * ldx + k + 1]);
                        F.axpyin(y[i * ldy + k + 2], dat[i * A.ld + j], x[col[i * A.ld + j] * ldx + k + 2]);
                        F.axpyin(y[i * ldy + k + 3], dat[i * A.ld + j], x[col[i * A.ld + j] * ldx + k + 3]);
                    }
                    for (; k < blockSize; ++k)
                        F.axpyin(y[i * ldy + k], dat[i * A.ld + j], x[col[i * A.ld + j] * ldx + k]);
                }
            }
        }

        template <class Field>
        inline void fspmm(const Field &F, const Sparse<Field, SparseMatrix_t::ELL> &A, size_t blockSize,
                          typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                          FieldCategories::UnparametricTag) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);

            for (index_t i = 0; i < A.m; ++i) {
                for (index_t j = 0; j < A.ld; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        y[i * ldy + k] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k];
                        y[i * ldy + k + 1] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k + 1];
                        y[i * ldy + k + 2] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k + 2];
                        y[i * ldy + k + 3] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k + 3];
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k];
                }
            }
        }

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

        template <class Field>
        inline void fspmm_simd_aligned(const Field &F, const Sparse<Field, SparseMatrix_t::ELL> &A, size_t blockSize,
                                       typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                                       FieldCategories::UnparametricTag) {
            using simd = Simd<typename Field::Element>;
            using vect_t = typename simd::vect_t;
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);

            for (index_t i = 0; i < A.m; ++i) {
                for (index_t j = 0; j < A.ld; ++j) {
                    vect_t vx1, vx2, vy1, vy2, vdat;
                    size_t k = 0;
                    vdat = simd::set1(dat[i * A.ld + j]);
                    for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                        vy1 = simd::load(y + i * ldy + k);
                        vy2 = simd::load(y + i * ldy + k + simd::vect_size);
                        vx1 = simd::load(x + col[i * A.ld + j] * ldx + k);
                        vx2 = simd::load(x + col[i * A.ld + j] * ldx + k + simd::vect_size);
                        simd::store(y + i * ldy + k, simd::fmadd(vy1, vx1, vdat));
                        simd::store(y + i * ldy + k + simd::vect_size, simd::fmadd(vy2, vx2, vdat));
                    }
                    for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                        vy1 = simd::load(y + i * ldy + k);
                        vx1 = simd::load(x + col[i * A.ld + j] * ldx + k);
                        simd::store(y + i * ldy + k, simd::fmadd(vy1, vx1, vdat));
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k];
                }
            }
        }

        template <class Field>
        inline void fspmm_simd_unaligned(const Field &F, const Sparse<Field, SparseMatrix_t::ELL> &A, size_t blockSize,
                                         typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                                         FieldCategories::UnparametricTag) {
            using simd = Simd<typename Field::Element>;
            using vect_t = typename simd::vect_t;
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);

            for (index_t i = 0; i < A.m; ++i) {
                for (index_t j = 0; j < A.ld; ++j) {
                    vect_t vx1, vx2, vy1, vy2, vdat;
                    size_t k = 0;
                    vdat = simd::set1(dat[i * A.ld + j]);
                    for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                        vy1 = simd::loadu(y + i * ldy + k);
                        vy2 = simd::loadu(y + i * ldy + k + simd::vect_size);
                        vx1 = simd::loadu(x + col[i * A.ld + j] * ldx + k);
                        vx2 = simd::loadu(x + col[i * A.ld + j] * ldx + k + simd::vect_size);
                        simd::storeu(y + i * ldy + k, simd::fmadd(vy1, vx1, vdat));
                        simd::storeu(y + i * ldy + k + simd::vect_size, simd::fmadd(vy2, vx2, vdat));
                    }
                    for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                        vy1 = simd::loadu(y + i * ldy + k);
                        vx1 = simd::loadu(x + col[i * A.ld + j] * ldx + k);
                        simd::storeu(y + i * ldy + k, simd::fmadd(vy1, vx1, vdat));
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k];
                }
            }
        }

#endif

        template <class Field>
        inline void fspmm(const Field &F, const Sparse<Field, SparseMatrix_t::ELL> &A, size_t blockSize,
                          typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                          const int64_t kmax) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            index_t block = (A.ld) / kmax;
            for (index_t i = 0; i < A.m; ++i) {
                index_t j_loc = 0, j = 0;
                for (index_t l = 0; l < (index_t)block; ++l) {
                    j_loc += kmax;
                    for (; j < j_loc; ++j) {
                        size_t k = 0;
                        for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                            y[i * ldy + k] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k];
                            y[i * ldy + k + 1] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k + 1];
                            y[i * ldy + k + 2] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k + 2];
                            y[i * ldy + k + 3] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k + 3];
                        }
                        for (; k < blockSize; ++k) {
                            y[i * ldy + k] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k];
                        }
                    }
                    // TODO : replace with freduce
                    for (size_t k = 0; k < blockSize; ++k) {
                        F.reduce(y[i * ldy + k]);
                    }
                }
                for (; j < A.ld; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        y[i * ldy + k] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k];
                        y[i * ldy + k + 1] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k + 1];
                        y[i * ldy + k + 2] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k + 2];
                        y[i * ldy + k + 3] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k + 3];
                    }
                    for (; k < blockSize; ++k) {
                        y[i * ldy + k] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k];
                    }
                }
                // TODO : replace with freduce
                for (size_t k = 0; k < blockSize; ++k) {
                    F.reduce(y[i * ldy + k]);
                }
            }
        }

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

        template <class Field>
        inline void fspmm_simd_aligned(const Field &F, const Sparse<Field, SparseMatrix_t::ELL> &A, size_t blockSize,
                                       typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                                       const int64_t kmax) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            using simd = Simd<typename Field::Element>;
            using vect_t = typename simd::vect_t;
            index_t block = (A.ld) / kmax;
            for (index_t i = 0; i < A.m; ++i) {
                index_t j_loc = 0, j = 0;
                for (index_t l = 0; l < (index_t)block; ++l) {
                    j_loc += kmax;
                    for (; j < j_loc; ++j) {
                        vect_t vx1, vx2, vy1, vy2, vdat;
                        size_t k = 0;
                        vdat = simd::set1(dat[i * A.ld + j]);
                        for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                            vy1 = simd::load(y + i * ldy + k);
                            vy2 = simd::load(y + i * ldy + k + simd::vect_size);
                            vx1 = simd::load(x + col[i * A.ld + j] * ldx + k);
                            vx2 = simd::load(x + col[i * A.ld + j] * ldx + k + simd::vect_size);
                            simd::store(y + i * ldy + k, simd::fmadd(vy1, vx1, vdat));
                            simd::store(y + i * ldy + k + simd::vect_size, simd::fmadd(vy2, vx2, vdat));
                        }
                        for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                            vy1 = simd::load(y + i * ldy + k);
                            vx1 = simd::load(x + col[i * A.ld + j] * ldx + k);
                            simd::store(y + i * ldy + k, simd::fmadd(vy1, vx1, vdat));
                        }
                        for (; k < blockSize; ++k)
                            y[i * ldy + k] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k];
                    }
                    // TODO : replace with freduce
                    for (size_t k = 0; k < blockSize; ++k) {
                        F.reduce(y[i * ldy + k]);
                    }
                }
                for (; j < A.ld; ++j) {
                    vect_t vx1, vx2, vy1, vy2, vdat;
                    size_t k = 0;
                    vdat = simd::set1(dat[i * A.ld + j]);
                    for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                        vy1 = simd::load(y + i * ldy + k);
                        vy2 = simd::load(y + i * ldy + k + simd::vect_size);
                        vx1 = simd::load(x + col[i * A.ld + j] * ldx + k);
                        vx2 = simd::load(x + col[i * A.ld + j] * ldx + k + simd::vect_size);
                        simd::store(y + i * ldy + k, simd::fmadd(vy1, vx1, vdat));
                        simd::store(y + i * ldy + k + simd::vect_size, simd::fmadd(vy2, vx2, vdat));
                    }
                    for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                        vy1 = simd::load(y + i * ldy + k);
                        vx1 = simd::load(x + col[i * A.ld + j] * ldx + k);
                        simd::store(y + i * ldy + k, simd::fmadd(vy1, vx1, vdat));
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k];
                }
                // TODO : replace with freduce
                for (size_t k = 0; k < blockSize; ++k) {
                    F.reduce(y[i * ldy + k]);
                }
            }
        }

        template <class Field>
        inline void fspmm_simd_unaligned(const Field &F, const Sparse<Field, SparseMatrix_t::ELL> &A, size_t blockSize,
                                         typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                                         const int64_t kmax) {
            assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            using simd = Simd<typename Field::Element>;
            using vect_t = typename simd::vect_t;
            index_t block = (A.ld) / kmax;
            for (index_t i = 0; i < A.m; ++i) {
                index_t j_loc = 0, j = 0;
                for (index_t l = 0; l < (index_t)block; ++l) {
                    j_loc += kmax;
                    for (; j < j_loc; ++j) {
                        vect_t vx1, vx2, vy1, vy2, vdat;
                        size_t k = 0;
                        vdat = simd::set1(dat[i * A.ld + j]);
                        for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                            vy1 = simd::loadu(y + i * ldy + k);
                            vy2 = simd::loadu(y + i * ldy + k + simd::vect_size);
                            vx1 = simd::loadu(x + col[i * A.ld + j] * ldx + k);
                            vx2 = simd::loadu(x + col[i * A.ld + j] * ldx + k + simd::vect_size);
                            simd::storeu(y + i * ldy + k, simd::fmadd(vy1, vx1, vdat));
                            simd::storeu(y + i * ldy + k + simd::vect_size, simd::fmadd(vy2, vx2, vdat));
                        }
                        for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                            vy1 = simd::loadu(y + i * ldy + k);
                            vx1 = simd::loadu(x + col[i * A.ld + j] * ldx + k);
                            simd::storeu(y + i * ldy + k, simd::fmadd(vy1, vx1, vdat));
                        }
                        for (; k < blockSize; ++k)
                            y[i * ldy + k] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k];
                    }
                    // TODO : replace with freduce
                    for (size_t k = 0; k < blockSize; ++k) {
                        F.reduce(y[i * ldy + k]);
                    }
                }
                for (; j < A.ld; ++j) {
                    vect_t vx1, vx2, vy1, vy2, vdat;
                    size_t k = 0;
                    vdat = simd::set1(dat[i * A.ld + j]);
                    for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                        vy1 = simd::loadu(y + i * ldy + k);
                        vy2 = simd::loadu(y + i * ldy + k + simd::vect_size);
                        vx1 = simd::loadu(x + col[i * A.ld + j] * ldx + k);
                        vx2 = simd::loadu(x + col[i * A.ld + j] * ldx + k + simd::vect_size);
                        simd::storeu(y + i * ldy + k, simd::fmadd(vy1, vx1, vdat));
                        simd::storeu(y + i * ldy + k + simd::vect_size, simd::fmadd(vy2, vx2, vdat));
                    }
                    for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                        vy1 = simd::loadu(y + i * ldy + k);
                        vx1 = simd::loadu(x + col[i * A.ld + j] * ldx + k);
                        simd::storeu(y + i * ldy + k, simd::fmadd(vy1, vx1, vdat));
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] += dat[i * A.ld + j] * x[col[i * A.ld + j] * ldx + k];
                }
                // TODO : replace with freduce
                for (size_t k = 0; k < blockSize; ++k) {
                    F.reduce(y[i * ldy + k]);
                }
            }
        }

#endif // SIMD

        template <class Field>
        inline void fspmm_mone(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_ZO> &A, size_t blockSize,
                               typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                               FieldCategories::GenericTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.m; ++i) {
                for (index_t j = 0; j < A.ld; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        F.subin(y[i * ldy + k], x[col[i * A.ld + j] * ldx + k]);
                        F.subin(y[i * ldy + k + 1], x[col[i * A.ld + j] * ldx + k + 1]);
                        F.subin(y[i * ldy + k + 2], x[col[i * A.ld + j] * ldx + k + 2]);
                        F.subin(y[i * ldy + k + 3], x[col[i * A.ld + j] * ldx + k + 3]);
                    }
                    for (; k < blockSize; ++k)
                        F.subin(y[i * ldy + k], x[col[i * A.ld + j] * ldx + k]);
                }
            }
        }

        template <class Field>
        inline void fspmm_one(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_ZO> &A, size_t blockSize,
                              typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                              FieldCategories::GenericTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.m; ++i) {
                for (index_t j = 0; j < A.ld; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        F.addin(y[i * ldy + k], x[col[i * A.ld + j] * ldx + k]);
                        F.addin(y[i * ldy + k + 1], x[col[i * A.ld + j] * ldx + k + 1]);
                        F.addin(y[i * ldy + k + 2], x[col[i * A.ld + j] * ldx + k + 2]);
                        F.addin(y[i * ldy + k + 3], x[col[i * A.ld + j] * ldx + k + 3]);
                    }
                    for (; k < blockSize; ++k)
                        F.addin(y[i * ldy + k], x[col[i * A.ld + j] * ldx + k]);
                }
            }
        }

        template <class Field>
        inline void fspmm_mone(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_ZO> &A, size_t blockSize,
                               typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                               FieldCategories::UnparametricTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.m; ++i) {
                for (index_t j = 0; j < A.ld; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        y[i * ldy + k] -= x[col[i * A.ld + j] * ldx + k];
                        y[i * ldy + k + 1] -= x[col[i * A.ld + j] * ldx + k + 1];
                        y[i * ldy + k + 2] -= x[col[i * A.ld + j] * ldx + k + 2];
                        y[i * ldy + k + 3] -= x[col[i * A.ld + j] * ldx + k + 3];
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] -= x[col[i * A.ld + j] * ldx + k];
                }
            }
        }

        template <class Field>
        inline void fspmm_one(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_ZO> &A, size_t blockSize,
                              typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                              FieldCategories::UnparametricTag) {
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.m; ++i) {
                for (index_t j = 0; j < A.ld; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        y[i * ldy + k] += x[col[i * A.ld + j] * ldx + k];
                        y[i * ldy + k + 1] += x[col[i * A.ld + j] * ldx + k + 1];
                        y[i * ldy + k + 2] += x[col[i * A.ld + j] * ldx + k + 2];
                        y[i * ldy + k + 3] += x[col[i * A.ld + j] * ldx + k + 3];
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] += x[col[i * A.ld + j] * ldx + k];
                }
            }
        }

        // #ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

        template <class Field>
        inline void fspmm_one_simd_aligned(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_ZO> &A, size_t blockSize,
                                           typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_,
                                           int ldy, FieldCategories::UnparametricTag) {
            using simd = Simd<typename Field::Element>;
            using vect_t = typename simd::vect_t;
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.m; ++i) {
                for (index_t j = 0; j < A.ld; ++j) {
                    vect_t vx1, vx2, vy1, vy2;
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                        vy1 = simd::load(y + i * ldy + k);
                        vy2 = simd::load(y + i * ldy + k + simd::vect_size);
                        vx1 = simd::load(x + col[i * A.ld + j] * ldx + k);
                        vx2 = simd::load(x + col[i * A.ld + j] * ldx + k + simd::vect_size);
                        simd::store(y + i * ldx + k, simd::add(vy1, vx1));
                        simd::store(y + i * ldx + k + simd::vect_size, simd::add(vy2, vx2));
                    }
                    for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                        vy1 = simd::load(y + i * ldy + k);
                        vx1 = simd::load(x + col[i * A.ld + j] * ldy + k);
                        simd::store(y + i * ldx + k, simd::add(vy1, vx1));
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] += x[col[i * A.ld + j] * ldx + k];
                }
            }
        }

        template <class Field>
        inline void fspmm_one_simd_unaligned(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_ZO> &A, size_t blockSize,
                                             typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_,
                                             int ldy, FieldCategories::UnparametricTag) {
            using simd = Simd<typename Field::Element>;
            using vect_t = typename simd::vect_t;
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.m; ++i) {
                for (index_t j = 0; j < A.ld; ++j) {
                    vect_t vx1, vx2, vy1, vy2;
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                        vy1 = simd::loadu(y + i * ldy + k);
                        vy2 = simd::loadu(y + i * ldy + k + simd::vect_size);
                        vx1 = simd::loadu(x + col[i * A.ld + j] * ldx + k);
                        vx2 = simd::loadu(x + col[i * A.ld + j] * ldx + k + simd::vect_size);
                        simd::storeu(y + i * ldx + k, simd::add(vy1, vx1));
                        simd::storeu(y + i * ldx + k + simd::vect_size, simd::add(vy2, vx2));
                    }
                    for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                        vy1 = simd::loadu(y + i * ldy + k);
                        vx1 = simd::loadu(x + col[i * A.ld + j] * ldy + k);
                        simd::storeu(y + i * ldx + k, simd::add(vy1, vx1));
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] += x[col[i * A.ld + j] * ldx + k];
                }
            }
        }

        template <class Field>
        inline void fspmm_mone_simd_aligned(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_ZO> &A, size_t blockSize,
                                            typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_,
                                            int ldy, FieldCategories::UnparametricTag) {
            using simd = Simd<typename Field::Element>;
            using vect_t = typename simd::vect_t;
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.m; ++i) {
                for (index_t j = 0; j < A.ld; ++j) {
                    vect_t vx1, vx2, vy1, vy2;
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                        vy1 = simd::load(y + i * ldy + k);
                        vy2 = simd::load(y + i * ldy + k + simd::vect_size);
                        vx1 = simd::load(x + col[i * A.ld + j] * ldx + k);
                        vx2 = simd::load(x + col[i * A.ld + j] * ldx + k + simd::vect_size);
                        simd::store(y + i * ldx + k, simd::sub(vy1, vx1));
                        simd::store(y + i * ldx + k + simd::vect_size, simd::sub(vy2, vx2));
                    }
                    for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                        vy1 = simd::load(y + i * ldy + k);
                        vx1 = simd::load(x + col[i * A.ld + j] * ldy + k);
                        simd::store(y + i * ldx + k, simd::sub(vy1, vx1));
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] -= x[col[i * A.ld + j] * ldx + k];
                }
            }
        }

        template <class Field>
        inline void fspmm_mone_simd_unaligned(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_ZO> &A, size_t blockSize,
                                              typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_,
                                              int ldy, FieldCategories::UnparametricTag) {
            using simd = Simd<typename Field::Element>;
            using vect_t = typename simd::vect_t;
            assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
            assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
            assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
            for (index_t i = 0; i < A.m; ++i) {
                for (index_t j = 0; j < A.ld; ++j) {
                    vect_t vx1, vx2, vy1, vy2;
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                        vy1 = simd::loadu(y + i * ldy + k);
                        vy2 = simd::loadu(y + i * ldy + k + simd::vect_size);
                        vx1 = simd::loadu(x + col[i * A.ld + j] * ldx + k);
                        vx2 = simd::loadu(x + col[i * A.ld + j] * ldx + k + simd::vect_size);
                        simd::storeu(y + i * ldx + k, simd::sub(vy1, vx1));
                        simd::storeu(y + i * ldx + k + simd::vect_size, simd::sub(vy2, vx2));
                    }
                    for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                        vy1 = simd::loadu(y + i * ldy + k);
                        vx1 = simd::loadu(x + col[i * A.ld + j] * ldy + k);
                        simd::storeu(y + i * ldx + k, simd::sub(vy1, vx1));
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] -= x[col[i * A.ld + j] * ldx + k];
                }
            }
        }

        // #endif /*  __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS */

    } // ell_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_ELL_spmm_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
