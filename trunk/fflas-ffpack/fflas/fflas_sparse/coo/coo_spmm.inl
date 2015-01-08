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

#ifndef __FFLASFFPACK_fflas_sparse_coo_spmm_INL
#define __FFLASFFPACK_fflas_sparse_coo_spmm_INL

namespace FFLAS {
namespace sparse_details_impl {

// template <class Field>
// inline void fspmm(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A,
//                   int blockSize, typename Field::ConstElement_ptr x_,
//                   typename Field::Element_ptr y_, FieldCategories::GenericTag) {
//     assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
//     assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
//     for (index_t i = 0; i < A.nnz; ++i) {
//         int k = 0;
//         for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
//             F.axpyin(y[row[i] * blockSize + k], dat[i],
//                      x[col[i] * blockSize + k]);
//             F.axpyin(y[row[i] * blockSize + k + 1], dat[i],
//                      x[col[i] * blockSize + k + 1]);
//             F.axpyin(y[row[i] * blockSize + k + 2], dat[i],
//                      x[col[i] * blockSize + k + 2]);
//             F.axpyin(y[row[i] * blockSize + k + 3], dat[i],
//                      x[col[i] * blockSize + k + 3]);
//         }
//         for (; k < blockSize; ++k)
//             F.axpyin(y[row[i] * blockSize + k], dat[i],
//                      x[col[i] * blockSize + k]);
//     }
// }

template <class Field>
inline void fspmm(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A, int blockSize,
                  typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                  FieldCategories::GenericTag) {
    assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
    assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    for (index_t i = 0; i < A.nnz; ++i) {
        int k = 0;
        for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
            F.axpyin(y[row[i] * ldy + k], dat[i], x[col[i] * ldx + k]);
            F.axpyin(y[row[i] * ldy + k + 1], dat[i], x[col[i] * ldx + k + 1]);
            F.axpyin(y[row[i] * ldy + k + 2], dat[i], x[col[i] * ldx + k + 2]);
            F.axpyin(y[row[i] * ldy + k + 3], dat[i], x[col[i] * ldx + k + 3]);
        }
        for (; k < blockSize; ++k)
            F.axpyin(y[row[i] * ldy + k], dat[i], x[col[i] * ldx + k]);
    }
}

// template <class Field>
// inline void fspmm(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A,
//                   int blockSize, typename Field::ConstElement_ptr x_,
//                   typename Field::Element_ptr y_,
//                   FieldCategories::UnparametricTag) {
//     assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
//     assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
//     for (index_t i = 0; i < A.nnz; ++i) {
//         int k = 0;
//         for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
//             y[row[i] * blockSize + k] +=
//                 dat[i] * x[col[i] * blockSize + k];
//             y[row[i] * blockSize + k + 1] +=
//                 dat[i] * x[col[i] * blockSize + k + 1];
//             y[row[i] * blockSize + k + 2] +=
//                 dat[i] * x[col[i] * blockSize + k + 2];
//             y[row[i] * blockSize + k + 3] +=
//                 dat[i] * x[col[i] * blockSize + k + 3];
//         }
//         for (; k < blockSize; ++k)
//             y[row[i] * blockSize + k] +=
//                 dat[i] * x[col[i] * blockSize + k];
//     }
// }

template <class Field>
inline void fspmm(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A, int blockSize,
                  typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                  FieldCategories::UnparametricTag) {
    assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
    assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    for (index_t i = 0; i < A.nnz; ++i) {
        int k = 0;
        for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
            y[row[i] * ldy + k] += dat[i] * x[col[i] * ldx + k];
            y[row[i] * ldy + k + 1] += dat[i] * x[col[i] * ldx + k + 1];
            y[row[i] * ldy + k + 2] += dat[i] * x[col[i] * ldx + k + 2];
            y[row[i] * ldy + k + 3] += dat[i] * x[col[i] * ldx + k + 3];
        }
        for (; k < blockSize; ++k)
            y[row[i] * ldy + k] += dat[i] * x[col[i] * ldx + k];
    }
}

#ifdef __FFLASFFPACK_USE_SIMD

// template <class Field>
// inline void fspmm_simd_aligned(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A,
//                   int blockSize, typename Field::ConstElement_ptr x_,
//                   typename Field::Element_ptr y_,
//                   FieldCategories::UnparametricTag) {
//     assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
//     assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
//     using simd = Simd<typename Field::Element>;
//     using vect_t = typename simd::vect_t;

//     for (index_t i = 0; i < A.nnz; ++i) {
//         vect_t vy, vx, vdat;
//         vdat = simd::set1(dat[i]);
//         int k = 0;
//         for (; k < ROUND_DOWN(blockSize, simd::vect_size);
//              k += simd::vect_size) {
//             vy = simd::load(y + row[i] * blockSize + k);
//             vx = simd::load(x + col[i] * blockSize + k);
//             simd::store(y + row[i] * blockSize + k, simd::fmadd(vy, vdat, vx));
//         }
//         for (; k < blockSize; ++k)
//             y[row[i] * blockSize + k] +=
//                 dat[i] * x[col[i] * blockSize + k];
//     }
// }

// template <class Field>
// inline void fspmm_simd_unaligned(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A,
//                   int blockSize, typename Field::ConstElement_ptr x_,
//                   typename Field::Element_ptr y_,
//                   FieldCategories::UnparametricTag) {
//     assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
//     auto x = x_;
//     auto y = y_;
//     using simd = Simd<typename Field::Element>;
//     using vect_t = typename simd::vect_t;

//     for (index_t i = 0; i < A.nnz; ++i) {
//         vect_t vy, vx, vdat;
//         vdat = simd::set1(dat[i]);
//         int k = 0;
//         for (; k < ROUND_DOWN(blockSize, simd::vect_size);
//              k += simd::vect_size) {
//             vy = simd::loadu(y + row[i] * blockSize + k);
//             vx = simd::loadu(x + col[i] * blockSize + k);
//             simd::storeu(y + row[i] * blockSize + k, simd::fmadd(vy, vdat, vx));
//         }
//         for (; k < blockSize; ++k)
//             y[row[i] * blockSize + k] +=
//                 dat[i] * x[col[i] * blockSize + k];
//     }
// }

template <class Field>
inline void fspmm_simd_aligned(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A, int blockSize,
                               typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                               FieldCategories::UnparametricTag) {
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_t;
    assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
    assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    for (index_t i = 0; i < A.nnz; ++i) {
        vect_t vy, vx, vdat;
        vdat = simd::set1(dat[i]);
        int k = 0;
        for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
            vy = simd::load(y + row[i] * ldy + k);
            vx = simd::load(x + col[i] * ldx + k);
            simd::store(y + row[i] * ldy + k, simd::fmadd(vy, vdat, vx));
        }
        for (; k < blockSize; ++k)
            y[row[i] * ldy + k] += dat[i] * x[col[i] * ldx + k];
    }
}

template <class Field>
inline void fspmm_simd_unaligned(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A, int blockSize,
                                 typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                                 FieldCategories::UnparametricTag) {
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_t;
    assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
    assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    auto x = x_;
    auto y = y_;
    for (index_t i = 0; i < A.nnz; ++i) {
        vect_t vy, vx, vdat;
        vdat = simd::set1(dat[i]);
        int k = 0;
        for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
            vy = simd::loadu(y + row[i] * ldy + k);
            vx = simd::loadu(x + col[i] * ldx + k);
            simd::storeu(y + row[i] * ldy + k, simd::fmadd(vy, vdat, vx));
        }
        for (; k < blockSize; ++k)
            y[row[i] * ldy + k] += dat[i] * x[col[i] * ldx + k];
    }
}

#endif // SIMD

// template <class Field>
// inline void fspmm(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A,
//                   int blockSize, typename Field::ConstElement_ptr x_,
//                   typename Field::Element_ptr y_, const int64_t kmax) {
//     // TODO
// }

template <class Field>
inline void fspmm(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A, int blockSize,
                  typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                  const int64_t kmax) {
    // TODO
}

// template <class Field>
// inline void fspmm_one(const Field &F,
//                      const Sparse<Field, SparseMatrix_t::COO_ZO> &A,
//                      int blockSize, typename Field::ConstElement_ptr x_,
//                      typename Field::Element_ptr y_) {
//     assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
//     assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
//     for (index_t i = 0; i < A.nnz; ++i) {
//         int k = 0;
//         for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
//             F.addin(y[row[i] * blockSize + k], x[col[i] * blockSize + k]);
//             F.addin(y[row[i] * blockSize + k + 1],
//                  x[col[i] * blockSize + k + 1]);
//             F.addin(y[row[i] * blockSize + k + 2],
//                  x[col[i] * blockSize + k + 2]);
//             F.addin(y[row[i] * blockSize + k + 3],
//                  x[col[i] * blockSize + k + 3]);
//         }
//         for (; k < blockSize; ++k)
//             F.addin(y[row[i] * blockSize + k], x[col[i] * blockSize + k]);
//     }
// }

// template <class Field>
// inline void fspmm_mone(const Field &F,
//                      const Sparse<Field, SparseMatrix_t::COO_ZO> &A,
//                      int blockSize, typename Field::ConstElement_ptr x_,
//                      typename Field::Element_ptr y_) {
//     assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
//     assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
//     for (index_t i = 0; i < A.nnz; ++i) {
//         int k = 0;
//         for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
//             F.subin(y[row[i] * blockSize + k], x[col[i] * blockSize + k]);
//             F.subin(y[row[i] * blockSize + k + 1],
//                  x[col[i] * blockSize + k + 1]);
//             F.subin(y[row[i] * blockSize + k + 2],
//                  x[col[i] * blockSize + k + 2]);
//             F.subin(y[row[i] * blockSize + k + 3],
//                  x[col[i] * blockSize + k + 3]);
//         }
//         for (; k < blockSize; ++k)
//             F.subin(y[row[i] * blockSize + k], x[col[i] * blockSize + k]);
//     }
// }

template <class Field>
inline void fspmm_one(const Field &F, const Sparse<Field, SparseMatrix_t::COO_ZO> &A, int blockSize,
                      typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                      FieldCategories::GenericTag) {
    assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    for (index_t i = 0; i < A.nnz; ++i) {
        int k = 0;
        for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
            F.addin(y[row[i] * ldy + k], x[col[i] * ldx + k]);
            F.addin(y[row[i] * ldy + k + 1], x[col[i] * ldx + k + 1]);
            F.addin(y[row[i] * ldy + k + 2], x[col[i] * ldx + k + 2]);
            F.addin(y[row[i] * ldy + k + 3], x[col[i] * ldx + k + 3]);
        }
        for (; k < blockSize; ++k)
            F.addin(y[row[i] * ldy + k], x[col[i] * ldx + k]);
    }
}

template <class Field>
inline void fspmm_mone(const Field &F, const Sparse<Field, SparseMatrix_t::COO_ZO> &A, int blockSize,
                       typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                       FieldCategories::GenericTag) {
    assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    for (index_t i = 0; i < A.nnz; ++i) {
        int k = 0;
        for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
            F.subin(y[row[i] * ldy + k], x[col[i] * ldx + k]);
            F.subin(y[row[i] * ldy + k + 1], x[col[i] * ldx + k + 1]);
            F.subin(y[row[i] * ldy + k + 2], x[col[i] * ldx + k + 2]);
            F.subin(y[row[i] * ldy + k + 3], x[col[i] * ldx + k + 3]);
        }
        for (; k < blockSize; ++k)
            F.subin(y[row[i] * ldy + k], x[col[i] * ldx + k]);
    }
}

#ifdef __FFLASFFPACK_USE_SIMD

// template <class Field>
// inline void fspmm_one_simd_aligned(const Field &F,
//                      const Sparse<Field, SparseMatrix_t::COO> &A, int blockSize,
//                      typename Field::ConstElement_ptr x_,
//                      typename Field::Element_ptr y_, FieldCategories::UnparametricTag) {
//     using simd = Simd<typename Field::Element>;
//     using vect_t = typename simd::vect_size;
//     assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
//     assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
//     for (index_t i = 0; i < A.nnz; ++i) {
//         vect_t vy, vx;
//         int k = 0;
//         for (; k < ROUND_DOWN(blockSize, simd::vect_size);
//              k += simd::vect_size) {
//             vy = simd::load(y + row[i] * blockSize + k);
//             vx = simd::load(x + col[i] * blockSize + k);
//             simd::store(y + row[i] * blockSize + k, simd::add(vy, vx));
//         }
//         for (; k < blockSize; ++k)
//             y[row[i] * blockSize + k] += x[col[i] * blockSize + k];
//     }
// }

// template <class Field>
// inline void fspmm_one_simd_unaligned(const Field &F,
//                      const Sparse<Field, SparseMatrix_t::COO> &A, int blockSize,
//                      typename Field::ConstElement_ptr x_,
//                      typename Field::Element_ptr y_, FieldCategories::UnparametricTag) {
//     using simd = Simd<typename Field::Element>;
//     using vect_t = typename simd::vect_size;
//     assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
//     assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
//     for (index_t i = 0; i < A.nnz; ++i) {
//         vect_t vy, vx;
//         int k = 0;
//         for (; k < ROUND_DOWN(blockSize, simd::vect_size);
//              k += simd::vect_size) {
//             vy = simd::loadu(y + row[i] * blockSize + k);
//             vx = simd::loadu(x + col[i] * blockSize + k);
//             simd::storeu(y + row[i] * blockSize + k, simd::add(vy, vx));
//         }
//         for (; k < blockSize; ++k)
//             y[row[i] * blockSize + k] += x[col[i] * blockSize + k];
//     }
// }

// template <class Field>
// inline void fspmm_mone_simd_aligned(const Field &F,
//                      const Sparse<Field, SparseMatrix_t::COO> &A, int blockSize,
//                      typename Field::ConstElement_ptr x_,
//                      typename Field::Element_ptr y_, FieldCategories::UnparametricTag) {
//     using simd = Simd<typename Field::Element>;
//     using vect_t = typename simd::vect_size;
//     assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
//     assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
//     for (index_t i = 0; i < A.nnz; ++i) {
//         vect_t vy, vx;
//         int k = 0;
//         for (; k < ROUND_DOWN(blockSize, simd::vect_size);
//              k += simd::vect_size) {
//             vy = simd::load(y + row[i] * blockSize + k);
//             vx = simd::load(x + col[i] * blockSize + k);
//             simd::store(y + row[i] * blockSize + k, simd::sub(vy, vx));
//         }
//         for (; k < blockSize; ++k)
//             y[row[i] * blockSize + k] -= x[col[i] * blockSize + k];
//     }
// }

// template <class Field>
// inline void fspmm_mone_simd_unaligned(const Field &F,
//                      const Sparse<Field, SparseMatrix_t::COO> &A, int blockSize,
//                      typename Field::ConstElement_ptr x_,
//                      typename Field::Element_ptr y_, FieldCategories::UnparametricTag) {
//     using simd = Simd<typename Field::Element>;
//     using vect_t = typename simd::vect_size;
//     assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
//     assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
//     assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
//     for (index_t i = 0; i < A.nnz; ++i) {
//         vect_t vy, vx;
//         int k = 0;
//         for (; k < ROUND_DOWN(blockSize, simd::vect_size);
//              k += simd::vect_size) {
//             vy = simd::loadu(y + row[i] * blockSize + k);
//             vx = simd::loadu(x + col[i] * blockSize + k);
//             simd::storeu(y + row[i] * blockSize + k, simd::sub(vy, vx));
//         }
//         for (; k < blockSize; ++k)
//             y[row[i] * blockSize + k] -= x[col[i] * blockSize + k];
//     }
// }

template <class Field>
inline void fspmm_one_simd_aligned(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A, int blockSize,
                                   typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_,
                                   int ldy) {
    assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_size;
    for (index_t i = 0; i < A.nnz; ++i) {
        vect_t vy, vx;
        int k = 0;
        for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
            vy = simd::load(y + row[i] * ldy + k);
            vx = simd::load(x + col[i] * ldx + k);
            simd::store(y + row[i] * ldy + k, simd::add(vy, vx));
        }
        for (; k < blockSize; ++k)
            y[row[i] * ldy + k] += x[col[i] * ldx + k];
    }
}

template <class Field>
inline void fspmm_one_simd_unaligned(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A, int blockSize,
                                     typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_,
                                     int ldy) {
    assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_size;
    for (index_t i = 0; i < A.nnz; ++i) {
        vect_t vy, vx;
        int k = 0;
        for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
            vy = simd::loadu(y + row[i] * ldy + k);
            vx = simd::loadu(x + col[i] * ldx + k);
            simd::storeu(y + row[i] * ldy + k, simd::add(vy, vx));
        }
        for (; k < blockSize; ++k)
            y[row[i] * ldy + k] += x[col[i] * ldx + k];
    }
}

template <class Field>
inline void fspmm_mone_simd_aligned(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A, int blockSize,
                                    typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_,
                                    int ldy) {
    assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_size;
    for (index_t i = 0; i < A.nnz; ++i) {
        vect_t vy, vx;
        int k = 0;
        for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
            vy = simd::load(y + row[i] * ldy + k);
            vx = simd::load(x + col[i] * ldx + k);
            simd::store(y + row[i] * ldy + k, simd::sub(vy, vx));
        }
        for (; k < blockSize; ++k)
            y[row[i] * ldy + k] -= x[col[i] * ldx + k];
    }
}

template <class Field>
inline void fspmm_mone_simd_unaligned(const Field &F, const Sparse<Field, SparseMatrix_t::COO> &A, int blockSize,
                                      typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_,
                                      int ldy) {
    assume_aligned(row, A.row, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_size;
    for (index_t i = 0; i < A.nnz; ++i) {
        vect_t vy, vx;
        int k = 0;
        for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
            vy = simd::loadu(y + row[i] * ldy + k);
            vx = simd::loadu(x + col[i] * ldx + k);
            simd::storeu(y + row[i] * ldy + k, simd::sub(vy, vx));
        }
        for (; k < blockSize; ++k)
            y[row[i] * ldy + k] -= x[col[i] * ldx + k];
    }
}

#endif

} // coo_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_coo_spmm_INL