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

#ifndef __FFLASFFPACK_fflas_sparse_CSR_spmm_INL
#define __FFLASFFPACK_fflas_sparse_CSR_spmm_INL

namespace FFLAS {
namespace sparse_details_impl {

template <class Field>
inline void fspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A, int blockSize,
                  typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                  FieldCategories::GenericTag) {
    assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
    assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    for (index_t i = 0; i < A.m; ++i) {
        auto start = st[i], stop = st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                F.axpyin(y[i * ldy + k], dat[j], x[col[j] * ldx + k]);
                F.axpyin(y[i * ldy + k + 1], dat[j], x[col[j] * ldx + k + 1]);
                F.axpyin(y[i * ldy + k + 2], dat[j], x[col[j] * ldx + k + 2]);
                F.axpyin(y[i * ldy + k + 3], dat[j], x[col[j] * ldx + k + 3]);
            }
            for (; k < blockSize; ++k)
                F.axpyin(y[i * ldy + k], dat[j], x[col[j] * ldx + k]);
        }
    }
}

template <class Field>
inline void fspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A, int blockSize,
                  typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                  FieldCategories::UnparametricTag) {
    assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
    assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    for (index_t i = 0; i < A.m; ++i) {
        auto start = st[i], stop = st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                y[i * ldy + k] += dat[j] * x[col[j] * ldx + k];
                y[i * ldy + k + 1] += dat[j] * x[col[j] * ldx + k + 1];
                y[i * ldy + k + 2] += dat[j] * x[col[j] * ldx + k + 2];
                y[i * ldy + k + 3] += dat[j] * x[col[j] * ldx + k + 3];
            }
            for (; k < blockSize; ++k)
                y[i * ldy + k] += dat[j] * x[col[j] * ldx + k];
        }
    }
}

#ifdef __FFLASFFPACK_USE_SIMD

template <class Field>
inline void fspmm_simd_aligned(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A, int blockSize,
                               typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                               FieldCategories::UnparametricTag) {
    // std::cout << "spmm simd Unparam aligned" << std::endl;
    assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
    assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_t;
    for (index_t i = 0; i < A.m; ++i) {
        auto start = st[i], stop = st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            vect_t y1, x1, y2, x2, vdat;
            int k = 0;
            vdat = simd::set1(dat[j]);
            for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                y1 = simd::load(y+i*ldy+k);
                y2 = simd::load(y+i*ldy+k+simd::vect_size);
                x1 = simd::load(x + col[j] * ldx + k);
                x2 = simd::load(x + col[j] * ldx + k + simd::vect_size);
                y1 = simd::fmadd(y1, x1, vdat);
                y2 = simd::fmadd(y2, x2, vdat);
                simd::store(y + i * ldy + k, y1);
                simd::store(y + i * ldy + k + simd::vect_size, y2);
            }
            for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                y1 = simd::load(y+i*ldy+k);
                x1 = simd::load(x + col[j] * ldx + k);
                y1 = simd::fmadd(y1, x1, vdat);
                simd::store(y + i * ldy + k, y1);
            }
            for (; k < blockSize; ++k) {
                y[i * ldy + k] += dat[j] * x[col[j] * ldx + k];
            }
        }
    }
}

template <class Field>
inline void fspmm_simd_unaligned(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A, int blockSize,
                                typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                                FieldCategories::UnparametricTag) {
    // std::cout << "spmm simd Unparam unaligned" << std::endl;
    assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
    assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_t;
    for (index_t i = 0; i < A.m; ++i) {
        auto start = st[i], stop = st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            vect_t y1, x1, y2, x2, vdat;
            int k = 0;
            vdat = simd::set1(dat[j]);
            for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                y1 = simd::loadu(y+i*ldy+k);
                y2 = simd::loadu(y+i*ldy+k+simd::vect_size);
                x1 = simd::loadu(x + A.col[j] * ldx + k);
                x2 = simd::loadu(x + A.col[j] * ldx + k + simd::vect_size);
                y1 = simd::fmadd(y1, x1, vdat);
                y2 = simd::fmadd(y2, x2, vdat);
                simd::storeu(y + i * ldy + k, y1);
                simd::storeu(y + i * ldy + k + simd::vect_size, y2);
            }
            for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                y1 = simd::loadu(y+i*ldy+k);
                x1 = simd::loadu(x + col[j] * ldx + k);
                y1 = simd::fmadd(y1, x1, vdat);
                simd::storeu(y + i * ldy + k, y1);
            }
            for (; k < blockSize; ++k) {
                y[i * ldy + k] += dat[j] * x[col[j] * ldx + k];
            }
        }
    }
}
#endif

template <class Field>
inline void fspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A, int blockSize,
                  typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                  const int64_t kmax) {
    assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
    assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
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
                for (int k = 0; k < blockSize; ++k) {
                    y[i * ldy + k] += dat[j] * x[col[j] * ldx + k];
                }
            }
            // TODO : replace with freduce
	    FFLAS::freduce(F,blockSize,y+i*ldy,1);
	    // for (int k = 0; k < blockSize; ++k) {
		// F.reduce(y[i * ldy + k]);
	    // }
        }
        for (; j < j_end; ++j) {
            for (int k = 0; k < blockSize; ++k) {
                y[i * ldy + k] += dat[j] * x[col[j] * ldx + k];
            }
        }
	FFLAS::freduce(F,blockSize,y+i*ldy,1);
	// for (int k = 0; k < blockSize; ++k) {
	    // F.reduce(y[i * ldy + k]);
	// }
    }
}

#ifdef __FFLASFFPACK_USE_SIMD

template <class Field>
inline void fspmm_simd_unaligned(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A, int blockSize,
                                 typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                                 const int64_t kmax) {
    assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
    assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_t;
    for (index_t i = 0; i < A.m; ++i) {
        index_t j = st[i];
        index_t j_loc = j;
        index_t j_end = st[i + 1];
        index_t block = (j_end - j_loc) / kmax;
        for (index_t l = 0; l < (index_t)block; ++l) {
            j_loc += kmax;
            for (; j < j_loc; ++j) {
                vect_t y1, x1, y2, x2, vdat;
                int k = 0;
                vdat = simd::set1(dat[j]);
                for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                    y1 = simd::loadu(y+i*ldy+k);
                    y2 = simd::loadu(y+i*ldy+k+simd::vect_size);
                    x1 = simd::loadu(x + col[j] * ldx + k);
                    x2 = simd::loadu(x + col[j] * ldx + k + simd::vect_size);
                    y1 = simd::fmadd(y1, x1, vdat);
                    y2 = simd::fmadd(y2, x2, vdat);
                    simd::storeu(y + i * ldy + k, y1);
                    simd::storeu(y + i * ldy + k + simd::vect_size, y2);
                }
                for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                    y1 = simd::loadu(y+i*ldy+k);
                    x1 = simd::loadu(x + col[j] * ldx + k);
                    y1 = simd::fmadd(y1, x1, vdat);
                    simd::storeu(y + i * ldy + k, y1);
                }
                for (; k < blockSize; ++k) {
                    y[i * ldy + k] += dat[j] * x[col[j] * ldx + k];
                }
            }
            // TODO : replace with freduce
	    //
	    FFLAS::freduce(F,blockSize,y+i*ldy,1);
	    // for (int k = 0; k < blockSize; ++k) {
		// F.reduce(y[i * ldy + k]);
	    // }
        }
        for (; j < j_end; ++j) {
            vect_t y1, x1, y2, x2, vdat;
            y1 = simd::zero();
            y2 = simd::zero();
            int k = 0;
            vdat = simd::set1(dat[j]);
            for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                y1 = simd::loadu(y+i*ldy+k);
                y2 = simd::loadu(y+i*ldy+k+simd::vect_size);
                x1 = simd::loadu(x + col[j] * ldx + k);
                x2 = simd::loadu(x + col[j] * ldx + k + simd::vect_size);
                y1 = simd::fmadd(y1, x1, vdat);
                y2 = simd::fmadd(y2, x2, vdat);
                simd::storeu(y + i * ldy + k, y1);
                simd::storeu(y + i * ldy + k + simd::vect_size, y2);
            }
            for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                y1 = simd::loadu(y+i*ldy+k);
                x1 = simd::loadu(x + col[j] * ldx + k);
                y1 = simd::fmadd(y1, x1, vdat);
                simd::storeu(y + i * ldy + k, y1);
            }
            for (; k < blockSize; ++k) {
                y[i * ldy + k] += dat[j] * x[col[j] * ldx + k];
            }
        }
	    FFLAS::freduce(F,blockSize,y+i*ldy,1);
	// for (int k = 0; k < blockSize; ++k) {
	    // F.reduce(y[i * ldy + k]);
	// }
    }
}

template <class Field>
inline void fspmm_simd_aligned(const Field &F, const Sparse<Field, SparseMatrix_t::CSR> &A, int blockSize,
                               typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                               const int64_t kmax) {
    assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
    assume_aligned(dat, A.dat, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_t;
    for (index_t i = 0; i < A.m; ++i) {
        index_t j = st[i];
        index_t j_loc = j;
        index_t j_end = st[i + 1];
        index_t block = (j_end - j_loc) / kmax;
        for (index_t l = 0; l < (index_t)block; ++l) {
            j_loc += kmax;
            for (; j < j_loc; ++j) {
                vect_t y1, x1, y2, x2, vdat;
                int k = 0;
                vdat = simd::set1(dat[j]);
                for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                    y1 = simd::load(y+i*ldy+k);
                    y2 = simd::load(y+i*ldy+k+simd::vect_size);
                    x1 = simd::load(x + col[j] * ldx + k);
                    x2 = simd::load(x + col[j] * ldx + k + simd::vect_size);
                    y1 = simd::fmadd(y1, x1, vdat);
                    y2 = simd::fmadd(y2, x2, vdat);
                    simd::store(y + i * ldy + k, y1);
                    simd::store(y + i * ldy + k + simd::vect_size, y2);
                }
                for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                    y1 = simd::load(y+i*ldy+k);
                    x1 = simd::load(x + col[j] * ldx + k);
                    y1 = simd::fmadd(y1, x1, vdat);
                    simd::store(y + i * ldy + k, y1);
                }
                for (; k < blockSize; ++k) {
                    y[i * ldy + k] += dat[j] * x[col[j] * ldx + k];
                }
            }
            // TODO : replace with freduce
	    FFLAS::freduce(F,blockSize,y+i*ldy,1);
	    // for (int k = 0; k < blockSize; ++k) {
		// F.reduce(y[i * ldy + k]);
	    // }
        }
        for (; j < j_end; ++j) {
            vect_t y1, x1, y2, x2, vdat;
            int k = 0;
            vdat = simd::set1(dat[j]);
            for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                y1 = simd::load(y+i*ldy+k);
                y2 = simd::load(y+i*ldy+k+simd::vect_size);
                x1 = simd::load(x + col[j] * ldx + k);
                x2 = simd::load(x + col[j] * ldx + k + simd::vect_size);
                y1 = simd::fmadd(y1, x1, vdat);
                y2 = simd::fmadd(y2, x2, vdat);
                simd::store(y + i * ldy + k, y1);
                simd::store(y + i * ldy + k + simd::vect_size, y2);
            }
            for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                y1 = simd::load(y+i*ldy+k);
                x1 = simd::load(x + col[j] * ldx + k);
                y1 = simd::fmadd(y1, x1, vdat);
                simd::store(y + i * ldy + k, y1);
            }
            for (; k < blockSize; ++k) {
                y[i * ldy + k] += dat[j] * x[col[j] * ldx + k];
            }
        }
	FFLAS::freduce(F,blockSize,y+i*ldy,1);
	// for (int k = 0; k < blockSize; ++k) {
	    // F.reduce(y[i * ldy + k]);
	// }
    }
}

#endif // SIMD

template <class Field>
inline void fspmm_one(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A, int blockSize,
                      typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                      FieldCategories::GenericTag) {
    assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    for (index_t i = 0; i < A.m; ++i) {
        auto start = st[i], stop = st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                F.addin(y[i * ldy + k], x[col[j] * ldx + k]);
                F.addin(y[i * ldy + k + 1], x[col[j] * ldx + k + 1]);
                F.addin(y[i * ldy + k + 2], x[col[j] * ldx + k + 2]);
                F.addin(y[i * ldy + k + 3], x[col[j] * ldx + k + 3]);
            }
            for (; k < blockSize; ++k)
                F.addin(y[i * ldy + k], x[col[j] * ldx + k]);
        }
    }
}

template <class Field>
inline void fspmm_mone(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A, int blockSize,
                       typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_, int ldy,
                       FieldCategories::GenericTag) {
    assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    for (index_t i = 0; i < A.m; ++i) {
        auto start = st[i], stop = st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                F.subin(y[i * ldy + k], x[col[j] * ldx + k]);
                F.subin(y[i * ldy + k + 1], x[col[j] * ldx + k + 1]);
                F.subin(y[i * ldy + k + 2], x[col[j] * ldx + k + 2]);
                F.subin(y[i * ldy + k + 3], x[col[j] * ldx + k + 3]);
            }
            for (; k < blockSize; ++k)
                F.subin(y[i * ldy + k], x[col[j] * ldx + k]);
        }
    }
}

#ifdef __FFLASFFPACK_USE_SIMD

template <class Field>
inline void fspmm_one_simd_aligned(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A, int blockSize,
                                   typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_,
                                   int ldy, FieldCategories::UnparametricTag) {
    assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_t;
    for (index_t i = 0; i < A.m; ++i) {
        auto start = st[i], stop = st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            vect_t y1, x1, y2, x2;
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                y1 = simd::load(y+i*ldy+k);
                y2 = simd::load(y+i*ldy+k+simd::vect_size);
                x1 = simd::load(x + col[j] * ldx + k);
                x2 = simd::load(x + col[j] * ldx + k + simd::vect_size);
                simd::store(y + i * ldy + k, simd::add(y1, x1));
                simd::store(y + i * ldy + k + simd::vect_size, simd::add(y2, x2));
            }
            for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                y1 = simd::load(y+i*ldy+k);
                x1 = simd::load(x + col[j] * ldx + k);
                simd::store(y + i * ldy + k, simd::add(y1, x1));
            }
            for (; k < blockSize; ++k) {
                y[i * ldy + k] += x[col[j] * ldx + k];
            }
        }
    }
}

template <class Field>
inline void fspmm_one_simd_unaligned(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A, int blockSize,
                                     typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_,
                                     int ldy, FieldCategories::UnparametricTag) {
    assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_t;
    for (index_t i = 0; i < A.m; ++i) {
        auto start = st[i], stop = st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            vect_t y1, x1, y2, x2;
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                y1 = simd::loadu(y+i*ldy+k);
                y2 = simd::loadu(y+i*ldy+k+simd::vect_size);
                x1 = simd::loadu(x + col[j] * ldx + k);
                x2 = simd::loadu(x + col[j] * ldx + k + simd::vect_size);
                simd::storeu(y + i * ldy + k, simd::add(y1, x1));
                simd::storeu(y + i * ldy + k + simd::vect_size, simd::add(y2, x2));
            }
            for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                y1 = simd::loadu(y+i*ldy+k);
                x1 = simd::loadu(x + col[j] * ldx + k);
                simd::storeu(y + i * ldy + k, simd::add(y1, x1));
            }
            for (; k < blockSize; ++k) {
                y[i * ldy + k] += x[col[j] * ldx + k];
            }
        }
    }
}

template <class Field>
inline void fspmm_mone_simd_aligned(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A, int blockSize,
                                    typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_,
                                    int ldy, FieldCategories::UnparametricTag) {
    assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_t;
    for (index_t i = 0; i < A.m; ++i) {
        auto start = st[i], stop = st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            vect_t y1, x1, y2, x2;
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                y1 = simd::load(y+i*ldy+k);
                y2 = simd::load(y+i*ldy+k+simd::vect_size);
                x1 = simd::load(x + col[j] * ldx + k);
                x2 = simd::load(x + col[j] * ldx + k + simd::vect_size);
                simd::store(y + i * ldy + k, simd::sub(y1, x1));
                simd::store(y + i * ldy + k + simd::vect_size, simd::sub(y2, x2));
            }
            for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                y1 = simd::load(y+i*ldy+k);
                x1 = simd::load(x + col[j] * ldx + k);
                simd::store(y + i * ldy + k, simd::sub(y1, x1));
            }
            for (; k < blockSize; ++k) {
                y[i * ldy + k] -= x[col[j] * ldx + k];
            }
        }
    }
}

template <class Field>
inline void fspmm_mone_simd_unaligned(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_ZO> &A, int blockSize,
                                      typename Field::ConstElement_ptr x_, int ldx, typename Field::Element_ptr y_,
                                      int ldy, FieldCategories::UnparametricTag) {
    assume_aligned(st, A.st, (size_t)Alignment::CACHE_LINE);
    assume_aligned(col, A.col, (size_t)Alignment::CACHE_LINE);
    assume_aligned(x, x_, (size_t)Alignment::DEFAULT);
    assume_aligned(y, y_, (size_t)Alignment::DEFAULT);
    using simd = Simd<typename Field::Element>;
    using vect_t = typename simd::vect_t;
    for (index_t i = 0; i < A.m; ++i) {
        auto start = st[i], stop = st[i + 1];
        for (index_t j = start; j < stop; ++j) {
            vect_t y1, x1, y2, x2;
            int k = 0;
            for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                y1 = simd::load(y+i*ldy+k);
                y2 = simd::load(y+i*ldy+k+simd::vect_size);
                x1 = simd::loadu(x + col[j] * ldx + k);
                x2 = simd::loadu(x + col[j] * ldx + k + simd::vect_size);
                simd::storeu(y + i * ldy + k, simd::sub(y1, x1));
                simd::storeu(y + i * ldy + k + simd::vect_size, simd::sub(y2, x2));
            }
            for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                y1 = simd::load(y+i*ldy+k);
                x1 = simd::loadu(x + col[j] * ldx + k);
                simd::storeu(y + i * ldy + k, simd::sub(y1, x1));
            }
            for (; k < blockSize; ++k) {
                y[i * ldy + k] -= x[col[j] * ldx + k];
            }
        }
    }
}

#endif //__FFLASFFPACK_USE_SIMD

} // CSR_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_CSR_spmm_INL
