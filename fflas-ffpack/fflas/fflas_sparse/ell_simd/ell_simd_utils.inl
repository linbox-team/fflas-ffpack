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

#ifndef __FFLASFFPACK_fflas_sparse_ELL_simd_utils_INL
#define __FFLASFFPACK_fflas_sparse_ELL_simd_utils_INL

namespace FFLAS {

    template <class Field> inline void sparse_delete(const Sparse<Field, SparseMatrix_t::ELL_simd> &A) {
        fflas_delete(A.dat);
        fflas_delete(A.col);
    }

    template <class Field> inline void sparse_delete(const Sparse<Field, SparseMatrix_t::ELL_simd_ZO> &A) {
        fflas_delete(A.col);
    }

    template <class Field> inline void sparse_print(const Sparse<Field, SparseMatrix_t::ELL_simd> &A) {
        for (size_t i = 0; i < A.nChunks; ++i) {
            for (size_t k = 0; k < A.chunk; ++k) {
                std::cout << i *A.chunk + k << " : ";
                for (size_t j = 0; j < A.ld; ++j) {
                    std::cout << A.dat[i * A.ld * A.chunk + j * A.chunk + k] << " ";
                }
                std::cout << std::endl;
            }
        }
    }

    template <class Field, class IndexT>
    inline void sparse_init(const Field &F, Sparse<Field, SparseMatrix_t::ELL_simd> &A, const IndexT *row,
                            const IndexT *col, typename Field::ConstElement_ptr dat, uint64_t rowdim, uint64_t coldim,
                            uint64_t nnz) {
#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
        using simd = Simd<typename Field::Element>;
        A.chunk = simd::vect_size;
#else
        A.chunk = 8;
#endif
        A.kmax = Protected::DotProdBoundClassic(F, F.one);
        A.m = rowdim;
        A.n = coldim;
        A.nnz = nnz;
        std::vector<uint64_t> rows(A.m + 1, 0);
        for (uint64_t i = 0; i < A.nnz; ++i)
            rows[row[i] + 1]++;
        A.maxrow = *(std::max_element(rows.begin(), rows.end()));
        A.ld = A.maxrow;
        if (A.kmax > A.maxrow)
            A.delayed = true;
        for (size_t i = 1; i <= A.m; ++i) {
            rows[i] += rows[i - 1];
        }

        index_t m = (A.m % A.chunk == 0) ? A.m : ROUND_DOWN(A.m, A.chunk) + A.chunk;
        // cout << A.m << " " << ROUND_DOWN(A.m, simd::vect_size)+simd::vect_size <<
        // " " << m/A.chunk << endl;
        A.nChunks = m / A.chunk;

        A.col = fflas_new<index_t>(A.nChunks * A.chunk * A.ld, Alignment::CACHE_LINE);
        A.dat = fflas_new(F, A.nChunks * A.chunk * A.ld, Alignment::CACHE_LINE);

        A.nElements = A.nChunks * A.chunk * A.ld;

        for (size_t i = 0; i < A.nChunks * A.chunk * A.ld; ++i) {
            A.col[i] = 0;
            F.assign(A.dat[i], F.zero);
        }

        for (size_t i = 0; i < A.nChunks; ++i) {
            for (size_t k = 0; k < A.chunk; ++k) {
                if (i * A.chunk + k < rowdim) {
                    uint64_t start = rows[i * A.chunk + k], stop = rows[i * A.chunk + k + 1];
                    // cout << "start " << start << " stop " << stop << endl;
                    for (size_t j = 0; j < stop - start; ++j) {
                        A.dat[i * A.chunk * A.ld + j * A.chunk + k] = dat[start + j];
                        A.col[i * A.chunk * A.ld + j * A.chunk + k] = col[start + j];
                    }
                }
            }
        }
    }

    template <class Field, class IndexT>
    inline void sparse_init(const Field &F, Sparse<Field, SparseMatrix_t::ELL_simd_ZO> &A, const IndexT *row,
                            const IndexT *col, typename Field::ConstElement_ptr dat, uint64_t rowdim, uint64_t coldim,
                            uint64_t nnz) {
#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
        using simd = Simd<typename Field::Element>;
        A.chunk = simd::vect_size;
#else
        A.chunk = 8;
#endif
        A.kmax = Protected::DotProdBoundClassic(F, F.one);
        A.m = rowdim;
        A.n = coldim;
        A.nnz = nnz;
        std::vector<uint64_t> rows(A.m + 1, 0);
        for (uint64_t i = 0; i < A.nnz; ++i)
            rows[row[i] + 1]++;
        A.maxrow = *(std::max_element(rows.begin(), rows.end()));
        A.ld = A.maxrow;
        if (A.kmax > A.maxrow)
            A.delayed = true;
        for (size_t i = 1; i <= A.m; ++i) {
            rows[i] += rows[i - 1];
        }

        index_t m = (A.m % A.chunk == 0) ? A.m : ROUND_DOWN(A.m, A.chunk) + A.chunk;
        // cout << A.m << " " << ROUND_DOWN(A.m, simd::vect_size)+simd::vect_size <<
        // " " << m/A.chunk << endl;
        A.nChunks = m / A.chunk;

        A.col = fflas_new<index_t>(A.nChunks * A.chunk * A.ld, Alignment::CACHE_LINE);

        A.nElements = A.nChunks * A.chunk * A.ld;

        for (size_t i = 0; i < A.nChunks * A.chunk * A.ld; ++i) {
            A.col[i] = 0;
        }

        for (size_t i = 0; i < A.nChunks; ++i) {
            for (size_t k = 0; k < A.chunk; ++k) {
                if (i * A.chunk + k < rowdim) {
                    uint64_t start = rows[i * A.chunk + k], stop = rows[i * A.chunk + k + 1];
                    // cout << "start " << start << " stop " << stop << endl;
                    for (size_t j = 0; j < stop - start; ++j) {
                        A.col[i * A.chunk * A.ld + j * A.chunk + k] = col[start + j];
                    }
                }
            }
        }
    }
}
#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
