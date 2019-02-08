/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by  Bastien Vialla <bastien.vialla@lirmm.fr>
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

namespace FFLAS {

    template <class Field> inline void sparse_delete(const Sparse<Field, SparseMatrix_t::CSR> &A) {
        fflas_delete(A.dat);
        fflas_delete(A.col);
        fflas_delete(A.st);
    }

    template <class Field> inline void sparse_delete(const Sparse<Field, SparseMatrix_t::CSR_ZO> &A) {
        fflas_delete(A.col);
        fflas_delete(A.st);
    }

    template <class Field> inline std::ostream& sparse_print(std::ostream& os, const Sparse<Field, SparseMatrix_t::CSR> &A) {
        // for (size_t i = 0; i <= A.m; ++i)
        //     std::cout << A.st[i] << " ";
        // std::cout << std::endl;
        for (index_t i = 0; i < A.m; ++i) {
            auto start = A.st[i], stop = A.st[i + 1];
            index_t j = 0;
            index_t diff = stop - start;
            os << i << " : ";
            for (; j < diff; ++j) {
                os << '(' << A.col[start + j] << ',' << A.dat[start+j] << ") ";
            }
            os << std::endl;
        }
        return os;
    }

    template <class IndexT>
    inline void sparse_init(const Givaro::Modular<Givaro::Integer> &F, Sparse<Givaro::Modular<Givaro::Integer>, SparseMatrix_t::CSR> &A, const IndexT *row, const IndexT *col,
                            Givaro::Integer* dat, uint64_t rowdim, uint64_t coldim, uint64_t nnz) {
        A.m = rowdim;
        A.n = coldim;
        A.nnz = nnz;
        A.nElements = nnz;
        std::vector<uint64_t> rows(rowdim, 0);
        for (uint64_t i = 0; i < A.nnz; ++i)
            rows[row[i]]++;

        A.delayed = true;

        A.col = fflas_new<index_t>(nnz, Alignment::CACHE_LINE);
        A.st = fflas_new<index_t>(rowdim + 1, Alignment::CACHE_LINE);
        A.dat = fflas_new(F, nnz, Alignment::CACHE_LINE);

        for(size_t i = 0 ; i < nnz ; ++i){
            if(col[i] >= coldim){
                std::cout << "Error col index too big" << std::endl;
            }
        }

        for (size_t i = 0; i < nnz; ++i) {
            A.col[i] = static_cast<index_t>(col[i]);
            A.dat[i] = dat[i];
        }

        A.st[0] = 0;
        for (size_t i = 1; i <= rowdim; ++i) {
            A.st[i] = A.st[i - 1] + rows[i - 1];
        }
    }

    template <class IndexT>
    inline void sparse_init(const Givaro::ZRing<Givaro::Integer> &F, Sparse<Givaro::ZRing<Givaro::Integer>, SparseMatrix_t::CSR_ZO> &A, const IndexT *row, const IndexT *col,
                            Givaro::Integer* dat, uint64_t rowdim, uint64_t coldim, uint64_t nnz) {
        A.m = rowdim;
        A.n = coldim;
        A.nnz = nnz;
        A.nElements = nnz;
        std::vector<uint64_t> rows(rowdim, 0);
        for (uint64_t i = 0; i < A.nnz; ++i)
            rows[row[i]]++;

        A.delayed = true;

        A.col = fflas_new<index_t>(nnz, Alignment::CACHE_LINE);
        A.st = fflas_new<index_t>(rowdim + 1, Alignment::CACHE_LINE);

        for(size_t i = 0 ; i < nnz ; ++i){
            if(col[i] >= coldim){
                std::cout << "Error col index too big" << std::endl;
            }
        }

        for (size_t i = 0; i < nnz; ++i) {
            A.col[i] = static_cast<index_t>(col[i]);
        }

        A.st[0] = 0;
        for (size_t i = 1; i <= rowdim; ++i) {
            A.st[i] = A.st[i - 1] + rows[i - 1];
        }
    }

    template <class IndexT, size_t RECINT_SIZE>
    inline void sparse_init(const Givaro::ZRing<RecInt::rmint<RECINT_SIZE>> &F, Sparse<Givaro::ZRing<RecInt::rmint<RECINT_SIZE>>, SparseMatrix_t::CSR_ZO> &A, const IndexT *row, const IndexT *col,
                            typename Givaro::ZRing<RecInt::rmint<RECINT_SIZE>>::Element_ptr dat, uint64_t rowdim, uint64_t coldim, uint64_t nnz) {
        A.m = rowdim;
        A.n = coldim;
        A.nnz = nnz;
        A.nElements = nnz;
        std::vector<uint64_t> rows(rowdim, 0);
        for (uint64_t i = 0; i < A.nnz; ++i)
            rows[row[i]]++;

        A.delayed = true;

        A.col = fflas_new<index_t>(nnz, Alignment::CACHE_LINE);
        A.st = fflas_new<index_t>(rowdim + 1, Alignment::CACHE_LINE);

        for(size_t i = 0 ; i < nnz ; ++i){
            if(col[i] >= coldim){
                std::cout << "Error col index too big" << std::endl;
            }
        }

        for (size_t i = 0; i < nnz; ++i) {
            A.col[i] = static_cast<index_t>(col[i]);
        }

        A.st[0] = 0;
        for (size_t i = 1; i <= rowdim; ++i) {
            A.st[i] = A.st[i - 1] + rows[i - 1];
        }
    }

    template <class IndexT, size_t RECINT_SIZE>
    inline void sparse_init(const Givaro::ZRing<RecInt::rmint<RECINT_SIZE>> &F, Sparse<Givaro::ZRing<RecInt::rmint<RECINT_SIZE>>, SparseMatrix_t::CSR> &A, const IndexT *row, const IndexT *col,
                            typename Givaro::ZRing<RecInt::rmint<RECINT_SIZE>>::Element_ptr dat, uint64_t rowdim, uint64_t coldim, uint64_t nnz) {
        A.m = rowdim;
        A.n = coldim;
        A.nnz = nnz;
        A.nElements = nnz;
        std::vector<uint64_t> rows(rowdim, 0);
        for (uint64_t i = 0; i < A.nnz; ++i)
            rows[row[i]]++;

        A.delayed = true;

        A.col = fflas_new<index_t>(nnz, Alignment::CACHE_LINE);
        A.st = fflas_new<index_t>(rowdim + 1, Alignment::CACHE_LINE);

        for(size_t i = 0 ; i < nnz ; ++i){
            if(col[i] >= coldim){
                std::cout << "Error col index too big" << std::endl;
            }
        }

        for (size_t i = 0; i < nnz; ++i) {
            A.col[i] = static_cast<index_t>(col[i]);
        }

        A.st[0] = 0;
        for (size_t i = 1; i <= rowdim; ++i) {
            A.st[i] = A.st[i - 1] + rows[i - 1];
        }
    }

    template <class Field, class IndexT>
    inline void sparse_init(const Field &F, Sparse<Field, SparseMatrix_t::CSR> &A, const IndexT *row, const IndexT *col,
                            typename Field::ConstElement_ptr dat, uint64_t rowdim, uint64_t coldim, uint64_t nnz) {
        A.kmax = Protected::DotProdBoundClassic(F, F.one);
        A.m = rowdim;
        A.n = coldim;
        A.nnz = nnz;
        A.nElements = nnz;
        std::vector<uint64_t> rows(rowdim, 0);
        for (uint64_t i = 0; i < A.nnz; ++i)
            rows[row[i]]++;

        A.maxrow = *(std::max_element(rows.begin(), rows.end()));

        if (A.kmax > A.maxrow)
            A.delayed = true;

        A.col = fflas_new<index_t>(nnz, Alignment::CACHE_LINE);
        A.st = fflas_new<index_t>(rowdim + 1, Alignment::CACHE_LINE);
        A.stend = fflas_new<index_t>(rowdim + 1, Alignment::CACHE_LINE);
        A.dat = fflas_new(F, nnz, Alignment::CACHE_LINE);

        for (size_t i = 0; i < nnz; ++i) {
            A.col[i] = static_cast<index_t>(col[i]);
            A.dat[i] = dat[i];
        }
        A.st[0] = 0;
        for (size_t i = 1; i <= rowdim; ++i) {
            A.st[i] = A.st[i - 1] + rows[i - 1];
        }
        for(size_t i = 0 ; i < rowdim ; ++i){
            A.stend[i] = A.st[i+1];
        }
    }

    template <class Field, class IndexT>
    inline void sparse_init(const Field &F, Sparse<Field, SparseMatrix_t::CSR_ZO> &A, const IndexT *row, const IndexT *col,
                            typename Field::ConstElement_ptr dat, uint64_t rowdim, uint64_t coldim, uint64_t nnz) {
        A.delayed = true;
        A.m = rowdim;
        A.n = coldim;
        A.nnz = nnz;
        A.nElements = nnz;
        std::vector<uint64_t> rows(A.m, 0);
        for (uint64_t i = 0; i < A.nnz; ++i)
            rows[row[i]]++;
        A.maxrow = *(std::max_element(rows.begin(), rows.end()));
        A.col = fflas_new<index_t>(nnz, Alignment::CACHE_LINE);
        A.st = fflas_new<index_t>(rowdim + 1, Alignment::CACHE_LINE);
        for (size_t i = 0; i < nnz; ++i) {
            A.col[i] = static_cast<index_t>(col[i]);
        }
        for (size_t i = 0; i <= rowdim; ++i) {
            A.st[i] = 0;
        }
        for (size_t i = 0; i < nnz; ++i) {
            A.st[row[i] + 1]++;
        }
        for (size_t i = 1; i <= rowdim; ++i) {
            A.st[i] += A.st[i - 1];
        }
    }
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
