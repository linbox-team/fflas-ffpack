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

#ifndef __FFLASFFPACK_fflas_sparse_CSR_HYB_utils_INL
#define __FFLASFFPACK_fflas_sparse_CSR_HYB_utils_INL

// #define CSR_HYB_DEBUG 1

namespace FFLAS {

    template <class Field> inline void sparse_delete(const Sparse<Field, SparseMatrix_t::CSR_HYB> &A) {
        fflas_delete(A.dat);
        fflas_delete(A.col);
        fflas_delete(A.st);
    }

    namespace csr_hyb_details {

        struct Info {
            uint64_t size = 0;
            uint64_t perm = 0;
            uint64_t begin = 0;

            Info(uint64_t it, uint64_t s, uint64_t p) : size(s), perm(p), begin(it) {}
            Info() = default;
            Info(const Info &) = default;
            Info(Info &&) = default;

            Info &operator=(const Info &) = default;
            Info &operator=(Info &&) = default;
        };

        template <class ValT, class IdxT> struct Coo {
            using Self = Coo<ValT, IdxT>;

            ValT val = 0;
            IdxT row = 0;
            IdxT col = 0;

            Coo(ValT v, IdxT r, IdxT c) : val(v), row(r), col(c) {}
            Coo() = default;
            Coo(const Self &) = default;
            Coo(Self &&) = default;

            Self &operator=(const Self &) = default;
            Self &operator=(Self &&) = default;
        };
    }

    template <class Field, class IndexT>
    inline void sparse_init(const Field &F, Sparse<Field, SparseMatrix_t::CSR_HYB> &A, const IndexT *row, const IndexT *col,
                            typename Field::ConstElement_ptr dat, uint64_t rowdim, uint64_t coldim, uint64_t nnz) {
        using namespace csr_hyb_details;
        using coo = Coo<typename Field::Element, index_t>;

        A.kmax = Protected::DotProdBoundClassic(F, F.one);
        A.m = rowdim;
        A.n = coldim;
        A.nnz = nnz;
        A.nElements = nnz;
        std::vector<coo> data(nnz);
        for (uint64_t i = 0; i < nnz; ++i) {
            // data.emplace_back(dat[i], col[i], row[i]);
            data[i].val = dat[i];
            data[i].col = col[i];
            data[i].row = row[i];
        }

        std::vector<uint64_t> rows(rowdim, 0);

        for (uint64_t i = 0; i < A.nnz; ++i)
            rows[row[i]]++;

        A.maxrow = *(std::max_element(rows.begin(), rows.end()));

        if (A.kmax > A.maxrow)
            A.delayed = true;

        rows.resize(3 * (rowdim + 1));
        for (auto &x : rows)
            x = 0;

        for (uint64_t i = 0; i < data.size(); ++i) {
            auto x = data[i];
            if (F.isOne(x.val)) {
                rows[3 * x.row + 1]++;
                A.nOnes++;
            } else if (F.isMOne(x.val)) {
                rows[3 * x.row]++;
                A.nMOnes++;
            } else {
                rows[3 * x.row + 2]++;
                A.nOthers++;
            }
        }

        A.col = fflas_new<index_t>(nnz, Alignment::CACHE_LINE);
        A.st = fflas_new<index_t>(4 * (rowdim + 1), Alignment::CACHE_LINE);
        A.dat = fflas_new(F, A.nOthers, Alignment::CACHE_LINE);

        for (uint64_t i = 0; i < 4 * (rowdim + 1); ++i)
            A.st[i] = 0;

        for (size_t i = 0; i < nnz; ++i) {
            A.col[i] = static_cast<index_t>(data[i].col);
        }

        data.shrink_to_fit();

        // sort nnz by row with order -1 1 L

        std::sort(data.begin(), data.end(), [&F](const coo &a, const coo &b) {
                  return (a.row < b.row) || ((a.row == b.row) && (F.isMOne(a.val) && !F.isMOne(b.val))) ||
                  ((a.row == b.row) && (F.isMOne(a.val) && F.isMOne(b.val) && (a.col < b.col))) ||
                  ((a.row == b.row) && (F.isOne(a.val) && !F.isOne(b.val) && !F.isMOne(b.val))) ||
                  ((a.row == b.row) && (F.isOne(a.val) && F.isOne(b.val) && (a.col < b.col))) ||
                  ((a.row == b.row) && (!F.isOne(a.val) && !F.isMOne(a.val) && !F.isOne(b.val) && !F.isMOne(b.val)) &&
                   (a.col < b.col));
                  });

#ifdef CSR_HYB_DEBUG
        for (auto &x : data) {
            cout << "(" << x.row << "," << x.col << "," << x.val << ") ";
        }
        cout << endl;
#endif

        uint64_t it = 0;
        for (size_t i = 0; i < data.size(); ++i) {
            if (!F.isOne(data[i].val) && !F.isMOne(data[i].val)) {
                A.dat[it] = data[i].val;
                ++it;
            }
        }

        A.st[1] = rows[0];
        A.st[2] = rows[1] + A.st[1];
        A.st[3] = 0;
        A.st[4] = rows[2] + A.st[2];

        for (uint64_t i = 1; i < rowdim; ++i) {
            A.st[4 * i + 1] = rows[3 * i] + A.st[4 * i];
            A.st[4 * i + 2] = rows[3 * i + 1] + A.st[4 * i + 1];
            A.st[4 * i + 3] = rows[3 * (i - 1) + 2] + A.st[4 * (i - 1) + 3];
            A.st[4 * (i + 1)] = rows[3 * i + 2] + A.st[4 * i + 2];
        }

#ifdef CSR_HYB_DEBUG
        for (uint64_t i = 0; i < it; ++i)
            cout << A.dat[i] << " ";
        cout << endl;
        for (uint64_t i = 0; i < nnz; ++i)
            cout << A.col[i] << " ";
        cout << endl;

        for (uint64_t i = 0; i < rowdim; ++i)
            cout << "(" << A.st[4 * i] << " , " << A.st[4 * i + 1] << " , " << A.st[4 * i + 2] << " , " << A.st[4 * i + 3]
            << ") " << endl;
        cout << endl;
        cout << endl;
        for (uint64_t i = 0; i < rowdim; ++i) {
            index_t start = A.st[4 * i], stop = A.st[4 * i + 1];
            index_t diff = stop - start;
            cout << i << endl;
            cout << "  -1 : ";
            for (uint64_t j = 0; j < diff; ++j) {
                cout << A.col[start + j] << " ";
            }
            cout << endl;
            start = A.st[4 * i + 1], stop = A.st[4 * i + 2];
            diff = stop - start;
            cout << "  1 : ";
            for (uint64_t j = 0; j < diff; ++j) {
                cout << A.col[start + j] << " ";
            }
            cout << endl;
            start = A.st[4 * i + 2], stop = A.st[4 * (i + 1)];
            diff = stop - start;
            index_t startDat = A.st[4 * i + 3];
            cout << "  l : ";
            for (uint64_t j = 0; j < diff; ++j) {
                cout << "(" << A.col[start + j] << " , " << A.dat[startDat + j] << ") ";
            }
            cout << endl;
        }
#endif
    }
} // FFLAS
#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
