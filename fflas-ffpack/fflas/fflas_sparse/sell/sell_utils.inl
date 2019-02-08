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

#ifndef __FFLASFFPACK_fflas_sparse_sell_utils_INL
#define __FFLASFFPACK_fflas_sparse_sell_utils_INL

namespace FFLAS {

    template <class Field>
    inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::SELL_ZO> &A, typename Field::ConstElement_ptr x,
                      typename Field::Element_ptr y, FieldCategories::ModularTag) {
        fspmv(F, A, x, y, FieldCategories::UnparametricTag());
        freduce(F, A.m, y, 1);
    }

    template <class Field> inline void sparse_delete(const Sparse<Field, SparseMatrix_t::SELL> &A) {
        fflas_delete(A.dat);
        fflas_delete(A.col);
        fflas_delete(A.st);
        fflas_delete(A.chunkSize);
    }

    template <class Field> inline void sparse_delete(const Sparse<Field, SparseMatrix_t::SELL_ZO> &A) {
        fflas_delete(A.col);
        fflas_delete(A.st);
        fflas_delete(A.chunkSize);
    }

    namespace sell_details {

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

    template <class Field> inline void sparse_print(const Sparse<Field, SparseMatrix_t::SELL> &A) {
        uint64_t it = 0;
        for (size_t i = 0; i < A.nChunks; ++i) {
            for (size_t k = 0; k < A.chunk; ++k) {
                std::cout << i *A.chunk + k << " : ";
                for (size_t j = 0; j < A.chunkSize[i]; ++j) {
                    std::cout << A.dat[it + j * A.chunk + k] << " ";
                }
                std::cout << std::endl;
            }
            it += A.chunkSize[i] * A.chunk;
        }
    }

    template <class Field, class IndexT>
    inline void sparse_init(const Field &F, Sparse<Field, SparseMatrix_t::SELL> &A, const IndexT *row, const IndexT *col,
                            typename Field::ConstElement_ptr dat, uint64_t rowdim, uint64_t coldim, uint64_t nnz,
                            uint64_t sigma = 0) {
        using namespace sell_details;
        using coo = Coo<typename Field::Element, IndexT>;
        if (!sigma)
            sigma = rowdim;
        A.kmax = Protected::DotProdBoundClassic(F, F.one);
        A.m = rowdim;
        A.n = coldim;
        A.nnz = nnz;
#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
        using simd = Simd<typename Field::Element>;
        A.chunk = simd::vect_size;
#else
        A.chunk = 8;
#endif
        index_t m = (A.m % A.chunk == 0) ? A.m : ROUND_DOWN(A.m, A.chunk) + A.chunk;
        A.nChunks = (m / A.chunk);

        std::vector<coo> data;
        std::vector<Info> infos(A.nChunks * A.chunk);

        for (uint64_t i = 0; i < nnz; ++i) {
            data.emplace_back(dat[i], row[i], col[i]);
        }

        IndexT currow = row[0];

        infos[currow].begin = 0;
        for (uint64_t i = 0; i < nnz; ++i) {
            if (row[i] != currow) {
                currow = row[i];
                infos[currow].begin = i;
            }
            infos[row[i]].size++;
        }

        A.maxrow = (std::max_element(infos.begin(), infos.end(),
                                     [](const Info &a, const Info &b) { return a.size >= b.size; }))->size;

        // cout << "maxrow : " << A.maxrow << endl;

        if (A.maxrow < A.kmax)
            A.delayed = true;

        for (uint64_t i = 0; i < rowdim; ++i) {
            infos[i].perm = i;
        }

#ifdef SELL_DEBUG
        for (auto &x : infos) {
            cout << x.size << " ";
        }
        std::cout << std::endl;
#endif

        uint64_t it = 0;
        for (; it < ROUND_DOWN(rowdim, sigma); it += sigma) {
            std::sort(infos.begin() + it, infos.begin() + it + sigma,
                      [](const Info &a, const Info &b) { return a.size >= b.size; });
        }
        if (it != rowdim) {
            std::sort(infos.begin() + it, infos.end(), [](Info a, Info b) { return a.size >= b.size; });
        }

        // cout << "sorted : " << std::is_sorted(infos.begin(), infos.end(), [](Info
        // a, Info b){
        // 	return a.size >= b.size;
        // }) << endl;

        for (size_t i = 0; i < infos.size(); ++i) {
            if (infos[i].begin > nnz)
                std::cout << "ERROR sort " << i << " size : " << infos[i].size << " begin : " << infos[i].begin
                << " perm : " << infos[i].perm << std::endl;
        }

#ifdef SELL_DEBUG
        for (auto &x : infos) {
            cout << x.size << " ";
        }
        std::cout << std::endl;
#endif

        A.perm = fflas_new<index_t>(rowdim, Alignment::CACHE_LINE);

        // cout << "perm : ";
        for (uint64_t i = 0; i < rowdim; ++i) {
            // cout << "(" << i << " , " << infos[i].perm << ") ";
            A.perm[infos[i].perm] = i;
        }

        // for(size_t i = 0 ; i < A.m ; ++i)
        //           cout << A.perm[i] << " ";
        //       cout << endl;
        // cout << endl;

        // add info if rowdim%chunk != 0, with empty infos (size = 0, begin = 0)
        // infos.resize(A.nChunks*A.chunk);

        // for(auto & x:infos)
        // 	if(x.begin > nnz)
        // 		cout << "ERROR resize" << endl;

        A.chunkSize = fflas_new<index_t>(A.nChunks, Alignment::CACHE_LINE);

        for (uint64_t i = 0; i < A.nChunks; ++i)
            A.chunkSize[i] = 0;

        for (uint64_t i = 0; i < A.nChunks; ++i) {
            for (uint64_t j = 0; j < A.chunk; ++j) {
                if (infos[i * A.chunk + j].size >= A.chunkSize[i])
                    A.chunkSize[i] = infos[i * A.chunk + j].size;
            }
        }

#ifdef SELL_DEBUG
        for (uint64_t i = 0; i < A.nChunks; ++i)
            cout << "chunk " << i << " : " << A.chunkSize[i] << endl;
        ;
#endif
        uint64_t sum = 0;
        for (uint64_t i = 0; i < A.nChunks; ++i)
            sum += A.chunkSize[i];
#ifdef SELL_DEBUG
        cout << "sum :  " << sum << " chunk : " << A.chunk << endl;
#endif
        A.col = fflas_new<index_t>(sum * A.chunk, Alignment::CACHE_LINE);
        A.dat = fflas_new(F, sum * A.chunk, Alignment::CACHE_LINE);
        A.nElements = sum * A.chunk;

        for (uint64_t i = 0; i < sum * A.chunk; ++i) {
            A.col[i] = 0;
            F.assign(A.dat[i], F.zero);
        }

        it = 0;
        for (uint64_t i = 0; i < A.nChunks; ++i) {
            for (uint64_t k = 0; k < A.chunk; ++k) {
                uint64_t start = infos[i * A.chunk + k].begin;
#ifdef SELL_DEBUG
                cout << it << " " << start << " " << infos[i * A.chunk + k].size << endl;
                cout << "	";
#endif
                for (uint64_t j = 0; j < infos[i * A.chunk + k].size; ++j) {
                    if (it + k + j * A.chunk >= sum * A.chunk)
                        std::cout << "error : " << it + k + j *A.chunk << " " << sum *A.chunk << std::endl;
                    A.dat[it + k + j * A.chunk] = data[start + j].val;
                    A.col[it + k + j * A.chunk] = data[start + j].col;
#ifdef SELL_DEBUG
                    cout << data[start + j].val << " ";
#endif
                }
#ifdef SELL_DEBUG
                cout << endl;
#endif
            }
            it += A.chunkSize[i] * A.chunk;
        }
        A.st = fflas_new<uint64_t>(A.nChunks, Alignment::CACHE_LINE);
        A.st[0] = 0;
        for (uint64_t i = 1; i < A.nChunks; ++i) {
            A.st[i] = A.chunkSize[i - 1] * A.chunk;
        }
        for (uint64_t i = 1; i < A.nChunks; ++i) {
            A.st[i] += A.st[i - 1];
        }
#ifdef SELL_DEBUG
        cout << "st : ";
        for (uint64_t i = 0; i < A.nChunks; ++i)
            cout << A.st[i] << " ";
        cout << endl;
#endif
    }

    template <class Field, class IndexT>
    inline void sparse_init(const Field &F, Sparse<Field, SparseMatrix_t::SELL_ZO> &A, const IndexT *row, const IndexT *col,
                            typename Field::ConstElement_ptr dat, uint64_t rowdim, uint64_t coldim, uint64_t nnz) {}
} // FFLAS
#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
