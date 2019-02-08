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

/** @file fflas/fflas_fspmv_sell.inl
 * NO DOC
 */

#ifndef __FFLASFFPACK_fflas_sparse_sell_H
#define __FFLASFFPACK_fflas_sparse_sell_H

namespace FFLAS { /*  SELL */

    template <class _Field> struct Sparse<_Field, SparseMatrix_t::SELL> {
        using Field = _Field;
        bool delayed = false;
        int chunk = 0;
        index_t kmax = 0;
        index_t m = 0;
        index_t n = 0;
        index_t maxrow = 0;
        index_t sigma = 0;
        index_t nChunks = 0;
        uint64_t nnz = 0;
        uint64_t nElements = 0;
        index_t *perm = nullptr;
        uint64_t *st = nullptr;
        index_t *chunkSize = nullptr;
        index_t *col = nullptr;
        typename _Field::Element_ptr dat;
    };

    template <class _Field>
    struct Sparse<_Field, SparseMatrix_t::SELL_ZO>
    : public Sparse<_Field, SparseMatrix_t::SELL> {
        using Field = _Field;
        typename _Field::Element cst = 1;
    };

} // FFLAS

#include "fflas-ffpack/fflas/fflas_sparse/sell/sell_utils.inl"
#include "fflas-ffpack/fflas/fflas_sparse/sell/sell_spmv.inl"
#if defined(__FFLASFFPACK_USE_OPENMP)
#include "fflas-ffpack/fflas/fflas_sparse/sell/sell_pspmv.inl"
#endif
// #include "fflas-ffpack/fflas/fflas_sparse/sell/sell_spmm.inl"
// #include "fflas-ffpack/fflas/fflas_sparse/sell/sell_pspmm.inl"

#endif // __FFLASFFPACK_fflas_sparse_SELL_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
