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

/** @file fflas/fflas_fspmv_HYB_ZO.inl
 * NO DOC
 */

#ifndef __FFLASFFPACK_fflas_sparse_HYB_ZO_H
#define __FFLASFFPACK_fflas_sparse_HYB_ZO_H

namespace FFLAS { /*  HYB_ZO */

    template <class _Field> struct Sparse<_Field, SparseMatrix_t::HYB_ZO> {
        using Field = _Field;
        typedef Sparse<_Field, SparseMatrix_t::HYB_ZO> Self_t;
        bool delayed = false;
        uint64_t kmax = 0;
        index_t m = 0;
        index_t n = 0;
        uint64_t nnz = 0;
        uint64_t maxrow = 0;
        uint64_t nElements = 0;
        Sparse<_Field, SparseMatrix_t::CSR> *dat = nullptr;
        Sparse<_Field, SparseMatrix_t::CSR_ZO> *one = nullptr;
        Sparse<_Field, SparseMatrix_t::CSR_ZO> *mone = nullptr;

    };

} // FFLAS

#include "fflas-ffpack/fflas/fflas_sparse/hyb_zo/hyb_zo_utils.inl"
#include "fflas-ffpack/fflas/fflas_sparse/hyb_zo/hyb_zo_spmv.inl"
#include "fflas-ffpack/fflas/fflas_sparse/hyb_zo/hyb_zo_spmm.inl"
#if defined(__FFLASFFPACK_USE_OPENMP)
#include "fflas-ffpack/fflas/fflas_sparse/hyb_zo/hyb_zo_pspmv.inl"
#include "fflas-ffpack/fflas/fflas_sparse/hyb_zo/hyb_zo_pspmm.inl"
#endif


#endif // __FFLASFFPACK_fflas_sparse_HYB_ZO_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
