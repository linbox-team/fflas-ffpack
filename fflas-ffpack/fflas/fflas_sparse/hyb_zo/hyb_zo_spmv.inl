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

#ifndef __FFLASFFPACK_fflas_sparse_HYB_ZO_spmv_INL
#define __FFLASFFPACK_fflas_sparse_HYB_ZO_spmv_INL

namespace FFLAS {
    namespace sparse_details_impl {

        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::HYB_ZO> &A, typename Field::ConstElement_ptr x,
                          typename Field::Element_ptr y, FieldCategories::GenericTag) {
            if (A.one != nullptr)
                sparse_details_impl::fspmv_one(F, *(A.one), x, y, FieldCategories::GenericTag());
            if (A.mone != nullptr)
                sparse_details_impl::fspmv_mone(F, *(A.mone), x, y, FieldCategories::GenericTag());
            if (A.dat != nullptr)
                sparse_details_impl::fspmv(F, *(A.dat), x, y, FieldCategories::GenericTag());
        }

        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::HYB_ZO> &A, typename Field::ConstElement_ptr x,
                          typename Field::Element_ptr y, FieldCategories::UnparametricTag) {
            if (A.one != nullptr)
                sparse_details_impl::fspmv_one(F, *(A.one), x, y, FieldCategories::UnparametricTag());
            if (A.mone != nullptr)
                sparse_details_impl::fspmv_mone(F, *(A.mone), x, y, FieldCategories::UnparametricTag());
            if (A.dat != nullptr)
                sparse_details_impl::fspmv(F, *(A.dat), x, y, FieldCategories::UnparametricTag());
        }

        template <class Field>
        inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::HYB_ZO> &A, typename Field::ConstElement_ptr x,
                          typename Field::Element_ptr y, uint64_t kmax) {
            if (A.one != nullptr)
                sparse_details_impl::fspmv_one(F, *(A.one), x, y, FieldCategories::UnparametricTag());
            if (A.mone != nullptr)
                sparse_details_impl::fspmv_mone(F, *(A.mone), x, y, FieldCategories::UnparametricTag());
            if (A.dat != nullptr)
                sparse_details_impl::fspmv(F, *(A.dat), x, y, kmax);
        }

    } // HYB_ZO_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_HYB_ZO_spmv_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
