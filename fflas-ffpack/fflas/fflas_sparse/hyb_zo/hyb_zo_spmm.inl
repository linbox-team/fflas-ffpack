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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __FFLASFFPACK_fflas_sparse_HYB_ZO_spmm_INL
#define __FFLASFFPACK_fflas_sparse_HYB_ZO_spmm_INL

namespace FFLAS{
	namespace sparse_details_impl{
	
	template<class Field>
	inline void fspmm(const Field & F, const Sparse<Field, SparseMatrix_t::HYB_ZO> & A, int blockSize, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, FieldCategories::GenericTag){
		using Element = typename Field::Element;
		if(A.one != nullptr)
			sparse_details_impl::fspmm_zo(F, *(A.one), blockSize, x, y, [F](Element & a, const Element & b){F.addin(a, b);}, FieldCategories::GenericTag());
		if(A.mone != nullptr)
			sparse_details_impl::fspmm_zo(F, *(A.mone), blockSize, x, y, [F](Element & a, const Element & b){F.subin(a, b);}, FieldCategories::GenericTag());
		if(A.dat != nullptr)
			sparse_details_impl::fspmm(F, *(A.dat), blockSize, x, y, FieldCategories::GenericTag());
	}

	template<class Field>
	inline void fspmm(const Field & F, const Sparse<Field, SparseMatrix_t::HYB_ZO> & A, int blockSize, typename Field::ConstElement_ptr x, int ldx,
		      typename Field::Element_ptr y, int ldy, FieldCategories::GenericTag){
		using Element = typename Field::Element;
		if(A.one != nullptr)
			sparse_details_impl::fspmm_zo(F, *(A.one), blockSize, x, ldx, y, ldy, [F](Element & a, const Element & b){F.addin(a, b);}, FieldCategories::GenericTag());
		if(A.mone != nullptr)
			sparse_details_impl::fspmm_zo(F, *(A.mone), blockSize, x, ldx, y, ldy, [F](Element & a, const Element & b){F.subin(a, b);}, FieldCategories::GenericTag());
		if(A.dat != nullptr)
			sparse_details_impl::fspmm(F, *(A.dat), blockSize, x, ldx, y, ldy, FieldCategories::GenericTag());
	}

	template<class Field>
	inline void fspmm(const Field & F, const Sparse<Field, SparseMatrix_t::HYB_ZO> & A, int blockSize, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, FieldCategories::UnparametricTag){
		using Element = typename Field::Element;
		if(A.one != nullptr)
			sparse_details_impl::fspmm_zo(F, *(A.one), blockSize, x, y, [F](Element & a, const Element & b){F.addin(a, b);}, FieldCategories::UnparametricTag());
		if(A.mone != nullptr)
			sparse_details_impl::fspmm_zo(F, *(A.mone), blockSize, x, y, [F](Element & a, const Element & b){F.subin(a, b);}, FieldCategories::UnparametricTag());
		if(A.dat != nullptr)
			sparse_details_impl::fspmm(F, *(A.dat), blockSize, x, y, FieldCategories::UnparametricTag());
	}

	template<class Field>
	inline void fspmm(const Field & F, const Sparse<Field, SparseMatrix_t::HYB_ZO> & A, int blockSize, typename Field::ConstElement_ptr x, int ldx,
		      typename Field::Element_ptr y, int ldy, FieldCategories::UnparametricTag){
		using Element = typename Field::Element;
		if(A.one != nullptr)
			sparse_details_impl::fspmm_zo(F, *(A.one), blockSize, x, ldx, y, ldy, [F](Element & a, const Element & b){F.addin(a, b);}, FieldCategories::UnparametricTag());
		if(A.mone != nullptr)
			sparse_details_impl::fspmm_zo(F, *(A.mone), blockSize, x, ldx, y, ldy [F](Element & a, const Element & b){F.subin(a, b);}, FieldCategories::UnparametricTag());
		if(A.dat != nullptr)
			sparse_details_impl::fspmm(F, *(A.dat), blockSize, x, ldx, y, ldy, FieldCategories::UnparametricTag());
	}

	template<class Field>
	inline void fspmm(const Field & F, const Sparse<Field, SparseMatrix_t::HYB_ZO> & A, int blockSize, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, uint64_t kmax){
		using Element = typename Field::Element;
		if(A.one != nullptr)
			sparse_details_impl::fspmm_zo(F, *(A.one), blockSize, x, y, [F](Element & a, const Element & b){F.addin(a, b);}, FieldCategories::UnparametricTag());
		if(A.mone != nullptr)
			sparse_details_impl::fspmm_zo(F, *(A.mone), blockSize, x, y, [F](Element & a, const Element & b){F.subin(a, b);}, FieldCategories::UnparametricTag());
		if(A.dat != nullptr)
			sparse_details_impl::fspmm(F, *(A.dat), blockSize, x, y, kmax);
	}

	template<class Field>
	inline void fspmm(const Field & F, const Sparse<Field, SparseMatrix_t::HYB_ZO> & A, int blockSize, typename Field::ConstElement_ptr x, int ldx,
		      typename Field::Element_ptr y, int ldy, uint64_t kmax){
		using Element = typename Field::Element;
		if(A.one != nullptr)
			sparse_details_impl::fspmm_zo(F, *(A.one), blockSize, x, ldx, y, ldy, [F](Element & a, const Element & b){F.addin(a, b);}, FieldCategories::UnparametricTag());
		if(A.mone != nullptr)
			sparse_details_impl::fspmm_zo(F, *(A.mone), blockSize, x, ldx, y, ldy, [F](Element & a, const Element & b){F.subin(a, b);}, FieldCategories::UnparametricTag());
		if(A.dat != nullptr)
			sparse_details_impl::fspmm(F, *(A.dat), blockSize, x, ldx, y, ldy, kmax);
	}

}// HYB_ZO_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_HYB_ZO_spmm_INL