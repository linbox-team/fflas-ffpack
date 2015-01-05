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

#ifndef __FFLASFFPACK_fflas_sparse_coo_spmv_INL
#define __FFLASFFPACK_fflas_sparse_coo_spmv_INL

namespace FFLAS{
	namespace sparse_details_impl{
	template<class Field>
	inline void fspmv(const Field & F, const Sparse<Field, SparseMatrix_t::COO> & A, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, FieldCategories::GenericTag){
		index_t j = 0;
		for (; j < ROUND_DOWN(A.nnz, 4) ; j+=4){
			F.axpyin(y[A.row[j]],A.dat[j],x[A.col[j]]);
			F.axpyin(y[A.row[j+1]],A.dat[j+1],x[A.col[j+1]]);
			F.axpyin(y[A.row[j+2]],A.dat[j+2],x[A.col[j+2]]);
			F.axpyin(y[A.row[j+3]],A.dat[j+3],x[A.col[j+3]]);
		}
		for(; j < A.nnz ; ++j)
		{
			F.axpyin(y[A.row[j]],A.dat[j],x[A.col[j]]);
		}
	}

	template<class Field>
	inline void fspmv(const Field & F, const Sparse<Field, SparseMatrix_t::COO> & A, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, FieldCategories::UnparametricTag){
		index_t j = 0;
		for (; j < ROUND_DOWN(A.nnz, 4) ; j+=4){
			y[A.row[j]] += A.dat[j]*x[A.col[j]];
			y[A.row[j+1]] += A.dat[j+1]*x[A.col[j+1]];
			y[A.row[j+2]] += A.dat[j+2]*x[A.col[j+2]];
			y[A.row[j+3]] += A.dat[j+3]*x[A.col[j+3]];
		}
		for(; j < A.nnz ; ++j)
		{
			y[A.row[j]] += A.dat[j]*x[A.col[j]];
		}
	}

	template<class Field>
	inline void fspmv(const Field & F, const Sparse<Field, SparseMatrix_t::COO> & A, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, const int64_t kmax){
		size_t w = 0 ;
		index_t last_i = 0;
		typename Field::Element e ;
		F.init(e,y[last_i]);
		size_t accu = 0 ;

		while ( w < A.nnz) {
			if ( A.row[w] == last_i ) { // same line
				if (accu < (size_t)kmax) {
					e += A.dat[w] * x[A.col[w]] ;
					accu += 1 ;
				}
				else {
					F.axpyin(e,A.dat[w],x[A.col[w]]);
					accu = 0 ;
				}
			}
			else { // new line
				F.init(y[last_i],e);
				last_i = A.row[w] ;
				F.init(e,y[last_i]);
				e += A.dat[w] * x[A.col[w]];
				accu = 1 ;
			}
			++w ;
		}
		F.init(y[last_i],e);
	}

	template<class Field, class Func>
	inline void fspmv(const Field & F, const Sparse<Field, SparseMatrix_t::COO_ZO> & A, typename Field::ConstElement_ptr x,
			typename Field::Element_ptr y, Func && func){
		index_t j = 0;
		for (; j < ROUND_DOWN(A.nnz, 4) ; j+=4){
			func(y[A.row[j]], x[A.col[j]]);
			func(y[A.row[j+1]], x[A.col[j+1]]);
			func(y[A.row[j+2]], x[A.col[j+2]]);
			func(y[A.row[j+3]], x[A.col[j+3]]);
		}
		for(; j < A.nnz ; ++j)
		{
			func(y[A.row[j]], x[A.col[j]]);
		}
	}
}// coo_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_coo_spmv_INL