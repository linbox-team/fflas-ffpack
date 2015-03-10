/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   BB <bbboyer@ncsu.edu>
 *              Bastien Vialla <bastien.vialla@lirmm.fr>
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
 * Lesser General Public License for more ellr_details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/** @file fflas/fflas_fspmv_ellr.inl
 * NO DOC
 */

#ifndef __FFLASFFPACK_fflas_fflas_spmv_ellr_INL
#define __FFLASFFPACK_fflas_fflas_spmv_ellr_INL

namespace FFLAS { /*  ELLR */

	template<class Field>
	struct ELLR {
		size_t m ;
		size_t n ;
		size_t  ld ;
		index_t * row ;
		index_t * col ;
		typename Field::Element_ptr dat;
	};

	template<class Field>
	struct ELLR_sub : public ELLR<Field> {
	};

	template<class Field>
	struct ELLR_ZO  : ELLR<Field>{
		typename Field::Element cst ;
	};

	namespace ellr_details {

		// y = A x + b y ; (generic)
		template<class Field>
		void fspmv(
			      const Field& F,
			      const size_t m,
			      const size_t n,
			      const size_t ld,
			      const index_t * row,
			      const index_t * col,
			      const typename Field::Element_ptr  dat,
			      const typename Field::Element_ptr x ,
			      typename Field::Element_ptr y,
			      FieldCategories::GenericTag
			     )
		{
			for (size_t i = 0 ; i < m ; ++i) {
				for (index_t j = 0 ; j < row[i] ; ++j) {
					F.axpyin(y[i],dat[i*ld+j],x[col[i*ld+j]]);
				}
			}
		}

		// y = A x + b y ; (generic)
		template<class Field>
		void fspmv(
			      const Field& F,
			      const size_t m,
			      const size_t n,
			      const size_t ld,
			      const index_t * row,
			      const index_t * col,
			      const typename Field::Element_ptr  dat,
			      const typename Field::Element_ptr x ,
			      typename Field::Element_ptr y,
			      FieldCategories::UnparametricTag
			     )
		{
			for (size_t i = 0 ; i < m ; ++i) {
				for (index_t j = 0 ; j < row[i] ; ++j) {
					y[i] += dat[i*ld+j]*x[col[i*ld+j]];
				}
			}
		}

		// delayed by kmax
		template<class Field>
		void fspmv(
			      const Field& F,
			      const size_t m,
			      const size_t n,
			      const size_t ld,
			      const index_t * row,
			      const index_t * col,
			      const typename Field::Element_ptr  dat,
			      const typename Field::Element_ptr x ,
			      typename Field::Element_ptr y,
			      const index_t & kmax
			     )
		{

			for (size_t i = 0 ; i < m ; ++i) {
				index_t block = (row[i])/kmax ; // use DIVIDE_INTO from fspmvgpu
				index_t j = 0;
				index_t j_loc = 0 ;
				for (index_t l = 0 ; l < block ; ++l) {
					j_loc += kmax ;
					for ( ; j < j_loc ; ++j ) {
						y[i] += dat[i*ld+j] * x[col[i*ld+j]];
					}
					F.reduce(y[i]);
				}
				for ( ; j < row[i]  ; ++j) {
					y[i] += dat[i*ld+j] * x[col[i*ld+j]];
				}
				F.reduce (y[i]);
			}
		}

		// generic
		template<class Field, bool add>
		void fspmv_zo(
				 const Field & F,
				 const size_t m,
				 const size_t n,
				 const size_t ld,
				 const index_t * row,
				 const index_t * col,
				 const typename Field::Element_ptr x ,
				 typename Field::Element_ptr y,
				 FieldCategories::GenericTag
				)
		{
			if(add){
				for (size_t i = 0 ; i < m ; ++i) {
					for (index_t j = 0 ; j < row[i] ; ++j){
						F.addin(y[i], x[col[i*ld+j]]);
					}
				}
			}else{
				for (size_t i = 0 ; i < m ; ++i) {
					for (index_t j = 0 ; j < row[i] ; ++j){
						F.subin(y[i], x[col[i*ld+j]]);
					}
				}
			}
		}

		template<class Field, bool add>
		void fspmv_zo(
				 const Field & F,
				 const size_t m,
				 const size_t n,
				 const size_t ld,
				 const index_t * row,
				 const index_t * col,
				 const typename Field::Element_ptr x ,
				 typename Field::Element_ptr y,
				 FieldCategories::UnparametricTag
				)
		{
			if(add){
				for (size_t i = 0 ; i < m ; ++i) {
					for (index_t j = 0 ; j < row[i] ; ++j){
						y[i] += x[col[i*ld+j]];
					}
				}
			}else{
				for (size_t i = 0 ; i < m ; ++i) {
					for (index_t j = 0 ; j < row[i] ; ++j){
						y[i] -= x[col[i*ld+j]];
					}
				}
			}
		}

	} // ellr_details

	/* ******* */
	/* ELLR_sub */
	/* ******* */

	// y = A x + b y ; (generic)
	// it is supposed that no reduction is needed.
	template<class Field>
	void fspmv(
		      const Field& F,
		      const ELLR_sub<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     )
	{
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::category());
		fspmv(F, A, x, y, typename FieldTraits<Field>::category());
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const ELLR_sub<typename Field::Element> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::ModularTag)
	{
		ellr_details::fspmv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag ());
		freduce (F,A.m,y.dat,1);
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const ELLR_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::UnparametricTag)
	{
		ellr_details::fspmv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag() );
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const ELLR_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::GenericTag)
	{
		ellr_details::fspmv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,y.dat, FieldCategories::GenericTag() );
	}

	/* ***** */
	/* ELL_R */
	/* ***** */

	// y = A x + b y ; (generic)
	// reductions are delayed.
	template<class Field>
	void fspmv(
		      const Field& F,
		      const ELLR<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     )
	{
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::category());
		fspmv(F,A,x,y, typename FieldTraits<Field>::category());
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const ELLR<Field> & A,
		      const VECT<Field> & x,
		      VECT<Field> & y
		      , FieldCategories::GenericTag
		     )
	{
		ellr_details::fspmv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,y.dat, FieldCategories::GenericTag());
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const ELLR<Field> & A,
		      const VECT<Field> & x,
		      VECT<Field> & y
		      , FieldCategories::UnparametricTag
		     )
	{
		ellr_details::fspmv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag());
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const ELLR<Field> & A,
		      const VECT<Field> & x,
		      VECT<Field> & y
		      , FieldCategories::ModularTag
		     )
	{
		size_t kmax = Protected::DotProdBoundClassic(F,F.one) ;
		ellr_details::fspmv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,y.dat, (index_t) kmax);
	}

	// y = A x + b y ; (generic)
	// it is supposed that no reduction is needed.
	template<class Field>
	inline void fspmv(
			     const Field & F,
			     const ELLR_ZO<Field> & A,
			     const VECT<Field> & x,
			     const typename Field::Element & b,
			     VECT<Field> & y
			    )
	{
		details::init_y(F, y.m, b, y.dat, typename FieldTraits<Field>::category());
		fspmv(F, A, x, y, typename FieldTraits<Field>::category() );
	}

	template<class Field>
	inline void fspmv(
			     const Field & F,
			     const ELLR_ZO<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y
			     , FieldCategories::GenericTag
			    )
	{
		if (A.cst == F.one) {
			ellr_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.row,A.col,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else if (A.cst == F.mOne) {
			ellr_details::fspmv_zo<Field,false>(F,A.m,A.n,A.ld,A.row,A.col,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else {
			ellr_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.row,A.col,x.dat,y.dat, FieldCategories::GenericTag());
			fscalin(F,A.n,A.cst,y.dat,1);
		}
	}

	template<class Field>
	inline void fspmv(
			     const Field & F,
			     const ELLR_ZO<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y
			     , FieldCategories::UnparametricTag
			    )
	{
		if (A.cst == F.one) {
			ellr_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.row,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else if (A.cst == F.mOne) {
			ellr_details::fspmv_zo<Field,false>(F,A.m,A.n,A.ld,A.row,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else {
			ellr_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.row,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
			fscalin(F,A.n,A.cst,y.dat,1);
		}
	}

	template<class Field>
	inline void fspmv(
			     const Field & F,
			     const ELLR_ZO<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y
			     , FieldCategories::ModularTag
			    )
	{
		if (A.cst == F.one) {
			ellr_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.row,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else if (A.cst == F.mOne) {
			ellr_details::fspmv_zo<Field,false>(F,A.m,A.n,A.ld,A.row,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else {
			ellr_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.row,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
			fscalin(F,A.n,A.cst,y.dat,1);
		}
		freduce (F,y.m,y.dat,1);
	}
} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_spmv_ellr_INL

