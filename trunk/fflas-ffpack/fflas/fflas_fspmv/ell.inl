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
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/** @file fflas/fflas_fspmv_ell.inl
 * NO DOC
 */

#ifndef __FFLASFFPACK_fflas_fflas_spmv_ell_INL
#define __FFLASFFPACK_fflas_fflas_spmv_ell_INL

namespace FFLAS { namespace ell_details {

	// y = A x + b y ; (generic)
	template<class Field>
	inline void fspmv(
			     const Field& F,
			     const size_t m,
			     const size_t n,
			     const size_t ld,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const typename Field::Element_ptr x ,
			     typename Field::Element_ptr y,
			     FieldCategories::GenericTag
			    )
	{
		for (size_t i = 0 ; i < m ; ++i, dat+=ld, col+=ld) {
			for (index_t j = 0 ; j < ld ; ++j) {
				F.axpyin(y[i],dat[j],x[col[j]]);
			}
		}
	}


	template<class Field>
	inline void fspmv(
			     const Field & F,
			     const size_t m,
			     const size_t n,
			     const size_t ld,
			     const index_t * col,
			     typename Field::Element_ptr dat,
			     typename Field::Element_ptr x ,
			     typename Field::Element_ptr y,
			     FieldCategories::UnparametricTag
			    )
	{
		for (size_t i = 0 ;  i < m ; ++i) {
			for (index_t j = 0 ; j < ld ; ++j) {
				y[i] += dat[i*ld+j]*x[col[i*ld+j]];
			}
		}
	}


	// delayed by kmax
	template<class Field>
	inline void fspmv(
			     const Field& F,
			     const size_t m,
			     const size_t n,
			     const size_t ld,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const typename Field::Element_ptr x ,
			     typename Field::Element_ptr y,
			     const index_t & kmax
			    )
	{
		index_t block = (ld)/kmax ; // use DIVIDE_INTO from fspmvgpu
		for (size_t i = 0 ; i < m ; ++i) {
			index_t j = 0;
			index_t j_loc = 0 ;
			for (size_t l = 0 ; l < block ; ++l) {
				j_loc += kmax ;
				for ( ; j < j_loc ; ++j ) {
					y[i] += dat[i*ld+j] * x[col[i*ld+j]];
				}
				F.init(y[i],y[i]);
			}
			for ( ; j < ld  ; ++j) {
				y[i] += dat[i*ld+j] * x[col[i*ld+j]];
			}
			F.init(y[i],y[i]);
		}
	}

} // details
} // FFLAS

namespace FFLAS {

	/* ******* */
	/* ELL_sub */
	/* ******* */

	// y = A x + b y ; (generic)
	// it is supposed that no reduction is needed.
	template<class Field>
	inline void fspmv(
			     const Field& F,
			     const ELL_sub<Field> & A,
			     const VECT<Field> & x,
			     const typename Field::Element & b,
			     VECT<Field> & y
			    )
	{
		details::init_y(F, y.m, b, y.dat, typename FieldTraits<Field>::value());
		fspmv( F, A, x, y, typename FieldTraits<Field>::category());
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const ELL_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::ModularTag)
	{
		ell_details::fspmv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag ());
		finit(F,A.m,y.dat,1);
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const ELL_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::UnparametricTag)
	{
		ell_details::fspmv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag());
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const ELL_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::GenericTag)
	{
		ell_details::fspmv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,y.dat, FieldCategories::GenericTag());
	}

	/* *** */
	/* ELL */
	/* *** */

	// y = A x + b y ; (generic)
	// reductions are delayed.
	template<class Field>
	void fspmv(
		      const Field& F,
		      const ELL<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     )
	{
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::value());
		fspmv(F,A,x,y,typename FieldTraits<Field>::category());
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const ELL<Field> & A,
		      const VECT<Field> & x,
		      VECT<Field> & y,
		      FieldCategories::GenericTag
		     )
	{
		ell_details::fspmv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,y.dat, FieldCategories::GenericTag());
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const ELL<Field> & A,
		      const VECT<Field> & x,
		      VECT<Field> & y, 
		      FieldCategories::UnparametricTag
		     )
	{
		ell_details::fspmv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag());
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const ELL<Field> & A,
		      const VECT<Field> & x,
		      VECT<Field> & y,
		      FieldCategories::ModularTag
		     )
	{
		size_t kmax = Protected::DotProdBoundClassic(F,F.one) ;
		ell_details::fspmv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,y.dat, (index_t) kmax);
	}
} // FFLAS

namespace FFLAS { namespace ell_details { /*  ZO */

	// generic
	template<class Field, bool add>
	inline void fspmv_zo(
				const Field & F,
				const size_t m,
				const size_t n,
				const size_t ld,
				const index_t * col,
				const typename Field::Element_ptr x ,
				typename Field::Element_ptr y,
				FieldCategories::GenericTag
			    )
	{
		if(add){
			for (size_t i = 0 ; i < m ; ++i) {
				for (index_t j = 0 ; j < ld ; ++j) {
					F.addin(y[i],x[col[i*ld+j]]);
				}
			}	
		}else{
			for (size_t i = 0 ; i < m ; ++i) {
				for (index_t j = 0 ; j < ld ; ++j) {
					F.subin(y[i],x[col[i*ld+j]]);
				}
			}	
		}
	}


	template<class Field, bool add>
	inline void fspmv_zo(
				const Field & F,
				const size_t m,
				const size_t n,
				const size_t ld,
				const index_t * col,
				const typename Field::Element_ptr x ,
				typename Field::Element_ptr y,
				FieldCategories::UnparametricTag
			    )
	{
		if(add){
			for ( size_t i = 0 ;  i < m   ; ++i ) {
				for (index_t j = 0 ; j < ld ; ++j) {
					y[i] += x[col[i*ld+j]];
				}
			}
		}else{
			for ( size_t i = 0 ;  i < m   ; ++i ) {
				for (index_t j = 0 ; j < ld ; ++j) {
					y[i] -= x[col[i*ld+j]];
				}
			}
		}
	}


} // details
} // FFLAS

namespace FFLAS { /*  ZO */

	/* ******* */
	/* ELL_ZO  */
	/* ******* */

	// y = A x + b y ; (generic)
	// it is supposed that no reduction is needed.
	template<class Field>
	inline void fspmv(
			     const Field & F,
			     const ELL_ZO<Field> & A,
			     const VECT<Field> & x,
			     const typename Field::Element & b,
			     VECT<Field> & y
			    )
	{
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::value());
		fspmv(F, A, x, y, typename FieldTraits<Field>::category() );
	}

	template<class Field>
	inline void fspmv(
			     const Field & F,
			     const ELL_ZO<Field > & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::GenericTag
			    )
	{
		if (A.cst == F.one) {
			ell_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.col,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else if (A.cst == F.mOne) {
			ell_details::fspmv_zo<Field,false>(F,A.m,A.n,A.ld,A.col,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else {
			auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x.dat, 1, x1, 1);
			ell_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.col,x1,y.dat, FieldCategories::GenericTag());
			fflas_delete(x1);
		}
	}

	template<class Field>
	inline void fspmv(
			     const Field & F,
			     const ELL_ZO<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::UnparametricTag
			    )
	{
		if (A.cst == F.one) {
			ell_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else if (A.cst == F.mOne) {
			ell_details::fspmv_zo<Field,false>(F,A.m,A.n,A.ld,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else {
			auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x.dat, 1, x1, 1);
			ell_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.col,x1,y.dat, FieldCategories::UnparametricTag());
			fflas_delete(x1);
		}
	}

	template<class Field>
	inline void fspmv(
			     const Field & F,
			     const ELL_ZO<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::ModularTag
			    )
	{
		if (A.cst == F.one) {
			ell_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else if (A.cst == F.mOne) {
			ell_details::fspmv_zo<Field,false>(F,A.m,A.n,A.ld,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else {
			auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x.dat, 1, x1, 1);
			ell_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.col,x1,y.dat, FieldCategories::UnparametricTag());
			fflas_delete(x1);
		}
		finit(F,A.m,y.dat,1);
	}


} // FFLAS

namespace FFLAS { /*  conversions */

	template<class Field>
	inline void print_ell(const ELL<Field> & M)
	{
		for(size_t j = 0 ; j < M.m ; ++j){
			std::cout << j << " | ";
			for(size_t i = 0 ; i < M.ld ; ++i)
				std::cout << (int64_t)M.col[j*M.ld+i] << " ";
			std::cout << std::endl;
		}
	}

	template<class Field, class ColT, class RowT>
	inline void sp_ell_from_csr(
				    const Field & F,
				    const size_t CSR_m,
				    const size_t CSR_n,
				    const size_t nnz, 
				    const ColT * CSR_col,
				    const RowT * CSR_row,
				    const typename Field::Element_ptr CSR_dat,
				    size_t & ELL_m,
				    size_t & ELL_n,
				    size_t & ld,
				    size_t & chunk,
				    index_t *& ELL_col,
				    typename Field::Element_ptr& ELL_dat,
				    const bool ZO
				   )
	{
		ld = 0;
		ELL_m = CSR_m;
		ELL_n = CSR_n;
		for(size_t i = 0 ; i < CSR_m ; ++i){
			if(CSR_row[i+1]-CSR_row[i] > ld){
				ld = CSR_row[i+1]-CSR_row[i];
			}
		}
		std::cout << "ELL no simd " << std::endl; 
		chunk = 1;
		ELL_col = fflas_new<index_t >(ld*ELL_m, Alignment::CACHE_LINE);
		if(!ZO){
			ELL_dat = fflas_new<typename Field::Element >(ld*ELL_m, Alignment::CACHE_LINE);
		}
		size_t it = 0;
		for(size_t i = 0 ; i < CSR_m ; ++i){
			size_t start = CSR_row[i], stop = CSR_row[i+1];
			for(size_t j = 0 ; j < ld ; ++j){
				if(start + j < stop){
					if(!ZO){
						ELL_dat[i*ld+j] = CSR_dat[start+j];
					}
					ELL_col[i*ld+j] = CSR_col[start+j];
				}
				else{
					if(!ZO){
						F.init(ELL_dat[i*ld+j], F.zero);
					}
					ELL_col[i*ld+j] = 0;
				}
			}
		}
	}

	template<class Field, class ColT, class RowT>
	inline void sp_ell_from_coo(
				    const Field & F,
				    const size_t COO_m,
				    const size_t COO_n,
				    const size_t COO_nnz,
				    const ColT * COO_col,
				    const RowT * COO_row,
				    const typename Field::Element_ptr COO_dat,
				    size_t & ELL_m,
				    size_t & ELL_n,
				    size_t & ld,
				    size_t & chunk,
				    index_t * &ELL_col,
				    typename Field::Element_ptr &ELL_dat,
				    const bool ZO
				   )
	{
		index_t * row = fflas_new<index_t >(COO_m+1);
		for(size_t i = 0 ; i <= COO_m+1 ; ++i){
            row[i] = 0;
        }
        for(size_t i = 0 ; i < COO_nnz ; ++i){
            row[COO_row[i]+1]++;
        }
        for(size_t i = 1 ; i <= COO_m ; ++i){
            row[i] += row[i-1];
        }
		sp_ell_from_csr(F, COO_m, COO_n, COO_nnz, COO_col, row, COO_dat, ELL_m, ELL_n, ld, chunk, ELL_col, ELL_dat, ZO);
		fflas_delete(row);
	}

} // FFLAS

namespace FFLAS{ /* delete */

	template<class Field>
	inline void sp_delete(const ELL<Field> & m){
		fflas_delete(m.dat);
		fflas_delete(m.col);
	}

	template<class Field>
	inline void sp_delete(const ELL_sub<Field> & m){
		fflas_delete(m.dat);
		fflas_delete(m.col);
	}

	template<class Field>
	inline void sp_delete(const ELL_ZO<Field> & m){
		fflas_delete(m.col);
	}

}

#endif // __FFLASFFPACK_fflas_fflas_spmv_ell_INL
