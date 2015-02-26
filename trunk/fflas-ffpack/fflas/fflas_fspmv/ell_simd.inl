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

#ifndef __FFLASFFPACK_fflas_fflas_spmv_ell_simd_INL
#define __FFLASFFPACK_fflas_fflas_spmv_ell_simd_INL

#include "fflas-ffpack/fflas/fflas_simd.h"

namespace FFLAS { namespace ell_simd_details {

	// y = A x + b y ; (generic)
	template<class Field>
	inline void fspmv(
			  const Field& F,
			  const size_t m,
			  const size_t n,
			  const size_t ld,
			  const size_t chunk,
			  const index_t * col,
			  const typename Field::Element_ptr dat,
			  const typename Field::Element_ptr x ,
			  typename Field::Element_ptr y,
			  FieldCategories::GenericTag
			 )
	{
		size_t end = (m % chunk == 0) ? m : m + m%chunk;
		for ( size_t i = 0 ; i < end/chunk ; ++i, y+=chunk ) {
			for (index_t j = 0 ; j < ld ; ++j, dat+=chunk, col+=chunk) {
				for(size_t k = 0 ; k < chunk ; ++k)
					F.axpyin(y[k], dat[k], x[col[k]]);
			}
		}
	}


	template<class Field>
	inline void fspmv(
			  const Field & F,
			  const size_t m,
			  const size_t n,
			  const size_t ld,
			  const size_t chunk,
			  const index_t * col,
			  typename Field::Element_ptr dat,
			  typename Field::Element_ptr x ,
			  typename Field::Element_ptr y,
			  FieldCategories::UnparametricTag
			 )
	{
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;
		int end = (m%chunk == 0) ? m : m+(chunk-(m%chunk));
		std::cout << "end " << end << std::endl;
		vect_t X,Y,D,T ;
		for(int i = 0 ; i < end ; i+=chunk) {
			Y = simd::load(y+i);
			for (index_t j = 0 ; j < ld ; ++j) {
				D = simd::load(dat+i*ld+j*chunk);
				X = simd::gather(x,col+i*ld+j*chunk);
				Y = simd::fmadd(Y,D,X);
			}
			simd::stream(y+i,Y);
		}
	}

	/* NO SIMD VERSION */
	// 			size_t end = (m%chunk == 0)? m : m+m%chunk;
	// 			for( size_t i = 0 ; i < end ; i+=chunk) {
	// 				for (index_t j = 0 ; j < ld ; ++j, dat+=chunk, col+=chunk) {
	// 					for(size_t k = 0 ; k < chunk ; ++k)
	// 					{
	// 						y[i+k] += dat[k]*x[col[k]];
	// 					}
	// 				}
	// 			}


	// delayed by kmax
	template<class Field>
	inline void fspmv(
			  const Field& F,
			  const size_t m,
			  const size_t n,
			  const size_t ld,
			  const size_t chunk,
			  const index_t * col,
			  const typename Field::Element_ptr dat,
			  const typename Field::Element_ptr x ,
			  typename Field::Element_ptr y,
			  const index_t & kmax
			 )
	{
		index_t block = (ld)/kmax ; // use DIVIDE_INTO from fspmvgpu
		size_t end = (m%chunk == 0)? m : m+m%chunk;
		using simd = Simd<typename Field::Element >;
		using vect_t = typename simd::vect_t;

		vect_t X, Y, D, C, Q, TMP, NEGP, INVP, MIN, MAX, P;
		double p = (typename Field::Element)F.characteristic();

		P = simd::set1(p);
		NEGP = simd::set1(-p);
		INVP = simd::set1(1/p);
		MIN = simd::set1(F.minElement());
		MAX = simd::set1(F.maxElement());

		for ( size_t i = 0; i < end/chunk ; ++i ) {
			index_t j = 0 ;
			index_t j_loc = 0 ;
			Y = simd::load(y+i*chunk);
			for (size_t l = 0 ; l < block ; ++l) {
				j_loc += kmax ;

				for ( ; j < j_loc ; ++j) {
					D = simd::load(dat+i*chunk*ld+j*chunk);
					X = simd::gather(x,col+i*chunk*ld+j*chunk);
					Y = simd::fmadd(Y,D,X);
				}
				simd::mod(Y,P, INVP, NEGP, MIN, MAX, Q, TMP);
			}
			for ( ; j < ld ; ++j) {
				D = simd::load(dat+i*chunk*ld+j*chunk);
				X = simd::gather(x,col+i*chunk*ld+j*chunk);
				Y = simd::fmadd(Y,D,X);
			}
			simd::mod(Y,P, INVP, NEGP, MIN, MAX, Q, TMP);
			simd::store(y+i*chunk,Y);
		}
	}
	/* NO SIMD VERSION */
	// 	for ( size_t i = 0; i < end/chunk ; ++i ) {
	// 	index_t j = 0 ;
	// 	index_t j_loc = 0 ;
	// 	for (size_t l = 0 ; l < block ; ++l) {
	// 		j_loc += kmax ;

	// 		for ( ; j < j_loc ; ++j) {
	// 			for(size_t k = 0 ; k < chunk ; ++k)
	// 			{
	// 				y[i*chunk+k] += dat[i*chunk*ld+j*chunk+k]*x[col[i*chunk*ld+j*chunk+k]];
	// 			}
	// 		}
	// 		for(size_t k = 0 ; k < chunk ; ++k)
	// 		{
	// 			F.reduce (y[i*chunk+k]);
	// 		}
	// 	}
	// 	for ( ; j < ld ; ++j) {
	// 		for(size_t k = 0 ; k < chunk ; ++k)
	// 		{
	// 			y[i*chunk+k] += dat[i*chunk*ld+j*chunk+k]*x[col[i*chunk*ld+j*chunk+k]];
	// 		}
	// 	}
	// 	for(size_t k = 0 ; k < chunk ; ++k)
	// 	{
	// 		F.reduce (y[i*chunk+k]);
	// 	}
	// }

} // details
} // FFLAS

namespace FFLAS { namespace ell_simd_details { /*  ZO */

	// generic
	template<class Field, bool add>
	inline void fspmv_zo(
			     const Field & F,
			     const size_t m,
			     const size_t n,
			     const size_t ld,
			     const size_t chunk,
			     const index_t * col,
			     const typename Field::Element_ptr x ,
			     typename Field::Element_ptr y,
			     FieldCategories::GenericTag
			    )
	{
		size_t end = (m % chunk == 0) ? m : m + (chunk - m%chunk);
		if(add){
			for ( size_t i = 0 ; i < end/chunk ; ++i ) {
				for (index_t j = 0 ; j < ld ; ++j) {
					for(size_t k = 0 ; k < chunk ; ++k)
						F.addin(y[i*chunk+k], x[col[i*chunk*ld+j*chunk+k]]);
				}
			}
		}else{
			for ( size_t i = 0 ; i < end/chunk ; ++i ) {
				for (index_t j = 0 ; j < ld ; ++j) {
					for(size_t k = 0 ; k < chunk ; ++k)
						F.subin(y[i*chunk+k], x[col[i*chunk*ld+j*chunk+k]]);
				}
			}
		}
	}


	template<class Field, bool add>
	inline void fspmv_zo(const Field & F,
			     const size_t m,
			     const size_t n,
			     const size_t ld,
			     const size_t chunk,
			     const index_t * col,
			     const typename Field::Element_ptr x,
			     typename Field::Element_ptr y,
			     FieldCategories::UnparametricTag)
	{
		using simd = Simd<typename Field::Element >;
		using vect_t = typename simd::vect_t;
		size_t end = (m%chunk == 0)? m : m+(chunk-m%chunk);
		vect_t X,Y,D ;
		if(add){ /*  this is compile time decision */
			for( size_t i = 0 ; i < end/chunk ; ++i ) {
				Y = simd::load(y+i*simd::vect_size);
				for (index_t j = 0 ; j < ld ; ++j) {
					X = simd::gather(x,col+i*chunk*ld+j*chunk);
					Y = simd::add(Y,X);
				}
				simd::stream(y+i*chunk,Y);
			}
		}else{
			for( size_t i = 0 ; i < end/chunk ; ++i ) {
				Y = simd::load(y+i*simd::vect_size);
				for (index_t j = 0 ; j < ld ; ++j) {
					X = simd::gather(x,col+i*chunk*ld+j*chunk);
					Y = simd::sub(Y,X);
				}
				simd::stream(y+i*chunk,Y);
			}
		}
	}

	/* NO SIMD VERSION */

	// size_t end = (m%chunk == 0)? m : m+(chunk-m%chunk);
	// if(add){
	// 	for( size_t i = 0 ; i < end/chunk ; ++i ) {
	// 		for (index_t j = 0 ; j < ld ; ++j) {
	// 			for(size_t k = 0 ; k < chunk ; ++k) {
	// 				y[i*chunk+k] += x[col[i*chunk*ld+j*chunk+k]];
	// 			}
	// 		}
	// 	}
	// }else{
	// 	for( size_t i = 0 ; i < end/chunk ; ++i ) {
	// 		for (index_t j = 0 ; j < ld ; ++j) {
	// 			for(size_t k = 0 ; k < chunk ; ++k) {
	// 				y[i*chunk+k] -= x[col[i*chunk*ld+j*chunk+k]];
	// 			}
	// 		}
	// 	}
	// }

} // details
} // FFLAS

namespace FFLAS {

	/* ******* */
	/* ELL_simd_sub */
	/* ******* */

	// y = A x + b y ; (generic)
	// it is supposed that no reduction is needed.
	template<class Field>
	inline void fspmv(
			  const Field& F,
			  const ELL_simd_sub<Field> & A,
			  const VECT<Field> & x,
			  const typename Field::Element & b,
			  VECT<Field> & y
			 )
	{
		details::init_y(F, y.m, b, y.dat, typename FieldTraits<Field>::category());
		fspmv( F, A, x, y, typename FieldTraits<Field>::category());
	}

	template<class Field>
	inline void fspmv(const Field & F,
			  const ELL_simd_sub<Field> & A,
			  const VECT<Field> & x,
			  VECT<Field> & y,
			  FieldCategories::ModularTag)
	{
		ell_simd_details::fspmv(F,A.m,A.n,A.ld,A.chunk,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag ());
		freduce (F,A.m,y.dat,1);
	}

	template<class Field>
	inline void fspmv(const Field & F,
			  const ELL_simd_sub<Field> & A,
			  const VECT<Field> & x,
			  VECT<Field> & y,
			  FieldCategories::UnparametricTag)
	{
		ell_simd_details::fspmv(F,A.m,A.n,A.ld,A.chunk,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag());
	}

	template<class Field>
	inline void fspmv(const Field & F,
			  const ELL_simd_sub<Field> & A,
			  const VECT<Field> & x,
			  VECT<Field> & y,
			  FieldCategories::GenericTag)
	{
		ell_simd_details::fspmv(F,A.m,A.n,A.ld,A.chunk,A.col,A.dat,x.dat,y.dat, FieldCategories::GenericTag());
	}

	/* *** */
	/* ELL_simd */
	/* *** */

	// y = A x + b y ; (generic)
	// reductions are delayed.
	template<class Field>
	void fspmv(
		   const Field& F,
		   const ELL_simd<Field> & A,
		   const VECT<Field> & x,
		   const typename Field::Element & b,
		   VECT<Field> & y
		  )
	{
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::category());
		fspmv(F,A,x,y,typename FieldTraits<Field>::category());
	}

	template<class Field>
	void fspmv(
		   const Field& F,
		   const ELL_simd<Field> & A,
		   const VECT<Field> & x,
		   VECT<Field> & y,
		   FieldCategories::GenericTag
		  )
	{
		ell_simd_details::fspmv(F,A.m,A.n,A.ld,A.chunk,A.col,A.dat,x.dat,y.dat, FieldCategories::GenericTag());
	}

	template<class Field>
	void fspmv(
		   const Field& F,
		   const ELL_simd<Field> & A,
		   const VECT<Field> & x,
		   VECT<Field> & y,
		   FieldCategories::UnparametricTag
		  )
	{
		ell_simd_details::fspmv(F,A.m,A.n,A.ld,A.chunk,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag());
	}

	template<class Field>
	void fspmv(
		   const Field& F,
		   const ELL_simd<Field> & A,
		   const VECT<Field> & x,
		   VECT<Field> & y,
		   FieldCategories::ModularTag
		  )
	{
		size_t kmax = Protected::DotProdBoundClassic(F,F.one) ;
		ell_simd_details::fspmv(F,A.m,A.n,A.ld,A.chunk,A.col,A.dat,x.dat,y.dat, (index_t) kmax);
	}
} // FFLAS

namespace FFLAS { /*  ZO */

	/* ******* */
	/* ELL_simd_ZO  */
	/* ******* */

	// y = A x + b y ; (generic)
	// it is supposed that no reduction is needed.
	template<class Field>
	inline void fspmv(
			  const Field & F,
			  const ELL_simd_ZO<Field> & A,
			  const VECT<Field> & x,
			  const typename Field::Element & b,
			  VECT<Field> & y
			 )
	{
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::category());
		fspmv(F, A, x, y, typename FieldTraits<Field>::category() );
	}

	template<class Field>
	inline void fspmv(
			  const Field & F,
			  const ELL_simd_ZO<Field > & A,
			  const VECT<Field> & x,
			  VECT<Field> & y,
			  FieldCategories::GenericTag
			 )
	{
		if (A.cst == F.one) {
			ell_simd_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.chunk,A.col,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else if (A.cst == F.mOne) {
			ell_simd_details::fspmv_zo<Field,false>(F,A.m,A.n,A.ld,A.chunk,A.col,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else {
			auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x.dat, 1, x1, 1);
			ell_simd_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.chunk,A.col,x1,y.dat, FieldCategories::GenericTag());
			fflas_delete(x1);
		}
	}

	template<class Field>
	inline void fspmv(
			  const Field & F,
			  const ELL_simd_ZO<Field> & A,
			  const VECT<Field> & x,
			  VECT<Field> & y,
			  FieldCategories::UnparametricTag
			 )
	{
		if (A.cst == F.one) {
			ell_simd_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.chunk,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else if (A.cst == F.mOne) {
			ell_simd_details::fspmv_zo<Field,false>(F,A.m,A.n,A.ld,A.chunk,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else {
			auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x.dat, 1, x1, 1);
			ell_simd_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.chunk,A.col,x1,y.dat, FieldCategories::UnparametricTag());
			fflas_delete(x1);
		}
	}

	template<class Field>
	inline void fspmv(
			  const Field & F,
			  const ELL_simd_ZO<Field> & A,
			  const VECT<Field> & x,
			  VECT<Field> & y,
			  FieldCategories::ModularTag
			 )
	{
		if (A.cst == F.one) {
			ell_simd_details::fspmv_zo<Field,true>(F,A.m,A.n,A.ld,A.chunk,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else if (A.cst == F.mOne) {
			ell_simd_details::fspmv_zo<Field,false>(F,A.m,A.n,A.ld,A.chunk,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else {
			auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x.dat, 1, x1, 1);
			ell_simd_details::fspmv_zo<Field,true>(F,A.m,A.n,A.chunk,A.ld,A.col,x1,y.dat, FieldCategories::UnparametricTag());
			fflas_delete(x1);
		}
		freduce (F,A.m,y.dat,1);
	}
} // FFLAS

namespace FFLAS{

	template<class Field>
	inline void print_ell(const ELL_simd<Field> & M)
	{
		size_t m  = 0, ld = M.ld, chunk = M.chunk;
		m = (M.m % M.chunk == 0) ? M.m : M.m + (M.chunk-(M.m%M.chunk));
		// std::cout << "limit : " << m*ld << std::endl;
		for(size_t i = 0 ; i < m ; i+=chunk){
			for(size_t k = 0 ; k < chunk ; ++k){
				std::cout << i+k << " | ";
				for(size_t j = 0 ; j < ld ; ++j){
					// std::cout << i*ld+j*chunk+k << " " << std::endl;
					std::cout << M.dat[i*ld+j*chunk+k] << " ";
				}
				std::cout << std::endl;
			}
		}
	}

	template<class Field, class ColT, class RowT>
	inline void sp_ell_simd_from_csr(
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

		chunk = Simd<typename Field::Element >::vect_size;

		size_t m = (CSR_m%chunk == 0) ? CSR_m : CSR_m+(chunk-CSR_m%chunk);
		// std::cout << "ELL simd ; mod " << CSR_m%chunk << " ; m : " << m << " CSR_m : " << CSR_m << std::endl;

		ELL_col = fflas_new<index_t>(ld*m, Alignment::CACHE_LINE);
		if(!ZO){
			ELL_dat = fflas_new<typename Field::Element>(ld*m, Alignment::CACHE_LINE);
		}

		size_t i = 0;
		size_t end = CSR_m/chunk;
		for(; i < end ; ++i){
			for(size_t k = 0 ; k < chunk ; ++k){
				size_t start = CSR_row[i*chunk+k], stop = CSR_row[i*chunk+k+1];
				for(size_t j = 0 ; j < ld ; ++j){
					if(start + j < stop){
						if(!ZO){
							ELL_dat[i*chunk*ld+j*chunk+k] = CSR_dat[start+j];
						}
						ELL_col[i*ld*chunk+j*chunk+k] = CSR_col[start+j];
					}else{
						if(!ZO){
							ELL_dat[i*ld*chunk+j*chunk+k] = 0;
						}
						ELL_col[i*ld*chunk+j*chunk+k] = 0;
					}
				}
			}
		}
		if(CSR_m != m)
		{
			for(size_t j = 0 ; j < ld ; ++j){
				for(size_t k = 0 ; k < chunk ; ++k){
					if(i*chunk+k < CSR_m){
						size_t start = CSR_row[i*chunk+k], stop = CSR_row[i*chunk+k+1];
						if(start + j < stop){
							if(!ZO){
								ELL_dat[i*chunk*ld+j*chunk+k] = CSR_dat[start+j];
							}
							ELL_col[i*ld*chunk+j*chunk+k] = CSR_col[start+j];
						}else{
							if(!ZO){
								ELL_dat[i*ld*chunk+j*chunk+k] = 0;
							}
							ELL_col[i*ld*chunk+j*chunk+k] = 0;
						}
					}else{
						if(!ZO){

								ELL_dat[i*ld*chunk+j*chunk+k] = 0;
							}
							ELL_dat[i*ld*chunk+j*chunk+k] = 0;
					}
				}
			}
		}
	}

	template<class Field, class ColT, class RowT>
	inline void sp_ell_simd_from_coo(
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
		sp_ell_simd_from_csr(F, COO_m, COO_n, COO_nnz, COO_col, row, COO_dat, ELL_m, ELL_n, ld, chunk, ELL_col, ELL_dat, ZO);
		fflas_delete(row);
	}

	template<class Field>
	void sp_delete(ELL_simd<Field> & M){
		fflas_delete(M.dat);
		fflas_delete(M.col);
	}

	template<class Field>
	void sp_delete(ELL_simd_sub<Field> & M){
		fflas_delete(M.dat);
		fflas_delete(M.col);
	}

	template<class Field>
	void sp_delete(ELL_simd_ZO<Field> & M){
		fflas_delete(M.col);
	}

}// FFLAS

#endif // __FFLASFFPACK_fflas_fflas_spmv_ell_INL
