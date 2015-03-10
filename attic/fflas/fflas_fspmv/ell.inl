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

	// row major, lda == n
	template<class Field>
	inline void fspmm(
			     const Field& F,
			     const size_t m,
			     const size_t n,
			     const size_t ld,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const int blockSize,
			     const typename Field::Element_ptr x ,
			     typename Field::Element_ptr y,
			     FieldCategories::GenericTag
			    )
	{
		for (size_t i = 0 ; i < m ; ++i, dat+=ld, col+=ld) {
			for (index_t j = 0 ; j < ld ; ++j) {
				for(int k = 0 ; k < blockSize ; ++k)
					F.axpyin(y[i+k],dat[j],x[col[j]+k]);
			}
		}
	}

	// row major, lda != n
	template<class Field>
	inline void fspmm(
			     const Field& F,
			     const size_t m,
			     const size_t n,
			     const size_t ld,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const int lda,
			     const int blockSize,
			     const typename Field::Element_ptr x ,
			     typename Field::Element_ptr y,
			     FieldCategories::GenericTag
			    )
	{
		for (size_t i = 0 ; i < m ; ++i, dat+=ld, col+=ld) {
			for (index_t j = 0 ; j < ld ; ++j) {
				for(int k = 0 ; k < blockSize ; ++k)
					F.axpyin(y[i*lda+k],dat[j],x[col[j]*lda+k]);
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

	// row major, lda == n
	template<class Field, bool Aligned>
	inline void fspmm(
			     const Field & F,
			     const size_t m,
			     const size_t n,
			     const size_t ld,
			     const index_t * col,
			     typename Field::Element_ptr dat,
			     const int blockSize,
			     typename Field::Element_ptr x ,
			     typename Field::Element_ptr y,
			     FieldCategories::UnparametricTag
			    )
	{
#if defined(__FFLASFFPACK_USE_SIMD)
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;

		vect_t vx, vy, vdat;
		for(size_t i = 0 ; i <m ; ++i){
			for(size_t j = 0 ; j < ld ; ++j){
				int k = 0;
				vdat = simd::set1(dat+i*ld+j);
				for(; k < blockSize ; k+=simd::vect_size){
					if(Aligned){
						vy = simd::load(y+i+k);
						vx = simd::load(x+col[j+i*ld+j]+k);
						simd::store(y+i+k, simd::fmadd(vy, vx, vdat));
					}else{
						vy = simd::loadu(y+i+k);
						vx = simd::loadu(x+col[j+i*ld]+k);
						simd::storeu(y+i+k, simd::fmadd(vy, vx, vdat));
					}
				}
				for(; k < blockSize ; ++k){
					y[i+k] += dat[i*ld+j]*x[col[i*ld+j]+k];
				}
			}
		}
#else
		for (size_t i = 0 ;  i < m ; ++i) {
			for (index_t j = 0 ; j < ld ; ++j) {
				for(int k = 0 ; k < blockSize ; ++k)
					y[i+k] += dat[i*ld+j]*x[col[i*ld+j]+k];
			}
		}
#endif
	}

	// row major, lda != n
	template<class Field, bool Aligned>
	inline void fspmm(
			     const Field & F,
			     const size_t m,
			     const size_t n,
			     const size_t ld,
			     const index_t * col,
			     typename Field::Element_ptr dat,
			     const int lda,
			     const int blockSize,
			     typename Field::Element_ptr x ,
			     typename Field::Element_ptr y,
			     FieldCategories::UnparametricTag
			    )
	{
#if defined(__FFLASFFPACK_USE_SIMD)
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;

		vect_t vx, vy, vdat;
		for(size_t i = 0 ; i <m ; ++i){
			for(size_t j = 0 ; j < ld ; ++j){
				int k = 0;
				vdat = simd::set1(dat+i*ld+j);
				for(; k < blockSize ; k+=simd::vect_size){
					if(Aligned){
						vy = simd::load(y+i*lda+k);
						vx = simd::load(x+col[j+i*ld+j]*lda+k);
						simd::store(y+i*lda+k, simd::fmadd(vy, vx, vdat));
					}else{
						vy = simd::loadu(y+i*lda+k);
						vx = simd::loadu(x+col[j+i*ld]*lda+k);
						simd::storeu(y+i*lda+k, simd::fmadd(vy, vx, vdat));
					}
				}
				for(; k < blockSize ; ++k){
					y[i*lda+k] += dat[i*ld+j]*x[col[i*ld+j*lda]+k];
				}
			}
		}
#else
		for (size_t i = 0 ;  i < m ; ++i) {
			for (index_t j = 0 ; j < ld ; ++j) {
				for(int k = 0 ; k < blockSize ; ++k)
					y[i*lda+k] += dat[i*ld+j]*x[col[i*ld+j]*lda];
			}
		}
#endif
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
				F.reduce (y[i]);
			}
			for ( ; j < ld  ; ++j) {
				y[i] += dat[i*ld+j] * x[col[i*ld+j]];
			}
			F.reduce (y[i]);
		}
	}

	// delayed by kmax, row major, lda == n
	template<class Field, bool Aligned>
	inline void fspmm(
			     const Field& F,
			     const size_t m,
			     const size_t n,
			     const size_t ld,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const int blockSize,
			     const typename Field::Element_ptr x ,
			     typename Field::Element_ptr y,
			     const index_t & kmax
			    )
	{
#if defined(__FFLASFFPACK_USE_SIMD)
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;
		vect_t vdat, vx, vy;
		index_t block = (ld)/kmax ; // use DIVIDE_INTO from fspmvgpu
		for (size_t i = 0 ; i < m ; ++i) {
			index_t j = 0;
			index_t j_loc = 0 ;
			for (size_t l = 0 ; l < block ; ++l) {
				j_loc += kmax ;
				for ( ; j < j_loc ; ++j ) {
					int k = 0;
					vdat = simd::set1(dat+i*ld+j);
					for(; k < blockSize ; k += simd::vect_size){
						if(Aligned){
							vy = simd::load(y+i+k);
							vx = simd::load(x+col[i*ld+j]+k);
							simd::store(y+i+k, simd::fmadd(vy, vdat, vx));
						}else{
							vy = simd::loadu(y+i+k);
							vx = simd::loadu(x+col[i*ld+j]+k);
							simd::storeu(y+i+k, simd::fmadd(vy, vdat, vx));
						}
					}
					for(; k < blockSize ; ++k)
						y[i+k] += dat[i*ld+j] * x[col[i*ld+j]+k];
				}
				for(int k = 0 ; k < blockSize ; ++k)
					F.reduce (y[i+k]);
			}
			for ( ; j < ld  ; ++j) {
				int k = 0;
				vdat = simd::set1(dat[i*ld+j]);
				for( ; k < blockSize ; ++k){
					if(Aligned){
						vy = simd::load(y+i+k);
						vx = simd::load(x+col[i*ld+j]+k);
						simd::store(y+i+k, simd::fmadd(vy, vdat, vx));
					}else{
						vy = simd::loadu(y+i+k);
						vx = simd::loadu(x+col[i*ld+j]+k);
						simd::storeu(y+i+k, simd::fmadd(vy, vdat, vx));
					}
				}
			}
			for(int k = 0 ; k < blockSize ; ++k)
					F.reduce (y[i+k]);
		}
#else
		index_t block = (ld)/kmax ; // use DIVIDE_INTO from fspmvgpu
		for (size_t i = 0 ; i < m ; ++i) {
			index_t j = 0;
			index_t j_loc = 0 ;
			for (size_t l = 0 ; l < block ; ++l) {
				j_loc += kmax ;
				for ( ; j < j_loc ; ++j ) {
					for(int k = 0 ; k < blockSize ; ++k)
						y[i+k] += dat[i*ld+j] * x[col[i*ld+j]+k];
				}
				for(int k = 0 ; k < blockSize ; ++k)
					F.reduce (y[i+k]);
			}
			for ( ; j < ld  ; ++j) {
				for(int k = 0 ; k < blockSize ; ++k)
						y[i+k] += dat[i*ld+j] * x[col[i*ld+j]+k];
			}
			for(int k = 0 ; k < blockSize ; ++k)
					F.reduce (y[i+k]);
		}
#endif // __FFLASFFPACK_USE_SIMD
	}

	// delayed by kmax, row major, lda != n
	template<class Field, bool Aligned>
	inline void fspmm(
			     const Field& F,
			     const size_t m,
			     const size_t n,
			     const size_t ld,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const int lda,
			     const int blockSize,
			     const typename Field::Element_ptr x ,
			     typename Field::Element_ptr y,
			     const index_t & kmax
			    )
	{
#if defined(__FFLASFFPACK_USE_SIMD)
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;

		vect_t vdat, vx, vy;

		index_t block = (ld)/kmax ; // use DIVIDE_INTO from fspmvgpu
		for (size_t i = 0 ; i < m ; ++i) {
			index_t j = 0;
			index_t j_loc = 0 ;
			for (size_t l = 0 ; l < block ; ++l) {
				j_loc += kmax ;
				for ( ; j < j_loc ; ++j ) {
					int k = 0;
					vdat = simd::set1(dat+i*ld+j);
					for(; k < blockSize ; k += simd::vect_size){
						if(Aligned){
							vy = simd::load(y+i*lda+k);
							vx = simd::load(x+col[i*ld+j]*lda+k);
							simd::store(y+i*lda+k, simd::fmadd(vy, vdat, vx));
						}else{
							vy = simd::loadu(y+i*lda+k);
							vx = simd::loadu(x+col[i*ld+j]*lda+k);
							simd::storeu(y+i*lda+k, simd::fmadd(vy, vdat, vx));
						}
					}
					for(; k < blockSize ; ++k)
						y[i*lda+k] += dat[i*ld+j] * x[col[i*ld+j]*lda+k];
				}
				for(int k = 0 ; k < blockSize ; ++k)
					F.reduce (y[i*lda+k]);
			}
			for ( ; j < ld  ; ++j) {
				int k = 0;
				vdat = simd::set1(dat[i*ld+j]);
				for( ; k < blockSize ; ++k){
					if(Aligned){
						vy = simd::load(y+i*lda+k);
						vx = simd::load(x+col[i*ld+j]*lda+k);
						simd::store(y+i*lda+k, simd::fmadd(vy, vdat, vx));
					}else{
						vy = simd::loadu(y+i*lda+k);
						vx = simd::loadu(x+col[i*ld+j]*lda+k);
						simd::storeu(y+i*lda+k, simd::fmadd(vy, vdat, vx));
					}
				}
			}
			for(int k = 0 ; k < blockSize ; ++k)
					F.reduce(y[i*lda+k]);
		}
#else
		index_t block = (ld)/kmax ; // use DIVIDE_INTO from fspmvgpu
		for (size_t i = 0 ; i < m ; ++i) {
			index_t j = 0;
			index_t j_loc = 0 ;
			for (size_t l = 0 ; l < block ; ++l) {
				j_loc += kmax ;
				for ( ; j < j_loc ; ++j ) {
					for(int k = 0 ; k < blockSize ; ++k)
						y[i*lda+k] += dat[i*ld+j] * x[col[i*ld+j]*lda+k];
				}
				for(int k = 0 ; k < blockSize ; ++k)
					F.reduce (y[i*lda+k]);
			}
			for ( ; j < ld  ; ++j) {
				for(int k = 0 ; k < blockSize ; ++k)
						y[i*lda+k] += dat[i*ld+j] * x[col[i*ld+j]*lda+k];
			}
			for(int k = 0 ; k < blockSize ; ++k)
					F.reduce (y[i*lda+k]);
		}
#endif // __FFLASFFPACK_USE_SIMD
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
		details::init_y(F, y.m, b, y.dat, typename FieldTraits<Field>::category());
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
		freduce (F,A.m,y.dat,1);
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

	template<class Field>
	inline void fspmm(
			     const Field& F,
			     const ELL_sub<Field> & A,
			     const int lda,
			     const int blockSize,
			     const VECT<Field> & x,
			     const typename Field::Element & b,
			     VECT<Field> & y
			    )
	{
		details::init_y(F, lda, y.m, b, y.dat, typename FieldTraits<Field>::category());
		fspmm( F, A, lda, blockSize, x, y, typename FieldTraits<Field>::category());
	}

	template<class Field>
	inline void fspmm(const Field & F,
			     const ELL_sub<Field> & A,
			     const int lda,
			     const int blockSize,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::ModularTag)
	{
		using simd = Simd<typename Field::Element>;
		if(A.n == lda){
			if(blockSize % simd::vect_size == 0){
				ell_details::fspmm<Field, true>(F,A.m,A.n,A.ld,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}else{
				ell_details::fspmm<Field, false>(F,A.m,A.n,A.ld,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}
		}else{
			if(blockSize % simd::vect_size == 0){
				ell_details::fspmm<Field, true>(F,A.m,A.n,A.ld,A.col,A.dat,lda, blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}else{
				ell_details::fspmm<Field, false>(F,A.m,A.n,A.ld,A.col,A.dat,lda, blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}
		}
		freduce (F,blockSize,A.m,y.dat,lda);
	}

	template<class Field>
	inline void fspmm(const Field & F,
			     const ELL_sub<Field> & A,
			     const int lda,
			     const int blockSize,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::UnparametricTag)
	{
		using simd = Simd<typename Field::Element>;
		if(A.n == lda){
			if(blockSize % simd::vect_size == 0){
				ell_details::fspmm<Field, true>(F,A.m,A.n,A.ld,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}else{
				ell_details::fspmm<Field, false>(F,A.m,A.n,A.ld,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}
		}else{
			if(blockSize % simd::vect_size == 0){
				ell_details::fspmm<Field, true>(F,A.m,A.n,A.ld,A.col,A.dat,lda, blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}else{
				ell_details::fspmm<Field, false>(F,A.m,A.n,A.ld,A.col,A.dat,lda, blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}
		}
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const ELL_sub<Field> & A,
			     const int lda,
			     const int blockSize,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::GenericTag)
	{
		if(A.n == lda){
			ell_details::fspmm(F,A.m,A.n,A.ld,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::GenericTag ());
		}else{
			ell_details::fspmm(F,A.m,A.n,A.ld,A.col,A.dat,lda,blockSize,x.dat,y.dat, FieldCategories::GenericTag ());
		}
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
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::category());
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


	template<class Field>
	void fspmm(
		      const Field& F,
		      const ELL<Field> & A,
		      const int lda,
		      const int blockSize,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     )
	{
		details::init_y(F, lda, y.m, b, y.dat,  typename FieldTraits<Field>::category());
		fspmv(F,A,lda, blockSize,x,y,typename FieldTraits<Field>::category());
	}

	template<class Field>
	void fspmm(
		      const Field& F,
		      const ELL<Field> & A,
		      const int lda,
		      const int blockSize,
		      const VECT<Field> & x,
		      VECT<Field> & y,
		      FieldCategories::GenericTag
		     )
	{
		if(A.n == lda){
			if(A.n == lda){
				ell_details::fspmm(F,A.m,A.n,A.ld,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::GenericTag());
			}else{
				ell_details::fspmm(F,A.m,A.n,A.ld,A.col,A.dat,lda,blockSize,x.dat,y.dat, FieldCategories::GenericTag());
			}
		}
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const ELL<Field> & A,
		      const int lda,
		      const int blockSize,
		      const VECT<Field> & x,
		      VECT<Field> & y,
		      FieldCategories::UnparametricTag
		     )
	{
		using simd = Simd<typename Field::Element>;
		if(A.n == lda){
			if(blockSize % simd::vect_size == 0){
				ell_details::fspmm<Field, true>(F,A.m,A.n,A.ld,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}else{
				ell_details::fspmm<Field, false>(F,A.m,A.n,A.ld,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}
		}else{
			if(blockSize % simd::vect_size == 0){
				ell_details::fspmm<Field, true>(F,A.m,A.n,A.ld,A.col,A.dat,lda, blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}else{
				ell_details::fspmm<Field, false>(F,A.m,A.n,A.ld,A.col,A.dat,lda, blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}
		}
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const ELL<Field> & A,
		      const int lda,
		      const int blockSize,
		      const VECT<Field> & x,
		      VECT<Field> & y,
		      FieldCategories::ModularTag
		     )
	{
		size_t kmax = Protected::DotProdBoundClassic(F,F.one);
		using simd = Simd<typename Field::Element>;
		if(A.n == lda){
			if(blockSize % simd::vect_size == 0){
				ell_details::fspmm<Field, true>(F,A.m,A.n,A.ld,A.col,A.dat,blockSize,x.dat,y.dat,kmax);
			}else{
				ell_details::fspmm<Field, false>(F,A.m,A.n,A.ld,A.col,A.dat,blockSize,x.dat,y.dat,kmax);
			}
		}else{
			if(blockSize % simd::vect_size == 0){
				ell_details::fspmm<Field, true>(F,A.m,A.n,A.ld,A.col,A.dat,lda, blockSize,x.dat,y.dat,kmax);
			}else{
				ell_details::fspmm<Field, false>(F,A.m,A.n,A.ld,A.col,A.dat,lda, blockSize,x.dat,y.dat,kmax);
			}
		}
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
	inline void fspmm_zo(
				const Field & F,
				const size_t m,
				const size_t n,
				const size_t ld,
				const index_t * col,
				const int blockSize,
				const typename Field::Element_ptr x ,
				typename Field::Element_ptr y,
				FieldCategories::GenericTag
			    )
	{
		if(add){
			for (size_t i = 0 ; i < m ; ++i) {
				for (index_t j = 0 ; j < ld ; ++j) {
					for(int k = 0 ; k < blockSize ; ++k)
						F.addin(y[i+k],x[col[i*ld+j]+k]);
				}
			}
		}else{
			for (size_t i = 0 ; i < m ; ++i) {
				for (index_t j = 0 ; j < ld ; ++j) {
					for(int k = 0 ; k < blockSize ; ++k)
						F.subin(y[i+k],x[col[i*ld+j]+k]);
				}
			}
		}
	}

	template<class Field, bool add>
	inline void fspmm_zo(
				const Field & F,
				const size_t m,
				const size_t n,
				const size_t ld,
				const index_t * col,
				const int lda,
				const int blockSize,
				const typename Field::Element_ptr x ,
				typename Field::Element_ptr y,
				FieldCategories::GenericTag
			    )
	{
		if(add){
			for (size_t i = 0 ; i < m ; ++i) {
				for (index_t j = 0 ; j < ld ; ++j) {
					for(int k = 0 ; k < blockSize ; ++k)
						F.addin(y[i*lda+k],x[col[i*ld+j]*lda+k]);
				}
			}
		}else{
			for (size_t i = 0 ; i < m ; ++i) {
				for (index_t j = 0 ; j < ld ; ++j) {
					for(int k = 0 ; k < blockSize ; ++k)
						F.subin(y[i*lda+k],x[col[i*ld+j]*lda+k]);
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

	template<class Field, bool add, bool Aligned>
	inline void fspmm_zo(
				const Field & F,
				const size_t m,
				const size_t n,
				const size_t ld,
				const index_t * col,
				const int blockSize,
				const typename Field::Element_ptr x ,
				typename Field::Element_ptr y,
				FieldCategories::UnparametricTag
			    )
	{
#if defined(__FFLASFFPACK_USE_SIMD)
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;

		vect_t vy, vx;
		if(add){
			for(size_t i = 0 ; i < m ; ++i){
				for(size_t j = 0 ; j < ld ; ++j){
					int k = 0;
					for(; k < blockSize ; k+=simd::vect_size){
						if(Aligned){
							simd::store(y+i+k, simd::add(simd::load(y+i+k), simd::load(x+col[i*ld+j]+k)));
						}else{
							simd::storeu(y+i+k, simd::add(simd::loadu(y+i+k), simd::loadu(x+col[i*ld+j]+k)));
						}
					}
					for(; k < blockSize ; ++k){
						y[i+k] += x[col[i*ld+j]+k];
					}
				}
			}
		}else{
			for(size_t i = 0 ; i < m ; ++i){
				for(size_t j = 0 ; j < ld ; ++j){
					int k = 0;
					for(; k < blockSize ; k+=simd::vect_size){
						if(Aligned){
							simd::store(y+i+k, simd::sub(simd::load(y+i+k), simd::load(x+col[i*ld+j]+k)));
						}else{
							simd::storeu(y+i+k, simd::sub(simd::loadu(y+i+k), simd::loadu(x+col[i*ld+j]+k)));
						}
					}
					for(; k < blockSize ; ++k){
						y[i+k] -= x[col[i*ld+j]+k];
					}
				}
			}
		}
#else
		if(add){
			for ( size_t i = 0 ;  i < m   ; ++i ) {
				for (index_t j = 0 ; j < ld ; ++j) {
					for(int k = 0 ; k < blockSize ; ++k)
						y[i+k] += x[col[i*ld+j]];
				}
			}
		}else{
			for ( size_t i = 0 ;  i < m   ; ++i ) {
				for (index_t j = 0 ; j < ld ; ++j) {
					for(int k = 0 ; k < blockSize ; ++k)
						y[i+k] -= x[col[i*ld+j]];
				}
			}
		}
#endif
	}

	template<class Field, bool add, bool Aligned>
	inline void fspmm_zo(
				const Field & F,
				const size_t m,
				const size_t n,
				const size_t ld,
				const index_t * col,
				const int lda,
				const int blockSize,
				const typename Field::Element_ptr x ,
				typename Field::Element_ptr y,
				FieldCategories::UnparametricTag
			    )
	{
#if defined(__FFLASFFPACK_USE_SIMD)
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;

		vect_t vy, vx;
		if(add){
			for(size_t i = 0 ; i < m ; ++i){
				for(size_t j = 0 ; j < ld ; ++j){
					int k = 0;
					for(; k < blockSize ; k+=simd::vect_size){
						if(Aligned){
							simd::store(y+i*lda+k, simd::add(simd::load(y+i*lda+k), simd::load(x+col[i*ld+j]*lda+k)));
						}else{
							simd::storeu(y+i*lda+k, simd::add(simd::loadu(y+i*lda+k), simd::loadu(x+col[i*ld+j]*lda+k)));
						}
					}
					for(; k < blockSize ; ++k){
						y[i*lda+k] += x[col[i*ld+j]*lda+k];
					}
				}
			}
		}else{
			for(size_t i = 0 ; i < m ; ++i){
				for(size_t j = 0 ; j < ld ; ++j){
					int k = 0;
					for(; k < blockSize ; k+=simd::vect_size){
						if(Aligned){
							simd::store(y+i*lda+k, simd::sub(simd::load(y+i*lda+k), simd::load(x+col[i*ld+j]*lda+k)));
						}else{
							simd::storeu(y+i*lda+k, simd::sub(simd::loadu(y+i*lda+k), simd::loadu(x+col[i*ld+j]*lda+k)));
						}
					}
					for(; k < blockSize ; ++k){
						y[i*lda+k] -= x[col[i*ld+j]*lda+k];
					}
				}
			}
		}
#else
		if(add){
			for ( size_t i = 0 ;  i < m   ; ++i ) {
				for (index_t j = 0 ; j < ld ; ++j) {
					for(int k = 0 ; k < blockSize ; ++k)
						y[i*lda+k] += x[col[i*ld+j]]*lda;
				}
			}
		}else{
			for ( size_t i = 0 ;  i < m   ; ++i ) {
				for (index_t j = 0 ; j < ld ; ++j) {
					for(int k = 0 ; k < blockSize ; ++k)
						y[i*lda+k] -= x[col[i*ld+j]]*lda;
				}
			}
		}
#endif
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
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::category());
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
		freduce (F,A.m,y.dat,1);
	}

	template<class Field>
	inline void fspmm(
			     const Field & F,
			     ELL_ZO<Field> & A,
			     const int lda,
			     const int blockSize,
			     const VECT<Field> & x,
			     const typename Field::Element & b,
			     VECT<Field> & y
			    )
	{
		details::init_y(F, lda, A.m*blockSize, b, y.dat,  typename FieldTraits<Field>::category());
		fspmm(F, A, lda, blockSize, x, y, typename FieldTraits<Field>::category() );
	}

	template<class Field>
	inline void fspmm(
			     const Field & F,
			     ELL_ZO<Field> & A,
			     const int lda,
			     const int blockSize,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::GenericTag
			    )
	{
		if(A.n == lda){
			if (A.cst == F.one) {
				ell_details::fspmm_zo<Field,true>(F,A.m,A.n,A.ld,A.col,blockSize,x.dat,y.dat, FieldCategories::GenericTag());
			}
			else if (A.cst == F.mOne) {
				ell_details::fspmm_zo<Field,false>(F,A.m,A.n,A.ld,A.col,blockSize,x.dat,y.dat, FieldCategories::GenericTag());
			}
			else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, lda, x1, 1);
				ell_details::fspmm_zo<Field,true>(F,A.m,A.n,A.ld,A.col,blockSize,x1,y.dat, FieldCategories::GenericTag());
				fflas_delete(x1);
			}
		}else{
			if (A.cst == F.one) {
				ell_details::fspmm_zo<Field,true>(F,A.m,A.n,A.ld,A.col,lda,blockSize,x.dat,y.dat, FieldCategories::GenericTag());
			}
			else if (A.cst == F.mOne) {
				ell_details::fspmm_zo<Field,false>(F,A.m,A.n,A.ld,A.col,lda,blockSize,x.dat,y.dat, FieldCategories::GenericTag());
			}
			else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, lda, x1, 1);
				ell_details::fspmm_zo<Field,true>(F,A.m,A.n,A.ld,A.col,lda,blockSize,x1,y.dat, FieldCategories::GenericTag());
				fflas_delete(x1);
			}
		}

	}

	template<class Field>
	inline void fspmm(
			     const Field & F,
			     ELL_ZO<Field> & A,
			     const int lda,
			     const int blockSize,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::UnparametricTag
			    )
	{
		using simd=Simd<typename Field::Element>;
		if(A.n == lda){
			if (A.cst == F.one){
				if(blockSize % simd::vect_size == 0)
					ell_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.ld,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					ell_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.ld,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}else if (A.cst == F.mOne) {
				if(blockSize % simd::vect_size == 0)
					ell_details::fspmm_zo<Field, false, true>(F,A.m,A.n,A.ld,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					ell_details::fspmm_zo<Field, false, false>(F,A.m,A.n,A.ld,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, lda, x1, 1);
				if(blockSize % simd::vect_size == 0)
					ell_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.ld,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					ell_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.ld,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				fflas_delete(x1);
			}
		}else{
			if (A.cst == F.one){
				if(blockSize % simd::vect_size == 0)
					ell_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.ld,A.col,lda,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					ell_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.ld,A.col,lda,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}
			else if (A.cst == F.mOne) {
				if(blockSize % simd::vect_size == 0)
					ell_details::fspmm_zo<Field, false, true>(F,A.m,A.n,A.ld,A.col,lda,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					ell_details::fspmm_zo<Field, false, false>(F,A.m,A.n,A.ld,A.col,lda,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}
			else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, lda, x1, 1);
				if(blockSize % simd::vect_size == 0)
					ell_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.ld,A.col,lda,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					ell_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.ld,A.col,lda,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				fflas_delete(x1);
			}
		}
	}

	template<class Field>
	inline void fspmm(
			     const Field & F,
			     ELL_ZO<Field> & A,
			     const int lda,
			     const int blockSize,
			     const VECT<Field> & x,
			     VECT<Field> & y
			     , FieldCategories::ModularTag
			    )
	{
		using simd=Simd<typename Field::Element>;
		if(A.n == lda){
			if (A.cst == F.one){
				if(blockSize % simd::vect_size == 0)
					ell_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.ld,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					ell_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.ld,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}else if (A.cst == F.mOne) {
				if(blockSize % simd::vect_size == 0)
					ell_details::fspmm_zo<Field, false, true>(F,A.m,A.n,A.ld,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					ell_details::fspmm_zo<Field, false, false>(F,A.m,A.n,A.ld,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, lda, x1, 1);
				if(blockSize % simd::vect_size == 0)
					ell_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.ld,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					ell_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.ld,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				fflas_delete(x1);
			}
		}else{
			if (A.cst == F.one){
				if(blockSize % simd::vect_size == 0)
					ell_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.ld,A.col,lda,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					ell_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.ld,A.col,lda,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}
			else if (A.cst == F.mOne) {
				if(blockSize % simd::vect_size == 0)
					ell_details::fspmm_zo<Field, false, true>(F,A.m,A.n,A.ld,A.col,lda,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					ell_details::fspmm_zo<Field, false, false>(F,A.m,A.n,A.ld,A.col,lda,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}
			else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, lda, x1, 1);
				if(blockSize % simd::vect_size == 0)
					ell_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.ld,A.col,lda,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					ell_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.ld,A.col,lda,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				fflas_delete(x1);
			}
		}
		freduce (F,blockSize,A.m,y.dat,lda);
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
		ELL_col = fflas_new<index_t >(ld*ELL_m, Alignment::CACHE_LINE);
		if(!ZO){
			ELL_dat = fflas_new<typename Field::Element >(ld*ELL_m, Alignment::CACHE_LINE);
		}
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
		sp_ell_from_csr(F, COO_m, COO_n, COO_nnz, COO_col, row, COO_dat, ELL_m, ELL_n, ld, ELL_col, ELL_dat, ZO);
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
