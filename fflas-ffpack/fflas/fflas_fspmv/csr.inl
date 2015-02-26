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

/** @file fflas/fflas_fspmv_csr.inl
 * NO DOC
 */

#ifndef __FFLASFFPACK_fflas_fflas_spmv_csr_INL
#define __FFLASFFPACK_fflas_fflas_spmv_csr_INL

namespace FFLAS { namespace csr_details {

	// y = A x + b y ; (generic)
	template<class Field>
	inline void fspmv(
			     const Field& F,
			     const index_t m,
			     const index_t n,
			     const index_t * st,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const typename Field::Element_ptr x ,
			     typename Field::Element_ptr y
			     , FieldCategories::GenericTag
			    )
	{
		for (index_t i = 0 ; i < m ; ++i)
			for (index_t j = st[i] ; j < st[i+1] ; ++j)
				F.axpyin(y[i],dat[j],x[col[j]]);
	}

	// row major : lda == n
	template<class Field>
	inline void fspmm(
			     const Field& F,
			     const index_t m,
			     const index_t n,
			     const index_t * st,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const int blockSize,
			     const typename Field::Element_ptr x ,
			     typename Field::Element_ptr y
			     , FieldCategories::GenericTag
			    )
	{
		for (index_t i = 0 ; i < m ; ++i)
		{
			for (index_t j = st[i] ; j < st[i+1] ; ++j)
			{
				for(int k = 0 ; k < blockSize ; ++k)
				{
					F.axpyin(y[i*blockSize+k],dat[j],x[col[j]*blockSize+k]);
				}
			}
		}
	}

	// lda != n
	template<class Field>
	inline void fspmm(
			     const Field& F,
			     const index_t m,
			     const index_t n,
			     const index_t * st,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const int blockSize,
			     const typename Field::Element_ptr x ,
			     const int ldx,
			     typename Field::Element_ptr y,
			     const int ldy,
			     FieldCategories::GenericTag
			    )
	{
		for (index_t i = 0 ; i < m ; ++i)
		{
			for (index_t j = st[i] ; j < st[i+1] ; ++j)
			{
				for(int k = 0 ; k < blockSize ; ++k)
				{
					F.axpyin(y[i*ldy*blockSize+k],dat[j],x[col[j]*ldx*blockSize+k]);
				}
			}
		}
	}

	template<class Field>
	void fspmv(
		      const Field & F,
		      const index_t m,
		      const index_t n,
		      const index_t * st,
		      const index_t * col,
		      const typename Field::Element_ptr dat,
		      const typename Field::Element_ptr x ,
		      typename Field::Element_ptr y
		      , FieldCategories::UnparametricTag
		     )
	{
		for (index_t i = 0 ; i < m ; ++i)
			for (index_t j = st[i] ; j < st[i+1] ; ++j)
				y[i] += dat[j] * x[col[j]];
	}

	// Double
	template<>
	void fspmv(
		      const Givaro::DoubleDomain & F,
		      const index_t m,
		      const index_t n,
		      const index_t * st,
		      const index_t * col,
		      const typename Givaro::DoubleDomain::Element_ptr  dat,
		      const typename Givaro::DoubleDomain::Element_ptr x ,
		      typename Givaro::DoubleDomain::Element_ptr y
		      , FieldCategories::UnparametricTag
		     )
	{
#ifdef __FFLASFFPACK_HAVE_MKL
		// char * transa = (ta==FflasNoTrans)?'n':'t';
		char   transa = 'N';
		index_t m_ = (index_t) m ;
		double * yd =  FFLAS::fflas_new<double >(m);
		// mkl_dcsrgemv (bug, not zero based)
		mkl_cspblas_dcsrgemv
		(&transa, &m_, const_cast<double*>(dat), const_cast<index_t*>(st) , const_cast<index_t*>(col), const_cast<double*>(x), yd);
		faddin(F,m,yd,1,y,1);
		FFLAS::fflas_delete( yd );
#else
		for (index_t i = 0 ; i < m ; ++i)
			for (index_t j = st[i] ; j < st[i+1] ; ++j)
				y[i] += dat[j] * x[col[j]];
#endif // __FFLASFFPACK_HAVE_MKL
	}


	// Float
	template<>
	void fspmv(
		      const Givaro::FloatDomain & F,
		      const index_t m,
		      const index_t n,
		      const index_t * st,
		      const index_t * col,
		      const typename Givaro::FloatDomain::Element_ptr dat,
		      const typename Givaro::FloatDomain::Element_ptr x ,
		      typename Givaro::FloatDomain::Element_ptr y
		      , FieldCategories::UnparametricTag
		     )
	{
#ifdef __FFLASFFPACK_HAVE_MKL
		// char * transa = (ta==FflasNoTrans)?'n':'t';
		char   transa = 'n';
		index_t m_ = (index_t) m ;
		float * yd = FFLAS::fflas_new<float >(m);
		fscalin(F,m,n,y,1);
		// mkl_scsrgemv
		mkl_cspblas_scsrgemv
		(&transa, &m_, const_cast<float*>(dat), const_cast<index_t*>(st), const_cast<index_t*>(col), const_cast<float*>(x), yd);
		faddin(F,m,yd,1,y,1);
		FFLAS::fflas_delete( yd );
#else
		for (index_t i = 0 ; i < m ; ++i)
			for (index_t j = st[i] ; j < st[i+1] ; ++j)
				y[i] += dat[j] * x[col[j]];
#endif // __FFLASFFPACK_HAVE_MKL
	}

	// row major lda = n
	template<class Field, bool Aligned>
	inline void fspmm(
			     const Field& F,
			     const index_t m,
			     const index_t n,
			     const index_t * st,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const int blockSize,
			     const typename Field::Element_ptr x ,
			     typename Field::Element_ptr y
			     , FieldCategories::UnparametricTag
			    )
	{
#if defined(__FFLASFFPACK_USE_SIMD)
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;

		int end = (blockSize/simd::vect_size)*simd::vect_size;
		vect_t vdat, vx, vy;
		if(Aligned){
			for(index_t i = 0 ; i < m ; ++i){
				for(index_t j = st[i] ; j < st[i+1] ; ++j){
					int k = 0;
					vdat = simd::set1(dat+j);
					for(; k < end ; k+=simd::vect_size){
						vy = simd::load(y+i*blockSize+k);
						vx = simd::load(x+col[j]*blockSize+k);
						simd::store(y+i*blockSize+k, simd::fmadd(vy, vdat, vx));
					}
					for(; k < blockSize ; ++k){
						y[i*blockSize+k] += dat[j]*x[col[j]*blockSize+k];
					}
				}
			}
		}else{
			for(index_t i = 0 ; i < m ; ++i){
				for(index_t j = st[i] ; j < st[i+1] ; ++j){
					int k = 0;
					vdat = simd::set1(dat+j);
					for(; k < end ; k+=simd::vect_size){
						vy = simd::loadu(y+i*blockSize+k);
						vx = simd::loadu(x+col[j]*blockSize+k);
						vdat = simd::set1(dat+j);
						simd::storeu(y+i*blockSize+k, simd::fmadd(vy, vdat, vx));
					}
					for(; k < blockSize ; ++k){
						y[i*blockSize+k] += dat[j]*x[col[j]*blockSize+k];
					}
				}
			}
		}
#else
		for (index_t i = 0 ; i < m ; ++i)
		{
			for (index_t j = st[i] ; j < st[i+1] ; ++j)
			{
				for(int k = 0 ; k < blockSize ; ++k)
				{
					y[i*blockSize+k] += dat[j]*x[col[j]*blockSize+k];
				}
			}
		}
#endif // __FFLASFFPACK_USE_SIMD
	}

	// row major lda != n
	template<class Field, bool Aligned>
	inline void fspmm(
			     const Field& F,
			     const index_t m,
			     const index_t n,
			     const index_t * st,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const int blockSize,
			     const typename Field::Element_ptr x,
			     const int ldx,
			     typename Field::Element_ptr y,
			     const int ldy,
			     FieldCategories::UnparametricTag
			    )
	{
#if defined(__FFLASFFPACK_USE_SIMD)
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;

		int end = (blockSize/simd::vect_size)*simd::vect_size;

		vect_t vy, vx, vdat;
		if(Aligned){
			for(index_t i = 0 ; i < m ; ++i){
				for(index_t j = st[i] ; j < st[i+1] ; ++j){
					int k = 0;
					vdat = simd::set1(dat+j);
					for(; k < end ; k+=simd::vect_size){
						vy = simd::load(y+i*ldy*blockSize+k);
						vx = simd::load(x+col[j]*ldx*blockSize+k);
						simd::store(y+i*ldy*blockSize+k, simd::fmadd(vy, vdat, vx));
					}
					for(; k < blockSize ; ++k){
						y[i*ldy*blockSize+k] += dat[j]*x[col[j]*ldx*blockSize+k];
					}
				}
			}
		}else{
			for(index_t i = 0 ; i < m ; ++i){
				for(index_t j = st[i] ; j < st[i+1] ; ++j){
					int k = 0;
					vdat = simd::set1(dat+j);
					for(; k < end ; k+=simd::vect_size){
						vy = simd::loadu(y+i*ldy*blockSize+k);
						vx = simd::loadu(x+col[j]*ldx*blockSize+k);
						vdat = simd::set1(dat+j);
						simd::storeu(y+i*ldy*blockSize+k, simd::fmadd(vy, vdat, vx));
					}
					for(; k < blockSize ; ++k){
						y[i*ldy*blockSize+k] += dat[j]*x[col[j]*ldx*blockSize+k];
					}
				}
			}
		}
#else
		for (index_t i = 0 ; i < m ; ++i)
		{
			for (index_t j = st[i] ; j < st[i+1] ; ++j)
			{
				for(int k = 0 ; k < blockSize ; ++k)
				{
					y[i*ldy*blockSize+k] += dat[j]*x[col[j]*ldx*blockSize+k];
				}
			}
		}
#endif // __FFLASFFPACK_USE_SIMD
	}

	// delayed by kmax
	template<class Field>
	inline void fspmv(
			     const Field& F,
			     const index_t m,
			     const index_t n,
			     const index_t * st,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const typename Field::Element_ptr x ,
			     typename Field::Element * y,
			     const index_t & kmax
			    )
	{
		for (index_t i = 0 ; i < m ; ++i) {
			index_t j = st[i];
			index_t j_loc = j;
			index_t j_end = st[i+1];
			index_t block = (j_end - j_loc)/kmax ;
			for (index_t l = 0 ; l < (index_t) block ; ++l) {
				j_loc += kmax ;
				for ( ; j < j_loc ; ++j) {
					y[i] += dat[j] * x[col[j]];
				}
				F.reduce (y[i]);
			}
			for ( ; j < j_end ; ++j) {
				y[i] += dat[j] * x[col[j]];
			}
			F.reduce (y[i]);
		}
	}

	// row major lda = n
	template<class Field, bool Aligned>
	inline void fspmm(
			     const Field& F,
			     const index_t m,
			     const index_t n,
			     const index_t * st,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const int blockSize,
			     const typename Field::Element_ptr x ,
			     typename Field::Element * y,
			     const index_t & kmax
			    )
	{
#if defined(__FFLASFFPACK_USE_SIMD)
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vec_t;

		int end = (blockSize/simd::vect_size)*simd::vect_size;
		vect_t vx, vy, vdat;
		for(index_t i = 0 ; i < m ; ++i){
			index_t j = st[i];
			index_t j_loc = j;
			index_t j_end = st[i+1];
			index_t block = (j_end - j_loc)/kmax ;
			for (index_t l = 0 ; l < (index_t) block ; ++l) {
				j_loc += kmax ;
				for ( ; j < j_loc ; ++j) {
					int k = 0;
					vdat = simd::set1(dat+j);
					for(; k < end ; k+=simd::vect_size){
						if(Aligned){
							vy = simd::load(y+i*blockSize+k);
							vx = simd::load(x+col[j]*blockSize+k);
							simd::store(y+i*blockSize+k, simd::fmadd(vy, vdat, vx));
						}else{
							vy = simd::loadu(y+i*blockSize+k);
							vx = simd::loadu(x+col[j]*blockSize+k);
							simd::storeu(y+i*blockSize+k, simd::fmadd(vy, vdat, vx));
						}
					}
					for(; k < blockSize ; ++k)
					{
						y[i*blockSize+k] += dat[j]*x[col[j]*blockSize+k];
					}
				}
				freduce (F, blockSize, y+i*blockSize, 1);
			}
			for ( ; j < j_end ; ++j) {
				int k = 0;
				for(; k < end ; k+=simd::vect_size){
					vdat = simd::set1(dat+j);
					if(Aligned){
						vy = simd::load(y+i*blockSize+k);
						vx = simd::load(x+col[j]*blockSize+k);
						simd::store(y+i*blockSize+k, simd::fmadd(vy, vdat, vx));
					}else{
						vy = simd::loadu(y+i*blockSize+k);
						vx = simd::loadu(x+col[j]*blockSize+k);
						simd::storeu(y+i*blockSize+k, simd::fmadd(vy, vdat, vx));
					}
				}
				for(; k < blockSize ; ++k)
				{
					y[i*blockSize+k] += dat[j]*x[col[j]*blockSize+k];
				}
			}
			freduce (F, blockSize, y+i*blockSize, 1);
		}
#else
		for (index_t i = 0 ; i < m ; ++i) {
			index_t j = st[i];
			index_t j_loc = j;
			index_t j_end = st[i+1];
			index_t block = (j_end - j_loc)/kmax ;
			for (index_t l = 0 ; l < (index_t) block ; ++l) {
				j_loc += kmax ;
				for ( ; j < j_loc ; ++j) {
					for(int k = 0 ; k < blockSize ; ++k){
						y[i*blockSize+k] += dat[j] * x[col[j]*blockSize+k];
					}
				}
				freduce (F, blockSize, y+i*blockSize, 1);
			}
			for ( ; j < j_end ; ++j) {
				for(int k = 0 ; k < blockSize ; ++k){
					y[i*blockSize+k] += dat[j] * x[col[j]*blockSize+k];
				}
			}
			freduce (F, blockSize, y+i*blockSize, 1);
		}
#endif // __FFLASFFPACK_USE_SIMD
	}

	// row major lda != n
	template<class Field, bool Aligned>
	inline void fspmm(
			     const Field& F,
			     const index_t m,
			     const index_t n,
			     const index_t * st,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const int blockSize,
			     const typename Field::Element_ptr x ,
			     const int ldx,
			     typename Field::Element * y,
			     const int ldy,
			     const index_t & kmax
			    )
	{
#if defined(__FFLASFFPACK_USE_SIMD)
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vec_t;

		int end = (blockSize/simd::vect_size)*simd::vect_size;
		vect_t vx, vy, vdat;
		for(index_t i = 0 ; i < m ; ++i){
			index_t j = st[i];
			index_t j_loc = j;
			index_t j_end = st[i+1];
			index_t block = (j_end - j_loc)/kmax ;
			for (index_t l = 0 ; l < (index_t) block ; ++l) {
				j_loc += kmax ;
				for ( ; j < j_loc ; ++j) {
					int k = 0;
					vdat = simd::set1(dat+j);
					for(; k < end ; k+=simd::vect_size){
						if(Aligned){
							vy = simd::load(y+i*blockSize*ldy+k);
							vx = simd::load(x+col[j]*blockSize*ldx+k);
							simd::store(y+i*blockSize*ldy+k, simd::fmadd(vy, vdat, vx));
						}else{
							vy = simd::loadu(y+i*blockSize*ldy+k);
							vx = simd::loadu(x+col[j]*blockSize*ldx+k);
							simd::storeu(y+i*blockSize*ldy+k, simd::fmadd(vy, vdat, vx));
						}
					}
					for(; k < blockSize ; ++k)
					{
						y[i*blockSize*ldy+k] += dat[j]*x[col[j]*blockSize*ldx+k];
					}
				}
				freduce (F, blockSize, y+i*blockSize*ldy, 1);
			}
			for ( ; j < j_end ; ++j) {
				int k = 0;
				for(; k < end ; k+=simd::vect_size){
					vdat = simd::set1(dat+j);
					if(Aligned){
						vy = simd::load(y+i*blockSize*ldy+k);
						vx = simd::load(x+col[j]*blockSize*ldx+k);
						simd::store(y+i*blockSize*ldy+k, simd::fmadd(vy, vdat, vx));
					}else{
						vy = simd::loadu(y+i*blockSize*ldy+k);
						vx = simd::loadu(x+col[j]*blockSize*ldx+k);
						simd::storeu(y+i*blockSize*ldy+k, simd::fmadd(vy, vdat, vx));
					}
				}
				for(; k < blockSize ; ++k)
				{
					y[i*blockSize*ldy+k] += dat[j]*x[col[j]*blockSize*ldx+k];
				}
			}
			freduce (F, blockSize, y+i*blockSize*ldy, 1);
		}
#else
		for (index_t i = 0 ; i < m ; ++i) {
			index_t j = st[i];
			index_t j_loc = j;
			index_t j_end = st[i+1];
			index_t block = (j_end - j_loc)/kmax ;
			for (index_t l = 0 ; l < (index_t) block ; ++l) {
				j_loc += kmax ;
				for ( ; j < j_loc ; ++j) {
					for(int k = 0 ; k < blockSize ; ++k){
						y[i*blockSize*ldy+k] += dat[j] * x[col[j]*blockSize*ldx+k];
					}
				}
				freduce (F, blockSize, y+i*blockSize*ldy, 1);
			}
			for ( ; j < j_end ; ++j) {
				for(int k = 0 ; k < blockSize ; ++k){
					y[i*blockSize*ldy+k] += dat[j] * x[col[j]*blockSize*ldx+k];
				}
			}
			freduce (F, blockSize, y+i*blockSize*ldy, 1);
		}
#endif // __FFLASFFPACK_USE_SIMD
	}
} // details
} // FFLAS

namespace FFLAS {

	/* ******* */
	/* CSR_sub */
	/* ******* */

	// y = A x + b y ; (generic)
	// it is supposed that no reduction is needed.
	template<class Field>
	inline void fspmv(
			     const Field& F,
			     const CSR_sub<Field> & A,
			     const VECT<Field> & x,
			     const typename Field::Element & b,
			     VECT<Field> & y
			    )
	{
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::category());
		fspmv( F, A, x, y, typename FieldTraits<Field>::category());
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const CSR_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::ModularTag)
	{
		csr_details::fspmv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag ());
		freduce (F,A.m,y.dat,1);
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const CSR_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::UnparametricTag)
	{
		csr_details::fspmv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag ());
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const CSR_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::GenericTag)
	{
		csr_details::fspmv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat, FieldCategories::GenericTag ());
	}


	template<class Field>
	inline void fspmm(
			     const Field& F,
			     const CSR_sub<Field> & A,
			     const int blockSize,
			     const VECT<Field> & x,
			     const int ldx,
			     const typename Field::Element & b,
			     VECT<Field> & y,
			     const int ldy
			    )
	{
		details::init_y(F, A.m, blockSize, b, y.dat, ldy,  typename FieldTraits<Field>::category());
		fspmm(F, A, blockSize, x, ldx, y, ldy, typename FieldTraits<Field>::category());
	}

	template<class Field>
	inline void fspmm(const Field & F,
			     const CSR_sub<Field> & A,
			     const int blockSize,
			     const VECT<Field> & x,
			     const int ldx,
			     VECT<Field> & y,
			     const int ldy,
			     FieldCategories::ModularTag)
	{
		using simd = Simd<typename Field::Element>;
		if(ldx == A.n && ldy == blockSize){
			if(blockSize % simd::vectize == 0)
				csr_details::fspmm<Field, true>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			else
				csr_details::fspmm<Field, false>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
		}else{
			if(blockSize % simd::vectize == 0)
				csr_details::fspmm<Field, true>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
			else
				csr_details::fspmm<Field, false>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
		}
		freduce (F,blockSize,A.m,y.dat,ldy);
	}

	template<class Field>
	inline void fspmm(const Field & F,
			     const CSR_sub<Field> & A,
			     const int blockSize,
			     const VECT<Field> & x,
			     const int ldx,
			     VECT<Field> & y,
			     const int ldy,
			     FieldCategories::UnparametricTag)
	{
		using simd = Simd<typename Field::Element>;
		if(ldx == A.n && ldy == blockSize){
			if(blockSize % simd::vectize == 0)
				csr_details::fspmm<Field, true>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			else
				csr_details::fspmm<Field, false>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
		}else{
			if(blockSize % simd::vectize == 0)
				csr_details::fspmm<Field, true>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
			else
				csr_details::fspmm<Field, false>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
		}
	}

	template<class Field>
	inline void fspmm(const Field & F,
			     const CSR_sub<Field> & A,
			     const int blockSize,
			     const VECT<Field> & x,
			     const int ldx,
			     VECT<Field> & y,
			     const int ldy,
			     FieldCategories::GenericTag)
	{
		if(ldx == A.n && ldy == blockSize){
			csr_details::fspmm(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::GenericTag ());
		}else{
			csr_details::fspmm(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::GenericTag ());
		}
	}

	/* *** */
	/* CSR */
	/* *** */

	// y = A x + b y ; (generic)
	// reductions are delayed.
	template<class Field>
	void fspmv(
		      const Field& F,
		      const CSR<Field> & A,
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
		      const CSR<Field> & A,
		      const VECT<Field> & x,
		      VECT<Field> & y
		      , FieldCategories::GenericTag
		     )
	{
		csr_details::fspmv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat, FieldCategories::GenericTag());
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const CSR<Field> & A,
		      const VECT<Field> & x,
		      VECT<Field> & y
		      , FieldCategories::UnparametricTag
		     )
	{
		csr_details::fspmv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag());
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const CSR<Field> & A,
		      const VECT<Field> & x,
		      VECT<Field> & y
		      , FieldCategories::ModularTag
		     )
	{
		index_t kmax = Protected::DotProdBoundClassic(F,F.one);
		if(kmax > A.maxrow)
		{
			csr_details::fspmv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag());
			freduce (F, A.m, y.dat, 1);
		}else{
			csr_details::fspmv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat,(index_t) kmax);
		}
	}

	template<class Field>
	inline void fspmm(
			     const Field& F,
			     const CSR<Field> & A,
			     const int blockSize,
			     const VECT<Field> & x,
			     const int ldx,
			     const typename Field::Element & b,
			     const int ldy,
			     VECT<Field> & y
			    )
	{
		details::init_y(F, A.m, blockSize, b, y.dat, ldy,  typename FieldTraits<Field>::category());
		fspmm(F, A, blockSize, x, ldx, y, ldy, typename FieldTraits<Field>::category());
	}

	template<class Field>
	inline void fspmm(const Field & F,
			     const CSR<Field> & A,
			     const int blockSize,
			     const VECT<Field> & x,
			     const int ldx,
			     VECT<Field> & y,
			     const int ldy,
			     FieldCategories::ModularTag)
	{
		using simd = Simd<typename Field::Element>;
		index_t kmax = Protected::DotProdBoundClassic(F,F.one);
		if(ldx == A.n && ldy == blockSize){
			if(blockSize % simd::vect_size == 0)
				csr_details::fspmm<Field, true>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,y.dat,kmax);
			else
				csr_details::fspmm<Field, false>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,y.dat,kmax);
		}else{
			if(blockSize % simd::vect_size == 0)
				csr_details::fspmm<Field, true>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy,kmax);
			else
				csr_details::fspmm<Field, false>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy,kmax);
		}
	}

	template<class Field>
	inline void fspmm(const Field & F,
			     const CSR<Field> & A,
			     const int blockSize,
			     const VECT<Field> & x,
			     const int ldx,
			     VECT<Field> & y,
			     const int ldy,
			     FieldCategories::UnparametricTag)
	{
		using simd = Simd<typename Field::Element>;
		if(ldx == A.n && ldy == blockSize){
			if(blockSize % simd::vect_size == 0)
				csr_details::fspmm<Field, true>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			else
				csr_details::fspmm<Field, false>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
		}else{
			if(blockSize % simd::vect_size == 0)
				csr_details::fspmm<Field, true>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
			else
				csr_details::fspmm<Field, false>(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat, ldx ,y.dat,ldy, FieldCategories::UnparametricTag ());
		}
	}

	template<class Field>
	void fspmm(
		      const Field& F,
		      const CSR<Field> & A,
		      const int blockSize,
		      const VECT<Field> & x,
		      const int ldx,
		      VECT<Field> & y,
		      const int ldy,
		      FieldCategories::GenericTag
		     )
	{
		if(ldx == A.n && ldy == blockSize){
			csr_details::fspmv(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::GenericTag());
		}else{
			csr_details::fspmv(F,A.m,A.n,A.st,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::GenericTag());
		}

	}

} // FFLAS

namespace FFLAS { namespace csr_details { /*  ZO */

	// generic
	template<class Field, bool add>
	inline void fspmv_zo(
				const Field & F,
				const index_t m,
				const index_t n,
				const index_t * st,
				const index_t * col,
				const typename Field::Element_ptr x ,
				typename Field::Element * y,
				FieldCategories::GenericTag
			       )
	{
		if(add){
			for (index_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					F.addin(y[i], x[col[j]]);
			}
		}else{
			for (index_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					F.subin(y[i], x[col[j]]);
			}
		}
	}

	//row major lda = n
	template<class Field, bool add>
	inline void fspmv_zo(
				const Field & F,
				const index_t m,
				const index_t n,
				const index_t * st,
				const index_t * col,
				const int blockSize,
				const typename Field::Element_ptr x ,
				typename Field::Element * y,
				FieldCategories::GenericTag
			    )
	{
		if(add){
			for (index_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					for(int k = 0 ; k < blockSize ; ++k)
						F.addin(y[i+k], x[col[j]+k]);
			}
		}else{
			for (index_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					for(int k = 0 ; k < blockSize ; ++k)
						F.subin(y[i+k], x[col[j]+k]);
			}
		}
	}

	//row major lda = n
	template<class Field, bool add>
	inline void fspmv_zo(
				const Field & F,
				const index_t m,
				const index_t n,
				const index_t * st,
				const index_t * col,
				const int lda,
				const int blockSize,
				const typename Field::Element_ptr x ,
				typename Field::Element * y,
				FieldCategories::GenericTag
			    )
	{
		if(add){
			for (index_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					for(int k = 0 ; k < blockSize ; ++k)
						F.addin(y[i*lda+k], x[col[j]*lda+k]);
			}
		}else{
			for (index_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					for(int k = 0 ; k < blockSize ; ++k)
						F.subin(y[i*lda+k], x[col[j]*lda+k]);
			}
		}
	}


	template<class Field, bool add>
	inline void fspmv_zo(
				const Field & F,
				const index_t m,
				const index_t n,
				const index_t * st,
				const index_t * col,
				const typename Field::Element_ptr x ,
				typename Field::Element_ptr y,
				FieldCategories::UnparametricTag
			       )
	{
		if(add){
			for (index_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					y[i] +=  x[col[j]];
			}
		}else{
			for (index_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					y[i] -=  x[col[j]];
			}
		}
	}

	//row major lda = n
	template<class Field, bool add, bool Aligned>
	inline void fspmm_zo(
				const Field & F,
				const index_t m,
				const index_t n,
				const index_t * st,
				const index_t * col,
				const int blockSize,
				const typename Field::Element_ptr x ,
				typename Field::Element * y,
				FieldCategories::UnparametricTag
			       )
	{
#if defined(__FFLASFFPACK_USE_SIMD)
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;

		vect_t vx, vy;
		if(add){
			for(index_t i = 0 ; i < m ; ++i){
				index_t start = st[i], stop = st[i+1];
				for(index_t j = start ; j < stop ; ++j){
					int k = 0;
					for(; k  < blockSize ; k+=simd::vect_size){
						if(Aligned){
							vy = simd::load(y+i+k);
							vx = simd::load(x+col[j]+k);
							simd::store(y+i+k, simd::add(vx, vy));
						}else{
							vy = simd::loadu(y+i+k);
							vx = simd::loadu(x+col[j]+k);
							simd::storeu(y+i+k, simd::add(vx, vy));
						}
					}
					for(; k < stop ; ++k){
						y[i] += x[col[j]];
					}
				}
			}
		}else{
			for(index_t i = 0 ; i < m ; ++i){
				index_t start = st[i], stop = st[i+1];
				for(index_t j = start ; j < stop ; ++j){
					int k = 0;
					for(; k  < blockSize ; k+=simd::vect_size){
						if(Aligned){
							vy = simd::load(y+i+k);
							vx = simd::load(x+col[j]+k);
							simd::store(y+i+k, simd::sub(vx, vy));
						}else{
							vy = simd::loadu(y+i+k);
							vx = simd::loadu(x+col[j]+k);
							simd::storeu(y+i+k, simd::sub(vx, vy));
						}
					}
					for(; k < stop ; ++k){
						y[i] -= x[col[j]];
					}
				}
			}
		}
#else
		if(add){
			for (index_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					for(int k = 0 ; k < blockSize ; ++k)
						y[i+k] += x[col[j]+k];
			}
		}else{
			for (index_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					for(int k = 0 ; k < blockSize ; ++k)
						y[i+k] -= x[col[j]+k];
			}
		}
#endif // __FFLASFFPACK_USE_SIMD
	}

	//row major lda != n
	template<class Field, bool add, bool Aligned>
	inline void fspmm_zo(
				const Field & F,
				const index_t m,
				const index_t n,
				const index_t * st,
				const index_t * col,
				const int lda,
				const int blockSize,
				const typename Field::Element_ptr x ,
				typename Field::Element * y,
				FieldCategories::UnparametricTag
			       )
	{
#if defined(__FFLASFFPACK_USE_SIMD)
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;

		vect_t vx, vy;
		if(add){
			for(index_t i = 0 ; i < m ; ++i){
				index_t start = st[i], stop = st[i+1];
				for(index_t j = start ; j < stop ; ++j){
					int k = 0;
					for(; k  < blockSize ; k+=simd::vect_size){
						if(Aligned){
							vy = simd::load(y+i*lda+k);
							vx = simd::load(x+col[j]*lda+k);
							simd::store(y+i*lda+k, simd::add(vx, vy));
						}else{
							vy = simd::loadu(y+i*lda+k);
							vx = simd::loadu(x+col[j]*lda+k);
							simd::storeu(y+i*lda+k, simd::add(vx, vy));
						}
					}
					for(; k < stop ; ++k){
						y[i*lda] += x[col[j]*lda];
					}
				}
			}
		}else{
			for(index_t i = 0 ; i < m ; ++i){
				index_t start = st[i], stop = st[i+1];
				for(index_t j = start ; j < stop ; ++j){
					int k = 0;
					for(; k  < blockSize ; k+=simd::vect_size){
						if(Aligned){
							vy = simd::load(y+i*lda+k);
							vx = simd::load(x+col[j]*lda+k);
							simd::store(y+i*lda+k, simd::sub(vx, vy));
						}else{
							vy = simd::loadu(y+i*lda+k);
							vx = simd::loadu(x+col[j]*lda+k);
							simd::storeu(y+i*lda+k, simd::sub(vx, vy));
						}
					}
					for(; k < stop ; ++k){
						y[i] -= x[col[j]*lda];
					}
				}
			}
		}
#else
		if(add){
			for (index_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					for(int k = 0 ; k < blockSize ; ++k)
						y[i*lda+k] += x[col[j]*lda+k];
			}
		}else{
			for (index_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					for(int k = 0 ; k < blockSize ; ++k)
						y[i*lda+k] -= x[col[j]*lda+k];
			}
		}
#endif // __FFLASFFPACK_USE_SIMD
	}


} // details
} // FFLAS

namespace FFLAS { /*  ZO */

	/* ******* */
	/* CSR_ZO  */
	/* ******* */

	// y = A x + b y ; (generic)
	// it is supposed that no reduction is needed.
	template<class Field>
	inline void fspmv(
			     const Field & F,
			     CSR_ZO<Field> & A,
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
			     CSR_ZO<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y
			     , FieldCategories::GenericTag
			    )
	{
		if (A.cst == F.one) {
			csr_details::fspmv_zo<Field,true>(F,A.m,A.n,A.st,A.col,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else if (A.cst == F.mOne) {
			csr_details::fspmv_zo<Field,false>(F,A.m,A.n,A.st,A.col,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else {
			auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x.dat, 1, x1, 1);
			csr_details::fspmv_zo<Field,true>(F,A.m,A.n,A.st,A.col,x1,y.dat, FieldCategories::GenericTag());
			fflas_delete(x1);
		}
	}

	template<class Field>
	inline void fspmv(
			     const Field & F,
			     CSR_ZO<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y
			     , FieldCategories::UnparametricTag
			    )
	{
		if (A.cst == F.one) {
			csr_details::fspmv_zo<Field,true>(F,A.m,A.n,A.st,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else if (A.cst == F.mOne) {
			csr_details::fspmv_zo<Field,false>(F,A.m,A.n,A.st,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else {
			auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x.dat, 1, x1, 1);
			csr_details::fspmv_zo<Field,true>(F,A.m,A.n,A.st,A.col,x1,y.dat, FieldCategories::UnparametricTag());
			fflas_delete(x1);
		}
	}

	template<class Field>
	inline void fspmv(
			     const Field & F,
			     CSR_ZO<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y
			     , FieldCategories::ModularTag
			    )
	{
		if (A.cst == F.one) {
			csr_details::fspmv_zo<Field,true>(F,A.m,A.n,A.st,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else if (A.cst == F.mOne) {
			csr_details::fspmv_zo<Field,false>(F,A.m,A.n,A.st,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else {
			auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x.dat, 1, x1, 1);
			csr_details::fspmv_zo<Field,true>(F,A.m,A.n,A.st,A.col,x1,y.dat, FieldCategories::UnparametricTag());
			fflas_delete(x1);
		}
		freduce (F,A.m,y.dat,1);
	}

	template<class Field>
	inline void fspmm(
			     const Field & F,
			     CSR_ZO<Field> & A,
			     const int blockSize,
			     const VECT<Field> & x,
			     const int ldx,
			     const typename Field::Element & b,
			     VECT<Field> & y,
			     const int ldy
			    )
	{
		details::init_y(F, A.m, blockSize, b, y.dat, ldy,  typename FieldTraits<Field>::category());
		fspmm(F, A, blockSize, x, ldx, y, ldy, typename FieldTraits<Field>::category() );
	}

	template<class Field>
	inline void fspmm(
			     const Field & F,
			     CSR_ZO<Field> & A,
			     const int blockSize,
			     const VECT<Field> & x,
			     const int ldx,
			     VECT<Field> & y,
			     const int ldy,
			     FieldCategories::GenericTag
			    )
	{
		using simd = Simd<typename Field::Element>;

		if(ldx == A.n && ldy == blockSize){
			if (A.cst == F.one) {
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,true, true>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,y.dat, FieldCategories::GenericTag());
				else
					csr_details::fspmm_zo<Field,true, false>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,y.dat, FieldCategories::GenericTag());
			}
			else if (A.cst == F.mOne) {
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,false, true>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,y.dat, FieldCategories::GenericTag());
				else
					csr_details::fspmm_zo<Field,false, false>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,y.dat, FieldCategories::GenericTag());
			}
			else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, ldx, x1, 1);
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,true, true>(F,A.m,A.n,A.st,A.col,blockSize,x1,y.dat, FieldCategories::GenericTag());
				else
					csr_details::fspmm_zo<Field,true, false>(F,A.m,A.n,A.st,A.col,blockSize,x1,y.dat, FieldCategories::GenericTag());
				fflas_delete(x1);
			}
		}else{
			if (A.cst == F.one) {
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,true, true>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::GenericTag());
				else
					csr_details::fspmm_zo<Field,true, false>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::GenericTag());
			}
			else if (A.cst == F.mOne) {
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,false, true>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::GenericTag());
				else
					csr_details::fspmm_zo<Field,false, false>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::GenericTag());
			}
			else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, ldx, x1, 1);
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,true, true>(F,A.m,A.n,A.st,A.col,blockSize,x1,blockSize,y.dat,ldy, FieldCategories::GenericTag());
				else
					csr_details::fspmm_zo<Field,true, false>(F,A.m,A.n,A.st,A.col,blockSize,x1,blockSize,y.dat,ldy, FieldCategories::GenericTag());
				fflas_delete(x1);
			}
		}
	}

	template<class Field>
	inline void fspmm(
			     const Field & F,
			     CSR_ZO<Field> & A,
			     const int blockSize,
			     const VECT<Field> & x,
			     const int ldx,
			     VECT<Field> & y,
			     const int ldy,
			     FieldCategories::UnparametricTag
			    )
	{
		using simd=Simd<typename Field::Element>;
		if(ldx == A.n && ldy == blockSize){
			if (A.cst == F.one) {
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,true, true>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag());
				else
					csr_details::fspmm_zo<Field,true, false>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag());
			}
			else if (A.cst == F.mOne) {
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,false, true>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag());
				else
					csr_details::fspmm_zo<Field,false, false>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag());
			}
			else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, ldx, x1, 1);
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,true, true>(F,A.m,A.n,A.st,A.col,blockSize,x1,y.dat, FieldCategories::UnparametricTag());
				else
					csr_details::fspmm_zo<Field,true, false>(F,A.m,A.n,A.st,A.col,blockSize,x1,y.dat, FieldCategories::UnparametricTag());
				fflas_delete(x1);
			}
		}else{
			if (A.cst == F.one) {
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,true, true>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag());
				else
					csr_details::fspmm_zo<Field,true, false>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag());
			}
			else if (A.cst == F.mOne) {
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,false, true>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag());
				else
					csr_details::fspmm_zo<Field,false, false>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag());
			}
			else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, ldx, x1, 1);
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,true, true>(F,A.m,A.n,A.st,A.col,blockSize,x1,blockSize,y.dat,ldy, FieldCategories::UnparametricTag());
				else
					csr_details::fspmm_zo<Field,true, false>(F,A.m,A.n,A.st,A.col,blockSize,x1,blockSize,y.dat,ldy, FieldCategories::UnparametricTag());
				fflas_delete(x1);
			}
		}
	}

	template<class Field>
	inline void fspmm(
			     const Field & F,
			     CSR_ZO<Field> & A,
			     const int blockSize,
			     const VECT<Field> & x,
			     const int ldx,
			     VECT<Field> & y,
			     const int ldy,
			     FieldCategories::ModularTag
			    )
	{
		using simd=Simd<typename Field::Element>;
		if(ldx == A.n && ldy == blockSize){
			if (A.cst == F.one) {
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,true, true>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag());
				else
					csr_details::fspmm_zo<Field,true, false>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag());
			}
			else if (A.cst == F.mOne) {
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,false, true>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag());
				else
					csr_details::fspmm_zo<Field,false, false>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag());
			}
			else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, ldx, x1, 1);
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,true, true>(F,A.m,A.n,A.st,A.col,blockSize,x1,y.dat, FieldCategories::UnparametricTag());
				else
					csr_details::fspmm_zo<Field,true, false>(F,A.m,A.n,A.st,A.col,blockSize,x1,y.dat, FieldCategories::UnparametricTag());
				fflas_delete(x1);
			}
		}else{
			if (A.cst == F.one) {
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,true, true>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag());
				else
					csr_details::fspmm_zo<Field,true, false>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag());
			}
			else if (A.cst == F.mOne) {
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,false, true>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag());
				else
					csr_details::fspmm_zo<Field,false, false>(F,A.m,A.n,A.st,A.col,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag());
			}
			else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, ldx, x1, 1);
				if(blockSize % simd::vect_size == 0)
					csr_details::fspmm_zo<Field,true, true>(F,A.m,A.n,A.st,A.col,blockSize,x1,blockSize,y.dat,ldy, FieldCategories::UnparametricTag());
				else
					csr_details::fspmm_zo<Field,true, false>(F,A.m,A.n,A.st,A.col,blockSize,x1,blockSize,y.dat,ldy, FieldCategories::UnparametricTag());
				fflas_delete(x1);
			}
		}
		freduce (F,blockSize,A.m,y.dat,ldy);
	}

} // FFLAS

namespace FFLAS { /*  conversions */

    template<class Field, class IdxT>
    void sp_csr_from_coo(const Field & F, const IdxT COO_rowdim, const IdxT COO_coldim,
                         const uint64_t COO_nnz, const IdxT * COO_row, const IdxT * COO_col,
                         const typename Field::Element_ptr COO_val, index_t & CSR_rowdim,
                         index_t & CSR_coldim, index_t & CSR_maxrow, index_t *& CSR_row, index_t *& CSR_col,
                         typename Field::Element_ptr & CSR_val, bool ZO){
        CSR_col = fflas_new<index_t>(COO_nnz, Alignment::CACHE_LINE);
        CSR_row = fflas_new<index_t>(COO_rowdim+1, Alignment::CACHE_LINE);
        if(!ZO){
            CSR_val = fflas_new(F, COO_nnz, 1, Alignment::CACHE_LINE);
        }
        for(size_t i = 0 ; i < COO_nnz ; ++i){
            CSR_col[i] = static_cast<index_t>(COO_col[i]);
            if(!ZO){
                CSR_val[i] = COO_val[i];
            }
        }
        CSR_rowdim = COO_rowdim;
        CSR_coldim = COO_coldim;
        for(size_t i = 0 ; i <= COO_rowdim ; ++i){
            CSR_row[i] = 0;
        }
        for(size_t i = 0 ; i < COO_nnz ; ++i){
            CSR_row[COO_row[i]+1]++;
        }
        CSR_maxrow = *(std::max_element(CSR_row, CSR_row+COO_rowdim+1));
        for(size_t i = 1 ; i <= COO_rowdim ; ++i){
            CSR_row[i] += CSR_row[i-1];
        }
    }

} // FFLAS


namespace FFLAS{ /* Delete matrix */
    template<class Field>
    void sp_delete(CSR<Field> & M){
        fflas_delete(M.dat);
        fflas_delete(M.st);
        fflas_delete(M.col);
    }

    template<class Field>
    void sp_delete(CSR_sub<Field> & M){
        fflas_delete(M.dat);
        fflas_delete(M.st);
        fflas_delete(M.col);
    }

    template<class Field>
    void sp_delete(CSR_ZO<Field> & M){
        fflas_delete(M.st);
        fflas_delete(M.col);
    }
}

namespace FFLAS { namespace csr_details {
#ifdef __FFLASFFPACK_HAVE_CUDA
	// Double
	template<>
	void fspmv(
		      const Givaro::DoubleDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      const index_t m,
		      const index_t n,
		      const index_t * st_d,
		      const index_t * col_d,
		      const double *  dat_d,
		      const double * x_d ,
		      const double & b,
		      double * y_d,
		      const cusparseMatDescr_t & descrA,
		      const cusparseHandle_t & handle
		     )
	{
		// char * transa = (ta==FflasNoTrans)?'n':'t';
		index_t nnz = st[m]-st[m-1];
		double one = 1.f;
		status= cusparseDcsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, nnz,
				       &one, descrA, dat_d, st_d, col_d,
				       x_d, &b, y_d);
	}

	// Float
	template<>
	void fspmv(
		      const Givaro::FloatDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      const index_t m,
		      const index_t n,
		      const index_t * st_d,
		      const index_t * col_d,
		      const float *  dat_d,
		      const float * x_d ,
		      const float & b,
		      float * y_d,
		      const cusparseMatDescr_t & descrA,
		      const cusparseHandle_t & handle
		     )
	{
		// char * transa = (ta==FflasNoTrans)?'n':'t';
		index_t nnz = st[m]-st[m-1];
		float one = 1.f;
		status= cusparseScsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, nnz,
				       &one, descrA, dat_d, st_d, col_d,
				       x_d, &b, y_d);

	}

	// need cuda freduce code (need nvcc)
#endif // __FFLASFFPACK_HAVE_CUDA

} // details
} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_spmv_csr_INL
