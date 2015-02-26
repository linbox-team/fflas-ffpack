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

/** @file fflas/fflas_fspmv_coo.inl
 * NO DOC
 */

#ifndef __FFLASFFPACK_fflas_fflas_spmv_coo_INL
#define __FFLASFFPACK_fflas_fflas_spmv_coo_INL

namespace FFLAS { namespace coo_details {

	// y = A x + b y ; (generic)
	template<class Field>
	inline void fspmv(
			     const Field& F,
			     const index_t m,
                             const index_t n,
                             const uint64_t z,
			     const index_t * row,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const typename Field::Element_ptr x ,
			     typename Field::Element_ptr y
			     , FieldCategories::GenericTag
			    )
	{
		for (index_t j = 0 ; j < z ; ++j)
			F.axpyin(y[row[j]],dat[j],x[col[j]]);
	}

	// row major lda = n
	template<class Field>
	inline void fspmm(const Field & F, const index_t m, const index_t n, const uint64_t z, const index_t * row,
					  const index_t * col, const typename Field::Element_ptr dat, const int blockSize, const typename Field::Element_ptr x,
					  typename Field::Element_ptr y, FieldCategories::GenericTag){
		for(uint64_t j = 0 ; j < z ; ++j){
			for(int k = 0 ; k < blockSize ; ++k){
				F.axpyin(y[row[j]*blockSize+k], dat[j], x[col[j]*blockSize+k]);
			}
		}
	}

	// lda != n
	template<class Field>
	inline void fspmm(const Field & F, const index_t m, const index_t n, const uint64_t z, const index_t * row,
					  const index_t * col, const typename Field::Element_ptr dat, const int blockSize,
					  const typename Field::Element_ptr x, const int ldx,
					  typename Field::Element_ptr y, const int ldy, FieldCategories::GenericTag){
		for(uint64_t j = 0 ; j < z ; ++j){
			for(int k = 0 ; k < blockSize ; ++k){
				F.axpyin(y[row[j]*ldy*blockSize+k], dat[j], x[col[j]*ldx*blockSize+k]);
			}
		}
	}


	template<class Field>
	void fspmv(
		      const Field & F,
		      const index_t m,
		      const index_t n,
		      const uint64_t z,
		      const index_t * row,
		      const index_t * col,
		      const typename Field::Element_ptr dat,
		      const typename Field::Element_ptr x ,
		      typename Field::Element_ptr y
		      , FieldCategories::UnparametricTag
		     )
	{
		for (size_t i = 0 ; i < z ; ++i) {
			y[row[i]] += dat[i] * x[col[i]];
		}
	}

	// row major lda = n
	template<class Field, bool Aligned>
	inline void fspmm(const Field & F, const index_t m, const index_t n, const uint64_t z, const index_t * row,
					  const index_t * col, const typename Field::Element_ptr dat, const int blockSize, const typename Field::Element_ptr x,
					  typename Field::Element_ptr y, FieldCategories::UnparametricTag){
#if defined(__FFLASFFPACK_USE_SIMD)
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;

		vect_t vy, vx, vdat;
		for(uint64_t j = 0 ; j < z ; ++j){
			int k = 0;
			for(; k < blockSize ; k += simd::vect_size){
				if(Aligned){
					vy =simd::load(y+row[j]*blockSize+k);
					vx =simd::load(x+col[j]*blockSize+k);
				}else{
					vy =simd::loadu(y+row[j]*blockSize+k);
					vx =simd::loadu(x+col[j]*blockSize+k);
				}
				vdat =simd::set1(dat[j]);
				if(Aligned){
					simd::store(y+row[j]*blockSize+k, simd::fmadd(vy, vdat, vx));
				}else{
					simd::storeu(y+row[j]*blockSize+k, simd::fmadd(vy, vdat, vx));
				}
			}
			for(; k < blockSize ; ++k){
				y[row[j]*blockSize+k] += dat[j] * x[col[j]+k];
			}
		}

#else
		for(index_t j = 0 ; j < z ; ++j){
			for(int k = 0 ; k < blockSize ; ++k){
				y[row[j]*blockSize+k] += dat[j] * x[col[j]*blockSize+k];
			}
		}
#endif // __FFLASFFPACK_USE_SIMD
	}


// row major lda != n
	template<class Field, bool Aligned>
	inline void fspmm(const Field & F, const index_t m, const index_t n, const uint64_t z, const index_t * row,
					  const index_t * col, const typename Field::Element_ptr dat, const int blockSize,
					  const typename Field::Element_ptr x, const int ldx,
					  typename Field::Element_ptr y, const int ldy, FieldCategories::UnparametricTag){
#if defined(__FFLASFFPACK_USE_SIMD)
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;

		vect_t vy, vx, vdat;
		for(uint64_t j = 0 ; j < z ; ++j){
			int k = 0;
			vdat =simd::set1(dat[j]);
			for(; k < blockSize ; k += simd::vect_size){
				if(Aligned){
					vy =simd::load(y+row[j]*ldy*blockSize+k);
					vx =simd::load(x+col[j]*ldx*blockSize+k);
					simd::store(y+row[j]*ldy*blockSize+k, simd::fmadd(vy, vdat, vx));
				}else{
					vy =simd::loadu(y+row[j]*ldy*blockSize+k);
					vx =simd::loadu(x+col[j]*ldx*blockSize+k);
					simd::storeu(y+row[j]*ldy*blockSize+k, simd::fmadd(vy, vdat, vx));
				}
			}
			for(; k < blockSize ; ++k){
				y[row[j]*ldy*blockSize+k] += dat[j] * x[col[j]*ldx*blockSize+k];
			}
		}
#else
		for(index_t j = 0 ; j < z ; ++j){
			for(int k = 0 ; k < blockSize ; ++k){
				y[row[j]*ldy*blockSize+k] += dat[j] * x[col[j]*ldx*blockSize+k];
			}
		}
#endif // __FFLASFFPACK_USE_SIMD
	}

	// Double
	template<>
	void fspmv(
		      const Givaro::DoubleDomain & F,
		      const index_t m,
		      const index_t n,
		      const uint64_t z,
		      const index_t * row,
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
		index_t z_ = (index_t) z ;
		double * yd =  FFLAS::fflas_new<double >(m);
		// mkl_dcoogemv (bug too ?)
		mkl_cspblas_dcoogemv
		(&transa, &m_, const_cast<double*>(dat), const_cast<index_t*>(row) , const_cast<index_t*>(col), &z_, const_cast<double*>(x), yd);
		faddin(F,m,yd,1,y,1);
		FFLAS::fflas_delete( yd );
#else

		for (size_t i = 0 ; i < z ; ++i) {
			y[row[i]] += dat[i] * x[col[i]];
		}
#endif // __FFLASFFPACK_HAVE_MKL
	}


	// Float
	template<>
	void fspmv(
		      const Givaro::FloatDomain & F,
		      const index_t m,
		      const index_t n,
		      const uint64_t z,
		      const index_t * row,
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
		index_t z_ = (index_t) z ;
		float * yd = FFLAS::fflas_new<float >(m);
		// mkl_scoogemv
		mkl_cspblas_scoogemv
		(&transa, &m_, const_cast<float*>(dat), const_cast<index_t*>(row), const_cast<index_t*>(col), &z_, const_cast<float*>(x), yd);
		faddin(F,m,yd,1,y,1);
		FFLAS::fflas_delete( yd );
#else
		for (index_t j = 0 ; j < (index_t) z ; ++j)
			y[row[j]] += dat[j] * x[col[j]];
#endif // __FFLASFFPACK_HAVE_MKL
	}


	// delayed by kmax
	template<class Field>
	inline void fspmv(
			     const Field& F,
			     const index_t m,
			     const index_t n,
			     const uint64_t z,
			     const index_t * row,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const typename Field::Element_ptr x ,
			     typename Field::Element_ptr y,
			     const index_t & kmax
			    )
	{
		size_t w = 0 ;
		index_t last_i = 0;
		typename Field::Element e ;
		F.init(e,y[last_i]);
		size_t accu = 0 ;

		while ( w < z) {
			if ( row[w] == last_i ) { // same line
				if (accu < (size_t)kmax) {
					e += dat[w] * x[col[w]] ;
					accu += 1 ;
				}
				else {
					F.axpyin(e,dat[w],x[col[w]]);
					accu = 0 ;
				}
			}
			else { // new line
				F.init(y[last_i],e);
				last_i = row[w] ;
				F.init(e,y[last_i]);
				e += dat[w] * x[col[w]];
				accu = 1 ;
			}
			++w ;
		}
		F.init(y[last_i],e);
	}

	// row major lda = n
	template<class Field>
	inline void fspmm(const Field & F, const index_t m, const index_t n, const uint64_t z, const index_t * row,
					  const index_t * col, const typename Field::Element_ptr dat, const int blockSize, const typename Field::Element_ptr x,
					  typename Field::Element_ptr y, const index_t & kmax){
		index_t currentRow = row[0];

		auto e = fflas_new(F, blockSize, 1);
		for(int i = 0 ; i < blockSize ; ++i){
			F.init(e[i], y[currentRow*blockSize+i]);
		}

		int accu = 0;

		for(uint64_t i = 0 ; i < z ; ++i){
			if(row[i] == currentRow){
				if(accu < kmax){
					for(int k = 0 ; k < blockSize ; ++k){
						e[k]+= dat[i]*x[col[i]*blockSize+k];
					}
					++accu;
				}else{
					for(int k = 0 ; k < blockSize ; ++k){
						F.axpyin(e[k], dat[i], x[col[i]*blockSize+k]);
					}
					accu = 0;
				}
			}else{
				for(int k = 0 ; k < blockSize ; ++k){
					F.init(y[currentRow*blockSize+k], e[k]);
				}
				currentRow = row[i];
				for(int k = 0 ; k < blockSize ; ++k){
					F.init(e[k], y[currentRow*blockSize+k]);
					e[k] += dat[i]*x[col[i]*blockSize+k];
				}
				accu = 1;
			}
			for(int k = 0 ; k < blockSize ; ++k){
					F.init(y[currentRow*blockSize+k], e[k]);
			}
		}
	}

	// row major lda != n
	template<class Field>
	inline void fspmm(const Field & F, const index_t m, const index_t n, const uint64_t z, const index_t * row,
					  const index_t * col, const typename Field::Element_ptr dat, const int blockSize,
					  const typename Field::Element_ptr x, const int ldx,
					  typename Field::Element_ptr y, const int ldy, const index_t & kmax){
		index_t currentRow = row[0];

		auto e = fflas_new(F, blockSize, 1);
		for(int i = 0 ; i < blockSize ; ++i){
			F.init(e[i], y[currentRow*ldy*blockSize+i]);
		}

		int accu = 0;

		for(uint64_t i = 0 ; i < z ; ++i){
			if(row[i] == currentRow){
				if(accu < kmax){
					for(int k = 0 ; k < blockSize ; ++k){
						e[k]+= dat[i]*x[col[i]*ldx*blockSize+k];
					}
					++accu;
				}else{
					for(int k = 0 ; k < blockSize ; ++k){
						F.axpyin(e[k], dat[i], x[col[i]*ldx*blockSize+k]);
					}
					accu = 0;
				}
			}else{
				for(int k = 0 ; k < blockSize ; ++k){
					F.init(y[currentRow*ldy*blockSize+k], e[k]);
				}
				currentRow = row[i];
				for(int k = 0 ; k < blockSize ; ++k){
					F.init(e[k], y[currentRow*ldy+k]);
					e[k] += dat[i]*x[col[i]*ldx*blockSize+k];
				}
				accu = 1;
			}
			for(int k = 0 ; k < blockSize ; ++k){
					F.init(y[currentRow*ldy*blockSize+k], e[k]);
			}
		}
	}

} // coo_details
} // FFLAS

namespace FFLAS {

	/* ******* */
	/* COO_sub */
	/* ******* */

	// y = A x + b y ; (generic)
	// it is supposed that no reduction is needed.
	template<class Field>
	inline void fspmv(
			     const Field& F,
			     const COO_sub<Field> & A,
			     const VECT<Field> & x,
			     const typename Field::Element & b,
			     VECT<Field> & y
			    )
	{
		details::init_y(F, A.m, b, y.dat,  typename FieldTraits<Field>::category());
		fspmv( F, A, x, y, typename FieldTraits<Field>::category());
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const COO_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::ModularTag)
	{
		coo_details::fspmv(F,A.m,A.n,A.z,A.row,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag ());
		freduce(F,A.m,y.dat,1);
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const COO_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::UnparametricTag)
	{
		coo_details::fspmv(F,A.m,A.n,A.z,A.row,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag ());
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const COO_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::GenericTag)
	{
		coo_details::fspmv(F,A.m,A.n,A.z,A.row,A.col,A.dat,x.dat,y.dat, FieldCategories::GenericTag());
	}

	template<class Field>
	inline void fspmm(
			     const Field& F,
			     const COO_sub<Field> & A,
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
			     const COO_sub<Field> & A,
			     const int blockSize,
			     const VECT<Field> & x,
			     const int ldx,
			     VECT<Field> & y,
			     const int ldy,
			     FieldCategories::ModularTag)
	{
		using simd = Simd<typename Field::Element>;
		if(ldx == A.n && ldy == blockSize)
		{
			if(blockSize % simd::vect_size == 0)
				coo_details::fspmm<Field, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			else
				coo_details::fspmm<Field, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
		}
		else{
			if(blockSize % simd::vect_size == 0)
				coo_details::fspmv<Field, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
			else
				coo_details::fspmm<Field, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
		}
		freduce(F,blockSize,A.m,y.dat,ldy);
	}

	template<class Field>
	inline void fspmm(const Field & F,
			     const COO_sub<Field> & A,
			     const int blockSize,
			     const VECT<Field> & x,
			     const int ldx,
			     VECT<Field> & y,
			     const int ldy,
			     FieldCategories::UnparametricTag)
	{
		using simd = Simd<typename Field::Element>;
		if(ldx == A.n && ldy == blockSize)
		{
			if(blockSize % simd::vect_size == 0)
				coo_details::fspmm<Field, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			else
				coo_details::fspmm<Field, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
		}
		else{
			if(blockSize % simd::vect_size == 0)
				coo_details::fspmv<Field, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
			else
				coo_details::fspmm<Field, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
		}
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const COO_sub<Field> & A,
			     const int blockSize,
			     const VECT<Field> & x,
			     const int ldx,
			     VECT<Field> & y,
			     const int ldy,
			     FieldCategories::GenericTag)
	{
		if(ldx == A.n && ldy == blockSize){
			coo_details::fspmm(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::GenericTag());
		}else{
			coo_details::fspmm(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::GenericTag());
		}
	}


	/* *** */
	/* COO */
	/* *** */

	// y = A x + b y ; (generic)
	// reductions are delayed.
	template<class Field>
	void fspmv(
		      const Field& F,
		      const COO<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     )
	{
		details::init_y(F, A.m, b, y.dat,  typename FieldTraits<Field>::category());
		fspmv(F,A,x,y, typename FieldTraits<Field>::category());
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const COO<Field> & A,
		      const VECT<Field> & x,
		      VECT<Field> & y
		      , FieldCategories::GenericTag
		     )
	{
		coo_details::fspmv(F,A.m,A.n,A.z,A.row,A.col,A.dat,x.dat,y.dat, FieldCategories::GenericTag());
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const COO<Field> & A,
		      const VECT<Field> & x,
		      VECT<Field> & y
		      , FieldCategories::UnparametricTag
		     )
	{
		coo_details::fspmv(F,A.m,A.n,A.z,A.row,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag());
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const COO<Field> & A,
		      const VECT<Field> & x,
		      VECT<Field> & y
		      , FieldCategories::ModularTag
		     )
	{
		index_t kmax = static_cast<index_t>(Protected::DotProdBoundClassic(F,F.one));
                coo_details::fspmv(F, A.m, A.n, A.z, A.row, A.col, A.dat, x.dat, y.dat, kmax);
	}

	template<class Field>
	void fspmm(
		      const Field& F,
		      const COO<Field> & A,
		      const int blockSize,
		      const VECT<Field> & x,
		      const int ldx,
		      const typename Field::Element & b,
		      VECT<Field> & y,
		      const int ldy
		     )
	{
		details::init_y(F, A.m, blockSize, b, y.dat, ldy, typename FieldTraits<Field>::category());
		fspmm(F, A, blockSize, x,ldx, y,ldy, typename FieldTraits<Field>::category());
	}

	template<class Field>
	void fspmm(
		      const Field& F,
		      const COO<Field> & A,
		      const int blockSize,
		      const VECT<Field> & x,
		      const int ldx,
		      VECT<Field> & y,
		      const int ldy,
		      FieldCategories::GenericTag
		     )
	{
		if(ldx == A.n && ldy == blockSize){
			coo_details::fspmm(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::GenericTag());
		}else{
			coo_details::fspmm(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::GenericTag());
		}
	}

	template<class Field>
	void fspmm(
		      const Field& F,
		      const COO<Field> & A,
		      const int blockSize,
		      const VECT<Field> & x,
		      const int ldx,
		      VECT<Field> & y,
		      const int ldy,
		      FieldCategories::UnparametricTag
		     )
	{
		using simd = Simd<typename Field::Element>;
		if(ldx == A.n && ldy == blockSize)
		{
			if(blockSize % simd::vect_size == 0)
				coo_details::fspmm<Field, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			else
				coo_details::fspmm<Field, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
		}
		else{
			if(blockSize % simd::vect_size == 0)
				coo_details::fspmm<Field, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
			else
				coo_details::fspmm<Field, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
		}
	}

	template<class Field>
	void fspmm(
		      const Field& F,
		      const COO<Field> & A,
		      const int blockSize,
		      const VECT<Field> & x,
		      const int ldx,
		      VECT<Field> & y,
		      const int ldy,
		      FieldCategories::ModularTag
		     )
	{
		index_t kmax = static_cast<index_t>(Protected::DotProdBoundClassic(F,F.one));
		if(ldx == A.n && ldy == blockSize)
        	coo_details::fspmm(F, A.m, A.n, A.z, A.row, A.col, A.dat, blockSize, x.dat, y.dat, kmax);
        else
        	coo_details::fspmm(F, A.m, A.n, A.z, A.row, A.col, A.dat, blockSize, x.dat, ldx, y.dat, ldy, kmax);
	}
} // FFLAS

namespace FFLAS { namespace coo_details { /*  ZO */

	// generic
	template<class Field, bool add>
	inline void fspmv_zo(
				const Field & F,
				const index_t m,
				const index_t n,
				const uint64_t z,
				const index_t * row,
				const index_t * col,
				const typename Field::Element_ptr x ,
				typename Field::Element_ptr y
				, FieldCategories::GenericTag
			       )
	{
		if (add) {
			for (index_t j = 0 ; j < z ; ++j)
				F.addin(y[row[j]], x[col[j]]);
		}
		else {
			for (index_t j = 0 ; j < z ; ++j)
				F.subin(y[row[j]], x[col[j]]);
		}
	}

	// row major lda = n
	template<class Field, bool add>
	inline void fspmm_zo(
				const Field & F,
				const index_t m,
				const index_t n,
				const uint64_t z,
				const index_t * row,
				const index_t * col,
				const int blockSize,
				const typename Field::Element_ptr x ,
				typename Field::Element_ptr y
				, FieldCategories::GenericTag
			       )
	{
		if (add) {
			for (index_t j = 0 ; j < z ; ++j)
				for(int k = 0 ; k < blockSize ; ++k)
					F.addin(y[row[j]+k], x[col[j]+k]);
		}
		else {
			for (index_t j = 0 ; j < z ; ++j)
				for(int k = 0 ; k < blockSize ; ++k)
					F.subin(y[row[j+k]], x[col[j]+k]);
		}
	}

	// row major lda != n
	template<class Field, bool add>
	inline void fspmm_zo(
				const Field & F,
				const index_t m,
				const index_t n,
				const uint64_t z,
				const index_t * row,
				const index_t * col,
				const int blockSize,
				const typename Field::Element_ptr x ,
				const int ldx,
				typename Field::Element_ptr y,
				const int ldy,
				FieldCategories::GenericTag
			    )
	{
		if (add) {
			for (index_t j = 0 ; j < z ; ++j)
				for(int k = 0 ; k < blockSize ; ++k)
					F.addin(y[row[j]*ldy+k], x[col[j]*ldx+k]);
		}
		else {
			for (index_t j = 0 ; j < z ; ++j)
				for(int k = 0 ; k < blockSize ; ++k)
					F.subin(y[row[j+k]*ldy], x[col[j]*ldx+k]);
		}
	}

	template<class Field, bool add>
	inline void fspmv_zo(
				const Field & F,
				const index_t m,
				const index_t n,
				const uint64_t z,
				const index_t * row,
				const index_t * col,
				const typename Field::Element_ptr x ,
				typename Field::Element_ptr y,
				FieldCategories::UnparametricTag
			       )
	{
		if (add){
			for (index_t j = 0 ; j < z ; ++j)
				y[row[j]] +=  x[col[j]];
		}else{
			for (index_t j = 0 ; j < z ; ++j)
				y[row[j]] -=  x[col[j]];
		}
	}

	// row major lda = n
	template<class Field, bool add, bool Aligned>
	inline void fspmm_zo(
				const Field & F,
				const index_t m,
				const index_t n,
				const uint64_t z,
				const index_t * row,
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

		int end = (blockSize/simd::vect_size)*simd::vect_size;

		vect_t vy, vx;
		if(add){
			for (index_t j = 0 ; j < z ; ++j)
			{
				int k = 0;
				for(; k < end ; k+=simd::vect_size){
					if(Aligned){
						vy = simd::load(y+row[j]+k);
						vx = simd::load(x+col[j]+k);
						simd::store(y+row[j]+k, simd::add(vy, vx));
					}else{
						vy = simd::loadu(y+row[j]+k);
						vx = simd::loadu(x+col[j]+k);
						simd::storeu(y+row[j]+k, simd::add(vy, vx));
					}
				}
				for(; k < blockSize ; ++k){
					y[row[j]+k] +=  x[col[j]+k];
				}
			}
		}else{
			for (index_t j = 0 ; j < z ; ++j)
			{
				int k = 0;
				for(; k < end ; k+=simd::vect_size){
					if(Aligned){
						vy = simd::load(y+row[j]+k);
						vx = simd::load(x+col[j]+k);
						simd::store(y+row[j]+k, simd::add(vy, vx));
					}else{
						vy = simd::loadu(y+row[j]+k);
						vx = simd::loadu(x+col[j]+k);
						simd::storeu(y+row[j]+k, simd::add(vy, vx));
					}
				}
				for(; k < blockSize ; ++k){
					y[row[j]+k] +=  x[col[j]+k];
				}
			}
		}
#else
		if (add){
			for (index_t j = 0 ; j < z ; ++j)
				for(int k = 0 ; k < blockSize ; ++k)
					y[row[j]+k] +=  x[col[j]+k];
		}else{
			for (index_t j = 0 ; j < z ; ++j)
				for(int k = 0 ; k < blockSize ; ++k)
					y[row[j]+k] -=  x[col[j]+k];
		}
#endif
	}

	// row major lda != n
	template<class Field, bool add, bool Aligned>
	inline void fspmm_zo(
				const Field & F,
				const index_t m,
				const index_t n,
				const uint64_t z,
				const index_t * row,
				const index_t * col,
				const int blockSize,
				const typename Field::Element_ptr x ,
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

		vect_t vy, vx;
		if(add){
			for (index_t j = 0 ; j < z ; ++j)
			{
				int k = 0;
				for(; k < end ; k+=simd::vect_size){
					if(Aligned){
						vy = simd::load(y+row[j]*ldy+k);
						vx = simd::load(x+col[j]*ldx+k);
						simd::store(y+row[j]*ldy+k, simd::add(vy, vx));
					}else{
						vy = simd::loadu(y+row[j]*ldy+k);
						vx = simd::loadu(x+col[j]*ldx+k);
						simd::storeu(y+row[j]*ldy+k, simd::add(vy, vx));
					}
				}
				for(; k < blockSize ; ++k){
					y[row[j]*ldy+k] +=  x[col[j]*ldx+k];
				}
			}
		}else{
			for (index_t j = 0 ; j < z ; ++j)
			{
				int k = 0;
				for(; k < end ; k+=simd::vect_size){
					if(Aligned){
						vy = simd::load(y+row[j]*ldy+k);
						vx = simd::load(x+col[j]*ldx+k);
						simd::store(y+row[j]*ldy+k, simd::sub(vy, vx));
					}else{
						vy = simd::loadu(y+row[j]*ldy+k);
						vx = simd::loadu(x+col[j]*ldx+k);
						simd::storeu(y+row[j]*ldy+k, simd::sub(vy, vx));
					}
				}
				for(; k < blockSize ; ++k){
					y[row[j]*ldy+k] +=  x[col[j]*ldx+k];
				}
			}
		}
#else
		if (add){
			for (index_t j = 0 ; j < z ; ++j)
				for(int k = 0 ; k < blockSize ; ++k)
					y[row[j]*ldy+k] +=  x[col[j]*ldx+k];
		}else{
			for (index_t j = 0 ; j < z ; ++j)
				for(int k = 0 ; k < blockSize ; ++k)
					y[row[j]*ldy+k] -=  x[col[j]*ldx+k];
		}
#endif
	}
} // details
} // FFLAS

namespace FFLAS { /*  ZO */

	/* ******* */
	/* COO_ZO  */
	/* ******* */

	// y = A x + b y ; (generic)
	// it is supposed that no reduction is needed.
	template<class Field>
	inline void fspmv(
			     const Field & F,
			     COO_ZO<Field> & A,
			     const VECT<Field> & x,
			     const typename Field::Element & b,
			     VECT<Field> & y
			    )
	{
		details::init_y(F, A.m, b, y.dat,  typename FieldTraits<Field>::category());
		fspmv(F, A, x, y, typename FieldTraits<Field>::category() );
	}

	template<class Field>
	inline void fspmv(
			     const Field & F,
			     COO_ZO<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y
			     , FieldCategories::GenericTag
			    )
	{
		if (A.cst == F.one) {
			coo_details::fspmv_zo<Field,true>(F,A.m,A.n,A.z,A.row,A.col,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else if (A.cst == F.mOne) {
			coo_details::fspmv_zo<Field,false>(F,A.m,A.n,A.z,A.row,A.col,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else {
			auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x.dat, 1, x1, 1);
			coo_details::fspmv_zo<Field,true>(F,A.m,A.n,A.z,A.row,A.col,x1,y.dat, FieldCategories::GenericTag());
			fflas_delete(x1);
		}
	}

	template<class Field>
	inline void fspmv(
			     const Field & F,
			     COO_ZO<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y
			     , FieldCategories::UnparametricTag
			    )
	{
		if (A.cst == F.one) {
			coo_details::fspmv_zo<Field,true>(F,A.m,A.n,A.z,A.row,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else if (A.cst == F.mOne) {
			coo_details::fspmv_zo<Field,false>(F,A.m,A.n,A.z,A.row,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else {
			auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x.dat, 1, x1, 1);
			coo_details::fspmv_zo<Field,true>(F,A.m,A.n,A.z,A.row,A.col,x1,y.dat, FieldCategories::UnparametricTag());
			fflas_delete(x1);
		}
	}

	template<class Field>
	inline void fspmv(
			     const Field & F,
			     COO_ZO<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y
			     , FieldCategories::ModularTag
			    )
	{
		if (A.cst == F.one) {
			coo_details::fspmv_zo<Field,true>(F,A.m,A.n,A.z,A.row,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else if (A.cst == F.mOne) {
			coo_details::fspmv_zo<Field,false>(F,A.m,A.n,A.z,A.row,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else {
			auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x.dat, 1, x1, 1);
			coo_details::fspmv_zo<Field,true>(F,A.m,A.n,A.z,A.row,A.col,x1,y.dat, FieldCategories::UnparametricTag());
			fflas_delete(x1);
		}
		freduce(F,A.m,y.dat,1);
	}

	template<class Field>
	inline void fspmm(
			     const Field & F,
			     COO_ZO<Field> & A,
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
			     COO_ZO<Field> & A,
			     const int blockSize,
			     const VECT<Field> & x,
			     const int ldx,
			     VECT<Field> & y,
			     const int ldy,
			     FieldCategories::GenericTag
			    )
	{
		if(ldx == A.n && ldy == blockSize){
			if (A.cst == F.one) {
				coo_details::fspmm_zo<Field,true>(F,A.m,A.n,A.z,A.row,A.col,blockSize,x.dat,y.dat, FieldCategories::GenericTag());
			}
			else if (A.cst == F.mOne) {
				coo_details::fspmm_zo<Field,false>(F,A.m,A.n,A.z,A.row,A.col,blockSize,x.dat,y.dat, FieldCategories::GenericTag());
			}
			else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_PAGESIZE);
				fscal(F, A.n, blockSize, A.cst, x.dat, ldx, x1, 1);
				coo_details::fspmm_zo<Field,true>(F,A.m,A.n,A.z,A.row,A.col,blockSize,x1,y.dat, FieldCategories::GenericTag());
				fflas_delete(x1);
			}
		}else{
			if (A.cst == F.one){
				coo_details::fspmm_zo<Field,true>(F,A.m,A.n,A.z,A.row,A.col,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::GenericTag());
			}
			else if (A.cst == F.mOne) {
				coo_details::fspmm_zo<Field,false>(F,A.m,A.n,A.z,A.row,A.col,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::GenericTag());
			}
			else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_PAGESIZE);
				fscal(F, A.n, blockSize, A.cst, x.dat, ldx, x1, 1);
				coo_details::fspmm_zo<Field,true>(F,A.m,A.n,A.z,A.row,A.col,blockSize,x1,1,y.dat,ldy, FieldCategories::GenericTag());
				fflas_delete(x1);
			}
		}
	}

	template<class Field>
	inline void fspmm(
			     const Field & F,
			     COO_ZO<Field> & A,
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
			if (A.cst == F.one){
				if(blockSize % simd::vect_size == 0)
					coo_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					coo_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}else if (A.cst == F.mOne) {
				if(blockSize % simd::vect_size == 0)
					coo_details::fspmm_zo<Field, false, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					coo_details::fspmm_zo<Field, false, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, ldx, x1, 1);
				if(blockSize % simd::vect_size == 0)
					coo_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					coo_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				fflas_delete(x1);
			}
		}else{
			if (A.cst == F.one){
				if(blockSize % simd::vect_size == 0)
					coo_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
				else
					coo_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
			}
			else if (A.cst == F.mOne) {
				if(blockSize % simd::vect_size == 0)
					coo_details::fspmm_zo<Field, false, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
				else
					coo_details::fspmm_zo<Field, false, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
			}
			else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, ldx, x1, 1);
				if(blockSize % simd::vect_size == 0)
					coo_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
				else
					coo_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
				fflas_delete(x1);
			}
		}
	}

	template<class Field>
	inline void fspmm(
			     const Field & F,
			     COO_ZO<Field> & A,
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
			if (A.cst == F.one){
				if(blockSize % simd::vect_size == 0)
					coo_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					coo_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}else if (A.cst == F.mOne) {
				if(blockSize % simd::vect_size == 0)
					coo_details::fspmm_zo<Field, false, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					coo_details::fspmm_zo<Field, false, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
			}else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, ldx, x1, 1);
				if(blockSize % simd::vect_size == 0)
					coo_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				else
					coo_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,y.dat, FieldCategories::UnparametricTag ());
				fflas_delete(x1);
			}
		}else{
			if (A.cst == F.one){
				if(blockSize % simd::vect_size == 0)
					coo_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
				else
					coo_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
			}
			else if (A.cst == F.mOne) {
				if(blockSize % simd::vect_size == 0)
					coo_details::fspmm_zo<Field, false, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
				else
					coo_details::fspmm_zo<Field, false, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
			}
			else {
				auto x1 = fflas_new(F, A.n, blockSize, Alignment::CACHE_LINE);
				fscal(F, A.n*blockSize, A.cst, x.dat, ldx, x1, 1);
				if(blockSize % simd::vect_size == 0)
					coo_details::fspmm_zo<Field, true, true>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
				else
					coo_details::fspmm_zo<Field, true, false>(F,A.m,A.n,A.z,A.row,A.col,A.dat,blockSize,x.dat,ldx,y.dat,ldy, FieldCategories::UnparametricTag ());
				fflas_delete(x1);
			}
		}
		freduce(F,blockSize,A.m,y.dat,ldy);
	}
} // FFLAS

namespace FFLAS{ /* delete Matrix */
    template<class Field>
    inline void sp_delete(COO<Field> & M){
        fflas_delete(M.col);
        fflas_delete(M.row);
        fflas_delete(M.dat);
    }

    template<class Field>
    inline void sp_delete(COO_sub<Field> & M){
        fflas_delete(M.col);
        fflas_delete(M.row);
        fflas_delete(M.dat);
    }

    template<class Field>
    inline void sp_delete(COO_ZO<Field> & M){
        fflas_delete(M.col);
        fflas_delete(M.row);
    }
}

#endif // __FFLASFFPACK_fflas_fflas_spmv_coo_INL
