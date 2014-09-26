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

#include "fflas-ffpack/fflas/fflas_simd.h"

namespace FFLAS { /*  ELL */

	/*
	 * When using SIMD, we suppose that the matrix is padded with rows of zeros if necessary.
	 */

	template<class _Element, bool _Simd =
#ifdef __FFLASFFPACK_USE_SIMD
	true
#else
	false
#endif // SiMD
	>
	struct ELL {
		size_t m = 0;
		size_t n = 0;
		size_t ld = 0;
#ifdef __FFLASFFPACK_USE_SIMD
		size_t chunk = Simd<_Element >::vect_size;
#else
		size_t chunk = 0;
#endif // SiMD
		index_t  * col = nullptr;
		_Element * dat = nullptr;
	};

	template<class _Element, bool _Simd =
#ifdef __FFLASFFPACK_USE_SIMD
	true
#else
	false
#endif // SiMD
	>
	struct ELL_sub : public ELL<_Element, _Simd> {
	};

	template<class _Element, bool _Simd =
#ifdef __FFLASFFPACK_USE_SIMD
	true
#else
	false
#endif // SiMD
	>
	struct ELL_ZO : public ELL<_Element, _Simd> {
		_Element cst = 1;
	};

} // FFLAS

namespace FFLAS { namespace details {

	// y = A x + b y ; (generic)
	template<class Field, bool simd_true>
	inline void fspmv(
			     const Field& F,
			     // const FFLAS_TRANSPOSE tA,
			     const size_t m,
			     const size_t n,
			     const size_t ld,
			     const size_t chunk,
			     const index_t * col,
			     const typename Field::Element *  dat,
			     const typename Field::Element * x ,
			     typename Field::Element * y
			     , FieldCategories::GenericTag
			    )
	{
		// std::cout << "details: generic " << std::endl;
		if(!simd_true)
		{
			for (size_t i = 0 ; i < m ; ++i) {
				// XXX can be delayed
				for (index_t j = 0 ; j < ld ; ++j) {
					F.axpyin(y[i],dat[i*ld+j],x[col[i*ld+j]]);
				}
			}
		}
		else
		{
			size_t end = (m % chunk == 0) ? m : m + m%chunk;
			for ( size_t i = 0 ; i < end/chunk ; ++i ) {
				for (index_t j = 0 ; j < ld ; ++j) {
					for(size_t k = 0 ; k < chunk ; ++k)
						F.axpyin(y[i*chunk+k], dat[i*chunk*ld+j*chunk+k], x[col[i*chunk*ld+j*chunk+k]]);
				}
			}
		}
	}


	template<class Field, bool simd_true>
	inline void fspmv(
			     const Field & F,
			     // const FFLAS_TRANSPOSE tA,
			     const size_t m,
			     const size_t n,
			     const size_t ld,
			     const size_t chunk,
			     const index_t * col,
			     const typename Field::Element *  dat,
			     const typename Field::Element * x ,
			     typename Field::Element * y
			     , FieldCategories::UnparametricTag
			    )
	{
		// std::cout << "details: double" << std::endl;
		if (!simd_true) {
			for ( size_t i = 0 ;  i < m   ; ++i ) {
				for (index_t j = 0 ; j < ld ; ++j) {
					y[i] += dat[i*ld+j]*x[col[i*ld+j]];
				}
			}
		}
		else {
#ifdef __FFLASFFPACK_USE_SIMD
			using simd = Simd<typename Field::Element >;
			using vect_t = typename simd::vect_t;
			size_t end = (m%chunk == 0)? m : m+m%chunk;
			vect_t X,Y,D ;
			for( size_t i = 0 ; i < end/chunk ; ++i ) {

				Y = simd::load(y+i*simd::vect_size);
				for (index_t j = 0 ; j < ld ; ++j) {
					D = simd::load(dat+i*chunk*ld+j*chunk);
					X = simd::gather(x,col+i*chunk*ld+j*chunk);
					Y = simd::fmadd(Y,D,X);
				}
				simd::store(y+i*chunk,Y);
			}
#else
			size_t end = (m%chunk == 0)? m : m+m%chunk;
			for( size_t i = 0 ; i < end/chunk ; ++i ) {
				for (index_t j = 0 ; j < ld ; ++j) {
					for(size_t k = 0 ; k < chunk ; ++k)
					{
						y[i*chunk+k] += dat[i*chunk*ld+j*chunk+k]*x[col[i*chunk*ld+j*chunk+k]];
					}
				}
			}
#endif
		}
	}


	// delayed by kmax
	template<class Field, bool simd_true >
	inline void fspmv(
			     const Field& F,
			     // const FFLAS_TRANSPOSE tA,
			     const size_t m,
			     const size_t n,
			     const size_t ld,
			     const size_t chunk,
			     const index_t * col,
			     const typename Field::Element *  dat,
			     const typename Field::Element * x ,
			     typename Field::Element * y,
			     const index_t & kmax
			     // , FieldCategories::ModularTag
			    )
	{
		index_t block = (ld)/kmax ; // use DIVIDE_INTO from fspmvgpu
		// std::cout << "Field delayed" << std::endl;
		if(!simd_true) {
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
		else {
			size_t end = (m%chunk == 0)? m : m+m%chunk;
#ifdef __FFLASFFPACK_USE_SIMD
			using simd = Simd<typename Field::Element >;
			using vect_t = typename simd::vect_t;
			vect_t X, Y, D, C, Q, TMP;
			double p = (typename Field::Element)F.characteristic();
			vect_t P = simd::set1(p);
			vect_t NEGP = simd::set1(-p);
			vect_t INVP = simd::set1(1./p);
			vect_t MIN = simd::set1(F.minElement());
			vect_t MAX = simd::set1(F.maxElement());


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
					vectorised::VEC_MOD(Y,Y,TMP, P, NEGP,INVP,MIN,MAX);
				}
				for ( ; j < ld ; ++j) {
					D = simd::load(dat+i*chunk*ld+j*chunk);
					X = simd::gather(x,col+i*chunk*ld+j*chunk);
					Y = simd::fmadd(Y,D,X);
				}
				vectorised::VEC_MOD(Y,Q,TMP, P, NEGP,INVP,MIN,MAX);
				simd::store(y+i*chunk,Y);
			}
#else
			for ( size_t i = 0; i < end/chunk ; ++i ) {
				index_t j = 0 ;
				index_t j_loc = 0 ;
				for (size_t l = 0 ; l < block ; ++l) {
					j_loc += kmax ;

					for ( ; j < j_loc ; ++j) {
						for(size_t k = 0 ; k < chunk ; ++k)
						{
							y[i*chunk+k] += dat[i*chunk*ld+j*chunk+k]*x[col[i*chunk*ld+j*chunk+k]];
						}
					}
					for(size_t k = 0 ; k < chunk ; ++k)
					{
						F.init(y[i*chunk+k], y[i*chunk+k]);
					}
				}
				for ( ; j < ld ; ++j) {
					for(size_t k = 0 ; k < chunk ; ++k)
					{
						y[i*chunk+k] += dat[i*chunk*ld+j*chunk+k]*x[col[i*chunk*ld+j*chunk+k]];
					}
				}
				for(size_t k = 0 ; k < chunk ; ++k)
				{
					F.init(y[i*chunk+k], y[i*chunk+k]);
				}
			}
#endif

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
	template<class Field, bool simd_true >
	inline void fspmv(
			     const Field& F,
			     // const FFLAS_TRANSPOSE tA,
			     const ELL_sub<typename Field::Element, simd_true> & A,
			     const VECT<typename Field::Element > & x,
			     const typename Field::Element & b,
			     VECT<typename Field::Element > & y
			    )
	{
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::value());
		fspmv( F, A, x, y, typename FieldTraits<Field>::category());
	}

	template<class Field, bool simd_true>
	inline void fspmv(const Field & F,
			     const ELL_sub<typename Field::Element, simd_true> & A,
			     const VECT<typename Field::Element > & x,
			     VECT<typename Field::Element > & y,
			     FieldCategories::ModularTag)
	{
		details::fspmv<Field,simd_true>(F,A.m,A.n,A.ld,A.chunk,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag ());
		finit(F,A.m,y.dat,1);
	}

	template<class Field, bool simd_true>
	inline void fspmv(const Field & F,
			     const ELL_sub<typename Field::Element, simd_true> & A,
			     const VECT<typename Field::Element > & x,
			     VECT<typename Field::Element > & y,
			     FieldCategories::UnparametricTag)
	{
		details::fspmv<Field,simd_true>(F,A.m,A.n,A.ld,A.chunk,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag() );
	}

	template<class Field, bool simd_true>
	inline void fspmv(const Field & F,
			     const ELL_sub<typename Field::Element, simd_true> & A,
			     const VECT<typename Field::Element > & x,
			     VECT<typename Field::Element > & y,
			     FieldCategories::GenericTag)
	{
		details::fspmv<Field,simd_true>(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,y.dat, FieldCategories::GenericTag() );
	}

	/* *** */
	/* ELL */
	/* *** */

	// y = A x + b y ; (generic)
	// reductions are delayed.
	template<class Field, bool simd_true>
	void fspmv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<typename Field::Element, simd_true> & A,
		      const VECT<typename Field::Element > & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element > & y
		     )
	{
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::value());
		fspmv(F,A,x,y, typename FieldTraits<Field>::category());
	}

	template<class Field, bool simd_true>
	void fspmv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<typename Field::Element, simd_true> & A,
		      const VECT<typename Field::Element > & x,
		      VECT<typename Field::Element > & y
		      , FieldCategories::GenericTag
		     )
	{
		details::fspmv<simd_true>(F,A.m,A.n,A.ld,A.chunk,A.col,A.dat,x.dat,y.dat, FieldCategories::GenericTag());
	}

	template<class Field, bool simd_true>
	void fspmv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<typename Field::Element, simd_true> & A,
		      const VECT<typename Field::Element > & x,
		      VECT<typename Field::Element > & y
		      , FieldCategories::UnparametricTag
		     )
	{
		details::fspmv<simd_true>(F,A.m,A.n,A.ld,A.chunk,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag());
	}

	template<class Field, bool simd_true>
	void fspmv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<typename Field::Element, simd_true> & A,
		      const VECT<typename Field::Element > & x,
		      VECT<typename Field::Element > & y
		      , FieldCategories::ModularTag
		     )
	{
		size_t kmax = Protected::DotProdBoundClassic(F,F.one) ;
		details::fspmv<Field,simd_true>(F,A.m,A.n,A.ld,A.chunk,A.col,A.dat,x.dat,y.dat, (index_t) kmax);
	}


} // FFLAS

namespace FFLAS { namespace details { /*  ZO */

	// generic
	template<class Field, bool add, bool simd_true>
	inline void fspmv_zo(
				const Field & F,
				// const FFLAS_TRANSPOSE tA,
				const size_t m,
				const size_t n,
				const size_t ld,
				const size_t chunk,
				const index_t * col,
				const typename Field::Element * x ,
				typename Field::Element * y
				, FieldCategories::GenericTag
			       )
	{
		if(!simd_true)
		{
			for (size_t i = 0 ; i < m ; ++i) {
				if (add == true) {
					for (index_t j = 0 ; j < ld ; ++j) {
						F.addin(y[i],x[col[i*ld+j]]);
					}
				}
				else {
					for (index_t j = 0 ; j < ld ; ++j) {
						F.subin(y[i],x[col[i*ld+j]]);
					}
				}
			}
		}
		else
		{
			size_t end = (m % chunk == 0) ? m : m + m%chunk;
			for ( size_t i = 0 ; i < end/chunk ; ++i ) {
				for (index_t j = 0 ; j < ld ; ++j) {
					if (add == true) {
						for(size_t k = 0 ; k < chunk ; ++k)
							F.addin(y[i*chunk+k], x[col[i*chunk*ld+j*chunk+k]]);
					}
					else {
						for(size_t k = 0 ; k < chunk ; ++k)
							F.subin(y[i*chunk+k], x[col[i*chunk*ld+j*chunk+k]]);
					}
				}
			}
		}
	}


	template<class Field, bool add, bool simd_true>
	inline void fspmv_zo(
				const Field & F,
				// const FFLAS_TRANSPOSE tA,
				const size_t m,
				const size_t n,
				const size_t ld,
				const size_t chunk,
				const index_t * col,
				const typename Field::Element * x ,
				typename Field::Element * y
				, FieldCategories::UnparametricTag
			       )
	{
		if (!simd_true) {
			for ( size_t i = 0 ;  i < m   ; ++i ) {
				if (add == true) {
					for (index_t j = 0 ; j < ld ; ++j) {
						y[i] += x[col[i*ld+j]];
					}
				}
				else {
					for (index_t j = 0 ; j < ld ; ++j) {
						y[i] += x[col[i*ld+j]];
					}
				}
			}
		}
		else {
#ifdef __FFLASFFPACK_USE_SIMD
			using simd = Simd<typename Field::Element >;
			using vect_t = typename simd::vect_t;
			size_t end = (m%chunk == 0)? m : m+m%chunk;
			vect_t X,Y,D ;
			for( size_t i = 0 ; i < end/chunk ; ++i ) {

				Y = simd::load(y+i*simd::vect_size);
				if (add == true) {
					for (index_t j = 0 ; j < ld ; ++j) {
						X = simd::gather(x,col+i*chunk*ld+j*chunk);
						Y = simd::add(Y,X);
					}
				}
				else{
					for (index_t j = 0 ; j < ld ; ++j) {
						X = simd::gather(x,col+i*chunk*ld+j*chunk);
						Y = simd::sub(Y,X);
					}
				}
				simd::store(y+i*chunk,Y);
			}
#else
			size_t end = (m%chunk == 0)? m : m+m%chunk;
			for( size_t i = 0 ; i < end/chunk ; ++i ) {
				for (index_t j = 0 ; j < ld ; ++j) {
					if (add == true) {
						for(size_t k = 0 ; k < chunk ; ++k) {
							y[i*chunk+k] += x[col[i*chunk*ld+j*chunk+k]];
						}
					}
					else  {
						for(size_t k = 0 ; k < chunk ; ++k) {
							y[i*chunk+k] -= x[col[i*chunk*ld+j*chunk+k]];
						}
					}
				}
			}
#endif
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
	template<class Field, bool simd_true >
	inline void fspmv(
			     const Field & F,
			     // const FFLAS_TRANSPOSE tA,
			     const ELL_ZO<typename Field::Element, simd_true> & A,
			     const VECT<typename Field::Element > & x,
			     const typename Field::Element & b,
			     VECT<typename Field::Element > & y
			    )
	{
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::value());
		fspmv(F, A, x, y, typename FieldTraits<Field>::category() );
	}

	template<class Field, bool simd_true>
	inline void fspmv(
			     const Field & F,
			     // const FFLAS_TRANSPOSE tA,
			     const ELL_ZO<typename Field::Element, simd_true > & A,
			     const VECT<typename Field::Element > & x,
			     VECT<typename Field::Element > & y
			     , FieldCategories::GenericTag
			    )
	{
		if (A.cst == F.one) {
			details::fspmv<Field,true,simd_true>(F,A.m,A.n,A.ld,A.col,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else if (A.cst == F.mOne) {
			details::fspmv<Field,false,simd_true>(F,A.m,A.n,A.ld,A.col,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else {
			typename Field::Element * xd = fflas_new<typename Field::Element >(A.n) ;
			fscal(F,A.n,A.cst,x.dat,1,xd,1);
			details::fspmv<Field,true,simd_true>(F,A.m,A.n,A.ld,A.col,x.dat,y.dat, FieldCategories::GenericTag());
		}
	}

	template<class Field, bool simd_true>
	inline void fspmv(
			     const Field & F,
			     // const FFLAS_TRANSPOSE tA,
			     const ELL_ZO<typename Field::Element, simd_true > & A,
			     const VECT<typename Field::Element > & x,
			     VECT<typename Field::Element > & y
			     , FieldCategories::UnparametricTag
			    )
	{
		if (A.cst == F.one) {
			details::fspmv<Field,true,simd_true>(F,A.m,A.n,A.ld,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else if (A.cst == F.mOne) {
			details::fspmv<Field,false,simd_true>(F,A.m,A.n,A.ld,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else {
			typename Field::Element * xd = fflas_new<typename Field::Element >(A.n) ;
			fscal(F,A.n,A.cst,x.dat,1,xd,1);
			details::fspmv<Field,true,simd_true>(F,A.m,A.n,A.ld,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
	}

	template<class Field, bool simd_true>
	inline void fspmv(
			     const Field & F,
			     // const FFLAS_TRANSPOSE tA,
			     const ELL_ZO<typename Field::Element, simd_true > & A,
			     const VECT<typename Field::Element > & x,
			     VECT<typename Field::Element > & y
			     , FieldCategories::ModularTag
			    )
	{
		if (A.cst == F.one) {
			details::fspmv<Field,true,simd_true>(F,A.m,A.n,A.ld,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else if (A.cst == F.mOne) {
			details::fspmv<Field,false,simd_true>(F,A.m,A.n,A.ld,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else {
			typename Field::Element * xd = fflas_new<typename Field::Element >(A.n) ;
			fscal(F,A.n,A.cst,x.dat,1,xd,1);
			details::fspmv<Field,true,simd_true>(F,A.m,A.n,A.ld,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		finit(F,y.m,y.dat,1);
	}


} // FFLAS

namespace FFLAS { /*  conversions */

	template<class Field, class ColT, class RowT>
	inline void sp_ell_from_csr(
				    const Field & F,
				    const size_t CSR_m,
				    const size_t CSR_n,
				    const index_t * CSR_col,
				    const index_t * CSR_row,
				    const typename Field::Element * CSR_dat,
				    size_t & ELL_m,
				    size_t & ELL_n,
				    size_t & ld,
				    size_t & chunk,
				    index_t * ELL_col,
				    typename Field::Element * ELL_dat,
				    const bool bSimd,
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
		if(bSimd){
#ifdef __FFLASFFPACK_USE_SIMD
			chunk = Simd<typename Field::Element >::vect_size;
			size_t m = (CSR_m%chunk == 0) ? CSR_m : CSR_m+CSR_m%chunk;

			ELL_col = fflas_new<index_t >(ld*m, Alignment::CACHE_LINE);
			if(!ZO){
				ELL_dat = fflas_new<typename Field::Element >(ld*m, Alignment::CACHE_LINE);
			}

			size_t i = 0;
			for( ; i < m ; i+= chunk){
				for(size_t k = 0 ; k < chunk ; ++k){
					size_t start = CSR_row[i*chunk+k], stop = CSR_row[i*chunk+k+1];
					for(size_t j = 0 ; j < ld ; ++j){
						if(start + j  < stop){
							if(!ZO){
								ELL_dat[i*chunk*ld+j*chunk+k] = CSR_dat[start+j];
							}
							ELL_col[i*chunk*ld+j*chunk+k] = CSR_col[start+j];
						}
						else{
							if(!ZO){
								F.init(ELL_dat[i*chunk*ld+j*chunk+k]);
							}
							ELL_col[i*chunk*ld+j*chunk+k] = 0;
						}
					}
				}
			}

			for(size_t k = 0 ; k < chunk ; ++k){
				if(i + k < CSR_m){
					size_t start = CSR_row[i*chunk+k], stop = CSR_row[i*chunk+k+1];
					for(size_t j = 0 ; j < ld ; ++j){
						if(start + j  < stop){
							if(!ZO){
								ELL_dat[i*chunk*ld+j*chunk+k] = CSR_dat[start+j];
							}
							ELL_col[i*chunk*ld+j*chunk+k] = CSR_col[start+j];
						}
						else{
							if(!ZO){
								F.init(ELL_dat[i*chunk*ld+j*chunk+k]);
							}
							ELL_col[i*chunk*ld+j*chunk+k] = 0;
						}
					}
				}
				else{
					for(size_t j = 0 ; j < ld ; ++j){
						if(!ZO){
							F.init(ELL_dat[i*chunk*ld+j*chunk+k]);
						}
						ELL_col[i*chunk*ld+j*chunk+k] = 0;
					}
				}
			}
#else
			FFLASFFPACK_abort("you should have SIMD...");
#endif
		}else{

			ELL_col = fflas_new<index_t >(ld*ELL_m, Alignment::CACHE_LINE);
			if(!ZO){
				ELL_dat = fflas_new<typename Field::Element >(ld*ELL_m, Alignment::CACHE_LINE);
			}

			for(size_t i = 0 ; i < CSR_m ; ++i){
				size_t start = CSR_row[i], stop = CSR_row[i+1];
				for(size_t j = 0 ; j < ld ; ++j){
					if(start + j < stop){
						if(!ZO){
							ELL_dat[i*ld+j] = CSR_dat[i*ld+j];
						}
						ELL_col[i*ld+j] = CSR_col[i*ld+j];
					}
					else{
						if(!ZO){
							F.init(ELL_dat[i*ld+j]);
						}
						ELL_col[i*ld+j] = 0;
					}
				}
			}
		}
	}

	template<class Field, class ColT, class RowT>
	inline void sp_ell_from_coo(
				    const Field & F,
				    const size_t COO_m,
				    const size_t COO_n,
				    const index_t * COO_col,
				    const index_t * COO_row,
				    const typename Field::Element * COO_dat,
				    size_t & ELL_m,
				    size_t & ELL_n,
				    size_t & ld,
				    size_t & chunk,
				    index_t * ELL_col,
				    typename Field::Element * ELL_dat,
				    const bool bSimd,
				    const bool ZO
				   )
	{
		index_t * row = fflas_new<index_t >(COO_m+1);
		for(size_t i = 0 ; i < COO_m+1 ; ++i){
			row[i] = 0;
		}
		for(size_t i = 0 ; i < COO_m ; ++i){
			row[COO_row[i]]++;
		}
		for(size_t i = 1 ; i < COO_m+1 ; ++i){
			row[i+1]+=row[i];
		}
		sp_ell_from_csr(F, COO_m, COO_n, COO_col, row, COO_dat, ELL_m, ELL_n, ld, chunk, ELL_col, ELL_dat, bSimd, ZO);
		fflas_delete(row);
	}

} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_spmv_ell_INL
