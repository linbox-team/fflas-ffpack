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

namespace FFLAS { /*  CSR */

	template<class Field>
	struct CSR {
		size_t m ;
		size_t n ;
		index_t  * st  ;
		index_t  * col ;
		typename Field::Element_ptr dat ;
	};

	template<class Field>
	struct CSR_sub : public CSR<Field> {
	};

	template<class Field>
	struct CSR_ZO : public CSR<Field> {
		typename Field::Element cst = 1;
	};

} // FFLAS

namespace FFLAS { namespace details {

	// y = A x + b y ; (generic)
	template<class Field>
	inline void fspmv(
			     const Field& F,
			     const size_t m,
			     const size_t n,
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

	template<class Field>
	void fspmv(
		      const Field & F,
		      const size_t m,
		      const size_t n,
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
		      const DoubleDomain & F,
		      const size_t m,
		      const size_t n,
		      const index_t * st,
		      const index_t * col,
		      const typename DoubleDomain::Element_ptr  dat,
		      const typename DoubleDomain::Element_ptr x ,
		      typename DoubleDomain::Element_ptr y
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
		      const FloatDomain & F,
		      const size_t m,
		      const size_t n,
		      const index_t * st,
		      const index_t * col,
		      const typename FloatDomain::Element_ptr dat,
		      const typename FloatDomain::Element_ptr x ,
		      typename FloatDomain::Element_ptr y
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


	// delayed by kmax
	template<class Field>
	inline void fspmv(
			     const Field& F,
			     const size_t m,
			     const size_t n,
			     const index_t * st,
			     const index_t * col,
			     const typename Field::Element_ptr dat,
			     const typename Field::Element_ptr x ,
			     typename Field::Element * y,
			     const index_t & kmax
			    )
	{
		for (size_t i = 0 ; i < m ; ++i) {
			index_t j = st[i];
			index_t j_loc = j;
			index_t j_end = st[i+1];
			index_t block = (j_end - j_loc)/kmax ;
			for (size_t l = 0 ; l < (size_t) block ; ++l) {
				j_loc += kmax ;
				for ( ; j < j_loc ; ++j) {
					y[i] += dat[j] * x[col[j]];
				}
				F.init(y[i],y[i]);
			}
			for ( ; j < j_end ; ++j) {
				y[i] += dat[j] * x[col[j]];
			}
			F.init(y[i],y[i]);
		}
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
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::value());
		fspmv( F, A, x, y, typename FieldTraits<Field>::category());
	}

	template<class Field, bool simd_true>
	inline void fspmv(const Field & F,
			     const CSR_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::ModularTag)
	{
		details::fspmv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag ());
		finit(F,A.m,y.dat,1);
	}

	template<class Field, bool simd_true>
	inline void fspmv(const Field & F,
			     const CSR_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::UnparametricTag)
	{
		details::fspmv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag ());
	}

	template<class Field, bool simd_true>
	inline void fspmv(const Field & F,
			     const CSR_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::GenericTag)
	{
		details::fspmv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat, FieldCategories::GenericTag ());
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
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::value());
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
		details::fspmv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat, FieldCategories::GenericTag());
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
		details::fspmv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag());
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
		size_t kmax = Protected::DotProdBoundClassic(F,F.one) ;
		details::fspmv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat,(index_t) kmax);
	}


} // FFLAS

namespace FFLAS { namespace details { /*  ZO */

	// generic
	template<class Field, bool add>
	inline void fspmv_zo(
				const Field & F,
				const size_t m,
				const size_t n,
				const index_t * st,
				const index_t * col,
				const typename Field::Element_ptr x ,
				typename Field::Element * y
				, FieldCategories::GenericTag
			       )
	{
		if(add){
			for (size_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					F.addin(y[i], x[col[j]]);
			}
		}else{
			for (size_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					F.subin(y[i], x[col[j]]);
			}
		}
	}

	template<class Field, bool add>
	inline void fspmv_zo(
				const Field & F,
				const size_t m,
				const size_t n,
				const index_t * st,
				const index_t * col,
				const typename Field::Element_ptr x ,
				typename Field::Element_ptr y
				, FieldCategories::UnparametricTag
			       )
	{
		if(add){
			for (size_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					y[i] +=  x[col[j]];
			}
		}else{
			for (size_t i = 0 ; i < m ; ++i) {
				for (index_t j = st[i] ; j < st[i+1] ; ++j)
					y[i] -=  x[col[j]];
			}
		}
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
		details::init_y(F, y.m, b, y.dat, typename FieldTraits<Field>::value());
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
			details::fspmv_zo<Field,true>(F,A.m,A.n,A.st,A.col,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else if (A.cst == F.mOne) {
			details::fspmv_zo<Field,false>(F,A.m,A.n,A.st,A.col,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else {
			details::fspmv_zo<Field,true>(F,A.m,A.n,A.st,A.col,x.dat,y.dat, FieldCategories::GenericTag());
			fscalin(F,A.n,A.cst,y.dat,1);
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
			details::fspmv_zo<Field,true>(F,A.m,A.n,A.st,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else if (A.cst == F.mOne) {
			details::fspmv_zo<Field,false>(F,A.m,A.n,A.st,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else {
			details::fspmv_zo<Field,true>(F,A.m,A.n,A.st,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
			fscalin(F,A.n,A.cst,y.dat,1);
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
			details::fspmv_zo<Field,true>(F,A.m,A.n,A.st,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else if (A.cst == F.mOne) {
			details::fspmv_zo<Field,false>(F,A.m,A.n,A.st,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else {
			details::fspmv_zo<Field,true>(F,A.m,A.n,A.st,A.col,x.dat,y.dat, FieldCategories::UnparametricTag());
			fscalin(F,A.n,A.cst,y.dat,1);
		}
		finit(F,y.m,y.dat,1);
	}


} // FFLAS

namespace FFLAS { /*  conversions */
} // FFLAS


namespace FFLAS { namespace details {
#ifdef __FFLASFFPACK_HAVE_CUDA
	// Double
	template<>
	void fspmv(
		      const DoubleDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      const size_t m,
		      const size_t n,
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
		size_t nnz = st[m]-st[m-1];
		double one = 1.f;
		status= cusparseDcsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, nnz,
				       &one, descrA, dat_d, st_d, col_d,
				       x_d, &b, y_d);
	}

	// Float
	template<>
	void fspmv(
		      const FloatDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      const size_t m,
		      const size_t n,
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
		size_t nnz = st[m]-st[m-1];
		float one = 1.f;
		status= cusparseScsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, nnz,
				       &one, descrA, dat_d, st_d, col_d,
				       x_d, &b, y_d);

	}

	// need cuda finit code (need nvcc)
#endif // __FFLASFFPACK_HAVE_CUDA


} // details
} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_spmv_csr_INL
