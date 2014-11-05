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

namespace FFLAS { /*  COO */

	template<class Field>
	struct COO {
		index_t m  = 0;
		index_t n  = 0;
		uint64_t z = 0;
		index_t maxrow = 0;
		index_t  * row  ;
		index_t  * col ;
		typename Field::Element_ptr dat;
	};

	template<class Field>
	struct COO_sub : public COO<Field> {
	};

	template<class Field>
	struct COO_ZO : public COO<Field >{
		typename Field::Element cst = 1;
	};

} // FFLAS

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

	// Double
	template<>
	void fspmv(
		      const DoubleDomain & F,
		      const index_t m,
		      const index_t n,
		      const uint64_t z,
		      const index_t * row,
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
		      const FloatDomain & F,
		      const index_t m,
		      const index_t n,
		      const uint64_t z,
		      const index_t * row,
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
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::value());
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
		finit(F,A.m,y.dat,1);
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
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::value());
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
				, FieldCategories::UnparametricTag
			       )
	{
		if (add){
                    std::cout << "add" << std::endl;
			for (index_t j = 0 ; j < z ; ++j)
				y[row[j]] +=  x[col[j]];
		}else{
			for (index_t j = 0 ; j < z ; ++j)
				y[row[j]] -=  x[col[j]];
		}
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
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::value());
                std::cout << "init_y ok" << std::endl;
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
            std::cout << "Modular ok" << std::endl;
		if (A.cst == F.one) {
                    std::cout << "one ok" << std::endl;
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
                std::cout << "spmv ok" << std::endl;
		finit(F,y.m, y.dat,1);
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

namespace FFLAS { /*  conversions */
} // FFLAS





#endif // __FFLASFFPACK_fflas_fflas_spmv_coo_INL
