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

/** @file fflas/fflas_fspmv_ellr.inl
 * NO DOC
 */

#ifndef __FFLASFFPACK_fflas_fflas_spmv_ellr_INL
#define __FFLASFFPACK_fflas_fflas_spmv_ellr_INL

namespace FFLAS { /*  ELLR */

	template<class Element>
	struct ELLR {
		size_t m ;
		size_t n ;
		size_t  ld ;
		index_t  * row ;
		index_t  * col ;
		Element * dat ;
	};

	template<class Element>
	struct ELLR_sub : public ELLR<Element> {
	};

	template<class Element>
	struct ELLR {
		size_t m ;
		size_t n ;
		size_t  ld ;
		index_t  * row ;
		index_t  * col ;
		Element cst ;
	};


	namespace details {

		// y = A x + b y ; (generic)
		template<class Field>
		void sp_fgemv(
			      const Field& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t ld,
			      const index_t * row,
			      const index_t * col,
			      const typename Field::Element *  dat,
			      const typename Field::Element * x ,
			      const typename Field::Element & b,
			      typename Field::Element * y
			     )
		{
			for (size_t i = 0 ; i < m ; ++i) {
				if (! F.isOne(b)) {
					if (F.isZero(b)) {
						F.assign(y[i],F.zero);
					}
					else if (F.isMone(b)) {
						F.negin(y[i]);
					}
					else {
						F.mulin(y[i],b);
					}
				}
				// XXX can be delayed
				for (index_t j = 0 ; j < row[i] ; ++j) {
					F.axpyin(y[i],dat[i*ld+j],x[col[i*ld+j]]);
				}
			}
		}

		// double
		template<>
		void sp_fgemv(
			      const DoubleDomain& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t ld,
			      const index_t * row,
			      const index_t * col,
			      const double*  dat,
			      const double* x ,
			      const double& b,
			      double * y
			     )
		{
			for (size_t i = 0 ; i < m ; ++i) {
				if ( b != 1) {
					if ( b == 0.) {
						y[i] = 0;
					}
					else if ( b == -1 ) {
						y[i]= -y[i];
					}
					else {
						y[i] = y[i] * b;
					}
				}
				for (index_t j = 0 ; j < row[i] ; ++j) {
					y[i] += dat[i*ld+j]*x[col[i*ld+j]];
				}
			}
		}

		// float
		template<>
		void sp_fgemv(
			      const FloatDomain& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t ld,
			      const index_t * row,
			      const index_t * col,
			      const float*  dat,
			      const float* x ,
			      const float& b,
			      float * y
			     )
		{
			for (size_t i = 0 ; i < m ; ++i) {
				if ( b != 1) {
					if ( b == 0.) {
						y[i] = 0;
					}
					else if ( b == -1 ) {
						y[i]= -y[i];
					}
					else {
						y[i] = y[i] * b;
					}
				}
				for (index_t j = 0 ; j < row[i] ; ++j) {
					y[i] += dat[i*ld+j]*x[col[i*ld+j]];
				}
			}
		}


		// delayed by kmax
		template<class Field>
		void sp_fgemv(
			      const Field& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t ld,
			      const index_t * row,
			      const index_t * col,
			      const typename Field::Element *  dat,
			      const typename Field::Element * x ,
			      // const typename Field::Element & b,
			      typename Field::Element * y,
			      const index_t & kmax
			     )
		{

			for (size_t i = 0 ; i < m ; ++i) {
				index_t block = (row[i])/kmax ; // use DIVIDE_INTO from fspmvgpu
				// y[i] = 0;
				index_t j = 0;
				index_t j_loc = 0 ;
				bool term = false ;
				for (size_t l = 0 ; l < block ; ++l) {
					j_loc += block ;
					for ( ; j < j_loc ; ++j ) {
						y[i] += dat[i*ld+j] * x[col[i*ld+j]];
					}
					F.init(y[i],y[i]);
				}
				for ( ; j < row[i]  ; ++j) {
					y[i] += dat[i*ld+j] * x[i*ld+col[j]];
				}
				F.init(y[i],y[i]);
			}
		}

		// generic
		template<class Field, bool add>
		void sp_fgemv_zo(
				 const Field & F,
				 // const FFLAS_TRANSPOSE tA,
				 const size_t m,
				 const size_t n,
				 const size_t ld,
				 const index_t * row,
				 const index_t * col,
				 const typename Field::Element * x ,
				 // const typename Field::Element & b,
				 typename Field::Element * y
				)
		{
			for (size_t i = 0 ; i < m ; ++i) {
				if (add == true) {
					for (index_t j = 0 ; j < row[i] ; ++j)
						F.addin(y[i], x[col[i*ld+j]]);
				}
				else{
					for (index_t j = 0 ; j < row[i] ; ++j)
						F.subin(y[i], x[col[i*ld+j]]);
				}
				// F.init(y[i],y[i]);
			}
		}

		// Double
		template<bool add>
		void sp_fgemv_zo(
				 const DoubleDomain & ,
				 // const FFLAS_TRANSPOSE tA,
				 const size_t m,
				 const size_t n,
				 const size_t ld,
				 const index_t * row,
				 const index_t * col,
				 const double * x ,
				 // const double & b,
				 double * y
				)
		{
			for (size_t i = 0 ; i < m ; ++i) {
				if (add == true) {
					for (index_t j = 0 ; j < row[i] ; ++j)
						y[i] +=  x[col[i*ld+j]];
				}
				else
				{
					for (index_t j = 0 ; j < row[i] ; ++j)
						y[i] -=  x[col[i*ld+j]];
				}
			}
		}

		// Float
		template<bool add>
		void sp_fgemv_zo(
				 const FloatDomain & ,
				 // const FFLAS_TRANSPOSE tA,
				 const size_t m,
				 const size_t n,
				 const size_t ld,
				 const index_t * row,
				 const index_t * col,
				 const float * x ,
				 const float & b,
				 float * y
				)
		{
			for (size_t i = 0 ; i < m ; ++i) {
				if (add == true) {
					for (index_t j = 0 ; j < row[i] ; ++j)
						y[i] +=  x[col[i*ld+j]];
				}
				else
				{
					for (index_t j = 0 ; j < row[i] ; ++j)
						y[i] -=  x[col[i*ld+j]];
				}
			}
		}

	} // details

	// y = A x + b y ; (generic)
	// it is supposed that no reduction is needed.
	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELLR_sub<typename Field::Element> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const DoubleDomain& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELLR_sub<double> & A,
		      const VECT<double> & x,
		      const double& b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const FloatDomain& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELLR_sub<float> & A,
		      const VECT<float> & x,
		      const float& b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const Modular<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELLR_sub<double> & A,
		      const VECT<double> & x,
		      const double& b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv(DoubleDomain(),A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,b,y.dat);
		finit(F,A.m,y.dat,1);
	}

	template<>
	void sp_fgemv(
		      const Modular<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELLR_sub<float> & A,
		      const VECT<float> & x,
		      const float& b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv(FloatDomain(),A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,b,y.dat);
		finit(F,A.m,y.dat,1);
	}

	template<>
	void sp_fgemv(
		      const ModularBalanced<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELLR_sub<double> & A,
		      const VECT<double> & x,
		      const double& b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv(DoubleDomain(),A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,b,y.dat);
		finit(F,A.m,y.dat,1);
	}

	template<>
	void sp_fgemv(
		      const ModularBalanced<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELLR_sub<float> & A,
		      const VECT<float> & x,
		      const float& b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv(FloatDomain(),A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,b,y.dat);
		finit(F,A.m,y.dat,1);
	}


	// y = A x + b y ; (generic)
	// reductions are delayed.
	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELLR<typename Field::Element> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const DoubleDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELLR<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const FloatDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELLR<float> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,b,y.dat);
	}



	template<>
	void sp_fgemv(
		      const FFPACK::Modular<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELLR<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		// std::cout << "here" << std::endl;
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasDouble) ;

		details::sp_fgemv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,y.dat,(index_t) kmax);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::ModularBalanced<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELLR<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasDouble) ;

		details::sp_fgemv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,y.dat,(index_t) kmax);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::Modular<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELLR<float> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasFloat) ;

		details::sp_fgemv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,y.dat,(index_t)kmax);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::ModularBalanced<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELLR<float> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasFloat) ;

		details::sp_fgemv(F,A.m,A.n,A.ld,A.row,A.col,A.dat,x.dat,y.dat,(index_t) kmax);
	}

	// this is the cst data special case.
	// Viewed as a submatrix.
	// it is assumed that no reduction is needed while adding.
	template<class Field>
	void sp_fgemv(
		      const Field & F,
		      // const FFLAS_TRANSPOSE tA,
		      ELLR_ZO<typename Field::Element> & A,
		      const VECT<typename Field::Element > & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element > & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);

		FFLASFFPACK_check(!F.isZero(A.cst));

		if (A.cst == F.one) {
			details::sp_fgemv_zo<Field,true>(F,A.m,A.n,A.ld,A.row,A.col,x.dat,y.dat);
		}
		else if (A.cst == F.mOne) {
			details::sp_fgemv_zo<Field,false>(F,A.m,A.n,A.ld,A.row,A.col,x.dat,y.dat);
		}
		else {
			typename Field::Element * xd = FFLAS::fflas_new<typename Field::Element >(A.n) ;
			fscal(F,A.n,A.cst,x.dat,1,xd,1);
			details::sp_fgemv_zo<Field,true>(F,A.m,A.n,A.ld,A.row,A.col,xd,y.dat);
		}

		finit(F,A.m,y.dat,1);
	}


} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_spmv_ellr_INL

