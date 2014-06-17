/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   BB <bbboyer@ncsu.edu>
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

/** @file fflas/fflas_sparse_fgemv.inl
*/

#ifndef __FFLASFFPACK_fflas_fflas_sparse_fgemv_INL
#define __FFLASFFPACK_fflas_fflas_sparse_fgemv_INL

#include "fflas-ffpack/field/modular-positive.h"
#include "fflas-ffpack/field/modular-balanced.h"

namespace FFLAS {
	template<class Element>
	struct VECT {
		size_t m ;
		size_t inc ;
		Element * dat ;
	};

	template<class Element>
	struct DNS {

		size_t n ;
		size_t ld ;
		Element * dat ;
	};


	template<class Element>
	struct ELL {
		size_t m ;
		size_t n ;
		size_t  ld ;
		size_t  * col ;
		Element * dat ;
	};

	template<class Element>
	struct ELLR {
		size_t m ;
		size_t n ;
		size_t  ld ;
		size_t  * row ;
		size_t  * col ;
		Element * dat ;
	};

#if 0
	template<class Element>
	struct SPADD {
		size_t ncsr;
		CSR<Element> * csr;
		size_t ncoo;
		COO<Element> * coo;
		size_t ndns;
		DNS<Element> * dns;
		size_t nell;
		ELL<Element> * ell;
		size_t nellr ;
		ELLR<Element> * ellr ;
		size_t ndia ;
		DIA<Element> * dia;

		SPADD() :
			ncsr(0)  ,csr(NULL)
			,ncoo(0) ,coo(NULL)
			,ndns(0) ,dns(NULL)
			,ndia(0) ,dia(NULL)
			,nell(0) ,ell(NULL)
			,nellr(0),ellr(NULL)
		{}
	};
#endif

} // FFLAS

namespace FFLAS { /*  CSR */

	template<class Element>
	struct CSR {
		size_t m ;
		size_t n ;
		size_t  * st  ;
		size_t  * col ;
		Element * dat ;
		// int mc ;
		// int ml ;
	};

	template<class Element>
	struct CSR_sub : public CSR<Element> {
		size_t i0 ;
		size_t j0 ;
	};

	template<class Element>
	struct CSR_ZO {
		size_t m ;
		size_t n ;
		size_t  * st  ;
		size_t  * col ;
		// Element * dat ;
		Element cst ;
		size_t i0 ;
		size_t j0 ;
	};

	// y = A x + b y ; (generic)
	namespace details {
		template<class Field>
		void sp_fgemv(
			      const Field& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t * st,
			      const size_t * col,
			      typename Field::Element *  dat,
			      const typename Field::Element * x ,
			      const typename Field::Element & b,
			      typename Field::Element * y
			     )
		{
#if 0
			if (tA == FflasNoTrans)  {
#endif
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
					for (size_t j = st[i] ; j < st[i+1] ; ++j)
						F.axpyin(y[i],dat[j],x[col[j]]);
				}
#if 0
			}
			else {
				if (F.isZero(b)) {
					for (size_t i = 0 ; i < m ; ++i) {
						F.assign(y[i],F.zero);
					}
				}
				else if (F.isMone(b)) {
					for (size_t i = 0 ; i < m ; ++i) {
						F.negin(y[i]);
					}
				}
				for (size_t i = 0 ; i < m ; ++i) {
					for (size_t j = st[i] ; j < st[i+1] ; ++j)
						F.axpyin(y[col[j]],dat[j],x[j]);
				}
			}
#endif
		}

		// Double
		template<>
		void sp_fgemv(
			      const DoubleDomain & ,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t * st,
			      const size_t * col,
			      double *  dat,
			      const double * x ,
			      const double & b,
			      double * y
			     )
		{

#ifdef __FFLASFFPACK_HAVE_MKL
			fscalin(DoubleDomain(),m,b,y,1);

			// char * transa = (ta==FflasNoTrans)?'n':'t';
			char * transa = 'n';
			mkl_cspblas_dcsrgemv (transa, &m, dat, st , col, x, y);
#else
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
				for (size_t j = st[i] ; j < st[i+1] ; ++j)
					y[i] += dat[j] * x[col[j]];
			}
#endif // __FFLASFFPACK_HAVE_MKL
		}


		// Float
		template<>
		void sp_fgemv(
			      const FloatDomain & ,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t * st,
			      const size_t * col,
			      float *  dat,
			      const float * x ,
			      const float & b,
			      float * y
			     )
		{
#ifdef __FFLASFFPACK_HAVE_MKL
			fscalin(FloatDomain(),m,b,y,1);
			// char * transa = (ta==FflasNoTrans)?'n':'t';
			char * transa = 'n';
			mkl_cspblas_scsrgemv (transa, &m, dat, st , col, x, y);
#else
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
				for (size_t j = st[i] ; j < st[i+1] ; ++j)
					y[i] += dat[j] * x[col[j]];
			}
#endif // __FFLASFFPACK_HAVE_MKL
		}






		template<class Field>
		void sp_fgemv(
			      const Field& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t * st,
			      const size_t * col,
			      typename Field::Element *  dat,
			      const typename Field::Element * x ,
			      // const typename Field::Element & b,
			      typename Field::Element * y,
			      const size_t & kmax
			     )
		{
			for (size_t i = 0 ; i < m ; ++i) {
				// y[i] = 0;
				size_t j = st[i];
				size_t j_loc = j;
				size_t j_end = st[i+1];
				size_t block = (j_end - j_loc)/kmax ;
				for (size_t l = 0 ; l < block ; ++l) {
					j_loc += block ;
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

		// generic
		template<class Field, bool add>
		void sp_fgemv_zo(
				 const Field & F,
				 // const FFLAS_TRANSPOSE tA,
				 const size_t m,
				 const size_t n,
				 const size_t * st,
				 const size_t * col,
				 const typename Field::Element * x ,
				 // const typename Field::Element & b,
				 typename Field::Element * y
				)
		{
			for (size_t i = 0 ; i < m ; ++i) {
				// if ( b != 1) {
				// if ( b == 0.) {
				// F.assign(y[i], F.zero);
				// }
				// else if ( b == -1 ) {
				// F.negin(y[i]);
				// }
				// else {
				// F.mulin(y[i],b);
				// }
				// }
				if (add == true) {
					for (size_t j = st[i] ; j < st[i+1] ; ++j)
						F.addin(y[i], x[col[j]]);
				}
				else{
					for (size_t j = st[i] ; j < st[i+1] ; ++j)
						F.subin(y[i], x[col[j]]);
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
				 const size_t * st,
				 const size_t * col,
				 const double * x ,
				 // const double & b,
				 double * y
				)
		{
			for (size_t i = 0 ; i < m ; ++i) {
				// if ( b != 1) {
				// if ( b == 0.) {
				// y[i] = 0;
				// }
				// else if ( b == -1 ) {
				// y[i]= -y[i];
				// }
				// else {
				// y[i] = y[i] * b;
				// }
				// }
				if (add == true) {
					for (size_t j = st[i] ; j < st[i+1] ; ++j)
						y[i] +=  x[col[j]];
				}
				else
				{
					for (size_t j = st[i] ; j < st[i+1] ; ++j)
						y[i] -=  x[col[j]];
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
				 const size_t * st,
				 const size_t * col,
				 const float * x ,
				 const float & b,
				 float * y
				)
		{
			for (size_t i = 0 ; i < m ; ++i) {
				// if ( b != 1) {
				// if ( b == 0.) {
				// y[i] = 0;
				// }
				// else if ( b == -1 ) {
				// y[i]= -y[i];
				// }
				// else {
				// y[i] = y[i] * b;
				// }
				// }
				if (add == true) {
					for (size_t j = st[i] ; j < st[i+1] ; ++j)
						y[i] +=  x[col[j]];
				}
				else
				{
					for (size_t j = st[i] ; j < st[i+1] ; ++j)
						y[i] -=  x[col[j]];
				}
			}
		}



	} // details

	// y = A x + b y ; (generic)
	// in CSR_sub, i0, j0 is an offset in the original vectors x and y
	// it is supposed that no reduction is needed.
	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      CSR_sub<typename Field::Element> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.st,A.col,A.dat,x.dat+A.j0,b,y.dat+A.i0);
	}

	template<>
	void sp_fgemv(
		      const DoubleDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      CSR_sub<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.st,A.col,A.dat,x.dat+A.j0,b,y.dat+A.i0);
	}

	template<>
	void sp_fgemv(
		      const FloatDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      CSR_sub<float> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.st,A.col,A.dat,x.dat+A.j0,b,y.dat+A.i0);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::Modular<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      CSR_sub<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		sp_fgemv(DoubleDomain(),A,x,b,y);
		finit(F,A.m,y.dat+A.i0,1);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::ModularBalanced<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      CSR_sub<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		sp_fgemv(DoubleDomain(),A,x,b,y);
		finit(F,A.m,y.dat+A.i0,1);

	}

	template<>
	void sp_fgemv(
		      const FFPACK::Modular<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      CSR_sub<float> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		sp_fgemv(FloatDomain(),A,x,b,y);
		finit(F,A.m,y.dat+A.i0,1);


	}

	template<>
	void sp_fgemv(
		      const FFPACK::ModularBalanced<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      CSR_sub<float> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		sp_fgemv(FloatDomain(),A,x,b,y);
		finit(F,A.m,y.dat+A.i0,1);

	}


	// y = A x + b y ; (generic)
	// reductions are delayed.
	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      CSR<typename Field::Element> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const DoubleDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      CSR<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const FloatDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      CSR<float> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,b,y.dat);
	}



	template<>
	void sp_fgemv(
		      const FFPACK::Modular<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      CSR<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		std::cout << "here" << std::endl;
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasDouble) ;

		details::sp_fgemv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat,kmax);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::ModularBalanced<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      CSR<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasDouble) ;

		details::sp_fgemv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat,kmax);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::Modular<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      CSR<float> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasFloat) ;

		details::sp_fgemv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat,kmax);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::ModularBalanced<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      CSR<float> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasFloat) ;

		details::sp_fgemv(F,A.m,A.n,A.st,A.col,A.dat,x.dat,y.dat,kmax);
	}


	// this is the cst data special case.
	// Viewed as a submatrix.
	// it is assumed that no reduction is needed while adding.
	template<class Field>
	void sp_fgemv(
		      const Field & F,
		      // const FFLAS_TRANSPOSE tA,
		      CSR_ZO<typename Field::Element> & A,
		      const VECT<typename Field::Element > & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element > & y
		     )
	{
		fscalin(F,A.m,b,y.dat+A.i0,1);

		assert(!F.isZero(A.cst));

		if (A.cst == F.one) {
			details::sp_fgemv_zo<Field,true>(F,A.m,A.n,A.st,A.col,A.dat,x.dat+A.j0,y.dat+A.i0);
		}
		else if (A.cst == F.mOne) {
			details::sp_fgemv_zo<Field,false>(F,A.m,A.n,A.st,A.col,A.dat,x.dat+A.j0,y.dat+A.i0);
		}
		else {
			typename Field::Element * xd = new typename Field::Element[A.col] ;
			fscal(F,A.col,A.cst,x.dat+A.j0,1,xd,1);
			details::sp_fgemv_zo<Field,true>(F,A.m,A.n,A.st,A.col,A.dat,xd,y.dat+A.i0);
		}

		finit(F,A.m,y.dat+A.i0,1);
	}

} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_sparse_fgemv_INL
