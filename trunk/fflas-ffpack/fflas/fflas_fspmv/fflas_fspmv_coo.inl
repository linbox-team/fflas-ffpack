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

	template<class Element>
	struct COO {
		size_t m ;
		size_t n ;
		size_t z ;
		index_t  * row  ;
		index_t  * col ;
		Element * dat ;
		// int mc ;
		// int ml ;
	};

	template<class Element>
	struct COO_sub : public COO<Element> {
		// size_t i0 ;
		// size_t j0 ;
	};

	template<class Element>
	struct COO_ZO {
		size_t m ;
		size_t n ;
		size_t z ;
		index_t  * row  ;
		index_t  * col ;
		// Element * dat ;
		Element cst ;
		// size_t i0 ;
		// size_t j0 ;
	};

	// y = A x + b y ; (generic)
	namespace details {
		template<class Field>
		void sp_fgemv(
			      const Field& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t z,
			      const index_t * row,
			      const index_t * col,
			      const typename Field::Element *  dat,
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
					// XXX can be delayed
					for (index_t j = 0 ; j < z ; ++j)
						F.axpyin(y[row[i]],dat[j],x[col[j]]);
				}
#if 0
			}
			else {
			}
#endif
		}

		// Double
		template<>
		void sp_fgemv(
			      const DoubleDomain & F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t z,
			      const index_t * row,
			      const index_t * col,
			      const double *  dat,
			      const double * x ,
			      const double & b,
			      double * y
			     )
		{
			// std::cout << m << 'x' << n << std::endl;
			// std::cout << "y : " ; for (size_t i = 0 ; i < m ; ++i) std::cout << y[i] << ' '  ; std::cout << std::endl;
			// std::cout << "x : " ; for (size_t i = 0 ; i < n ; ++i) std::cout << x[i] << ' '  ; std::cout << std::endl;
			// std::cout << "row : " ; for (size_t i = 0 ; i < m+1 ; ++i) std::cout << row[i] << ' '  ; std::cout << std::endl;
			// std::cout << "col : " ; for (size_t i = 0 ; i < row[m] ; ++i) std::cout << col[i] << ' '  ; std::cout << std::endl;
			// std::cout << "dat : " ; for (size_t i = 0 ; i < row[m] ; ++i) std::cout << dat[i] << ' '  ; std::cout << std::endl;


			// std::cout << "MKL ?" << std::endl;
#ifdef __FFLASFFPACK_HAVE_MKL
			// std::cout << "MKL" << std::endl;
			// fscalin(F,m,b,y,1);

			// char * transa = (ta==FflasNoTrans)?'n':'t';
			char   transa = 'N';
			index_t m_ = (index_t) m ;
			index_t z_ = (index_t) z ;
			// index_t n_ = n ;
			double * yd ;
			if ( b == 0) {
				yd = y;
			}
			else {
				yd = new double [m];
				// std::cout << "yd : " ; for (size_t i = 0 ; i < m ; ++i) std::cout << yd[i] << ' '  ; std::cout << std::endl;
				fscalin(F,m,b,y,1);
			}
			// mkl_dcoogemv (bug too ?)
			mkl_cspblas_dcoogemv
			(&transa, &m_, const_cast<double*>(dat), const_cast<index_t*>(row) , const_cast<index_t*>(col),
			 &z_, const_cast<double*>(x), yd);
			// std::cout << "yd : " ; for (size_t i = 0 ; i < m ; ++i) std::cout << yd[i] << ' '  ; std::cout << std::endl;
			// std::cout << "y : " ; for (size_t i = 0 ; i < m ; ++i) std::cout << y[i] << ' '  ; std::cout << std::endl;
			// std::cout << "x : " ; for (size_t i = 0 ; i < n ; ++i) std::cout << x[i] << ' '  ; std::cout << std::endl;

			if ( b != 0) {
				faddin(F,m,yd,1,y,1);
				delete[] yd ;
			}
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
				for (index_t j = 0 ; j < (index_t)z ; ++j)
					y[row[i]] += dat[j] * x[col[j]];
			}
#endif // __FFLASFFPACK_HAVE_MKL
		}


		// Float
		template<>
		void sp_fgemv(
			      const FloatDomain & F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t z,
			      const index_t * row,
			      const index_t * col,
			      const float *  dat,
			      const float * x ,
			      const float & b,
			      float * y
			     )
		{
#ifdef __FFLASFFPACK_HAVE_MKL
			// char * transa = (ta==FflasNoTrans)?'n':'t';
			char   transa = 'n';
			index_t m_ = (index_t) m ;
			index_t z_ = (index_t) z ;
			float * yd ;
			if ( b == 0) {
				yd = y;
			}
			else {
				yd = new float [m];
				fscalin(F,m,b,y,1);
			}
			// mkl_scoogemv
			mkl_cspblas_scoogemv
			(&transa, &m_, const_cast<float*>(dat), const_cast<index_t*>(row), const_cast<index_t*>(col),
			 &z_, const_cast<float*>(x), yd);
			if ( b != 0) {
				faddin(F,m,yd,1,y,1);
				delete[] yd ;
			}
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
				for (index_t j = 0 ; j < (index_t) z ; ++j)
					y[row[i]] += dat[j] * x[col[j]];
			}
#endif // __FFLASFFPACK_HAVE_MKL
		}






		// delayed by kmax
		template<class Field>
		void sp_fgemv(
			      const Field& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t z,
			      const index_t * row,
			      const index_t * col,
			      const typename Field::Element *  dat,
			      const typename Field::Element * x ,
			      // const typename Field::Element & b,
			      typename Field::Element * y,
			      const index_t & kmax
			     )
		{
			// XXX bug, do it as in linbox.
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

		// generic
		template<class Field, bool add>
		void sp_fgemv_zo(
				 const Field & F,
				 // const FFLAS_TRANSPOSE tA,
				 const size_t m,
				 const size_t n,
				 const size_t z,
				 const index_t * row,
				 const index_t * col,
				 const typename Field::Element * x ,
				 // const typename Field::Element & b,
				 typename Field::Element * y
				)
		{
			if (add == true) {
				for (index_t j = 0 ; j < z ; ++j)
					F.addin(y[row[j]], x[col[j]]);
			}
			else {
				for (index_t j = 0 ; j < z ; ++j)
					F.subin(y[row[j]], x[col[j]]);
			}
		}

		// Double
		template<bool add>
		void sp_fgemv_zo(
				 const DoubleDomain & ,
				 // const FFLAS_TRANSPOSE tA,
				 const size_t m,
				 const size_t n,
				 const size_t z,
				 const index_t * row,
				 const index_t * col,
				 const double * x ,
				 // const double & b,
				 double * y
				)
		{
			if (add == true) {
				for (index_t j = 0 ; j < z ; ++j)
					y[row[j]] +=  x[col[j]];
			}
			else
			{
				for (index_t j = 0 ; j < z ; ++j)
					y[row[j]] -=  x[col[j]];
			}
		}

		// Float
		template<bool add>
		void sp_fgemv_zo(
				 const FloatDomain & ,
				 // const FFLAS_TRANSPOSE tA,
				 const size_t m,
				 const size_t n,
				 const size_t z,
				 const index_t * row,
				 const index_t * col,
				 const float * x ,
				 const float & b,
				 float * y
				)
		{
			if (add == true) {
				for (index_t j = 0 ; j < z ; ++j)
					y[row[j]] +=  x[col[j]];
			}
			else
			{
				for (index_t j = 0 ; j < z ; ++j)
					y[row[j]] -=  x[col[j]];
			}
		}



	} // details

	// y = A x + b y ; (generic)
	// in COO_sub, i0, j0 is an offset in the original vectors x and y
	// it is supposed that no reduction is needed.
	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO_sub<typename Field::Element> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.z,A.row,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const DoubleDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO_sub<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.z,A.row,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const FloatDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO_sub<float> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.z,A.row,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::Modular<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO_sub<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		// std::cout << "there" << std::endl;
		sp_fgemv(DoubleDomain(),A,x,b,y);
		finit(F,A.m,y.dat,1);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::ModularBalanced<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO_sub<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		sp_fgemv(DoubleDomain(),A,x,b,y);
		finit(F,A.m,y.dat,1);

	}

	template<>
	void sp_fgemv(
		      const FFPACK::Modular<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO_sub<float> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		sp_fgemv(FloatDomain(),A,x,b,y);
		finit(F,A.m,y.dat,1);


	}

	template<>
	void sp_fgemv(
		      const FFPACK::ModularBalanced<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO_sub<float> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		sp_fgemv(FloatDomain(),A,x,b,y);
		finit(F,A.m,y.dat,1);

	}


	// y = A x + b y ; (generic)
	// reductions are delayed.
	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO<typename Field::Element> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.z,A.row,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const DoubleDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.z,A.row,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const FloatDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO<float> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.z,A.row,A.col,A.dat,x.dat,b,y.dat);
	}



	template<>
	void sp_fgemv(
		      const FFPACK::Modular<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		// std::cout << "here" << std::endl;
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasDouble) ;

		details::sp_fgemv(F,A.m,A.n,A.z,A.row,A.col,A.dat,x.dat,y.dat,(index_t) kmax);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::ModularBalanced<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasDouble) ;

		details::sp_fgemv(F,A.m,A.n,A.z,A.row,A.col,A.dat,x.dat,y.dat,(index_t) kmax);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::Modular<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO<float> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasFloat) ;

		details::sp_fgemv(F,A.m,A.n,A.z,A.row,A.col,A.dat,x.dat,y.dat,(index_t)kmax);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::ModularBalanced<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO<float> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasFloat) ;

		details::sp_fgemv(F,A.m,A.n,A.z,A.row,A.col,A.dat,x.dat,y.dat,(index_t) kmax);
	}


	// this is the cst data special case.
	// Viewed as a submatrix.
	// it is assumed that no reduction is needed while adding.
	template<class Field>
	void sp_fgemv(
		      const Field & F,
		      // const FFLAS_TRANSPOSE tA,
		      COO_ZO<typename Field::Element> & A,
		      const VECT<typename Field::Element > & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element > & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);

		FFLASFFPACK_check(!F.isZero(A.cst));

		if (A.cst == F.one) {
			details::sp_fgemv_zo<Field,true>(F,A.m,A.n,A.z,A.row,A.col,x.dat,y.dat);
		}
		else if (A.cst == F.mOne) {
			details::sp_fgemv_zo<Field,false>(F,A.m,A.n,A.z,A.row,A.col,x.dat,y.dat);
		}
		else {
			typename Field::Element * xd = new typename Field::Element [A.n] ;
			fscal(F,A.n,A.cst,x.dat,1,xd,1);
			details::sp_fgemv_zo<Field,true>(F,A.m,A.n,A.z,A.row,A.col,xd,y.dat);
		}

		finit(F,A.m,y.dat,1);
	}

} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_spmv_coo_INL

