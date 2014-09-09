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

#include "fflas-ffpack/utils/simd.h"

namespace FFLAS { /*  ELL */

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
		size_t  ld = 0;
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

	namespace details {

		// y = A x + b y ; (generic)
		template<class Field, bool simd_true>
		void sp_fgemv(
			      const Field& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t ld,
			      const index_t * col,
			      const typename Field::Element *  dat,
			      const typename Field::Element * x ,
			      const typename Field::Element & b,
			      typename Field::Element * y
			     )
		{
			if(!simd_true)
			{
				for (size_t i = 0 ; i < m ; ++i) {
					if (! F.isOne(b)) {
						if (F.isZero(b)) {
							F.assign(y[i],F.zero);
						}
						else if (F.isMOne(b)) {
							F.negin(y[i]);
						}
						else {
							F.mulin(y[i],b);
						}
					}
					// XXX can be delayed
					for (index_t j = 0 ; j < ld ; ++j) {
						F.axpyin(y[i],dat[i*ld+j],x[col[i*ld+j]]);
					}
				}
			}
			else
			{
#ifdef __FFLASFFPACK_USE_SIMD
				std::cout << "a" << std::endl;
				size_t i = 0 ;
				using simd = Simd<typename Field::Element>;
				using vect_t = typename simd::vect_t;
				for ( ; i < m ; i += simd::vect_size) {
					if ( b != 1) {
						if ( b == 0.) {
							for(size_t k = 0 ; k < simd::vect_size ; ++k)
								F.assign(y[i*simd::vect_size+k], F.zero);
						}
						else if ( b == -1 ) {
							for(size_t k = 0 ; k < simd::vect_size ; ++k)
								F.negin(y[i*simd::vect_size+k]);
						}
						else {
							for(size_t k = 0 ; k < simd::vect_size ; ++k)
								F.mulin(y[i*simd::vect_size+k],b);
						}
					}
					for (index_t j = 0 ; j < ld ; ++j) {
						// simd part
						for(size_t k = 0 ; k < simd::vect_size ; ++k)
							F.axpyin(y[i*simd::vect_size+k], dat[i*simd::vect_size*ld+j*simd::vect_size+k], x[col[i*simd::vect_size*ld+j*simd::vect_size+k]]);
					}
				}
				if ( b != 1) {
					if ( b == 0.) {
						for (size_t ii = i*simd::vect_size ; ii < m ; ++ii) {
							F.assign(y[ii], F.zero);
						}
					}
					else if ( b == -1 ) {
						for (size_t ii = i*simd::vect_size ; ii < m ; ++ii) {
							F.negin(y[ii]);
						}
					}
					else {
						for (size_t ii = i*simd::vect_size ; ii < m ; ++ii) {
							F.mulin(y[ii], b);
						}
					}
				}
				size_t deplacement = m -i*simd::vect_size ;
				for (index_t j = 0 ; j < ld ; ++j) {
					for (size_t ii = 0 ; ii < deplacement ; ++ii) {
						F.axpyin(y[i*simd::vect_size+ii], dat[i*simd::vect_size*ld+j*deplacement + ii], x[col[i*simd::vect_size*ld+j*deplacement + ii]]);
					}
				}
#else
				// #error "simd must be enabled"
#endif
			}
		}

		// double
		template<bool simd_true>
		void sp_fgemv(
			      const DoubleDomain& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t ld,
			      const index_t * col,
			      const double*  dat,
			      const double* x ,
			      const double& b,
			      double * y
			     )
		{
			if (!simd_true) {
				for ( size_t i = 0 ;  i < m   ; ++i ) {
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
					for (index_t j = 0 ; j < ld ; ++j) {
						y[i] += dat[i*ld+j]*x[col[i*ld+j]];
					}
				}
			}
			else {
#ifdef __FFLASFFPACK_USE_SIMD
				std::cout << "b" << std::endl;
				size_t i = 0 ;
				using simd = Simd<double>;
				using vect_t = typename simd::vect_t;
				vect_t X,Y,D ;
				for ( ; i < m ; i += simd::vect_size) {

					if ( b != 1) {
						if ( b == 0.) {
							// y[i] = 0;
							Y = simd::zero();
						}
						else if ( b == -1 ) {
							// y[i]= -y[i];
							Y = simd::load(y+i);
							Y = simd::sub(simd::zero(),Y);
						}
						else {
							// y[i] = y[i] * b;
							Y = simd::load(y+i);
							Y = simd::mul(Y,simd::set1(b));
						}
					}
					for (index_t j = 0 ; j < ld ; ++j) {
						D = simd::load(dat+i*simd::vect_size*ld+j*simd::vect_size);
						X = simd::gather(x,col+i*simd::vect_size*ld+j*simd::vect_size);
						Y = simd::madd(Y,D,X);
					}
					simd::store(y+i,Y);
				}
				if ( b != 1) {
					if ( b == 0.) {
						for (size_t ii = i*simd::vect_size ; ii < m ; ++ii) {
							y[ii] = 0;
						}
					}
					else if ( b == -1 ) {
						for (size_t ii = i*simd::vect_size ; ii < m ; ++ii) {
							y[ii]= -y[ii];
						}
					}
					else {
						for (size_t ii = i*simd::vect_size ; ii < m ; ++ii) {
							y[ii] = y[ii] * b;
						}
					}
				}
				size_t deplacement = m -i*simd::vect_size ;
				for (index_t j = 0 ; j < ld ; ++j) {
					for (size_t ii = 0 ; ii < deplacement ; ++ii) {
						y[i*simd::vect_size+ii] += dat[i*simd::vect_size*ld+j*deplacement + ii]*x[col[i*simd::vect_size*ld+j*deplacement + ii]];
					}
				}
#else
				// #error "simd must be enabled"
#endif
			}

		}

		// float
		template<bool simd_true>
		void sp_fgemv(
			      const FloatDomain& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t ld,
			      const index_t * col,
			      const float*  dat,
			      const float* x ,
			      const float& b,
			      float * y
			     )
		{
			if (!simd_true) {
				for ( size_t i = 0 ;  i < m   ; ++i ) {
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
					for (index_t j = 0 ; j < ld ; ++j) {
						y[i] += dat[i*ld+j]*x[col[i*ld+j]];
					}
				}
			}
			else {
#ifdef __FFLASFFPACK_USE_SIMD
				std::cout << "c" << std::endl;
				size_t i = 0 ;
				using simd = Simd<float>;
				using vect_t = typename simd::vect_t;
				vect_t X,Y,D ;
				for ( ; i < m ; i += simd::vect_size) {

					if ( b != 1) {
						if ( b == 0.) {
							// y[i] = 0;
							Y = simd::zero();
						}
						else if ( b == -1 ) {
							// y[i]= -y[i];
							Y = simd::load(y+i);
							Y = simd::sub(simd::zero(),Y);
						}
						else {
							// y[i] = y[i] * b;
							Y = simd::load(y+i);
							Y = simd::mul(Y,simd::set1(b));
						}
					}
					for (index_t j = 0 ; j < ld ; ++j) {
						D = simd::load(dat+i*simd::vect_size*ld+j*simd::vect_size);
						X = simd::gather(x,col+i*simd::vect_size*ld+j*simd::vect_size);
						Y = simd::madd(Y,D,X);
					}
					simd::store(y+i,Y);
				}
				if ( b != 1) {
					if ( b == 0.) {
						for (size_t ii = i*simd::vect_size ; ii < m ; ++ii) {
							y[ii] = 0;
						}
					}
					else if ( b == -1 ) {
						for (size_t ii = i*simd::vect_size ; ii < m ; ++ii) {
							y[ii]= -y[ii];
						}
					}
					else {
						for (size_t ii = i*simd::vect_size ; ii < m ; ++ii) {
							y[ii] = y[ii] * b;
						}
					}
				}
				size_t deplacement = m -i*simd::vect_size ;
				for (index_t j = 0 ; j < ld ; ++j) {
					for (size_t ii = 0 ; ii < deplacement ; ++ii) {
						y[i*simd::vect_size+ii] += dat[i*simd::vect_size*ld+j*deplacement + ii]*x[col[i*simd::vect_size*ld+j*deplacement + ii]];
					}
				}
#else
				// #error "simd must be enabled"
#endif
			}

		}

		// delayed by kmax
		//! @bug check field is M(B)<f|d>
		template<class Field, bool simd_true >
		void sp_fgemv(
			      const Field& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const size_t ld,
			      const index_t * col,
			      const typename Field::Element *  dat,
			      const typename Field::Element * x ,
			      // const typename Field::Element & b,
			      typename Field::Element * y,
			      const index_t & kmax			     )
		{
			std::cout << "1" << std::endl;
			if(!simd_true) {
				index_t block = (ld)/kmax ; // use DIVIDE_INTO from fspmvgpu
				for (size_t i = 0 ; i < m ; ++i) {
					index_t j = 0;
					index_t j_loc = 0 ;
					for (size_t l = 0 ; l < block ; ++l) {
						j_loc += block ;
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
#ifdef __FFLASFFPACK_USE_SIMD
#else
				// #error "simd must be enabled"
#endif

			}
		}
	} // details

	/* ******* */
	/* ELL_sub */
	/* ******* */


	// y = A x + b y ; (generic)
	// it is supposed that no reduction is needed.
	template<class Field, bool simd_true >
	void sp_fgemv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL_sub<typename Field::Element, 1> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     )
	{
		details::sp_fgemv<Field,simd_true>(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
	}

	template< bool simd_true >
	void sp_fgemv(
		      const DoubleDomain& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL_sub<double, simd_true> & A,
		      const VECT<double> & x,
		      const double& b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv<DoubleDomain,simd_true>(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
	}

	template< bool simd_true >
	void sp_fgemv(
		      const FloatDomain& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL_sub<float, simd_true> & A,
		      const VECT<float> & x,
		      const float& b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv<FloatDomain,simd_true>(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
	}

	template< bool simd_true >
	void sp_fgemv(
		      const FFPACK::Modular<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL_sub<double, simd_true> & A,
		      const VECT<double> & x,
		      const double& b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv<DoubleDomain,simd_true>(DoubleDomain(),A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
		finit(F,A.m,y.dat,1);
	}

	template< bool simd_true >
	void sp_fgemv(
		      const FFPACK::Modular<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL_sub<float, simd_true> & A,
		      const VECT<float> & x,
		      const float& b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv<FloatDomain,simd_true>(FloatDomain(),A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
		finit(F,A.m,y.dat,1);
	}

	template< bool simd_true >
	void sp_fgemv(
		      const FFPACK::ModularBalanced<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL_sub<double, simd_true> & A,
		      const VECT<double> & x,
		      const double& b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv<DoubleDomain,simd_true>(DoubleDomain(),A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
		finit(F,A.m,y.dat,1);
	}

	template<bool simd_true>
	void sp_fgemv(
		      const FFPACK::ModularBalanced<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL_sub<float, simd_true> & A,
		      const VECT<float> & x,
		      const float& b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv<FloatDomain,simd_true>(FloatDomain(),A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
		finit(F,A.m,y.dat,1);
	}

	/* *** */
	/* ELL */
	/* *** */

	// y = A x + b y ; (generic)
	// reductions are delayed.
	template<class Field, bool simd_true>
	void sp_fgemv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<typename Field::Element, simd_true> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     )
	{
		details::sp_fgemv<Field,simd_true>(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
	}

	template<bool simd_true>
	void sp_fgemv(
		      const DoubleDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<double,simd_true> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv<DoubleDomain,simd_true>(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
	}

	template<bool simd_true>
	void sp_fgemv(
		      const FloatDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<float, simd_true> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv<FloatDomain,simd_true>(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
	}



	template<bool simd_true>
	void sp_fgemv(
		      const FFPACK::Modular<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<double, simd_true> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		std::cout << "1" << std::endl;
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasDouble) ;

		details::sp_fgemv<FFPACK::Modular<double>,simd_true>(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,y.dat,(index_t) kmax);
	}

	template<bool simd_true>
	void sp_fgemv(
		      const FFPACK::ModularBalanced<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<double, simd_true> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasDouble) ;

		details::sp_fgemv<FFPACK::ModularBalanced<double>,simd_true>(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,y.dat,(index_t) kmax);
	}

	template<bool simd_true>
	void sp_fgemv(
		      const FFPACK::Modular<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<float, simd_true> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasFloat) ;

		details::sp_fgemv<FFPACK::Modular<float>,simd_true>(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,y.dat,(index_t)kmax);
	}

	template<bool simd_true>
	void sp_fgemv(
		      const FFPACK::ModularBalanced<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<float, simd_true> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasFloat) ;

		details::sp_fgemv<FFPACK::ModularBalanced<float>,simd_true>(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,y.dat,(index_t) kmax);
	}


} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_spmv_ell_INL

