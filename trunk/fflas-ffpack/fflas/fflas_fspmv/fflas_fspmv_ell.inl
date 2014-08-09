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

	template<class Element, size_t Simd>
	struct ELL {
		size_t m ;
		size_t n ;
		size_t  ld ;
		index_t  * col ;
		Element * dat ;
	};

	template<class Element, size_t Simd>
	struct ELL_sub : public ELL<Element, Simd> {
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
				for (index_t j = 0 ; j < ld ; ++j) {
					if (F.isZero(dat[i*ld+j]))
						break;
					F.axpyin(y[i],dat[i*ld+j],x[col[i*ld+j]]);
				}
			}
		}

		template<class Field>
		void sp_fgemv_vec(
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
			using simd = Simd<typename Field::Element>;
			using vect_t = typename simd::vect_t;

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
				for (index_t j = 0 ; j < ld ; ++j) {
					if (F.isZero(dat[i*ld+j]))
						break;
					F.axpyin(y[i],dat[i*ld+j],x[col[i*ld+j]]);
				}
			}
		}

#define BLOCKSIZE 4
		// double
		template<>
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
			      ,bool simd_true
			     )
		{
			if (simd_true) {
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
						// if (dat[i*ld+j] == 0) { break; }
						y[i] += dat[i*ld+j]*x[col[i*ld+j]];
					}
				}
			}
			else {
				size_t i = 0 ;
				using simd = Simd<double>;
				using vect_t = typename simd::vect_t;
				vect_t X,Y,D ;
				for ( ; i < m ; i += BLOCKSIZE) {

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
						D = simd::load(dat+i*BLOCKSIZE*ld+j*BLOCKSIZE);
						X = simd::gather(x,col+i*BLOCKSIZE*ld+j*BLOCKSIZE);
						Y = simd::madd(Y,D,X);
					}
					simd::store(y+i,Y);
				}
				if ( b != 1) {
					if ( b == 0.) {
						for (size_t ii = i*BLOCKSIZE ; ii < m ; ++ii) {
							y[ii] = 0;
						}
					}
					else if ( b == -1 ) {
						for (size_t ii = i*BLOCKSIZE ; ii < m ; ++ii) {
							y[ii]= -y[ii];
						}
					}
					else {
						for (size_t ii = i*BLOCKSIZE ; ii < m ; ++ii) {
							y[ii] = y[ii] * b;
						}
					}
				}
				size_t deplacement = m -i*BLOCKSIZE ;
				for (index_t j = 0 ; j < ld ; ++j) {
					for (size_t ii = 0 ; ii < deplacement ; ++ii) {
						y[i*BLOCKSIZE+ii] += dat[i*BLOCKSIZE*ld+j*deplacement + ii]*x[col[i*BLOCKSIZE*ld+j*deplacement + ii]];
					}
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
				for (index_t j = 0 ; j < ld ; ++j) {
					if (dat[i*ld+j] == 0)
						break;
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
			      const index_t * col,
			      const typename Field::Element *  dat,
			      const typename Field::Element * x ,
			      // const typename Field::Element & b,
			      typename Field::Element * y,
			      const index_t & kmax
			     )
		{

			index_t block = (ld)/kmax ; // use DIVIDE_INTO from fspmvgpu
			for (size_t i = 0 ; i < m ; ++i) {
				// y[i] = 0;
				index_t j = 0;
				index_t j_loc = 0 ;
				bool term = false ;
				for (size_t l = 0 ; l < block ; ++l) {
					j_loc += block ;
					for ( ; j < j_loc ; ++j ) {
						if (dat[i*ld+j] = 0) { break; }
						y[i] += dat[i*ld+j] * x[col[i*ld+j]];
					}
					F.init(y[i],y[i]);
				}
				if (! term ) {
					for ( ; j < ld  ; ++j) {
						if (dat[i*ld+j] = 0) { break; }
						y[i] += dat[i*ld+j] * x[i*ld+col[j]];
					}
					F.init(y[i],y[i]);
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
		      const ELL_sub<typename Field::Element, 1> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const DoubleDomain& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL_sub<double, 1> & A,
		      const VECT<double> & x,
		      const double& b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const FloatDomain& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL_sub<float, 1> & A,
		      const VECT<float> & x,
		      const float& b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const Modular<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL_sub<double, 1> & A,
		      const VECT<double> & x,
		      const double& b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv(DoubleDomain(),A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
		finit(F,A.m,y.dat,1);
	}

	template<>
	void sp_fgemv(
		      const Modular<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL_sub<float, 1> & A,
		      const VECT<float> & x,
		      const float& b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv(FloatDomain(),A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
		finit(F,A.m,y.dat,1);
	}

	template<>
	void sp_fgemv(
		      const ModularBalanced<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL_sub<double, 1> & A,
		      const VECT<double> & x,
		      const double& b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv(DoubleDomain(),A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
		finit(F,A.m,y.dat,1);
	}

	template<>
	void sp_fgemv(
		      const ModularBalanced<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL_sub<float, 1> & A,
		      const VECT<float> & x,
		      const float& b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv(FloatDomain(),A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
		finit(F,A.m,y.dat,1);
	}


	// y = A x + b y ; (generic)
	// reductions are delayed.
	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<typename Field::Element, 1> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const DoubleDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<double> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
	}

	template<>
	void sp_fgemv(
		      const FloatDomain & F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<float, 1> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		details::sp_fgemv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,b,y.dat);
	}



	template<>
	void sp_fgemv(
		      const FFPACK::Modular<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<double, 1> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		// std::cout << "here" << std::endl;
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasDouble) ;

		details::sp_fgemv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,y.dat,(index_t) kmax);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::ModularBalanced<double>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<double, 1> & A,
		      const VECT<double> & x,
		      const double & b,
		      VECT<double> & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasDouble) ;

		details::sp_fgemv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,y.dat,(index_t) kmax);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::Modular<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<float, 1> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasFloat) ;

		details::sp_fgemv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,y.dat,(index_t)kmax);
	}

	template<>
	void sp_fgemv(
		      const FFPACK::ModularBalanced<float>& F,
		      // const FFLAS_TRANSPOSE tA,
		      const ELL<float, 1> & A,
		      const VECT<float> & x,
		      const float & b,
		      VECT<float> & y
		     )
	{
		fscalin(F,A.m,b,y.dat,1);
		size_t kmax = Protected::DotProdBoundClassic(F,F.one,FflasFloat) ;

		details::sp_fgemv(F,A.m,A.n,A.ld,A.col,A.dat,x.dat,y.dat,(index_t) kmax);
	}


} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_spmv_ell_INL

