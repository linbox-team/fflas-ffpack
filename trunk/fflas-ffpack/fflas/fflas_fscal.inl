/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_faxpy.inl
 * Copyright (C) 2014 FFLAS-FFPACK group
 *
 * Written by BB <bbboyer@ncsu.edu>
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

#ifndef __FFLASFFPACK_fscal_INL
#define __FFLASFFPACK_fscal_INL


namespace FFLAS {

	/***************************/
	/*         LEVEL 1         */
	/***************************/


	template<class Field>
	inline void
	fscal( const Field& F, const size_t N,
	       const typename Field::Element a,
	       typename Field::ConstElement_ptr X, const size_t incX,
	       typename Field::Element_ptr Y, const size_t incY )
	{
		if (F.isOne(a)) {
			fcopy(F,N,X,incX,Y,incY);
			return ;
		}

		typename Field::ConstElement_ptr Xi = X;
		typename Field::Element_ptr Yi = Y;
		if (F.areEqual(a,F.mOne)){
			fneg(F,N,X,incX,Y,incY);
			return;
		}

		if (F.isZero(a)){
			fzero(F,N,Y,incY);
			return;
		}

		if (incX == 1 && incY == 1)
			for (size_t i = 0 ; i < N ; ++i)
				F.mul( Y[i], a, X[i] );
		else
		for (; Xi < X+N*incX; Xi+=incX, Yi+=incY )
			F.mul( *Yi, a, *Xi );
	}

	template<class Field>
	inline void
	fscalin (const Field& F, const size_t n, const typename Field::Element a,
	       typename Field::Element_ptr X, const size_t incX)
	{
		if (F.isOne(a))
			return ;

		if (F.isMOne(a)){
			fnegin(F,n,X,incX);
			return;
		}

		if (F.isZero(a)){
			fzero(F,n,X,incX);
			return;
		}

		typename Field::Element_ptr Xi = X ;

		if ( incX == 1)
			for (size_t i = 0 ; i < n ; ++i)
				F.mulin( X[i], a);
		else

			for (; Xi < X+n*incX; Xi+=incX )
				F.mulin( *Xi, a);
	}

	template<>
	inline void
	fscal( const DoubleDomain& , const size_t N,
	       const DoubleDomain::Element a,
	       DoubleDomain::ConstElement_ptr x, const size_t incx,
	       DoubleDomain::Element_ptr y, const size_t incy )
	{
		cblas_dcopy( (int)N, x, (int)incy, y, (int)incy);
		cblas_dscal( (int)N, a, y, (int)incy);
	}

	template<>
	inline void
	fscal( const FloatDomain& , const size_t N,
	       const FloatDomain::Element a,
	       FloatDomain::ConstElement_ptr x, const size_t incx,
	       FloatDomain::Element_ptr y, const size_t incy )
	{
		cblas_scopy( (int)N, x, (int)incy, y, (int)incy);
		cblas_sscal( (int)N, a, y, (int)incy);
	}

	template<>
	inline void
	fscalin( const DoubleDomain& , const size_t N,
	       const DoubleDomain::Element a,
	       DoubleDomain::Element_ptr y, const size_t incy )
	{

		cblas_dscal( (int)N, a, y, (int)incy);
	}

	template<>
	inline void
	fscalin( const FloatDomain& , const size_t N,
	       const FloatDomain::Element a,
	       FloatDomain::Element_ptr y, const size_t incy )
	{

		cblas_sscal( (int)N, a, y, (int)incy);
	}


	template<>
	inline void
	fscalin( const FFPACK:: Modular<float>& F , const size_t N,
	       const float a,
	       float * X, const size_t incX )
	{
		if(incX == 1) {
			float p, invp;
			p=(float)F.cardinality();
			invp=a/p;
			vectorised::scalp(X,a,X,N,p,invp,0,p-1);
		}
		else {
			float * Xi = X ;
			for (; Xi < X+N*incX; Xi+=incX )
				F.mulin( *Xi , a);

		}
	}

	template<>
	inline void
	fscalin( const FFPACK:: ModularBalanced<float>& F , const size_t N,
	       const float a,
	       float * X, const size_t incX )
	{
		if(incX == 1) {
			float p, invp;
			p=(float)F.cardinality();
			invp=a/p;
			float pmax = (p-1)/2 ;
			float pmin = pmax-p+1;
			vectorised::scalp(X,a,X,N,p,invp,pmin,pmax);
		}
		else {
			float * Xi = X ;
			for (; Xi < X+N*incX; Xi+=incX )
				F.mulin( *Xi , a);

		}
	}

	template<>
	inline void
	fscalin( const FFPACK:: Modular<double>& F , const size_t N,
	       const double a,
	       double * X, const size_t incX )
	{
		if(incX == 1) {
			double p, invp;
			p=(double)F.cardinality();
			invp=a/p;
			vectorised::scalp(X,a,X,N,p,invp,0,p-1);
		}
		else {
			double * Xi = X ;
			for (; Xi < X+N*incX; Xi+=incX )
				F.mulin( *Xi , a);

		}
	}

	template<>
	inline void
	fscal( const FFPACK:: Modular<double>& F , const size_t N,
	       const double a,
		 const double * X, const size_t incX,
		 double * Y, const size_t incY )
	{
		if(incX == 1 && incY==1) {
			double p, invp;
			p=(double)F.cardinality();
			invp=a/p;
			vectorised::scalp(Y,a,X,N,p,invp,0,p-1);
		}
		else {
			const double * Xi = X ;
			double * Yi = Y ;
			for (; Xi < X+N*incX; Xi+=incX,Yi+=incY )
				F.mul(*Yi, *Xi , a);

		}
	}

	template<>
	inline void
	fscalin( const FFPACK:: ModularBalanced<double>& F , const size_t N,
	       const double a,
	       double * X, const size_t incX )
	{
		if(incX == 1) {
			double p, invp;
			p=(double)F.cardinality();
			invp=a/p;
			double pmax = (p-1)/2 ;
			double pmin = pmax-p+1;
			vectorised::scalp(X,a,X,N,p,invp,pmin,pmax);
		}
		else {
			double * Xi = X ;
			for (; Xi < X+N*incX; Xi+=incX )
				F.mulin( *Xi , a);

		}
	}


	/***************************/
	/*         LEVEL 2         */
	/***************************/



	template<class Field>
	void
	fscalin (const Field& F, const size_t m , const size_t n,
	       const typename Field::Element a,
	       typename Field::Element_ptr A, const size_t lda)
	{
		if (F.isOne(a)) {
			return ;
		}
		else if (F.isZero(a)) {
			fzero(F,m,n,A,lda);
		}
		else if (F.isMOne(a)) {
			fnegin(F,m,n,A,lda);
		}
		else {
			if (lda == n) {
				fscalin(F,n*m,a,A,1);
			}
			else {
				for (size_t i = 0 ; i < m ; ++i)
					fscalin(F,n,a,A+i*lda,1);
			}

			return;
		}
	}

	template<class Field>
	void
	fscal (const Field& F, const size_t m , const size_t n,
	       const typename Field::Element a,
	       typename Field::ConstElement_ptr A, const size_t lda,
	       typename Field::Element_ptr B, const size_t ldb)
	{
		if (F.isOne(a)) {
			fcopy(F,m,n,A,lda,B,ldb) ;
		}
		else if (F.isZero(a)) {
			fzero(F,m,n,B,ldb);
		}
		else if (F.isMOne(a)) {
			fneg(F,m,n,A,lda,B,ldb);
		}
		else {
			if (n == lda && m == lda)
				fscal(F,m*n,a,A,lda,B,ldb);
			else {
				for (size_t i = 0; i < m ; ++i)
					fscal(F,n,a,A+i*lda,1,B+i*ldb,1);
			}
		}

		return;
	}

} // FFLAS

#endif // __FFLASFFPACK_fscal_INL
