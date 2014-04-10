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

	template<class Field>
	inline void
	fscal( const Field& F, const size_t N,
	       const typename Field::Element alpha,
	       const typename Field::Element * X, const size_t incX,
	       typename Field::Element * Y, const size_t incY )
	{
		typedef typename Field::Element Element ;

		if (F.isOne(alpha)) {
			fcopy(F,N,X,incX,Y,incY);
			return ;
		}

		const Element * Xi = X;
		Element * Yi = Y;
		if (F.areEqual(alpha,F.mOne)){
			fneg(F,N,X,incX,Y,incY);
			return;
		}

		if (F.isZero(alpha)){
			fzero(F,N,Y,incY);
			return;
		}

		if (incX == 1 && incY == 1)
			for (size_t i = 0 ; i < N ; ++i)
				F.mul( Y[i], alpha, X[i] );
		else
		for (; Xi < X+N*incX; Xi+=incX, Yi+=incY )
			F.mul( *Yi, alpha, *Xi );
	}

	template<class Field>
	void
	fscalin (const Field& F, const size_t n, const typename Field::Element alpha,
	       typename Field::Element * X, const size_t incX)
	{
		typedef typename Field::Element Element ;

		if (F.isOne(alpha))
			return ;

		if (F.isMOne(alpha)){
			fnegin(F,n,X,incX);
			return;
		}

		if (F.isZero(alpha)){
			fzero(F,n,X,incX);
			return;
		}

		Element * Xi = X ;

		if ( incX == 1)
			for (size_t i = 0 ; i < n ; ++i)
				F.mulin( X[i], alpha );
		else

			for (; Xi < X+n*incX; Xi+=incX )
				F.mulin( *Xi, alpha );
	}

	template<>
	inline void
	fscal( const DoubleDomain& , const size_t N,
	       const DoubleDomain::Element a,
	       const DoubleDomain::Element * x, const size_t incx,
	       DoubleDomain::Element * y, const size_t incy )
	{
		cblas_dcopy( (int)N, x, (int)incy, y, (int)incy);
		cblas_dscal( (int)N, a, y, (int)incy);
	}

	template<>
	inline void
	fscal( const FloatDomain& , const size_t N,
	       const FloatDomain::Element a,
	       const FloatDomain::Element * x, const size_t incx,
	       FloatDomain::Element * y, const size_t incy )
	{
		cblas_scopy( (int)N, x, (int)incy, y, (int)incy);
		cblas_sscal( (int)N, a, y, (int)incy);
	}

	template<>
	inline void
	fscalin( const DoubleDomain& , const size_t N,
	       const DoubleDomain::Element a,
	       DoubleDomain::Element * y, const size_t incy )
	{

		cblas_dscal( (int)N, a, y, (int)incy);
	}

	template<>
	inline void
	fscalin( const FloatDomain& , const size_t N,
	       const FloatDomain::Element a,
	       FloatDomain::Element * y, const size_t incy )
	{

		cblas_sscal( (int)N, a, y, (int)incy);
	}

	template<class Field>
	void
	fscalin (const Field& F, const size_t m , const size_t n,
	       const typename Field::Element alpha,
	       typename Field::Element * A, const size_t lda)
	{
		if (F.isOne(alpha)) {
			return ;
		}
		else if (F.isZero(alpha)) {
			fzero(F,m,n,A,lda);
		}
		else if (F.isMOne(alpha)) {
			fnegin(F,m,n,A,lda);
		}
		else {
			if (lda == n) {
				fscalin(F,n*m,alpha,A,1);
			}
			else {
				for (size_t i = 0 ; i < m ; ++i)
					fscalin(F,n,alpha,A+i*lda,1);
			}

			return;
		}
	}

	template<class Field>
	void
	fscal (const Field& F, const size_t m , const size_t n,
	       const typename Field::Element alpha,
	       const typename Field::Element * A, const size_t lda,
	       typename Field::Element * B, const size_t ldb)
	{
		if (F.isOne(alpha)) {
			fcopy(F,m,n,A,lda,B,ldb) ;
		}
		else if (F.isZero(alpha)) {
			fzero(F,m,n,B,ldb);
		}
		else if (F.isMOne(alpha)) {
			fneg(F,m,n,A,lda,B,ldb);
		}
		else {
			for (size_t i = 0; i < m ; ++i)
				fscal(F,n,alpha,A+i*lda,1,B+i*ldb,1);
		}

		return;
	}

} // FFLAS

#endif // __FFLASFFPACK_fscal_INL
