/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_faxpy.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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

#ifndef __FFLASFFPACK_faxpy_INL
#define __FFLASFFPACK_faxpy_INL



namespace FFLAS {

template<class Field>
inline void
faxpy( const Field& F, const size_t N,
		      const typename Field::Element a,
		      typename Field::ConstElement_ptr X, const size_t incX,
		      typename Field::Element_ptr Y, const size_t incY )
{

	if (F.isZero(a))
		return ;

	if (F.isOne(a))
		return fcopy(F,N,X,incX,Y,incY);

	if (F.isMOne(a))
		return fneg(F,N,X,incX,Y,incY);

	typename Field::ConstElement_ptr Xi = X;
	typename Field::Element_ptr Yi=Y;
	for (; Xi < X+N*incX; Xi+=incX, Yi+=incY )
		F.axpyin( *Yi, a, *Xi );
}

template<>
inline void
faxpy( const DoubleDomain& , const size_t N,
       const DoubleDomain::Element a,
       DoubleDomain::ConstElement_ptr x, const size_t incx,
       DoubleDomain::Element_ptr y, const size_t incy )
{

	cblas_daxpy( (int)N, a, x, (int)incx, y, (int)incy);
}

template<>
inline void
faxpy( const FloatDomain& , const size_t N,
       const FloatDomain::Element a,
       FloatDomain::ConstElement_ptr x, const size_t incx,
       FloatDomain::Element_ptr y, const size_t incy )
{

	cblas_saxpy( (int)N, a, x, (int)incx, y, (int)incy);
}

template<class Field>
inline void
faxpy( const Field& F, const size_t m, const size_t n,
		      const typename Field::Element a,
		      typename Field::ConstElement_ptr X, const size_t ldX,
		      typename Field::Element_ptr Y, const size_t ldY )
{

	if (F.isZero(a))
		return ;

	if (F.isOne(a))
		return fcopy(F,m,n,X,ldX,Y,ldY);

	if (F.isMOne(a))
		return fneg(F,m,n,X,ldX,Y,ldY);

	if (n == ldX && n == ldY)
		return faxpy(F,m*n,a,X,1,Y,1);

	typename Field::ConstElement_ptr Xi = X;
	typename Field::Element_ptr Yi=Y;
	for (; Xi < X+m*ldX; Xi+=ldX, Yi+=ldY )
		faxpy(F,n,a,Xi,1,Yi,1);
}

} // FFLAS

#endif // __FFLASFFPACK_faxpy_INL
