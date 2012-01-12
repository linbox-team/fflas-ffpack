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
		      const typename Field::Element * X, const size_t incX,
		      typename Field::Element * Y, const size_t incY )
{

	const typename Field::Element * Xi = X;
	typename Field::Element * Yi=Y;
	for (; Xi < X+N*incX; Xi+=incX, Yi+=incY )
		F.axpyin( *Yi, a, *Xi );
}

template<>
inline void
faxpy( const DoubleDomain& , const size_t N,
		      const DoubleDomain::Element a,
		      const DoubleDomain::Element * x, const size_t incx,
		      DoubleDomain::Element * y, const size_t incy )
{

	cblas_daxpy( (int)N, a, x, (int)incx, y, (int)incy);
}

} // FFLAS

#endif // __FFLASFFPACK_faxpy_INL
