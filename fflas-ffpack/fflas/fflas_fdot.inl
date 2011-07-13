/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas_fdot.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __FFLASFFPACK_fdot_INL
#define __FFLASFFPACK_fdot_INL


// Default implementation
// Specializations should be written
// to increase efficiency


namespace FFLAS {

	template<class Field>
	inline typename Field::Element
	fdot( const Field& F, const size_t N,
	      const typename Field::Element * x, const size_t incx,
	      const typename Field::Element * y, const size_t incy )
	{

		typename Field::Element d;
		const typename Field::Element* xi = x;
		const typename Field::Element* yi = y;
		F.init( d, 0 );
		for ( ; xi < x+N*incx; xi+=incx, yi+=incy )
			F.axpyin( d, *xi, *yi );
		return d;
	}

	template<>
	inline DoubleDomain::Element
	fdot( const DoubleDomain& , const size_t N,
	      const DoubleDomain::Element * x, const size_t incx,
	      const DoubleDomain::Element * y, const size_t incy )
	{

		return cblas_ddot( (int)N, x, (int)incx, y, (int)incy );
	}

} // FFLAS

#endif // __FFLASFFPACK_fdot_INL
