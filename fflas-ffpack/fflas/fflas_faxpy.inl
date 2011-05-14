/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_faxpy.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __FFLAFLAS_faxpy_INL
#define __FFLAFLAS_faxpy_INL



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

#endif // __FFLAFLAS_faxpy_INL
