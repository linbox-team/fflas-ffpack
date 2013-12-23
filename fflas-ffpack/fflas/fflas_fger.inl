/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_fger.inl
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

#ifndef __FFLASFFPACK_fger_INL
#define __FFLASFFPACK_fger_INL
namespace FFLAS {

	template<class Field>
	inline void
	fger (const Field& F, const size_t M, const size_t N,
	      const typename Field::Element alpha,
	      const typename Field::Element * x, const size_t incx,
	      const typename Field::Element * y, const size_t incy,
	      typename Field::Element * A, const size_t lda)
	{

		typename Field::Element tmp;
		const typename Field::Element* xi=x, *yj=y;
		typename Field::Element* Ai=A;

		if ( M < N ){
			if ( F.areEqual( alpha, F.one ) )
				for ( ; Ai < A+M*lda; Ai+=lda, xi+=incx ){
					yj = y;
					for (size_t j = 0; j < N; ++j, yj+=incy )
						F.axpyin( *(Ai+j), *xi, *yj );
				}
			else if ( F.areEqual( alpha, F.mOne ) )
				for ( ; Ai < A+M*lda; Ai+=lda, xi+=incx ){
					F.neg( tmp, *xi );
					yj = y;
					for (size_t j = 0; j < N; ++j, yj+=incy )
						F.axpyin( *(Ai+j), tmp, *yj );
				}
			else
				for ( ; Ai < A+M*lda; Ai+=lda, xi+=incx ){
					F.mul( tmp, alpha, *xi );
					yj = y;
					for (size_t j = 0; j < N; ++j, yj+=incy )
						F.axpyin( *(Ai+j), tmp, *yj );
				}
		} else {
			if ( F.areEqual( alpha, F.one ) ){
				for ( ; Ai < A+N; ++Ai, yj+=incy ){
					xi = x;
					for (size_t i = 0; i < M; ++i, xi+=incx )
						F.axpyin( *(Ai+i*lda), *xi, *yj );
				}
			}
			else if ( F.areEqual( alpha, F.mOne ) )
				for ( ; Ai < A+N; ++Ai, yj+=incy ){
					F.neg( tmp, *yj );
					xi = x;
					for (size_t i = 0; i < M; ++i, xi+=incx )
						F.axpyin( *(Ai+i*lda), *xi, tmp );
				}
			else
				for ( ; Ai < A+N; ++Ai, yj+=incy ){
					F.mul( tmp, alpha, *yj );
					xi = x;
					for (size_t i = 0; i < M; ++i, xi+=incx )
						F.axpyin( *(Ai+i*lda), *xi, tmp );
				}
		}

	}

	template<>
	inline void
	fger( const DoubleDomain& , const size_t M, const size_t N,
	      const DoubleDomain::Element alpha,
	      const DoubleDomain::Element * x, const size_t incx,
	      const DoubleDomain::Element * y, const size_t incy,
	      DoubleDomain::Element * A, const size_t lda)
	{

		FFLASFFPACK_check(lda);
		cblas_dger( CblasRowMajor, (int)M, (int)N, alpha, x, (int)incx, y, (int)incy, A, (int)lda );
	}

	template<>
	inline void
	fger( const FloatDomain& , const size_t M, const size_t N,
	      const FloatDomain::Element alpha,
	      const FloatDomain::Element * x, const size_t incx,
	      const FloatDomain::Element * y, const size_t incy,
	      FloatDomain::Element * A, const size_t lda)
	{

			FFLASFFPACK_check(lda);
		cblas_sger( CblasRowMajor, (int)M, (int)N, alpha, x, (int)incx, y, (int)incy, A, (int)lda );
	}

} // FFLAS
#endif // __FFLASFFPACK_fger_INL
