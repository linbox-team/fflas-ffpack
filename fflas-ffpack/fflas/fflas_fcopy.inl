/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_fcopy.inl
 * Copyright (C) 2007 Clement Pernet
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

#ifndef __FFLASFFPACK_fcopy_INL
#define __FFLASFFPACK_fcopy_INL

#include <string.h>

namespace FFLAS {

	template<class Field>
	inline void
	fcopy (const Field& F, const size_t N,
	       typename Field::Element * X, const size_t incX,
	       const typename Field::Element * Y, const size_t incY )
	{
		typename Field::Element * Xi = X;
		const typename Field::Element * Yi=Y;

		if (incX == 1 && incY == 1) {
			// memcpy(X,Y,N*sizeof(typename Field::Element)); // much faster (hopefully)
			for (; Xi < X+N; ++Xi, ++Yi)
				F.assign(*Xi,*Yi);

			return;
		}
		for (; Xi < X+N*incX; Xi+=incX, Yi+=incY )
			F.assign(*Xi,*Yi);
		return;
	}

	template<class Field>
	void fcopy (const Field& F, const size_t m, const size_t n,
		    typename Field::Element * A, const size_t lda,
		    const typename Field::Element * B, const size_t ldb )
	{
		FFLASFFPACK_check(n<=lda);
		FFLASFFPACK_check(n<=ldb);
		// if possible, copy one big block
		if (lda == n && ldb == n) {
			fcopy(F,m*n,A,1,B,1);
			return ;
		}
		// else, copy row after row
		for (size_t i = 0 ; i < m ; ++i) {
			fcopy(F,n,A+i*lda,1,B+i*ldb,1);
		}
		return;

	}


}


#endif // __FFLASFFPACK_fcopy_INL
