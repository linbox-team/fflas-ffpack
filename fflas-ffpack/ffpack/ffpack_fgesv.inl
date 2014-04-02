/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 FFLAS-FFACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 * BB <bbboyer@ncsu.edu>
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

#ifndef __FFLASFFPACK_ffpack_fgesv_INL
#define __FFLASFFPACK_ffpack_fgesv_INL


namespace FFPACK {



template <class Field>
	size_t
	fgesv (const Field& F,
	       const FFLAS::FFLAS_SIDE Side,
	       const size_t M, const size_t N,
	       typename Field::Element *A, const size_t lda,
	       typename Field::Element *B, const size_t ldb,
	       int * info)
	{

		size_t Na;
		if (Side == FFLAS::FflasLeft)
			Na = M;
		else
			Na = N;

		size_t* P = new size_t[Na];
		size_t* Q = new size_t[Na];

		size_t R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, Na, Na, A, lda, P, Q, FfpackLQUP);

		fgetrs (F, Side, M, N, R, A, lda, P, Q, B, ldb, info);

		delete[] P;
		delete[] Q;

		return R;
	}

	template <class Field>
	size_t
	fgesv (const Field& F,
	       const FFLAS::FFLAS_SIDE Side,
	       const size_t M, const size_t N, const size_t NRHS,
	       typename Field::Element *A, const size_t lda,
	       typename Field::Element *X, const size_t ldx,
	       const typename Field::Element *B, const size_t ldb,
	       int * info)
	{

		size_t* P = new size_t[N];
		size_t* Q = new size_t[M];

		size_t R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, N, A, lda, P, Q, FfpackLQUP);

		fgetrs (F, Side, M, N, NRHS, R, A, lda, P, Q, X, ldx, B, ldb, info);

		delete[] P;
		delete[] Q;

		return R;
	}

} //FFPACK

#endif // __FFLASFFPACK_ffpack_fgesv_INL
