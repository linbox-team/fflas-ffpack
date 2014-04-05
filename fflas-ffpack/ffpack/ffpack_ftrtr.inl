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

#ifndef __FFLASFFPACK_ffpack_ftrtr_INL
#define __FFLASFFPACK_ffpack_ftrtr_INL


namespace FFPACK {


template<class Field>
	void
	ftrtri (const Field& F, const FFLAS::FFLAS_UPLO Uplo, const FFLAS::FFLAS_DIAG Diag,
		const size_t N, typename Field::Element * A, const size_t lda)
	{
		if (N == 1){
			if (Diag == FFLAS::FflasNonUnit)
				F.invin (*A);
		}
		else {
			size_t N1 = N/2;
			size_t N2 = N - N1;
			ftrtri (F, Uplo, Diag, N1, A, lda);
			ftrtri (F, Uplo, Diag, N2, A + N1*(lda+1), lda);
			if (Uplo == FFLAS::FflasUpper){
				ftrmm (F, FFLAS::FflasLeft, Uplo, FFLAS::FflasNoTrans, Diag, N1, N2,
				       F.one, A, lda, A + N1, lda);
				ftrmm (F, FFLAS::FflasRight, Uplo, FFLAS::FflasNoTrans, Diag, N1, N2,
				       F.mOne, A + N1*(lda+1), lda, A + N1, lda);
			}
			else {
				ftrmm (F, FFLAS::FflasLeft, Uplo, FFLAS::FflasNoTrans, Diag, N2, N1,
				       F.one, A + N1*(lda+1), lda, A + N1*lda, lda);
				ftrmm (F, FFLAS::FflasRight, Uplo, FFLAS::FflasNoTrans, Diag, N2, N1,
				       F.mOne, A, lda, A + N1*lda, lda);
			}
		}
	}


	template<class Field>
	void
	ftrtrm (const Field& F, const FFLAS::FFLAS_DIAG diag, const size_t N,
		typename Field::Element * A, const size_t lda)
	{

		if (N == 1)
			return;
		size_t N1 = N/2;
		size_t N2 = N-N1;

		ftrtrm (F, diag, N1, A, lda);

		fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, N1, N1, N2, F.one,
		       A+N1, lda, A+N1*lda, lda, F.one, A, lda);

		ftrmm (F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans,
		       (diag == FFLAS::FflasUnit) ? FFLAS::FflasNonUnit : FFLAS::FflasUnit,
		       N1, N2, F.one, A + N1*(lda+1), lda, A + N1, lda);

		ftrmm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, diag, N2, N1,
		       F.one, A + N1*(lda+1), lda, A + N1*lda, lda);

		ftrtrm (F, diag, N2, A + N1*(lda+1), lda);

	}

	template<class Field>
	void trinv_left( const Field& F, const size_t N, const typename Field::Element * L, const size_t ldl,
			 typename Field::Element * X, const size_t ldx )
	{
		FFLAS::fcopy(F,N,N,L,ldl,X,ldx);
		ftrtri (F, FFLAS::FflasLower, FFLAS::FflasUnit, N, X, ldx);
		//invL(F,N,L,ldl,X,ldx);
	}

} // FFPACK

#endif // __FFLASFFPACK_ffpack_ftrtr_INL
