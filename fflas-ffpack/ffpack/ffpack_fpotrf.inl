/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2016 FFLAS-FFACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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

#ifndef __FFLASFFPACK_ffpack_fpotrf_INL
#define __FFLASFFPACK_ffpack_fpotrf_INL

namespace FFPACK {

	template <class Field>
	inline bool fpotrf (const Field& F, const FFLAS::FFLAS_UPLO UpLo, const size_t N,
						typename Field::Element_ptr A, const size_t lda){

		if (N==1){
		} else {
			size_t N1 = N>>1;
			size_t N2 = N-N1;

				// Comments written for the UpLo = FflasUpper case
			typename Field::Element_ptr A21 = A + N1*lda;
			typename Field::Element_ptr A12 = A + N1;
			typename Field::Element_ptr A22 = A21 + N1;

				// A1 = U1^T x D1 x U1
			if (!fpotrf (F, UpLo, N1, A, lda)) return false;

				// A12 <- U1^-T x A12
			FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasTrans, FFLAS::FflasUnit,
				   N1, N2, F.one, A, lda, A+N1, lda);

				// A12 <- D1^-1 x A12
			typename Field::Element_ptr Ai = A, A12i = A12;
			for (size_t i=0; i<N1; i++, Ai+=(lda+1), A12i+=lda){
				typename Field::Element inv;
				F.init (inv);
				F.inv (inv, Ai);
				fscal (F, N2, inv, A12i, 1);
			}
				// A22 <- A22 - A12^T x A12
			FFLAS::fsyrk (F, UpLo, FFLAS::FflasNoTrans, N2, N1, F.mOne, A12, lda, F.one, A22, lda);

				// A22 = U2^T x D2 x U2
			if (!fpotrf (F, Uplo, N2, A22, lda)) return false;

			return true;
		}
	}
	
} //FFPACK

#endif // __FFLASFFPACK_ffpack_fpotrf_INL
