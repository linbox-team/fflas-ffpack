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

#ifndef __FFLASFFPACK_ffpack_fsytrf_INL
#define __FFLASFFPACK_ffpack_fsytrf_INL

namespace FFPACK {

	template <class Field>
	inline bool fsytrf_BC (const Field& F, const FFLAS::FFLAS_UPLO UpLo, const size_t N,
						   typename Field::Element_ptr A, const size_t lda){

		typename Field::Element_ptr Ai = A, An = A;
		if (UpLo==FFLAS::FflasUpper){
			for (size_t i = 0; i<N; i++, Ai+=lda+1, An++){
				typename Field::Element_ptr tmp = FFLAS::fflas_new(F, 1, i);
				for (size_t j=0; j<i; ++j)
					F.div (tmp[j], An[j*lda], A[j*(lda+1)]);
				FFLAS::fgemv (F, FFLAS::FflasTrans, i, N-i, F.mOne, An, lda, tmp, 1, F.one, Ai, 1);
				FFLAS::fflas_delete(tmp);
				if (F.isZero(*Ai))
					return false;
				F.invin(*Ai);
				FFLAS::fscalin (F, N-i-1, *Ai, Ai+1, 1);
			}
		} else {
			for (size_t i = 0; i<N; i++, Ai+=lda+1, An+=lda){
				typename Field::Element_ptr tmp = FFLAS::fflas_new(F, 1, i);
				for (size_t j=0; j<i; ++j)
					F.div (tmp[j], An[j], A[j*(lda+1)]);
				FFLAS::fgemv (F, FFLAS::FflasNoTrans, N-i, i, F.mOne, An, lda, tmp, 1, F.one, Ai, lda);
				FFLAS::fflas_delete(tmp);
				if (F.isZero(*Ai))
					return false;
				F.invin(*Ai);
				FFLAS::fscalin (F, N-i-1, *Ai, Ai+lda, lda);
			}
		}
		return true;
	}

	template <class Field>
	inline bool fsytrf (const Field& F, const FFLAS::FFLAS_UPLO UpLo, const size_t N,
						typename Field::Element_ptr A, const size_t lda, size_t threshold){

		// if (N==1){
		// 	if (F.isZero(*A))
		// 		return false;
		// 	else{
		// 		F.invin(*A);
		// 		return true;
		// 	}
		if (N <= threshold){
				//std::cerr<<"Seuil"<<std::endl;
			return fsytrf_BC (F, UpLo, N, A, lda);
		}
		else {
				//std::cerr<<"Rec N = "<<N<<std::endl;

			size_t N1 = N>>1;
			size_t N2 = N-N1;
			size_t Arows, Acols;
			FFLAS::FFLAS_TRANSPOSE trans;
			FFLAS::FFLAS_SIDE side;
			if (UpLo==FFLAS::FflasUpper){side = FFLAS::FflasLeft; Arows = N1; Acols = N2;trans=FFLAS::FflasTrans;}
			else{side = FFLAS::FflasRight; Arows = N2; Acols = N1;trans=FFLAS::FflasNoTrans;}
				// Comments written for the UpLo = FflasUpper case
			typename Field::Element_ptr A12 = A + N1*((UpLo==FFLAS::FflasUpper)?1:lda);
			typename Field::Element_ptr A22 = A + N1*(lda+1);

				// A1 = U1^T x D1^-1 x U1
			if (!fsytrf (F, UpLo, N1, A, lda, threshold)) return false;

				// A12 <- U1^-T x A12
			FFLAS::ftrsm (F, side, UpLo, FFLAS::FflasTrans, FFLAS::FflasUnit, Arows, Acols, F.one, A, lda, A12, lda);

				// A22 <- A22 - A12^T x A12
			FFLAS::fsyrk (F, UpLo, trans, N2, N1, F.mOne, A12, lda, A, lda+1, F.one, A22, lda);

				// A22 = U2^T x D2^-1 x U2
			if (!fsytrf (F, UpLo, N2, A22, lda, threshold)) return false;
			return true;
		}
	}
	
} //FFPACK

#endif // __FFLASFFPACK_ffpack_fsytrf_INL
