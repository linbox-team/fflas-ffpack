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
#include "fflas-ffpack/utils/fflas_io.h"
namespace FFPACK {

	template <class Field>
	inline bool fsytrf_BC_Crout (const Field& F, const FFLAS::FFLAS_UPLO UpLo, const size_t N,
						   typename Field::Element_ptr A, const size_t lda,
						   typename Field::Element_ptr Dinv, const size_t incDinv){

		typename Field::Element_ptr Ai = A, An = A;
		if (UpLo==FFLAS::FflasUpper){
			for (size_t i = 0; i<N; i++, Ai+=lda+1, An++){
				typename Field::Element_ptr tmp = FFLAS::fflas_new(F, i);
				typename Field::Element_ptr Dinvj = Dinv;
				typename Field::Element_ptr Anj = An;
				for (size_t j=0; j<i; ++j, Anj+=lda,Dinvj+=incDinv)
					F.mul (tmp[j], *Anj, *Dinvj);
				FFLAS::fgemv (F, FFLAS::FflasTrans, i, N-i, F.mOne, An, lda, tmp, 1, F.one, Ai, 1);
				FFLAS::fflas_delete(tmp);
				if (F.isZero(*Ai))
					return false;
				F.inv (Dinv[i*incDinv], *Ai);
			}
		} else {
			for (size_t i = 0; i<N; i++, Ai+=lda+1, An+=lda){
				typename Field::Element_ptr tmp = FFLAS::fflas_new(F, i);
				typename Field::Element_ptr Dinvj = Dinv;
				typename Field::Element_ptr Anj = An;
				for (size_t j=0; j<i; ++j,Anj++,Dinvj+=incDinv)
					F.mul (tmp[j], *Anj, *Dinvj);
				FFLAS::fgemv (F, FFLAS::FflasNoTrans, N-i, i, F.mOne, An, lda, tmp, 1, F.one, Ai, lda);
				FFLAS::fflas_delete(tmp);
				if (F.isZero(*Ai))
					return false;
				F.inv (Dinv[i*incDinv], *Ai);
			}
		}
		return true;
	}

	template <class Field>
	inline size_t fsytrf_BC_RL (const Field& F, const FFLAS::FFLAS_UPLO UpLo, const size_t N,
								typename Field::Element_ptr A, const size_t lda,
								typename Field::Element_ptr Dinv, const size_t incDinv
								){
		typename Field::Element_ptr Ai = A, Dinvi = Dinv;
		if (UpLo==FFLAS::FflasUpper){
			for (size_t i=0; i<N; i++, Ai+=(lda+1), Dinvi+=incDinv){
				if (F.isZero(*Ai)) return false;
				F.inv (*Dinvi, *Ai);
				typename Field::Element mDinvi;
				F.init(mDinvi);
				F.neg(mDinvi,*Dinvi);
				typename Field::Element_ptr tmp = FFLAS::fflas_new(F, N-i-1);
					//FFLAS::fassign (F, N-i-1, Ai+1, 1, tmp, 1);
				FFLAS::fscal (F, N-i-1, mDinvi, Ai+1, 1, tmp, 1);
				typename Field::Element_ptr tmpi=tmp, Aj=Ai+lda+1;
				for (size_t j=i+1; j<N; j++,tmpi++, Aj+=lda+1)
					FFLAS::faxpy (F, N-j, *tmpi, Ai+j-i, 1, Aj, 1);
					//FFLAS::WriteMatrix(std::cerr<<"Fin iteration "<<i<< " tmp = "<<std::endl,F,1,N-i-1,tmp,N-i-1);
					FFLAS::fflas_delete(tmp);
						//FFLAS::WriteMatrix(std::cerr<<"Fin iteration "<<i<<"A = "<<std::endl,F,N,N,A,lda);
			}
		}else { // FflasLower
			for (size_t i=0; i<N; i++, Ai+=(lda+1), Dinvi+=incDinv){
				if (F.isZero(*Ai)) return false;
				F.inv (*Dinvi, *Ai);
				typename Field::Element mDinvi;
				F.init(mDinvi);
				F.neg(mDinvi,*Dinvi);
				typename Field::Element_ptr tmp = FFLAS::fflas_new(F, N-i-1);
					//FFLAS::fassign (F, N-i-1, Ai+lda, lda, tmp, 1);
				FFLAS::fscal (F, N-i-1, mDinvi, Ai+lda, lda, tmp, 1);
				typename Field::Element_ptr tmpi=tmp, Aj=Ai+lda+1;
				for (size_t j=i+1; j<N; j++,tmpi++, Aj+=lda+1)
					FFLAS::faxpy (F, N-j, *tmpi, Ai+(j-i)*lda, lda, Aj, lda);
				FFLAS::fflas_delete(tmp);
			}
		}
		return true;				
	}
		template <class Field>
		inline size_t fsytrf_RPM_BC_Crout (const Field& F, const FFLAS::FFLAS_UPLO UpLo, const size_t N,
										typename Field::Element_ptr A, const size_t lda,
										typename Field::Element_ptr Dinv, const size_t incDinv,
										size_t * P){
			
		typename Field::Element_ptr Ai = A, An = A;
		if (UpLo==FFLAS::FflasUpper){
			for (size_t i = 0; i<N; i++, Ai+=lda+1, An++){
				typename Field::Element_ptr tmp = FFLAS::fflas_new(F, i);
				typename Field::Element_ptr Dinvj = Dinv;
				typename Field::Element_ptr Anj = An;
				for (size_t j=0; j<i; ++j, Anj+=lda,Dinvj+=incDinv)
					F.mul (tmp[j], *Anj, *Dinvj);
				FFLAS::fgemv (F, FFLAS::FflasTrans, i, N-i, F.mOne, An, lda, tmp, 1, F.one, Ai, 1);
				FFLAS::fflas_delete(tmp);
				typename Field::Element Apiv = Ai;
				while (F.isZero(*Apiv) && Apiv!=Ai+N-i) Apiv++;
				size_t j=Apiv-Ai;
				if (j){ // No pivot found at Aii
					if (i+j<N) {// found a pivot somewhere else
						
					} else {// zero row

					}
				}
				
				F.inv (Dinv[i*incDinv], *Ai);
			}
		} else {
			for (size_t i = 0; i<N; i++, Ai+=lda+1, An+=lda){
				typename Field::Element_ptr tmp = FFLAS::fflas_new(F, i);
				typename Field::Element_ptr Dinvj = Dinv;
				typename Field::Element_ptr Anj = An;
				for (size_t j=0; j<i; ++j,Anj++,Dinvj+=incDinv)
					F.mul (tmp[j], *Anj, *Dinvj);
				FFLAS::fgemv (F, FFLAS::FflasNoTrans, N-i, i, F.mOne, An, lda, tmp, 1, F.one, Ai, lda);
				FFLAS::fflas_delete(tmp);
				if (F.isZero(*Ai))
					return false;
				F.inv (Dinv[i*incDinv], *Ai);
			}
		}
		return true;
	}

	template <class Field>
	inline bool fsytrf_nonunit (const Field& F, const FFLAS::FFLAS_UPLO UpLo, const size_t N,
								typename Field::Element_ptr A, const size_t lda,
								typename Field::Element_ptr Dinv, const size_t incDinv,
								size_t threshold){

		// if (N==1){
		// 	if (F.isZero(*A))
		// 		return false;
		// 	else{
		// 		return true;
		// 	}
		if (N <= threshold)
#ifdef FSYTRF_BC_RL
			return fsytrf_BC_RL (F, UpLo, N, A, lda, Dinv, incDinv);
#elif defined FSYTRF_BC_CROUT
			return fsytrf_BC_CROUT (F, UpLo, N, A, lda, Dinv, incDinv);
#else
			return fsytrf_BC_RL (F, UpLo, N, A, lda, Dinv, incDinv);
#endif
		else {
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
			if (!fsytrf_nonunit (F, UpLo, N1, A, lda, Dinv, incDinv, threshold)) return false;

				// A12 <- U1^-T x A12
			FFLAS::ftrsm (F, side, UpLo, FFLAS::FflasTrans, FFLAS::FflasNonUnit, Arows, Acols, F.one, A, lda, A12, lda);

				// A22 <- A22 - A12^T x D1 x A12 and A12 <- A12
			FFLAS::fsyrk (F, UpLo, trans, N2, N1, F.mOne, A12, lda, A, lda+1, F.one, A22, lda);

				// A22 = U2^T x D2^-1 x U2
			if (!fsytrf_nonunit (F, UpLo, N2, A22, lda, Dinv+N1, incDinv, threshold)) return false;
			return true;
		}
	}
	
	template <class Field>
	inline bool fsytrf (const Field& F, const FFLAS::FFLAS_UPLO UpLo, const size_t N,
						typename Field::Element_ptr A, const size_t lda,
						size_t threshold){
		typename Field::Element_ptr Dinv = FFLAS::fflas_new(F,N);
		bool success = fsytrf_nonunit (F, UpLo, N, A, lda, Dinv, 1, threshold);
		// FFLAS::WriteMatrix(std::cerr<<"After fsytrf_nonunit A = "<<std::endl,F,N,N,A, lda);
		// FFLAS::WriteMatrix(std::cerr<<"After fsytrf_nonunit Dinv = "<<std::endl,F,1,N,Dinv, 1);
		if (!success) return false;
		size_t incA = (UpLo==FFLAS::FflasUpper) ? 1 : lda;
		for (size_t i=0; i<N; i++)
			FFLAS::fscalin (F, N-i-1, Dinv[i], A+i*(lda+1)+incA, incA);
		FFLAS::fflas_delete(Dinv);
		return true;
	}
} //FFPACK

#endif // __FFLASFFPACK_ffpack_fsytrf_INL
