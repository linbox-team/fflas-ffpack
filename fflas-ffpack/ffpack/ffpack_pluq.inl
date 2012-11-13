/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* ffpack/ffpack_pluq.inl
 * Copyright (C) 2012 Clement Pernet
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

#ifndef __FFLASFFPACK_ffpack_pluq_INL
#define __FFLASFFPACK_ffpack_pluq_INL

#ifndef MIN
#define MIN(a,b) (a<b)?a:b
#endif
#ifndef MAX
#define MAX(a,b) (a<b)?b:a
#endif

using namespace FFLAS;

namespace FFPACK {
	template<class Field>
	inline size_t
	PLUQ (const Field& Fi, const FFLAS_DIAG Diag,
	      const size_t M, const size_t N,
	      typename Field::Element * A, const size_t lda, size_t*P, size_t *Q){
		
		if (M = 1){
			size_t piv = 0;
			while ((piv < N) && F.isZero (A[piv])) piv++;
			if (piv == N)
				return 0;
			F.assign (*A, A[piv]);
			F.assign (A[piv], F.Zero);
			if (Diag== FflasUnit){
				typename Field::Element invpivot;
				F.inv(invpivot, *A);
				for (size_t i=piv+1; i<N; ++i)
					F.mulin (A[i], invpivot);
			}
			P[0] = 0;
			Q[0] = piv;
			return 1;
		}
		if (N = 1){
			size_t piv = 0;
			while ((piv < M) && F.isZero (A +piv*lda)) piv++;
			if (piv == M)
				return 0;
			F.assign (*A, *(A+piv*lda));
			F.assign (*(A+piv*lda), F.Zero);
			if (Diag== FflasNonUnit){
				typename Field::Element invpivot;
				F.inv(invpivot, *A);
				for (size_t i=piv+1; i<M; ++i)
					F.mulin (*(A+i*lda), invpivot);
			}
			P[0] = piv;
			Q[0] = 0;
			return 1;
		}
		FFLAS_DIAG OppDiag = (Diag == FflasUnit)? FflasNonUnit : FflasUnit;
		size_t M2 = M >> 1;
		size_t N2 = N >> 1;
		size_t * P1 = new size_t [M2];
		size_t * P2 = new size_t [N2];
		size_t R1,R2,R3,R4;

		    // A1 = P1 [ L1 ] [ U1 V1 ] Q1
		    //        [ M1 ]
		R1 = PLUQ (Fi, Diag, M2, N2, A, lda, P1, Q1);
		typename Field::Element * A2 = A + N2;
		typename Field::Element * A3 = A + M2*lda;
		typename Field::Element * A3 = A3 + N2;
		typename Field::Element * F = A2 + R1*lda;
		typename Field::Element * G = A3 + R1;
		    // [ B1 ] <- P1^T A2
		    // [ B2 ]
		applyP (Fi, FflasLeft, FflasTrans, N-N2, 0, R1, A2, lda, P1);
		    // [ C1 C2 ] <- A3 Q1^T
		applyP (Fi, FflasRight, FflasTrans, M-M2, 0, R1, A3, lda, Q1);
		    // D <- L1^-1 B1
		ftrsm (Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R1, N-N2, F.one, A, lda, A2, lda);
		    // E <- C1 U1^-1 
		ftrsm (Fi, FflasRight, FflasUpper, FflasNoTrans, Diag, M-M2, R1, F.one, A, lda, A3, lda);
		    // F <- B2 - M1 D
		fgemm (Fi, FflasNoTrans, FflasNoTrans, M2-R1, N-N2, R1, F.mOne, A + R1*lda, lda, A2, lda, F.one, A2+R1*lda, lda);
		    // G <- C2 - E V1
		fgemm (Fi, FflasNoTrans, FflasNoTrans, M-M2, N2-R1, R1, F.mOne, A3, lda, A+R1, lda, F.one, A3+R1, lda);
		    // F = P2 [ L2 ] [ U2 V2 ] Q2
		    //        [ M2 ]
		size_t P2 = new size_t [M2-R1];
		size_t Q2 = new size_t [N-N2];
		R2 = PLUQ (Fi, Diag, M2-R1, N-N2, F, lda, P2, Q2);
		    // G = P3 [ L3 ] [ U3 V3 ] Q3
		    //        [ M3 ]
		size_t P3 = new size_t [M-M2];
		size_t Q3 = new size_t [N2-R1];
		R3 = PLUQ (Fi, Diag, M-M2, N2-R1, G, lda, P3, Q3);
		    // [ H1 H2 ] <- P3^T H Q2^T
		    // [ H3 H4 ]
		applyP (Fi, FflasRight, FflasTrans, M-M2, 0, R2, A4, lda, Q2);
		applyP (Fi, FflasLeft, FflasTrans, N-N2, 0, R3, A4, lda, P3);
		    // [ E1 ] <- P3^T E
		    // [ E2 ]
		applyP (Fi, FflasLeft, FflasTrans, N-N2, 0, R3, A3, lda, P3);
		    // [ M11 ] <- P2^T M1
		    // [ M12 ]
		applyP (Fi, FflasLeft, FflasTrans, R1, 0, R2, A+R1*lda, lda, P2);
		    // [ D1 D2 ] <- D Q2^T
		applyP (Fi, FflasRight, FflasTrans, M-M2, 0, R2, A2, lda, Q2);
		    // [ V1 V2 ] <- V1 Q3^T
		applyP (Fi, FflasRight, FflasTrans, R1, 0, R3, A+R1, lda, Q3);
		    // I <- H U2^-1
		ftrsm (Fi, FflasRight, FflasUpper, FflasNoTrans, Diag, M-M2, R2, F.one, F, lda, A4, lda);
		    // J <- L3^-1 I (in a temp)
		typename Field::Element temp = new typename Field::Element [R3*R2];
		for (size_t i=0; i<R3; ++i)
			fcopy (F, R2, temp + i*R2, 1, A4 + i*lda, 1);
		ftrsm (Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R3, R2, F.one, G, lda, temp, R2);
		
		    // K <- H3 U2^-1
		ftrsm (Fi, FflasRight, FflasUpper, FflasNoTrans, Diag, M-M2-R3, R2, F.one, F, lda, A4+R3*lda, lda);
		    // N <- L3^-1 H2
		ftrsm (Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R3, N-N2-R2, F.one, G, lda, A4+R2, R2);
		    // O <- N - J V2
		fgemm (Fi, FflasNoTrans, FflasNoTrans, R3, N-N2-R2, R2, F.mOne, temp, R2, F+R2, lda, F.one, A4+R2, lda);
		delete[] temp;
		    // R <- H4 - K V2 - M3 O
		fgemm (Fi, FflasNoTrans, FflasNoTrans, M-M2-R3, N-N2-R2, R2, F.mOne, A4+R3*lda, lda, F+R2, lda, F.one, A4+R2+R3*lda, lda);
		fgemm (Fi, FflasNoTrans, FflasNoTrans, M-M2-R3, N-N2-R2, R3, F.mOne, G+R3*lda, lda, A4+R2, lda, F.one, A4+R2+R3*lda, lda);
		    // H4 = P4 [ L4 ] [ U4 V4 ] Q4
		    //         [ M4 ]
		size_t P4 = new size_t [M-M2-R3];
		size_t Q4 = new size_t [N-N2-R2];
		R4 = PLUQ (Fi, Diag, M-M2-R3, N-N2-R2, R, lda, P4, Q4);
		    // [ E21 M31 0 K1 ] <- P4^T [ E2 M3 0 K ] 
		    // [ E22 M32 0 K2 ]
		applyP (Fi, FflasLeft, FflasTrans, N2+R2, 0, R4, A3+R3*lda, lda, P4);
		    // [ D21 D22 ]     [ D2 ]
		    // [ V21 V22 ]  <- [ V2 ] Q4^T
		    // [  0   0  ]     [  0 ]
		    // [ O1   O2 ]     [  O ]
		applyP (Fi, FflasRight, FflasTrans, M2+R3, 0, R4, A2+R2, lda, Q4);

	}
} // namespace FFPACK
