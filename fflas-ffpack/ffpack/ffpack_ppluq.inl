/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* ffpack/ffpack_ppluq.inl
 * Copyright (C) 2014 Ziad Sultan
 *
 * Written by Ziad.Sultan@imag.fr
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
 * License along with this library; if not, WRITE to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __FFLASFFPACK_ffpack_ppluq_INL
#define __FFLASFFPACK_ffpack_ppluq_INL


// #define __FFLASFFPACK_USE_OPENMP // you can't define this
#ifdef __FFLASFFPACK_USE_OPENMP
#define __FFLAS__TRSM_READONLY
//#include "parallel.h"

// do not use namespace in the library, especially std.
// using namespace std;
// using namespace FFLAS;
// using namespace FFPACK;



#define MEMCOPY
#define LEFTLOOKING
#define CROUT
#define BASECASE_K 256

namespace FFPACK {

	template<class Field>
	inline size_t
	pPLUQ(const Field& Fi, const FFLAS::FFLAS_DIAG Diag,
	      const size_t M, const size_t N,
	      typename Field::Element_ptr A, const size_t lda,
	      size_t* P, size_t* Q)
	  {

		  const FFLAS::CuttingStrategy method = FFLAS::BLOCK_THREADS;
		  size_t NUM = NUM_THREADS;
		  for (size_t i=0; i<M; ++i) P[i] = i;
		  for (size_t i=0; i<N; ++i) Q[i] = i;
		  if (std::min(M,N) == 0) return 0;
		  if (std::max(M,N) == 1) return (Fi.isZero(*A))? 0 : 1;
		  if (M == 1){
			  size_t piv = 0;
			  while ((piv < N) && Fi.isZero (A[piv])) piv++;
			  if (piv == N)
				  return 0;
			  if (piv){
				  Q[0] = piv;
				  Fi.assign (*A, A[piv]);
				  Fi.assign (A[piv], Fi.zero);
			  }
			  if (Diag== FFLAS::FflasUnit){
				  typename Field::Element invpivot;
				  Fi.inv(invpivot, *A);
				  for (size_t i=piv+1; i<N; ++i)
					  Fi.mulin (A[i], invpivot);
			  }
			  return 1;
		  }
		  if (N == 1){
			  size_t piv = 0;
			  while ((piv < M) && Fi.isZero (A[piv*lda])) piv++;
			  if (piv == M)
				  return 0;
			  if (piv){
				  P[0] = piv;
				  Fi.assign (*A, *(A+piv*lda));
				  Fi.assign (*(A+piv*lda), Fi.zero);
			  }
			  if (Diag== FFLAS::FflasNonUnit){
				  typename Field::Element invpivot;
				  Fi.inv(invpivot, *A);
				  for (size_t i=piv+1; i<M; ++i)
					  Fi.mulin (*(A+i*lda), invpivot);
			  }
			  return 1;
		  }

#ifdef BASECASE_K
		  if (std::min(M,N) < BASECASE_K)
			  return PLUQ_basecaseCrout(Fi, Diag, M, N, A, lda, P, Q);
#endif

		  FFLAS::FFLAS_DIAG OppDiag = (Diag == FFLAS::FflasUnit)? FFLAS::FflasNonUnit : FFLAS::FflasUnit;
		  size_t M2 = M >> 1;
		  size_t N2 = N >> 1;
		  //    size_t M2 = size_t(M/3);
		  //    size_t N2 = size_t(N/3);
		  size_t * P1 = new size_t [M2];
		  size_t * Q1 = new size_t [N2];
		  size_t R1,R2,R3,R4;

		  // A1 = P1 [ L1 ] [ U1 V1 ] Q1
		  //        [ M1 ]
		  //#pragma omp task shared(R1) firstprivate(M2,N2,lda) shared(A,P1,Q1)
		  R1 = pPLUQ (Fi, Diag, M2, N2, A, lda, P1, Q1);
		  typename Field::Element_ptr A2 = A + N2;
		  typename Field::Element_ptr A3 = A + M2*lda;
		  typename Field::Element_ptr A4 = A3 + N2;
		  typename Field::Element_ptr F = A2 + R1*lda;
		  typename Field::Element_ptr G = A3 + R1;

		  // [ B1 ] <- P1^T A2
		  // [ B2 ]
		  TASK(READ(Fi, P1), NOWRITE(), READWRITE(A2), papplyP, Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N-N2, 0, M2, A2, lda, P1);
		  //applyP( Fi, FflasLeft, FflasNoTrans, N-N2, 0, M2, A2, lda, P1);

		  // [ C1 C2 ] <- A3 Q1^T
		  TASK(READ(Fi, Q1), NOWRITE(), READWRITE(A3), papplyP, Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M-M2, 0, N2, A3, lda, Q1);
		  //papplyP( Fi, FflasRight, FflasTrans, M-M2, 0, N2, A3, lda, Q1);

		  WAIT;
		  // D <- L1^-1 B1
		  TASK(READ(Fi, A, R1), NOWRITE(), READWRITE(A2), pftrsm, Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R1, N-N2, Fi.one, A, lda, A2, lda,  method, NUM);
		  //    pftrsm( Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R1, N-N2, Fi.one, A, lda, A2 , lda,  method, NUM);
		  //ftrsm( Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R1, N-N2, Fi.one, A, lda, A2 , lda);

		  // E <- C1 U1^-1
		  TASK(READ(Fi, R1, A), NOWRITE(), READWRITE(A3), pftrsm,Fi, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag, M-M2, R1, Fi.one, A, lda, A3, lda,  method, NUM);
		  //pftrsm(Fi, FflasRight, FflasUpper, FflasNoTrans, Diag, M-M2, R1, Fi.one, A, lda, A3, lda,  method, NUM);
		  //ftrsm(Fi, FflasRight, FflasUpper, FflasNoTrans, Diag, M-M2, R1, Fi.one, A, lda, A3, lda);

		  WAIT;
		  // F <- B2 - M1 D
		  TASK(READ(Fi, A), NOWRITE(), READWRITE(A2), pfgemm, Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M2-R1, N-N2, R1, Fi.mOne, A + R1*lda, lda, A2, lda, Fi.one, A2+R1*lda, lda, method, NUM_THREADS);
		  //fgemm( Fi, FflasNoTrans, FflasNoTrans, M2-R1, N-N2, R1, Fi.mOne, A + R1*lda, lda, A2, lda, Fi.one, A2+R1*lda, lda);

		  // G <- C2 - E V1
		  TASK(READ(Fi, R1, A), NOWRITE(), READWRITE(A3), pfgemm, Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2, N2-R1, R1, Fi.mOne, A3, lda, A+R1, lda, Fi.one, A3+R1, lda,  method, NUM_THREADS);
		  //fgemm( Fi, FflasNoTrans, FflasNoTrans, M-M2, N2-R1, R1, Fi.mOne, A3, lda, A+R1, lda, Fi.one, A3+R1, lda);

		  WAIT;

		  size_t * P2 = new size_t [M2-R1];
		  size_t * Q2 = new size_t [N-N2];
		  // F = P2 [ L2 ] [ U2 V2 ] Q2
		  //        [ M2 ]
		  TASK(READ(Fi), WRITE(R2), READWRITE(F, P2, Q2), PPLUQ, R2, Fi, Diag, M2-R1, N-N2, F, lda, P2, Q2);
		  //R2 = PLUQ (Fi, Diag, M2-R1, N-N2, F, lda, P2, Q2);

		  size_t * P3 = new size_t [M-M2];
		  size_t * Q3 = new size_t [N2-R1];
		  // G = P3 [ L3 ] [ U3 V3 ] Q3
		  //        [ M3 ]
		  TASK(READ(Fi), WRITE(R3), READWRITE(G, P3, Q3), PPLUQ, R3, Fi, Diag, M-M2, N2-R1, G, lda, P3, Q3);
		  //		  R3 = PLUQ (Fi, Diag, M-M2, N2-R1, G, lda, P3, Q3);

		  //		  WAIT;


		  // H <- A4 - ED
		  TASK(READ(Fi, M2, N2, R1, A3, A2), NOWRITE(), READWRITE(A4), pfgemm, Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2, N-N2, R1, Fi.mOne, A3, lda, A2, lda, Fi.one, A4, lda, 0,  method, NUM_THREADS);
		  //fgemm( Fi, FflasNoTrans, FflasNoTrans, M-M2, N-N2, R1, Fi.mOne, A3, lda, A2, lda, Fi.one, A4, lda);
		  //		  std::cout<<"NUM "<<NUM_THREADS<<std::endl;

		  WAIT;

		  // [ H1 H2 ] <- P3^T H Q2^T
		  // [ H3 H4 ]
		  TASK(READ(Fi,Q2), NOWRITE(), READWRITE(A4), papplyP, Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M-M2, 0, N-N2, A4, lda, Q2);
		  //applyP( Fi, FflasRight, FflasTrans, M-M2, 0, N-N2, A4, lda, Q2);
		  WAIT;

		  TASK(READ(Fi,P3), NOWRITE(), READWRITE(A4), papplyP, Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N-N2, 0, M-M2, A4, lda, P3);
		  //applyP( Fi, FflasLeft, FflasNoTrans, N-N2, 0, M-M2, A4, lda, P3);
		  WAIT;

      // [ E1 ] <- P3^T E
      // [ E2 ]
		  TASK(READ(Fi,P3), NOWRITE(), READWRITE(A3), papplyP,  Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, R1, 0, M-M2, A3, lda, P3);
		  //applyP(  Fi, FflasLeft, FflasNoTrans, R1, 0, M-M2, A3, lda, P3);

      // [ M11 ] <- P2^T M1
      // [ M12 ]
		  TASK(READ(Fi,P2), NOWRITE(), READWRITE(A), papplyP,Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, R1, 0, M2-R1, A+R1*lda, lda, P2);
		  //applyP(Fi, FflasLeft, FflasNoTrans, R1, 0, M2-R1, A+R1*lda, lda, P2);

      // [ D1 D2 ] <- D Q2^T
		  TASK(READ(Fi, Q2), NOWRITE(), READWRITE(A2), papplyP, Fi, FFLAS::FflasRight, FFLAS::FflasTrans, R1, 0, N-N2, A2, lda, Q2);
		  //papplyP( Fi, FflasRight, FflasTrans, R1, 0, N-N2, A2, lda, Q2);

      // [ V1 V2 ] <- V1 Q3^T
		  TASK(READ(Fi, Q3), NOWRITE(), READWRITE(A), papplyP, Fi, FFLAS::FflasRight, FFLAS::FflasTrans, R1, 0, N2-R1, A+R1, lda, Q3);
		  //applyP( Fi, FflasRight, FflasTrans, R1, 0, N2-R1, A+R1, lda, Q3);

      // I <- H1 U2^-1
      // K <- H3 U2^-1
		  TASK(READ(Fi, R2, F), NOWRITE(), READWRITE(A4), pftrsm, Fi, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag, M-M2, R2, Fi.one, F, lda, A4, lda,  method, NUM);
		  //pftrsm( Fi, FflasRight, FflasUpper, FflasNoTrans, Diag, M-M2, R2, Fi.one, F, lda, A4, lda,  method, NUM);
		  //ftrsm( Fi, FflasRight, FflasUpper, FflasNoTrans, Diag, M-M2, R2, Fi.one, F, lda, A4, lda);
		  WAIT;

		  typename Field::Element_ptr temp = fflas_new (F, R3, R2);
		  /*    for (size_t i=0; i<R3; ++i)
			fcopy (Fi, R2, temp + i*R2, 1, A4 + i*lda, 1);
		  */
		  FFLAS::fcopy (Fi, R3, R2, A4 , lda, temp , R2);

    // J <- L3^-1 I (in a temp)
		  TASK(READ(Fi, R2, G), NOWRITE(), READWRITE(temp), pftrsm,  Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R3, R2, Fi.one, G, lda, temp, R2,  method, NUM);
		  //pftrsm( Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R3, R2, Fi.one, G, lda, temp, R2,  method, NUM);
		  //ftrsm( Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R3, R2, Fi.one, G, lda, temp, R2);

   // N <- L3^-1 H2
		  TASK(READ(Fi, R3, G), NOWRITE(), READWRITE(A4), pftrsm,Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R3, N-N2-R2, Fi.one, G, lda, A4+R2, lda,  method, NUM);
		  //    pftrsm(Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R3, N-N2-R2, Fi.one, G, lda, A4+R2, lda,  method, NUM);
		  //ftrsm(Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R3, N-N2-R2, Fi.one, G, lda, A4+R2, lda);
		  WAIT;

   // O <- N - J V2
		  TASK(READ(Fi, R2, temp, F), NOWRITE(), READWRITE(A4), pfgemm, Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, R3, N-N2-R2, R2, Fi.mOne, temp, R2, F+R2, lda, Fi.one, A4+R2, lda,  method, NUM_THREADS);
		  //fgemm( Fi, FflasNoTrans, FflasNoTrans, R3, N-N2-R2, R2, Fi.mOne, temp, R2, F+R2, lda, Fi.one, A4+R2, lda);

		  typename Field::Element_ptr R = A4 + R2 + R3*lda;
  // R <- H4 - K V2
		  TASK(READ(Fi, R2, A4, F), NOWRITE(), READWRITE(R), pfgemm, Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2-R3, N-N2-R2, R2, Fi.mOne, A4+R3*lda, lda, F+R2, lda, Fi.one, R, lda,  method, NUM_THREADS);
		  //fgemm( Fi, FflasNoTrans, FflasNoTrans, M-M2-R3, N-N2-R2, R2, Fi.mOne, A4+R3*lda, lda, F+R2, lda, Fi.one, R, lda);
		  WAIT;

		  fflas_delete (temp);
    // R <- R - M3 O
		  TASK(READ(Fi, R3, A4, G), NOWRITE(), READWRITE(R), pfgemm, Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2-R3, N-N2-R2, R3, Fi.mOne, G+R3*lda, lda, A4+R2, lda, Fi.one, R, lda,  method, NUM_THREADS);
		  //fgemm( Fi, FflasNoTrans, FflasNoTrans, M-M2-R3, N-N2-R2, R3, Fi.mOne, G+R3*lda, lda, A4+R2, lda, Fi.one, R, lda);

		  WAIT;

		  size_t * P4 = new size_t [M-M2-R3];
		  size_t * Q4 = new size_t [N-N2-R2];

    // H4 = P4 [ L4 ] [ U4 V4 ] Q4
    //         [ M4 ]
		  //TASK(READ(Fi), NOWRITE(R4), READWRITE(R, P4, Q4), PPLUQ, R4, Fi, Diag, M-M2-R3, N-N2-R2, R, lda, P4, Q4);
		  R4 = PLUQ (Fi, Diag, M-M2-R3, N-N2-R2, R, lda, P4, Q4);
		  //WAIT;

    // [ E21 M31 0 K1 ] <- P4^T [ E2 M3 0 K ]
    // [ E22 M32 0 K2 ]
		  TASK(READ(Fi, P4), NOWRITE(), READWRITE(A3), papplyP, Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N2+R2, 0, M-M2-R3, A3+R3*lda, lda, P4);
		  //applyP( Fi, FflasLeft, FflasNoTrans, N2+R2, 0, M-M2-R3, A3+R3*lda, lda, P4);

    // [ D21 D22 ]     [ D2 ]
    // [ V21 V22 ]  <- [ V2 ] Q4^T
    // [  0   0  ]     [  0 ]
    // [ O1   O2 ]     [  O ]
		  TASK(READ(Fi, Q4), NOWRITE(), READWRITE(A2), papplyP, Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M2+R3, 0, N-N2-R2, A2+R2, lda, Q4);
		  //applyP( Fi, FflasRight, FflasTrans, M2+R3, 0, N-N2-R2, A2+R2, lda, Q4);

		  size_t* MathP = new size_t[M];
    // P <- Diag (P1 [ I_R1    ] , P3 [ I_R3    ])
    //               [      P2 ]      [      P4 ]
		  composePermutationsP (MathP, P1, P2, R1, M2);
		  composePermutationsP (MathP+M2, P3, P4, R3, M-M2);
		  for (size_t i=M2; i<M; ++i)
			  MathP[i] += M2;

		  WAIT;
		  if (R1+R2 < M2){
			  // A <-  S^T A
			  //TASK(READ(R1, R2, R3, R4), NOWRITE(), READWRITE(A), papplyS, F, A, lda, N, M2, R1, R2, R3, R4);
			  MatrixApplyS(F, A, lda, N, M2, R1, R2, R3, R4);

			  // P <- P S
			  PermApplyS( MathP, 1,1, M2, R1, R2, R3, R4);
		  }

		  // Q<- Diag ( [ I_R1    ] Q1,  [ I_R2    ] Q2 )
		  //            [      Q3 ]      [      P4 ]
		  size_t * MathQ = new size_t [N];
		  composePermutationsQ (MathQ, Q1, Q3, R1, N2);
		  composePermutationsQ (MathQ+N2, Q2, Q4, R2, N-N2);
		  for (size_t i=N2; i<N; ++i)
			  MathQ[i] += N2;

		  if (R1 < N2){
			  // Q <- T Q
			  PermApplyT (MathQ, 1,1,N2, R1, R2, R3, R4);
			  // A <-   A T^T
			  //TASK(READ(R1, R2, R3, R4), NOWRITE(), READWRITE(A), papplyT, F, A, lda, M, N2, R1, R2, R3, R4);
			  MatrixApplyT(F, A, lda, M, N2, R1, R2, R3, R4);
		  }
		  MathPerm2LAPACKPerm (Q, MathQ, N);
		  MathPerm2LAPACKPerm (P, MathP, M);
		  WAIT;

    delete[] MathQ;
    delete[] MathP;
    delete[] P1;
    delete[] P2;
    delete[] P3;
    delete[] P4;
    delete[] Q1;
    delete[] Q2;
    delete[] Q3;
    delete[] Q4;
    return R1+R2+R3+R4;

	  }

}// namespace FFPACK

#endif
#endif // __FFLASFFPACK_ffpack_ppluq_INL
