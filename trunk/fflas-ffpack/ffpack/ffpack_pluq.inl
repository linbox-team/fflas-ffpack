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

#define CROUT
#define BASECASE_K 256

namespace FFPACK {

        // Base Case based on a CUP decomp with rotations
	template<class Field>
	inline size_t
	PLUQ_basecaseCrout (const Field& Fi, const FFLAS::FFLAS_DIAG Diag,
			    const size_t M, const size_t N,
			    typename Field::Element_ptr A, const size_t lda, size_t*P, size_t *Q)
	{
		int row = 0;
		int rank = 0;
		typename Field::Element_ptr CurrRow=A;
		size_t * MathP = new size_t[M];
		size_t * MathQ = new size_t[N];
		for (size_t i=0; i<M; ++i) MathP[i] = i;
		for (size_t i=0; i<N; ++i) MathQ[i] = i;
		while (((size_t)row<M) && ((size_t)rank<N)){
			    // Updating row where pivot will be searched for	
			fgemv(Fi, FFLAS::FflasTrans, rank, N-rank, Fi.mOne, A+rank, lda, CurrRow, 1, Fi.one, CurrRow+rank, 1);
			int i = rank-1;
			while(Fi.isZero (CurrRow[++i]) && (i<(int)N-1));
			if (!Fi.isZero (CurrRow[i])){
				    // found pivot
				    // Updating column below pivot
				// Q [rank] = i;
				// pivotRows [row] = true;
				// P [rank] = row;
				fgemv(Fi, FFLAS::FflasNoTrans, M-row-1, rank, Fi.mOne, CurrRow+lda, lda, A+i, lda, Fi.one, CurrRow+lda+i, lda);
				    // Normalization
				typename Field::Element invpiv;
				Fi.inv (invpiv, CurrRow[i]);
				if (Diag == FFLAS::FflasUnit)
					FFLAS::fscalin (Fi, N-i-1, invpiv, CurrRow+i+1,1);
				else
					FFLAS::fscalin (Fi, M-row-1, invpiv, CurrRow+i+lda,lda);

				if (i > rank){
					    // Column rotation to move pivot on the diagonal
					    // on U
					cyclic_shift_col(Fi, A+rank, rank, i-rank+1, lda);
					cyclic_shift_mathPerm(MathQ+rank, (size_t)(i-rank+1));
					    // on A
					cyclic_shift_col(Fi, CurrRow+lda+rank, M-row-1, i-rank+1, lda);
					Fi.assign(A[rank*(lda+1)], CurrRow[i]);
					FFLAS::fzero (Fi, i-rank, A+rank*(lda+1)+1, 1);
				}
				if (row > rank){
					    // Row rotation for L
					    // Optimization: delay this to the end
					cyclic_shift_row(Fi, A+rank*lda, row-rank+1, rank, lda);
					cyclic_shift_mathPerm(MathP+rank, (size_t) (row-rank+1) );
					    // Row rotation for U (not moving the 0 block)
					FFLAS::fcopy (Fi, N-i-1, CurrRow+i+1, 1, A+rank*lda+i+1, 1);
					Fi.assign(A[rank*(lda+1)], CurrRow[i]);
					FFLAS::fzero (Fi, row-rank, A+rank*(lda+1)+lda, lda);
					Fi.assign(CurrRow[i],Fi.zero); // only needed once here
				}
				rank++;
			}
			CurrRow+=lda;
			row++;
		}

		// size_t nonpiv = rank;
		// for (size_t i = 0; i<M; ++i)
		// 	if (!pivotRows[i])
		// 		MathP[nonpiv++] = i;
		// nonpiv = rank;
		// for (size_t i = 0; i<N; ++i)
		// 	if (!pivotCols[i])
		// 		MathQ[nonpiv++] = i;

		// std::cerr<<"MathP = ";
		// for (int i=0; i<M; ++i)
		// 	std::cerr<<MathP[i]<<" ";
		// std::cerr<<std::endl;
		// std::cerr<<"MathQ = ";
		// for (int i=0; i<N; ++i)
		// 	std::cerr<<MathQ[i]<<" ";
		// std::cerr<<std::endl;

		MathPerm2LAPACKPerm (Q, MathQ, N);
		delete[] MathQ;
		MathPerm2LAPACKPerm (P, MathP, M);
		delete[] MathP;
		FFLAS::fzero (Fi, M-rank, N-rank, A+rank*(1+lda), lda);

		return (size_t) rank;
	}


	template<class Field>
	inline size_t
	PLUQ (const Field& Fi, const FFLAS::FFLAS_DIAG Diag,
	      const size_t M, const size_t N,
	      typename Field::Element_ptr A, const size_t lda, size_t*P, size_t *Q)
	{

		for (size_t i=0; i<M; ++i) P[i] = i;
		for (size_t i=0; i<N; ++i) Q[i] = i;
		if (std::min (M,N) == 0) return 0;
		if (std::max (M,N) == 1) return (Fi.isZero(*A))? 0 : 1;
#ifndef BASECASE_K
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
				// for (size_t i=piv+1; i<N; ++i)
					// Fi.mulin (A[i], invpivot);
				FFLAS::fscalin(Fi,N-piv-1,invpivot,A+piv+1,1)
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
				// for (size_t i=piv+1; i<M; ++i)
					// Fi.mulin (*(A+i*lda), invpivot);
				FFLAS::fscalin(Fi,M-piv-1,invpivot,A+(piv+1)*lda,lda);
			}
			return 1;
		}
#endif
#ifdef BASECASE_K
		if (std::min(M,N) < BASECASE_K)
			return PLUQ_basecaseCrout (Fi, Diag, M, N, A, lda, P, Q);
#endif
		FFLAS::FFLAS_DIAG OppDiag = (Diag == FFLAS::FflasUnit)? FFLAS::FflasNonUnit : FFLAS::FflasUnit;
		size_t M2 = M >> 1;
		size_t N2 = N >> 1;
		size_t * P1 = new size_t [M2];
		size_t * Q1 = new size_t [N2];
		size_t R1,R2,R3,R4;

		    // A1 = P1 [ L1 ] [ U1 V1 ] Q1
		    //         [ M1 ]
		R1 = PLUQ (Fi, Diag, M2, N2, A, lda, P1, Q1);
		typename Field::Element_ptr A2 = A + N2;
		typename Field::Element_ptr A3 = A + M2*lda;
		typename Field::Element_ptr A4 = A3 + N2;
		typename Field::Element_ptr F = A2 + R1*lda;
		typename Field::Element_ptr G = A3 + R1;
		    // [ B1 ] <- P1^T A2
		    // [ B2 ]
		applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N-N2, 0, M2, A2, lda, P1);
		    // [ C1 C2 ] <- A3 Q1^T
		applyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M-M2, 0, N2, A3, lda, Q1);
		    // D <- L1^-1 B1
		ftrsm (Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R1, N-N2, Fi.one, A, lda, A2, lda);
		    // E <- C1 U1^-1
		ftrsm (Fi, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag, M-M2, R1, Fi.one, A, lda, A3, lda);
		    // F <- B2 - M1 D
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M2-R1, N-N2, R1, Fi.mOne, A + R1*lda, lda, A2, lda, Fi.one, A2+R1*lda, lda);
		    // G <- C2 - E V1
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2, N2-R1, R1, Fi.mOne, A3, lda, A+R1, lda, Fi.one, A3+R1, lda);
		    // H <- A4 - ED
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2, N-N2, R1, Fi.mOne, A3, lda, A2, lda, Fi.one, A4, lda);
                    // F = P2 [ L2 ] [ U2 V2 ] Q2
		    //        [ M2 ]
		size_t * P2 = new size_t [M2-R1];
		size_t * Q2 = new size_t [N-N2];
		R2 = PLUQ (Fi, Diag, M2-R1, N-N2, F, lda, P2, Q2);
		    // G = P3 [ L3 ] [ U3 V3 ] Q3
		    //        [ M3 ]
		size_t * P3 = new size_t [M-M2];
		size_t * Q3 = new size_t [N2-R1];
		R3 = PLUQ (Fi, Diag, M-M2, N2-R1, G, lda, P3, Q3);
		    // [ H1 H2 ] <- P3^T H Q2^T
		    // [ H3 H4 ]
		applyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M-M2, 0, N-N2, A4, lda, Q2);
		applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N-N2, 0, M-M2, A4, lda, P3);
		    // [ E1 ] <- P3^T E
		    // [ E2 ]
		applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, R1, 0, M-M2, A3, lda, P3);
		    // [ M11 ] <- P2^T M1
		    // [ M12 ]
		applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, R1, 0, M2-R1, A+R1*lda, lda, P2);
		    // [ D1 D2 ] <- D Q2^T
		applyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, R1, 0, N-N2, A2, lda, Q2);
		    // [ V1 V2 ] <- V1 Q3^T
		applyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, R1, 0, N2-R1, A+R1, lda, Q3);
		    // I <- H U2^-1
		    // K <- H3 U2^-1
		ftrsm (Fi, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag, M-M2, R2, Fi.one, F, lda, A4, lda);
		    // J <- L3^-1 I (in a temp)
		typename Field::Element_ptr temp = FFLAS::fflas_new (Fi, R3, R2);
		// for (size_t i=0; i<R3; ++i)
			// FFLAS::fcopy (Fi, R2, A4 + i*lda, 1, temp + i*R2, 1);
		FFLAS::fcopy (Fi, R3, R2, A4 , lda, temp , R2);
		ftrsm (Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R3, R2, Fi.one, G, lda, temp, R2);
		    // N <- L3^-1 H2
		ftrsm (Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R3, N-N2-R2, Fi.one, G, lda, A4+R2, lda);
		    // O <- N - J V2
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, R3, N-N2-R2, R2, Fi.mOne, temp, R2, F+R2, lda, Fi.one, A4+R2, lda);
		FFLAS::fflas_delete (temp);
		    // R <- H4 - K V2 - M3 O
		typename Field::Element_ptr R = A4 + R2 + R3*lda;
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2-R3, N-N2-R2, R2, Fi.mOne, A4+R3*lda, lda, F+R2, lda, Fi.one, R, lda);
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2-R3, N-N2-R2, R3, Fi.mOne, G+R3*lda, lda, A4+R2, lda, Fi.one, R, lda);
		    // H4 = P4 [ L4 ] [ U4 V4 ] Q4
		    //         [ M4 ]
		size_t * P4 = new size_t [M-M2-R3];
		size_t * Q4 = new size_t [N-N2-R2];
		R4 = PLUQ (Fi, Diag, M-M2-R3, N-N2-R2, R, lda, P4, Q4);
		    // [ E21 M31 0 K1 ] <- P4^T [ E2 M3 0 K ]
		    // [ E22 M32 0 K2 ]
		applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N2+R2, 0, M-M2-R3, A3+R3*lda, lda, P4);
		    // [ D21 D22 ]     [ D2 ]
		    // [ V21 V22 ]  <- [ V2 ] Q4^T
		    // [  0   0  ]     [  0 ]
		    // [ O1   O2 ]     [  O ]
		applyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M2+R3, 0, N-N2-R2, A2+R2, lda, Q4);

		    // P <- Diag (P1 [ I_R1    ] , P3 [ I_R3    ])
		    //               [      P2 ]      [      P4 ]
		size_t* MathP = new size_t[M];
		composePermutationsP (MathP, P1, P2, R1, M2);
		composePermutationsP (MathP+M2, P3, P4, R3, M-M2);
		delete[] P1;
		delete[] P2;
		delete[] P3;
		delete[] P4;
		for (size_t i=M2; i<M; ++i)
			MathP[i] += M2;
		if (R1+R2 < M2){
			    // P <- P S
			PermApplyS (MathP, 1,1,M2, R1, R2, R3, R4);
			    // A <-  S^T A
			MatrixApplyS (Fi, A, lda, N, M2, R1, R2, R3, R4);
		}
		MathPerm2LAPACKPerm (P, MathP, M);
		delete[] MathP;

		    // Q<- Diag ( [ I_R1    ] Q1,  [ I_R2    ] Q2 )
		    //            [      Q3 ]      [      P4 ]
		size_t * MathQ = new size_t [N];
		composePermutationsQ (MathQ, Q1, Q3, R1, N2);
		composePermutationsQ (MathQ+N2, Q2, Q4, R2, N-N2);
		delete[] Q1;
		delete[] Q2;
		delete[] Q3;
		delete[] Q4;
		for (size_t i=N2; i<N; ++i)
			MathQ[i] += N2;

		if (R1 < N2){
			    // Q <- T Q
			PermApplyT (MathQ, 1,1,N2, R1, R2, R3, R4);
			    // A <-   A T^T
			MatrixApplyT (Fi, A, lda, M, N2, R1, R2, R3, R4);
		}
		MathPerm2LAPACKPerm (Q, MathQ, N);
		delete[] MathQ;

		return R1+R2+R3+R4;
	}

	inline void
	RankProfilesFromPLUQ (size_t* RowRankProfile, size_t* ColumnRankProfile,
			      const size_t * P, const size_t * Q,
			      const size_t M, const size_t N, const size_t R)
	{
		size_t * RRP=new size_t[M];
		size_t * CRP=new size_t[N];
		for (size_t i=0;i < M; ++i)
			RRP [i]=i;
		for (size_t i=0;i < N; ++i)
			CRP [i]=i;

		// std::cerr<<"CRP = ";
		// for (size_t i=0; i<N; ++i)
		// 	std::cerr<<CRP[i]<<" ";
		// std::cerr<<std::endl;
		for (size_t i=0; i<M; ++i)
			if (P[i] != i){
				size_t tmp = RRP [i];
				RRP [i] = RRP [P [i]];
				RRP [P [i]] = tmp;
			}
		for (size_t i=0; i<R; ++i)
			RowRankProfile[i] = RRP[i];
		std::sort(RowRankProfile,RowRankProfile+R);
		for (size_t i=0; i<N; ++i){
			// std::cerr<<"Q["<<i<<"] ="<<Q[i]<<std::endl;
			if (Q[i] != i){
				size_t tmp = CRP [i];
				CRP [i] = CRP [Q [i]];
				CRP [Q [i]] = tmp;
			}
		}
		for (size_t i=0; i<R; ++i)
			ColumnRankProfile [i] = CRP [i];
		std::sort(ColumnRankProfile,ColumnRankProfile+R);
		// std::cerr<<"CRP = ";
		// for (size_t i=0; i<N; ++i)
		// 	std::cerr<<CRP[i]<<" ";
		// std::cerr<<std::endl;
		delete[] RRP;
		delete[] CRP;
	}


} // namespace FFPACK
#endif // __FFLASFFPACK_ffpack_pluq_INL
