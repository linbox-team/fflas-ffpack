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
#define MEMCOPY
#ifndef MIN
#define MIN(a,b) (a<b)?a:b
#endif
#ifndef MAX
#define MAX(a,b) (a<b)?b:a
#endif

//#define LEFTLOOKING
#define BASECASE_K 128
using namespace std;
namespace FFPACK {
    using namespace FFLAS;

	template<class Field>
	inline size_t
	PLUQ_basecaseV3 (const Field& Fi, const FFLAS_DIAG Diag,
		       const size_t M, const size_t N,
		       typename Field::Element * A, const size_t lda, size_t*P, size_t *Q){
		size_t row = 0;
		size_t col = 0;
		size_t rank = 0;
		size_t * MathP = new size_t[M];
		size_t * MathQ = new size_t[N];

		for (size_t i=0; i<M; ++i) MathP[i] = i;
		for (size_t i=0; i<N; ++i) MathQ[i] = i;
		for (size_t i=0; i<M; ++i) P[i] = i;
		for (size_t i=0; i<N; ++i) Q[i] = i;
		while ((col < N)||(row < M)){
			size_t piv2 = rank;
			size_t piv3 = rank;
			typename Field::Element * A1 = A + rank*lda;
			typename Field::Element * A2 = A + col;
			typename Field::Element * A3 = A + row*lda;
			    // search for pivot in A2
			if (row==M){
				piv3=col;
			}else
				while ((piv3 < col) && Fi.isZero (A3 [piv3])) piv3++;
			if (piv3 == col){
				if (col==N){
					row++;
					continue;
				}
#ifdef LEFTLOOKING
				    // Left looking style update 
				ftrsv (Fi, FflasLower, FflasNoTrans, 
				       (Diag==FflasUnit)?FflasNonUnit:FflasUnit,
				       rank, A, lda, A2, lda);
				fgemv (Fi, FflasNoTrans, M-rank, rank, Fi.mOne, 
				       A1,lda, A2, lda,
				       Fi.one, A2+rank*lda, lda);
#endif
				while ((piv2 < row) && Fi.isZero (A2 [piv2*lda])) piv2++;
				if (col<N) col++;
				if (piv2==M)
					continue;
			} else 
				piv2 = row;
			if (row<M)  row++;
			if (Fi.isZero (A [piv2*lda+piv3])){
				    // no pivot found
				    //cerr<<endl;
				continue;
			}
			    // At this point the pivot is located at x=piv2 y = piv3
//			P [rank] = piv2;
//			Q [rank] = piv3;			
			
			A2 = A+piv3;
			A3 = A+piv2*lda;
			
			    // update permutations (cyclic shift)


			if(piv2 > rank)
				cyclic_shift_mathPerm(MathP+rank, piv2-rank+1);
			
			if(piv3 > rank)
				cyclic_shift_mathPerm(MathQ+rank, piv3-rank+1);

			typename Field::Element invpiv;
			Fi.inv (invpiv, A3[piv3]);
			if (Diag==FflasUnit){
#ifdef LEFTLOOKING
				    // Normalizing the pivot row
				for (size_t i=piv3+1; i<N; ++i)
					Fi.mulin (A3[i], invpiv);
#endif
			}
			else
				    // Normalizing the pivot column
				for (size_t i=piv2+1; i<M; ++i)
					Fi.mulin (A2 [i*lda], invpiv);
			    // Update
#ifndef LEFTLOOKING
			for (size_t i=piv2+1; i<M; ++i)
			 	for (size_t j=piv3+1; j<N; ++j)
					Fi.maxpyin (A[i*lda+j], A2[i*lda], A3[j]);
#endif


			    // cyclic shift pivot column and row
			if(piv3 > rank || piv2 > rank)
			{

				cyclic_shift_row_col(A+rank*(1+lda), piv2-rank+1, piv3-rank+1, lda);

				cyclic_shift_row(A+rank*lda, piv2-rank+1, rank, lda);
				cyclic_shift_row(A+rank*lda+piv3+1, piv2-rank+1, N-1-piv3, lda);
				cyclic_shift_col(A+rank, rank, piv3-rank+1, lda);
				cyclic_shift_col(A+rank+(piv2+1)*lda, M-1-piv2, piv3-rank+1, lda);
			}


/*

			if(piv2 > rank)
				cyclic_shift_row(A+rank*lda, piv2-rank+1, N, lda);
			if(piv3 > rank)
				cyclic_shift_col(A+rank, M, piv3-rank+1, lda);
*/			

#ifdef LEFTLOOKING
			    // Need to update the cols already updated
			for (size_t i=piv2+1; i<M; ++i)
				for (size_t j=piv3+1; j<col; ++j)
					Fi.maxpyin (A[i*lda+j], 
						    A[i*lda+rank], A[rank*lda+j]);				
#endif
			rank++;
		}
		MathPerm2LAPACKPerm(P, MathP, M);
		MathPerm2LAPACKPerm(Q, MathQ, N);
			delete[] MathP;
			delete[] MathQ;

		return rank;
	}

	template<class Field>
	inline size_t
	PLUQ_basecaseV2 (const Field& Fi, const FFLAS_DIAG Diag,
			 const size_t M, const size_t N,
			 typename Field::Element * A, const size_t lda, size_t*P, size_t *Q){
		size_t row = 0;
		size_t col = 0;
		size_t rank = 0;
		vector<bool> pivotRows(M,false);
		vector<bool> pivotCols(N,false);
		size_t MathP[M];
		size_t MathQ[N];
		size_t npp=0;
		size_t npq=0;

#ifdef LEFTLOOKING
		typename Field::Element* Ltemp = new typename Field::Element[M*N];
		for (size_t i=0; i<M*N; ++i)
			Fi.assign(Ltemp[i],Fi.zero);
		typename Field::Element vtemp[M];
#endif				
		while ((col < N)||(row < M)){
			size_t piv2 = 0;
			size_t piv3 = 0;
			typename Field::Element * A1 = A;
			typename Field::Element * A2 = A + col;
			typename Field::Element * A3 = A + row*lda;
			    // search for pivot in A2
			if (row==M){
				piv3=col;
			}else
				while ((piv3 < col) && (pivotCols[piv3] || Fi.isZero (A3 [piv3]))) piv3++;
			if (piv3 == col){
				if (col==N){
					row++;
					continue;
				}
#ifdef LEFTLOOKING 
				for (size_t i=0; i<rank; ++i)
					Fi.assign (vtemp[i], A2 [MathP[i]*lda]);
				typename Field::Element * vtemp_it = vtemp +rank;
				for (size_t i=0; i<M; ++i)
					if (!pivotRows[i])
						Fi.assign (*(vtemp_it++), A2[i*lda]);	
				    // Left looking update 
				ftrsv (Fi, FflasLower, FflasNoTrans, 
				       (Diag==FflasUnit)?FflasNonUnit:FflasUnit,
				       rank, Ltemp, N, vtemp, 1);
				fgemv (Fi, FflasNoTrans, M-rank, rank, Fi.mOne, 
				       Ltemp + rank*N, N, 
				       vtemp, 1, Fi.one, vtemp + rank, 1);
				for (size_t i=0; i<rank; ++i)
					Fi.assign (A2 [MathP[i]*lda], vtemp[i]);
				vtemp_it = vtemp +rank;
				for (size_t i=0; i<M; ++i)
					if (!pivotRows[i])
						Fi.assign (A2[i*lda], *(vtemp_it++));
#endif
				while ((piv2 < row) && (pivotRows[piv2] || Fi.isZero (A2 [piv2*lda]))) piv2++;
				if (col<N) col++;
				if (piv2==M)
					continue;
			} else 
				piv2 = row;
				       
			if (row<M)  row++;
			if (Fi.isZero (A [piv2*lda+piv3])){
				    // no pivot found
				continue;
			}
			    // At this point the pivot is located at x=piv2 y = piv3
			A2 = A+piv3;
			A3 = A+piv2*lda;
			MathQ[rank] = piv3;
			MathP[rank] = piv2;
			pivotCols[piv3] = true;	       
			pivotRows[piv2] = true;	       
			typename Field::Element invpiv;
			Fi.inv (invpiv, A3[piv3]);
			if (Diag==FflasUnit){
#ifndef LEFTLOOKING
				    // Normalizing the pivot row
				for (size_t i=piv3+1; i<N; ++i)
					Fi.mulin (A3[i], invpiv);
#endif
			}
			else{
#ifdef LEFTLOOKING
				    // finding the row idx of row piv2 in Ltemp
				size_t Lpiv2 = rank;
				for (size_t i=0; i<piv2; ++i)
					if (!pivotRows[i])
						Lpiv2 ++;
				if (Lpiv2 != rank)
					cyclic_shift_row(Ltemp+rank*N, Lpiv2-rank+1, rank, N);
				    // Normalizing the pivot column
				typename Field::Element * Lt_it = Ltemp + rank*(N+1) + N;
				for (size_t i=0; i<row; ++i)
					if (!pivotRows[i]){
						Fi.assign (*Lt_it, Fi.mulin (A2 [i*lda], invpiv));
						Lt_it+= N;
					}
				
				for (size_t i=row; i<M; ++i){
					Fi.assign (*Lt_it,Fi.mulin (A2 [i*lda], invpiv));
					Lt_it+=N;
				}
#else
				      // Normalizing the pivot column 
				for (size_t i=piv2+1; i<row; ++i) 
					if (!pivotRows[i]) 
						Fi.mulin (A2 [i*lda], invpiv); 
				for (size_t i=row; i<M; ++i) 
					Fi.mulin (A2 [i*lda], invpiv); 
#endif
			}
			    // Update
#ifndef LEFTLOOKING
			for (size_t i=piv2+1; i<row; ++i)
				if (!pivotRows[i]){
					for (size_t j=piv3+1; j<col; ++j)
						if (!pivotCols[j])
							Fi.maxpyin (A[i*lda+j], A2[i*lda], A3[j]);
					for (size_t j=col; j<N; ++j)
						Fi.maxpyin (A[i*lda+j], A2[i*lda], A3[j]);
				}
			for (size_t i=row; i<M; ++i){
				for (size_t j=piv3+1; j<col; ++j)
					if (!pivotCols[j])
						Fi.maxpyin (A[i*lda+j], A2[i*lda], A3[j]);
				for (size_t j=col; j<N; ++j)
					Fi.maxpyin (A[i*lda+j], A2[i*lda], A3[j]);
			}
#endif
#ifdef LEFTLOOKING
			    // Need to update the cols already updated
			if (piv3<col)
				for (size_t i=piv2+1; i<M; ++i)
					for (size_t j=piv3+1; j<col; ++j)
						if (!pivotCols[j])
							Fi.maxpyin (A[i*lda+j], A2[i*lda], A3[j]);
#endif
			rank++; 
		}
#ifdef LEFTLOOKING
		delete[] Ltemp;
#endif
		    // Building permutations
		 size_t nonpiv = rank;
		 for (size_t i = 0; i<M; ++i)
			 if (!pivotRows[i])
				 MathP[nonpiv++] = i;
		 nonpiv = rank;
		 for (size_t i = 0; i<N; ++i)
			 if (!pivotCols[i])
				 MathQ[nonpiv++] = i;
		MathPerm2LAPACKPerm (Q, MathQ, N);
		applyP (Fi, FflasRight, FflasTrans, M, 0, N, A, lda, Q);
		MathPerm2LAPACKPerm (P, MathP, M);
		applyP (Fi, FflasLeft, FflasNoTrans, N, 0, M, A, lda, P);
		return rank;
	}

	template<class Field>
	inline size_t
	PLUQ_basecase (const Field& Fi, const FFLAS_DIAG Diag,
		       const size_t M, const size_t N,
		       typename Field::Element * A, const size_t lda, size_t*P, size_t *Q){
		size_t row = 0;
		size_t col = 0;
		size_t rank = 0;
		for (size_t i=0; i<M; ++i) P[i] = i;
		for (size_t i=0; i<N; ++i) Q[i] = i;
		while ((col < N)||(row < M)){
			size_t piv2 = rank;
			size_t piv3 = rank;
			typename Field::Element * A1 = A + rank*lda;
			typename Field::Element * A2 = A + col;
			typename Field::Element * A3 = A + row*lda;
			    // search for pivot in A2
			if (row==M){
				piv3=col;
			}else
				while ((piv3 < col) && Fi.isZero (A3 [piv3])) piv3++;
			if (piv3 == col){
				if (col==N){
					row++;
					continue;
				}
#ifdef LEFTLOOKING
				    // Left looking style update 
				ftrsv (Fi, FflasLower, FflasNoTrans, 
				       (Diag==FflasUnit)?FflasNonUnit:FflasUnit,
				       rank, A, lda, A2, lda);
				fgemv (Fi, FflasNoTrans, M-rank, rank, Fi.mOne, 
				       A1,lda, A2, lda,
				       Fi.one, A2+rank*lda, lda);
#endif
				while ((piv2 < row) && Fi.isZero (A2 [piv2*lda])) piv2++;
				if (col<N) col++;
				if (piv2==M)
					continue;
			} else 
				piv2 = row;
			if (row<M)  row++;
			if (Fi.isZero (A [piv2*lda+piv3])){
				    // no pivot found
				    //cerr<<endl;
				continue;
			}
			    // At this point the pivot is located at x=piv2 y = piv3
			P [rank] = piv2;
			Q [rank] = piv3;
			A2 = A+piv3;
			A3 = A+piv2*lda;
			typename Field::Element invpiv;
			Fi.inv (invpiv, A3[piv3]);
			if (Diag==FflasUnit){
#ifdef LEFTLOOKING
				    // Normalizing the pivot row
				for (size_t i=piv3+1; i<N; ++i)
					Fi.mulin (A3[i], invpiv);
#endif
			}
			else
				    // Normalizing the pivot column
				for (size_t i=piv2+1; i<M; ++i)
					Fi.mulin (A2 [i*lda], invpiv);
			    // Update
#ifndef LEFTLOOKING
			for (size_t i=piv2+1; i<M; ++i)
			 	for (size_t j=piv3+1; j<N; ++j)
					Fi.maxpyin (A[i*lda+j], A2[i*lda], A3[j]);
#endif
			    //Swapping pivot column
			if (piv3 > rank)
				for (size_t i=0; i<M; ++i){
					typename Field::Element tmp;
					Fi.assign (tmp, A[i*lda + rank]);
					Fi.assign (A[i*lda + rank], A2[i*lda]);
					Fi.assign (A2[i*lda], tmp);
				}
				    // Updating cols 

			    //Swapping pivot row
			if (piv2 > rank)
				for (size_t i=0; i<N; ++i){
					typename Field::Element tmp;
					Fi.assign (tmp, A1[i]);
					Fi.assign (A1[i], A3[i]);
					Fi.assign (A3[i], tmp);
				}
#ifdef LEFTLOOKING
			    // Need to update the cols already updated
			for (size_t i=piv2+1; i<M; ++i)
				for (size_t j=piv3+1; j<col; ++j)
					Fi.maxpyin (A[i*lda+j], 
						    A[i*lda+rank], A[rank*lda+j]);				
#endif
			rank++; 
		}
		return rank;
	}
	
	template<class Field>
	inline size_t
	PLUQ (const Field& Fi, const FFLAS_DIAG Diag,
	      const size_t M, const size_t N,
	      typename Field::Element * A, const size_t lda, size_t*P, size_t *Q){
		
		for (size_t i=0; i<M; ++i) P[i] = i;
		for (size_t i=0; i<N; ++i) Q[i] = i;
		if (MIN (M,N) == 0) return 0;
		if (MAX (M,N) == 1) return (Fi.isZero(*A))? 0 : 1;
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
			if (Diag== FflasUnit){
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
			if (Diag== FflasNonUnit){
				typename Field::Element invpivot;
				Fi.inv(invpivot, *A);
				for (size_t i=piv+1; i<M; ++i)
					Fi.mulin (*(A+i*lda), invpivot);
			}
			return 1;
		}

#ifdef BASECASE_K
		if (MIN(M,N) < BASECASE_K)
			return PLUQ_basecaseV2 (Fi, Diag, M, N, A, lda, P, Q);
#endif
		FFLAS_DIAG OppDiag = (Diag == FflasUnit)? FflasNonUnit : FflasUnit;
		size_t M2 = M >> 1;
		size_t N2 = N >> 1;
		size_t * P1 = new size_t [M2];
		size_t * Q1 = new size_t [N2];
		size_t R1,R2,R3,R4;

		    // A1 = P1 [ L1 ] [ U1 V1 ] Q1
		    //        [ M1 ]
		R1 = PLUQ (Fi, Diag, M2, N2, A, lda, P1, Q1);
		typename Field::Element * A2 = A + N2;
		typename Field::Element * A3 = A + M2*lda;
		typename Field::Element * A4 = A3 + N2;
		typename Field::Element * F = A2 + R1*lda;
		typename Field::Element * G = A3 + R1;
		    // [ B1 ] <- P1^T A2
		    // [ B2 ]
		applyP (Fi, FflasLeft, FflasNoTrans, N-N2, 0, M2, A2, lda, P1);
		    // [ C1 C2 ] <- A3 Q1^T
		applyP (Fi, FflasRight, FflasTrans, M-M2, 0, N2, A3, lda, Q1);
		    // D <- L1^-1 B1
		ftrsm (Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R1, N-N2, Fi.one, A, lda, A2, lda);
		    // E <- C1 U1^-1 
		ftrsm (Fi, FflasRight, FflasUpper, FflasNoTrans, Diag, M-M2, R1, Fi.one, A, lda, A3, lda);
		    // F <- B2 - M1 D
		fgemm (Fi, FflasNoTrans, FflasNoTrans, M2-R1, N-N2, R1, Fi.mOne, A + R1*lda, lda, A2, lda, Fi.one, A2+R1*lda, lda);
		    // G <- C2 - E V1
		fgemm (Fi, FflasNoTrans, FflasNoTrans, M-M2, N2-R1, R1, Fi.mOne, A3, lda, A+R1, lda, Fi.one, A3+R1, lda);
		    // H <- A4 - ED
		fgemm (Fi, FflasNoTrans, FflasNoTrans, M-M2, N-N2, R1, Fi.mOne, A3, lda, A2, lda, Fi.one, A4, lda);
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
		applyP (Fi, FflasRight, FflasTrans, M-M2, 0, N-N2, A4, lda, Q2);
		applyP (Fi, FflasLeft, FflasNoTrans, N-N2, 0, M-M2, A4, lda, P3);
		    // [ E1 ] <- P3^T E
		    // [ E2 ]
		applyP (Fi, FflasLeft, FflasNoTrans, R1, 0, M-M2, A3, lda, P3);
		    // [ M11 ] <- P2^T M1
		    // [ M12 ]
		applyP (Fi, FflasLeft, FflasNoTrans, R1, 0, M2-R1, A+R1*lda, lda, P2);
		    // [ D1 D2 ] <- D Q2^T
		applyP (Fi, FflasRight, FflasTrans, R1, 0, N-N2, A2, lda, Q2);
		    // [ V1 V2 ] <- V1 Q3^T
		applyP (Fi, FflasRight, FflasTrans, R1, 0, N2-R1, A+R1, lda, Q3);
		    // I <- H U2^-1
		    // K <- H3 U2^-1
		ftrsm (Fi, FflasRight, FflasUpper, FflasNoTrans, Diag, M-M2, R2, Fi.one, F, lda, A4, lda);
		    // J <- L3^-1 I (in a temp)
		typename Field::Element * temp = new typename Field::Element [R3*R2];
		for (size_t i=0; i<R3; ++i)
			fcopy (Fi, R2, temp + i*R2, 1, A4 + i*lda, 1);
		ftrsm (Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R3, R2, Fi.one, G, lda, temp, R2);
		    // N <- L3^-1 H2
		ftrsm (Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R3, N-N2-R2, Fi.one, G, lda, A4+R2, lda);
		    // O <- N - J V2
		fgemm (Fi, FflasNoTrans, FflasNoTrans, R3, N-N2-R2, R2, Fi.mOne, temp, R2, F+R2, lda, Fi.one, A4+R2, lda);
		delete[] temp;
		    // R <- H4 - K V2 - M3 O
		typename Field::Element * R = A4 + R2 + R3*lda;
		fgemm (Fi, FflasNoTrans, FflasNoTrans, M-M2-R3, N-N2-R2, R2, Fi.mOne, A4+R3*lda, lda, F+R2, lda, Fi.one, R, lda);
		fgemm (Fi, FflasNoTrans, FflasNoTrans, M-M2-R3, N-N2-R2, R3, Fi.mOne, G+R3*lda, lda, A4+R2, lda, Fi.one, R, lda);
		    // H4 = P4 [ L4 ] [ U4 V4 ] Q4
		    //         [ M4 ]
		size_t * P4 = new size_t [M-M2-R3];
		size_t * Q4 = new size_t [N-N2-R2];
		R4 = PLUQ (Fi, Diag, M-M2-R3, N-N2-R2, R, lda, P4, Q4);
		    // [ E21 M31 0 K1 ] <- P4^T [ E2 M3 0 K ] 
		    // [ E22 M32 0 K2 ]
		applyP (Fi, FflasLeft, FflasNoTrans, N2+R2, 0, M-M2-R3, A3+R3*lda, lda, P4);
		    // [ D21 D22 ]     [ D2 ]
		    // [ V21 V22 ]  <- [ V2 ] Q4^T
		    // [  0   0  ]     [  0 ]
		    // [ O1   O2 ]     [  O ]
		applyP (Fi, FflasRight, FflasTrans, M2+R3, 0, N-N2-R2, A2+R2, lda, Q4);
		
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
			applyS (MathP, 1,1,M2, R1, R2, R3, R4);
			    // A <-  S^T A
			applyS (A, lda, N, M2, R1, R2, R3, R4);
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
			applyT (MathQ, 1,1,N2, R1, R2, R3, R4);
			    // A <-   A T^T
			applyT (A, lda, M, N2, R1, R2, R3, R4);
		}
		MathPerm2LAPACKPerm (Q, MathQ, N);
		delete[] MathQ;
		
		return R1+R2+R3+R4;

	
	}

	template <class Element>
	inline void applyS (Element* A, const size_t lda, const size_t width, 
			    const size_t M2,
			    const size_t R1, const size_t R2, 
			    const size_t R3, const size_t R4){
		Element * tmp = new Element [(M2-R1-R2)*width];
		// std::cerr<<"ici"<<std::endl;
		for (size_t i = 0, j = R1+R2; j < M2; ++i, ++j)
			for (size_t k = 0; k<width; ++k)
				tmp [i*width + k] = A [j*lda + k];
		for (size_t i = M2, j = R1+R2; i < M2+R3+R4; ++i, ++j)
			for (size_t k = 0; k<width; ++k)
				A [j*lda + k] = A [i*lda +k];
		for (size_t i = 0, j = R1+R2+R3+R4; i < M2-R1-R2; ++i, ++j)
			for (size_t k = 0; k<width; ++k)
				A [j*lda + k] = tmp [i*width + k];
		delete[] tmp;
	}
	

	template <class Element>
	inline void applyT (Element* A, const size_t lda, const size_t width, 
			    const size_t N2,
			    const size_t R1, const size_t R2, 
			    const size_t R3, const size_t R4){
		Element * tmp = new Element[(N2-R1)*width];
		for (size_t k = 0; k < width; ++k){
			for (size_t i = 0, j = R1; j < N2; ++i, ++j){
				tmp [i + k*(N2-R1)] = A [k*lda + j];
			}
			
			for (size_t i = N2, j = R1; i < N2+R2; ++i, ++j)
				A [k*lda + j] = A [k*lda + i];
			
			for (size_t i = 0, j = R1+R2; i < R3; ++i, ++j)
				A [k*lda + j] = tmp [k*(N2-R1) + i];
			
			for (size_t i = N2+R2, j = R1+R2+R3; i < N2+R2+R4; ++i,++j)
				A [k*lda + j] = A [k*lda + i];
			for (size_t i = R3, j = R1+R2+R3+R4; i < N2-R1; ++i,++j)
				A [k*lda + j] = tmp [k*(N2-R1) + i];
		}
		delete[] tmp;
	}


            /**
	     * Conversion of a permutation from LAPACK format to Math format
	     */
	inline void LAPACKPerm2MathPerm (size_t * MathP, const size_t * LapackP, 
				  const size_t N){
		for (size_t i=0; i<N; i++)
			MathP[i] = i;
		for (size_t i=0; i<N; i++){
			if (LapackP[i] != i){
				size_t tmp = MathP[i];
				MathP[i] = MathP[LapackP[i]];
				MathP[LapackP[i]] = tmp;
			}
		}
	}
	    /**
	     * Conversion of a permutation from Maths format to LAPACK format
	     */
	inline void MathPerm2LAPACKPerm (size_t * LapackP, const size_t * MathP, 
					 const size_t N){
		size_t * T = new size_t[N];
		size_t * Tinv = new size_t[N]; 
		for (size_t i=0; i<N; i++){
			T[i] =i;
			Tinv[i] = i;
		}
		for (size_t i=0; i<N; i++){
			size_t j = Tinv [MathP [i]];
			LapackP [i] = j;
			size_t tmp = T[j];
			T[j] = T[i];
			Tinv[T[i]] = j;
			T[i] = tmp;
			Tinv[tmp] = i;
		}
		delete[] T;
		delete[] Tinv;
	}
	    /**
	     * Computes P1 [ I_R     ] stored in MathPermutation format
	     *             [     P_2 ]
	     */
	inline void composePermutationsP (size_t * MathP, 
				  const size_t * P1, 
				  const size_t * P2, 
				  const size_t R, const size_t N){
		for (size_t i=0; i<N; ++i)
			MathP[i] = i;
		LAPACKPerm2MathPerm (MathP, P1, N);

		for (size_t i=R; i<N; i++){
			if (P2[i-R] != i-R){
				size_t tmp = MathP[i];
				MathP[i] = MathP[P2[i-R]+R];
				MathP[P2[i-R]+R] = tmp;
			}
		}
	}
	inline void composePermutationsQ (size_t * MathP, 
				  const size_t * Q1, 
				  const size_t * Q2, 
				  const size_t R, const size_t N){
		for (size_t i=0; i<N; ++i)
			MathP[i] = i;
		LAPACKPerm2MathPerm (MathP, Q1, N);
		
		for (size_t i=R; i<N; i++){
			if (Q2[i-R] != i-R){
				size_t tmp = MathP[i];
				MathP[i] = MathP[Q2[i-R]+R];
				MathP[Q2[i-R]+R] = tmp;
			}
		}
	}
	inline void
	RankProfilesFromPLUQ (size_t* RowRankProfile, size_t* ColumnRankProfile,
			      const size_t * P, const size_t * Q, 
			      const size_t M, const size_t N, const size_t R){
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

	inline void 
	cyclic_shift_mathPerm (size_t * P,  const size_t s){
                size_t tmp;
                tmp = *(P+s-1);
		    //memmove(P+1, P, (s)*sizeof(size_t));	
		size_t * Pi = P;
		std::copy(Pi, Pi+s-1, Pi+1);

                *(P)=tmp;
	}

	template<typename Base_t>
	inline void cyclic_shift_row_col(Base_t * A, size_t m, size_t n, size_t lda) {
		    //	std::cerr << "BEF m: " << m << ", n: " << n << std::endl;
		if (m > 1) {
			const size_t mun(m-1);
			if (n > 1) {
				const size_t nun(n-1);
				
				Base_t * b = new Base_t[mun];
				Base_t * Ainun = A+nun;
				for(size_t i=0; i<mun; ++i, Ainun+=lda) b[i] = *Ainun;
				
				    // dc = [ d c ]
				Base_t * dc = new Base_t[n];            
				std::copy(Ainun-nun, Ainun, dc+1);
				
				    // this is d
				*dc = *Ainun;
				
				Base_t * Ai=A+(mun-1)*lda;
				for(size_t i=mun; i>0; --i, Ai-=lda)
					std::copy(Ai, Ai+nun, Ai+1+lda);
				
				std::copy(dc, dc+n, A);
				
				Base_t * Aipo = A+lda;
				for(size_t i=0; i<mun; ++i, Aipo+=lda) *Aipo = b[i];
				
				delete [] dc;
				delete [] b;
			} else if (n != 0) {
				Base_t * Ai=A+mun*lda;
				Base_t d = *Ai;
				for(; Ai != A; Ai-=lda) *Ai= *(Ai-lda);
				*A=d;
			}
		} else {
			if ((m!=0) && (n > 1)) {
				const size_t nun(n-1);
				Base_t d = A[nun];
				std::copy(A,A+nun,A+1);
				*A=d;
			}
		}
		
	}
	template<typename Base_t>
	inline void cyclic_shift_row(Base_t * A, size_t m, size_t n, size_t lda)
	{
		if (m > 1) {
			const size_t mun(m-1);
			const size_t nun(n-1);
			
			Base_t * b = new Base_t[n];
			Base_t * Ai = A+mun*lda;
			for(size_t i=0; i<n; ++i, Ai+=1) b[i] = *Ai;
			
			    // dc = [ d c ]
//			Base_t * dc = new Base_t[n];

			Base_t * Ac = A;
			for(int i=mun-1; i>=0; --i)
				std::copy(Ac+i*lda, Ac+i*lda+n, Ac+(i+1)*lda);
			
			Base_t * Aii = A;
			for(size_t i=0; i<n; ++i, Aii++) *Aii = b[i];
				
//			delete [] dc;
			delete [] b;
		}
	}

	template<typename Base_t>
	inline void cyclic_shift_col(Base_t * A, size_t m, size_t n, size_t lda)
	{
		if (n > 1) {
			const size_t mun(m-1);
			const size_t nun(n-1);		
			Base_t * b = new Base_t[m];
			Base_t * Ainun = A+nun;
			for(size_t i=0; i<m; ++i, Ainun+=lda) b[i] = *Ainun;
			
			    // dc = [ d c ]
			Base_t * tmp = new Base_t[nun]; 
			Base_t * Ac = A;
			for(int i=mun; i>=0; --i)
			{				
//				std::copy(A+i*lda, A+i*lda+1, A+1);
				std::copy(Ac+i*lda, Ac+i*lda+nun, tmp);
				std::copy(tmp, tmp+nun, Ac+i*lda+1);
			}
			
			Base_t * Aipo = A;
			for(size_t i=0; i<m; ++i, Aipo+=lda) *Aipo = b[i];
			delete [] tmp;
			delete [] b;
		}
	}
	
	    /*
	inline void cyclic_shift_permut (size_t * P,  const size_t N){
                size_t tmp;
                tmp = *P;
      		memmove(P, P+1, (N-1)*sizeof(size_t));
                *(P+N-1)=tmp;
	}
	    */
	
} // namespace FFPACK
#endif // __FFLASFFPACK_ffpack_pluq_INL
