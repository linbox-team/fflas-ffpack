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
#define MEMCOPY
#define LEFTLOOKING
	template<class Field>
	inline size_t
	PLUQ_basecaseV2 (const Field& Fi, const FFLAS::FFLAS_DIAG Diag,
			 const size_t M, const size_t N,
			 typename Field::Element * A, const size_t lda, size_t*P, size_t *Q)
	{
		typedef typename Field::Element Element;
		size_t row = 0;
		size_t col = 0;
		size_t rank = 0;
		std::vector<bool> pivotRows(M,false);
		std::vector<bool> pivotCols(N,false);
		size_t * MathP = new size_t[M];
		size_t * MathQ = new size_t[N];
		// size_t npp=0;
		// size_t npq=0;

#ifdef LEFTLOOKING
		Element* Ltemp = new Element[M*N];
		for (size_t i=0; i<M*N; ++i)
			Fi.assign(Ltemp[i],Fi.zero);
		// this is C99 (-Wno-vla)
		Element vtemp[M];
#endif
		while ((col < N)||(row < M)){
			size_t piv2 = 0;
			size_t piv3 = 0;
// 			Element * A1 = A;
			Element * A2 = A + col;
			Element * A3 = A + row*lda;
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
				Element * vtemp_it = vtemp +rank;
				for (size_t i=0; i<M; ++i)
					if (!pivotRows[i])
						Fi.assign (*(vtemp_it++), A2[i*lda]);
				    // Left looking update
				ftrsv (Fi, FFLAS::FflasLower, FFLAS::FflasNoTrans,
				       (Diag==FFLAS::FflasUnit)?FFLAS::FflasNonUnit:FFLAS::FflasUnit,
				       rank, Ltemp, N, vtemp, 1);
				fgemv (Fi, FFLAS::FflasNoTrans, M-rank, rank, Fi.mOne,
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
			Element invpiv;
			Fi.inv (invpiv, A3[piv3]);
			if (Diag==FFLAS::FflasUnit){
#ifndef LEFTLOOKING
				    // Normalizing the pivot row
				// for (size_t i=piv3+1; i<N; ++i)
					// Fi.mulin (A3[i], invpiv);
				FFLAS::fscalin(Fi,N-piv3-1,invpiv,A3+piv3+1,1);
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
				Element * Lt_it = Ltemp + rank*(N+1) + N;
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
				FFLAS::fscalin(Fi,M-row,invpiv,A2+row*lda,lda);
				// for (size_t i=row; i<M; ++i)
					// Fi.mulin (A2 [i*lda], invpiv);
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
		delete[] MathQ;
		applyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M, 0, N, A, lda, Q);

		MathPerm2LAPACKPerm (P, MathP, M);
		delete[] MathP;
		applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N, 0, M, A, lda, P);

		return rank;
	}
// Premiere tentative de Crout avortee: trop de copies compact<->disperse et 2 cyclic-shift
#if 0
	template<class Field>
	inline size_t
	PLUQ_basecaseCrout (const Field& Fi, const FFLAS::FFLAS_DIAG Diag,
			    const size_t M, const size_t N,
			    typename Field::Element * A, const size_t lda, size_t*P, size_t *Q)
	{
		typedef typename Field::Element Element;
		size_t row = 0;
		size_t col = 0;
		size_t rank = 0;
		std::vector<bool> pivotRows(M,false);
		std::vector<bool> pivotCols(N,false);
		size_t * MathP = new size_t[M];
		size_t * MathQ = new size_t[N];
		    // Elimination takes place on Acop
		    // The Compact L\U output will be progressively stored in A
		Element* Acop = new Element[M*N];
		size_t ldac=N;
		Element* Aci=Acop;
		for (size_t i=0; i<M; ++i)
			for (size_t j=0; i<N; ++j, ++Aci)
				Fi.assign (*Aci, A[i*lda+j]);
		// this is C99 (-Wno-vla)
		while ((col < N)||(row < M)){
			size_t piv2 = 0;
			size_t piv3 = 0;
			Element * Ac2 = Acop + col;
			Element * Ac3 = Acop + row*ldac;
			Element * A12 = A+rank;
			Element * A21 = A+rank*lda;
			Element * A22 = A21+rank;
			if (row==M){
				piv3=col;
			}else
				while ((piv3 < col) && (pivotCols[piv3] || Fi.isZero (Ac3 [piv3]))) piv3++;
			if (piv3 == col){
				    // No pivot found in bottom row -> need to update one column and search in it
				if (col==N){
					row++;
					continue;
				}
				    // Copying the U part of the column
				for (size_t i=0; i<rank; ++i)
					Fi.assign (A12 [i*lda], Ac2 [MathP[i]*ldac]);
				    // Copying the lower part to be updated
				Element * A_it = A22;
				for (size_t i=0; i<row; ++i)
					if (!pivotRows[i])
						Fi.assign (*(A_it++), Ac2[i*ldac]);
				for (size_t i=row; i<M; ++i)
					Fi.assign (*(A_it++), Ac2[i*ldac]);

				    //Update
				fgemv (Fi, FFLAS::FflasNoTrans, M-rank, rank, Fi.mOne,
				       A21, lda, A12, lda, Fi.one, A22, lda);

				// Copying back the updated column to Acop
				// could be avoided: search the min index of pivot directly on A22

				A_it = A22;
				for (size_t i=0; i<M; ++i)
					if (!pivotRows[i]){
						Fi.assign (Ac2[i*ldac], *(A_it));
						A_it+=lda;
					}
				while ((piv2 < row) && (pivotRows[piv2] || Fi.isZero (Ac2 [piv2*ldac]))) piv2++;
				if (col<N) col++;
				if (piv2==M)
					continue;
			} else
				piv2 = row;

			if (row<M)  row++;
			if (Fi.isZero (Acop [piv2*lda+piv3])){
				    // no pivot found
				continue;
			}
			    // At this point the pivot is located at x=piv2 y = piv3
			Ac2 = Acop+piv3;
			Ac3 = Acop+piv2*ldac;
			MathQ[rank] = piv3;
			MathP[rank] = piv2;
			pivotCols[piv3] = true;
			pivotRows[piv2] = true;
			    //Applying the permutation on L
			    // finding the row idx of row piv2 in L stored in A
			size_t Lpiv2=1;
			    // This value should be maintained, not computed every time!!!
			if (piv2==row)
				Lpiv2 = row-rank+1;
			else
				for (size_t i=0; i<piv2; ++i)
					if (!pivotRows[i])
						Lpiv2 ++;
			cyclic_shift_row(A21, Lpiv2, rank+1, lda);

			Element invpiv;
			Fi.inv (invpiv, Ac3[piv3]);
			if (Diag==FFLAS::FflasNonUnit){
				    // Normalizing the pivot column
				Element * L_it = A22 + lda;
				for (size_t i=0; i<row+1; ++i)
					if (!pivotRows[i]){
						Fi.assign (*L_it, Fi.mulin (Ac2 [i*ldac], invpiv));
						L_it+= lda;
					}

				for (size_t i=row+1; i<M; ++i){
					Fi.assign (*L_it,Fi.mulin (Ac2 [i*ldac], invpiv));
					Ldt_it+=lda;
				}
			}
			    //Copying the new row of U to the compact storage in A
			Element* U_it = A22+1;
			for (size_t i=0; i<col+1; ++i)
				if (!pivotCols[i])
					Fi.assign (*(U_it++), Ac3 [i]);
			for (size_t i=col+1; i<N; ++i)
				Fi.assign (*(U_it++), Ac3 [i]);


		      !!!!
			    //Then Updating it
			fgemv (Fi, FFLAS::FflasNoTrans, rank, N-rank-1, Fi.mOne, A12, lda, A21, 1, Fi.one, A22+1, 1);


			    // Divinding this new row by L
			if (Diag==FFLAS::FflasUnit){
				    // Normalizing the pivot row
				for (size_t i=piv3+1; i<N; ++i)
					Fi.assign (A21[i], Fi.mulin (Ac3[i], invpiv));
			}

			    // Need to update the cols already updated
			if (piv3<col)
				for (size_t i=piv2+1; i<M; ++i)
					for (size_t j=piv3+1; j<col; ++j)
						if (!pivotCols[j])
							Fi.assign(!!!!, Fi.maxpyin (Acop[i*lda+j], Ac2[i*ldac], Ac3[j]));
// A faire aussi dans A
			rank++;
		}
		delete[] Acop;
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
		delete[] MathQ;

		MathPerm2LAPACKPerm (P, MathP, M);
		delete[] MathP;

		return rank;
	}
#endif

    // First implem of base case
	template<class Field>
	inline size_t
	PLUQ_basecase (const Field& Fi, const FFLAS::FFLAS_DIAG Diag,
		       const size_t M, const size_t N,
		       typename Field::Element * A, const size_t lda, size_t*P, size_t *Q)
	{
		typedef typename Field::Element Element;
		size_t row = 0;
		size_t col = 0;
		size_t rank = 0;
		for (size_t i=0; i<M; ++i) P[i] = i;
		for (size_t i=0; i<N; ++i) Q[i] = i;
		while ((col < N)||(row < M)){
			size_t piv2 = rank;
			size_t piv3 = rank;
			Element * A1 = A + rank*lda;
			Element * A2 = A + col;
			Element * A3 = A + row*lda;
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
				ftrsv (Fi, FFLAS::FflasLower, FFLAS::FflasNoTrans,
				       (Diag==FFLAS::FflasUnit)?FFLAS::FflasNonUnit:FFLAS::FflasUnit,
				       rank, A, lda, A2, lda);
				fgemv (Fi, FFLAS::FflasNoTrans, M-rank, rank, Fi.mOne,
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
			Element invpiv;
			Fi.inv (invpiv, A3[piv3]);
			if (Diag==FFLAS::FflasUnit){
#ifdef LEFTLOOKING
				    // Normalizing the pivot row
				// for (size_t i=piv3+1; i<N; ++i)
					// Fi.mulin (A3[i], invpiv);
				FFLAS::fscalin(Fi,N-piv3-1,invpiv,A3+piv3+1,1);
#endif
			}
			else
				    // Normalizing the pivot column
				// for (size_t i=piv2+1; i<M; ++i)
					// Fi.mulin (A2 [i*lda], invpiv);
				FFLAS::fscalin(Fi,M-piv2-1,invpiv,A2+(piv2+1)*lda,lda);
			    // Update
#ifndef LEFTLOOKING
			for (size_t i=piv2+1; i<M; ++i)
			 	for (size_t j=piv3+1; j<N; ++j)
					Fi.maxpyin (A[i*lda+j], A2[i*lda], A3[j]);
#endif
			    //Swapping pivot column
			if (piv3 > rank)
				for (size_t i=0; i<M; ++i){
					Element tmp;
					Fi.assign (tmp, A[i*lda + rank]);
					Fi.assign (A[i*lda + rank], A2[i*lda]);
					Fi.assign (A2[i*lda], tmp);
				}
				    // Updating cols

			    //Swapping pivot row
			if (piv2 > rank)
				for (size_t i=0; i<N; ++i){
					Element tmp;
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
