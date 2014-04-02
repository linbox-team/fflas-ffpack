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
#define LEFTLOOKING
#define CROUT
#define BASECASE_K 256

namespace FFPACK {


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
		delete[] MathQ;
		applyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M, 0, N, A, lda, Q);

		MathPerm2LAPACKPerm (P, MathP, M);
		delete[] MathP;
		applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N, 0, M, A, lda, P);

		return rank;
	}

        // Base Case based on a CUP decomp with rotations
	template<class Field>
	inline size_t
	PLUQ_basecaseCrout (const Field& Fi, const FFLAS::FFLAS_DIAG Diag,
			    const size_t M, const size_t N,
			    typename Field::Element * A, const size_t lda, size_t*P, size_t *Q)
	{
		typedef typename Field::Element Element;
		int row = 0;
		int rank = 0;
		Element* CurrRow=A;
		size_t * MathP = new size_t[M];
		size_t * MathQ = new size_t[N];
		for (size_t i=0; i<M; ++i) MathP[i] = i;
		for (size_t i=0; i<N; ++i) MathQ[i] = i;
		while (((size_t)row<M) && ((size_t)rank<N)){
			    // Updating row where pivot will be searched for
			fgemv(Fi, FFLAS::FflasTrans, rank, N-rank, Fi.mOne, A+rank, lda, CurrRow, 1, Fi.one, CurrRow+rank, 1);
//			write_field(Fi,std::cerr<<"Avant la cherche"<<std::endl,A,M,N,lda);
			int i = rank-1;
			while(Fi.isZero (CurrRow[++i]) && (i<N-1));
//			std::cout<<"row = "<<row<<" i = "<<i<<" -> A[row,i] = "<<CurrRow[i]<<std::endl;
			if (!Fi.isZero (CurrRow[i])){
				    // found pivot
				    // Updating column below pivot
				// Q [rank] = i;
				// pivotRows [row] = true;
				// P [rank] = row;
				fgemv(Fi, FFLAS::FflasNoTrans, M-row-1, rank, Fi.mOne, CurrRow+lda, lda, A+i, lda, Fi.one, CurrRow+lda+i, lda);
				    // Normalization
				Element invpiv;
				Fi.inv (invpiv, CurrRow[i]);
				if (Diag == FFLAS::FflasUnit)
					fscalin (Fi, N-i-1, invpiv, CurrRow+i+1,1);
				else
					fscalin (Fi, M-row-1, invpiv, CurrRow+i+lda,lda);
				if (i > rank){
					    // Column rotation to move pivot on the diagonal
					    // on U
					cyclic_shift_col(A+rank, rank, i-rank+1, lda);
					cyclic_shift_mathPerm(MathQ+rank, (size_t)(i-rank+1));
					    // on A
					cyclic_shift_col(CurrRow+lda+rank, M-row-1, i-rank+1, lda);
					Fi.assign(A[rank*(lda+1)], CurrRow[i]);
					fzero (Fi, i-rank, A+rank*(lda+1)+1, 1);
				}
				if (row > rank){
					    // Row rotation for L
					    // Optimization: delay this to the end
					cyclic_shift_row(A+rank*lda, row-rank+1, rank, lda);
					cyclic_shift_mathPerm(MathP+rank, (size_t) (row-rank+1) );
					    // Row rotation for U (not moving the 0 block)
					fcopy (Fi, N-i-1, A+rank*lda+i+1, 1, CurrRow+i+1, 1);
					Fi.assign(A[rank*(lda+1)], CurrRow[i]);
					fzero (Fi, row-rank, A+rank*(lda+1)+lda, lda);
					Fi.assign(CurrRow[i],Fi.zero); // only needed once here
				}
//				write_field(Fi,std::cerr<<"Apres cylicRow"<<std::endl,A,M,N,lda);
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
		for (size_t i=(size_t)rank; i<M; ++i)
			fzero (Fi, N-rank, A+i*lda+rank, 1);
//		write_field(Fi,std::cerr<<"Fini"<<std::endl,A,M,N,lda);


		return (size_t) rank;
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

	template<class Field>
	inline size_t
	PLUQ (const Field& Fi, const FFLAS::FFLAS_DIAG Diag,
	      const size_t M, const size_t N,
	      typename Field::Element * A, const size_t lda, size_t*P, size_t *Q)
	{
		typedef typename Field::Element Element;

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
				Element invpivot;
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
				Element invpivot;
				Fi.inv(invpivot, *A);
				for (size_t i=piv+1; i<M; ++i)
					Fi.mulin (*(A+i*lda), invpivot);
			}
			return 1;
		}
#endif
#ifdef BASECASE_K
		if (std::min(M,N) < BASECASE_K)
			return
#ifdef CROUT
				PLUQ_basecaseCrout
#else
				PLUQ_basecaseV2
#endif
				(Fi, Diag, M, N, A, lda, P, Q);
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
		Element * A2 = A + N2;
		Element * A3 = A + M2*lda;
		Element * A4 = A3 + N2;
		Element * F = A2 + R1*lda;
		Element * G = A3 + R1;
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
		Element * temp = new Element [R3*R2];
		for (size_t i=0; i<R3; ++i)
			fcopy (Fi, R2, temp + i*R2, 1, A4 + i*lda, 1);
		ftrsm (Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R3, R2, Fi.one, G, lda, temp, R2);
		    // N <- L3^-1 H2
		ftrsm (Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R3, N-N2-R2, Fi.one, G, lda, A4+R2, lda);
		    // O <- N - J V2
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, R3, N-N2-R2, R2, Fi.mOne, temp, R2, F+R2, lda, Fi.one, A4+R2, lda);
		delete[] temp;
		    // R <- H4 - K V2 - M3 O
		Element * R = A4 + R2 + R3*lda;
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
