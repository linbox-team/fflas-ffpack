/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
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
	inline size_t fsytrf_UP_RPM_BC_RL (const Field& F, const size_t N,
									   typename Field::Element_ptr A, const size_t lda,
									   typename Field::Element_ptr Dinv, const size_t incDinv,
									   size_t * P){
		
	}

	template <class Field>
	inline size_t fsytrf_UP_RPM_BC_Crout (const Field& F, const size_t N,
										  typename Field::Element_ptr A, const size_t lda,
										  typename Field::Element_ptr Dinv, const size_t incDinv,
										  size_t * P){
		
		size_t * MathP = FFLAS::fflas_new<size_t>(N);
		size_t * twoBlocks = FFLAS::fflas_new<size_t>(N);
		for (size_t i=0; i<N; ++i){
			twoBlocks[i]=0;
			MathP[i] = i;
		}
	
		for (size_t row = 0, size_t rank = 0, typename Field::Element_ptr CurrRow=A;
			 row < N;
			 row++, CurrRow += lda+1){

				/* A =  [   U  | * | b | B ]  where U is rank x rank, b rank x 1
				 *      [------------------]
				 *      [      | 0 |   0   ]
				 *      [------------------]
				 *      [      | 0 |   c   ]
				 *      [      | 0 |   C   ]
				 */

				// tmp <-  D^-1 . b
			typename Field::Element_ptr tmp = FFLAS::fflas_new(F, rank);
			for (size_t j=0, typename Field::Element_ptr bj=A+row; j<rank; ++j, bj+=lda, Dinvj+=incDinv)
				F.mul (tmp[j], *Anj, *Dinvj);

				// c <- c - tmp^T . [ b | B ]
			fgemv (F, FFLAS::FflasTrans, rank, N-rank, Fi.mOne, A+row, lda, tmp, 1, F.one, CurrRow, 1);
			size_t i = 0;
			while (Fi.isZero (*(CurrRow + i++)) && i<N);
			i--;
			if (!Fi.isZero (CurrRow[i])){
				    // found pivot
				if (!i){ // pivot is on the leading diagonal -> 1x1 diagonal block
					F.inv (*Dinvj, CurrRow[i]);
					if (row>rank){ // some zero blocks are in the way
							// column rotations
							// On U
						cyclic_shift_col(F, A+rank, rank, row+1, lda);
							// On P
						cyclic_shift_mathPerm(MathP+rank, row+1);
						
							// Moving pivot row
						F.assign (A[rank*(lda+1)], CurrRow[i]);
						FFLAS::fzero (F, row-rank, A+rank*(lda+1)+1, 1);
						FFLAS::fassign (F, N-row-1, CurrRow+1, 1, A+rank*lda+row+1, 1);
							//Fi.assign(*(CurrRow+i),Fi.zero); // CP : not sure we need it
					}
					rank++;
				} else { // off diagonal pivot -> forming a  2x2 diagonal block
						
						/* Changing 
						 * A =  [   U  |      B         ]  where U is rank x rank
						 *      [-----------------------]
						 *rank->[      | 0 |   0        ]
						 *      [-----------------------]
						 * row->[      | 0 |   0 0 x  C ]
						 *      [      | 0 |   0 D ET F ]
						 *      [      | 0 |   x E y  G ]
						 *      [      | 0 |   * * *  H ]
						 * into 
						/* A =  [   U  |        B'      ]  where U is rank x rank, b rank x 1
						 *      [-----------------------]
						 *      [      | x y/2| 0  E' G']
						 *      [      | *  x | 0  0  C']
						 *      [-----------------------]
						 *      [             | 0  0  0 ]
						 *      [             | 0  D' F']
						 *      [             | 0  *  H']
						 */
						// A[rank,rank] <- x
					F.assign (A[rank*(lda+1)], CurrRow[i]);
						// A[rank,rank+1] <- y/2
					F.div (A[rank*(lda+1)+1], A[(row+i)*(lda+1)], two);
						// A[rank+1,rank+1] <- x
					F.assign (A[(rank+1)*(lda+1)], A[rank*(lda+1)]);
						// 0 E' G' <- 0 E G
					FFLAS::fzero (F, row-rank, A+rank*(lda+1)+2, 1);
					FFLAS::fassign(F, i-1, CurrRow+i*lda+1,1, A+rank*lda+row+2, 1);
					FFLAS::fassign(F, N-row-i-1, CurrRow+i*(lda+1)+1,1, A+rank*lda+row+i+1, 1);
						// 0 0 C' <- 0 0 C
					FFLAS::fzero (F, row-rank+i-1, A+(rank+1)*(lda+1)+1, 1);
					FFLAS::fassign (F, N-row-i-1, CurrRow+i+1,1, A+(rank+1)*lda+row+i+1,1);
						// Moving D 1 row down and 1 col right
					FFLAS::fassign (F, i-1, i-1, CurrRow+lda+1,lda,CurrRow+2*(lda+1),lda);
						// Moving F 1 row down
					FFLAS::fassign (F, i-1, N-row-i-1, CurrRow+lda+i+1,lda,CurrRow+lda+i+1,lda,CurrRow+2*lda+i+1,lda);
						// G'<- G' -y/2.x^-1.C'
					typename Field::Element x;
					F.init(x);
					F.div (x, A[rank*(lda+1)+1], A[rank*(lda+1)]);
					F.negin(x);
					FFLAS::faxpy (F, N-row-i-1, x, A+(rank+1)*lda+row+i+1, 1, A+rank*lda+row+i+1, 1);
						// Update permutation
					size_t Prow = MathP[row];
					size_t Pi = MathP[row+i];
					std::copy(MathP+row+1,MathP+row+i,Math+row+2);
					std::copy(MathP+rank,MathP+row,MathP+rank+2);
					MathP[rank]=Prow;
					MathP[rank+1]=Pi;

					twoBlocks[rank]=row;
					twoBlocks[rank+1]=row+i;
					row++;
					CurrRow+=lda+1;
					rank+=2;
				}
			} // if no pivot found then keep going
		}

		MathPerm2LAPACKPerm (P, MathP, N);
		for (size_t i=0; i<rank; i++){
			if (twoBlocks[i+1]>0){
				P[i] = -twoBlocks[i];
				P[i+1] = -twoBlocks[i+1];
				i++;
			}
		}
		FFLAS::fflas_delete (MathP, twoBlocks);
		FFLAS::fzero (Fi, N-rank, N-rank, A+rank*(1+lda), lda);
		return (size_t) rank;
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
