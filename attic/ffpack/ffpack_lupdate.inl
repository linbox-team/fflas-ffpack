/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* ffpack/ffpack_ludivine.inl
 * Copyright (C) 2005 Brice Boyer
 *
 * Written by Brice Boyer
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
	/*!
	 * @brief Updates an existing LU factorisation with more rows.
	 *
	 * @param F Field on which arithmetic is done
	 * @param Diag Is \p L unit ? (if so, \c FFLAS::FflasUnit)
	 * @param trans Not used yet, should be \c FFLAS::FflasNoTrans
	 * @param M rows in \p A
	 * @param N cols in \p A
	 * @param A \p A is already in \c LU factorisation
	 * @param lda leading dimension of \p A
	 * @param R rank of \p A
	 * @param K rows in \p B
	 * @param B more rows to append to \p A
	 * @param ldb leading dimension of \p B (not tested if != lda)
	 * @param P permutation for \c LU in \p A. Should be big enough so it can store the permutation for \c LU of \p A and \p B
	 * @param Q same as \p P
	 * @param LuTag see \c LUdivine.
	 * @param cutoff see \c LUdivine.
	 *
	 * @return rank of <code>A.append(B)</code>
	 * @bug may be bogus.
	 */
	template <class Field>
	size_t LUpdate (const Field& F,
			const FFLAS::FFLAS_DIAG Diag, const FFLAS::FFLAS_TRANSPOSE trans,
			const size_t M, const size_t N,
			typename Field::Element * A, const size_t lda,
			const size_t R,
			const size_t K,
			typename Field::Element * B, const size_t ldb,
			size_t*P, size_t *Q , const FFPACK::FFPACK_LUDIVINE_TAG LuTag , const size_t cutoff)
	{
		if (trans == FFLAS::FflasTrans)
			throw Failure(__func__,__FILE__,__LINE__,"Transposed version is not implemented yet");
		// size_t MN = std::min(M,N);

		size_t incRow, incCol, rowDim, colDim;
#if 0 /*  not working */
		if (trans == FFLAS::FflasTrans){
			incRow = 1;
			incCol = lda;
			colDim = M;
			rowDim = N;
		}
		else
#endif
		{ // trans == FFLAS::FflasNoTrans
			incRow = lda;
			incCol = 1;
			colDim = N;
			rowDim = M;
		}

		size_t Nup   = rowDim;
		size_t Ndown = K;


		if ( !K ) {  // no line to append
			std::cout << "no row to append" << std::endl;
			return R;
		}
		if ( !R ) { // A was null
			// std::cout << "A was 0" << std::endl;
			// LU on B
			size_t R2 = LUdivine(F,Diag,trans,K,N,B,lda,P,Q+Nup,LuTag,cutoff);
			if (!R2)
				return 0 ;
			// Non zero row permutations
			for (size_t i = Nup; i < Nup + R2; ++i)
				Q[i] += Nup;
			{ // Permutation of the 0 rows Could probably be improved !
				if (Diag == FFLAS::FflasNonUnit){
					for ( size_t i = 0, j = R ; i < R2; ++i, ++j){
						FFLAS::fcopy( F, colDim - j, B + i*incRow + j*incCol, incCol, A + j * (lda + 1),
							      incCol);
						typename Field::Element *Ai = B + i*incRow + j*incCol ;
						typename Field::Element *Aend = B + colDim*incCol ;
						for (; Ai != Aend + i*incRow ; Ai+=incCol)
							F.assign (*Ai, F.zero);
						///@todo std::swap ?
						size_t t = Q[j];
						Q[j]=Q[Nup+i];
						Q[Nup+i] = t;
					}
				}
				else { // Diag == FFLAS::FflasUnit
					for ( size_t i = 0, ii = R+1 ; i < R2; ++i, ++ii){
						if (ii < M)
							FFLAS::fcopy( F, colDim - ii,
								      B + i*incRow + ii*incCol, incCol,
								      A + (ii-1)*incRow + ii*incCol, incCol);
						else {
							std::cout << "dangerous zone" << std::endl;
							FFLAS::fcopy( F, colDim - ii,
								      B + i*incRow + ii*incCol, incCol,
								      B + (ii-M-1)*incRow + ii*incCol, incCol);
						}

						typename Field::Element *Ai   = B + i*incRow + ii*incCol ;
						typename Field::Element *Aend = B + colDim*incCol ;
						for (; Ai != Aend + i*incRow ; Ai+=incCol)
							F.assign (*Ai, F.zero);
						size_t t = Q[ii-1];
						Q[ii-1]=Q[Nup+i];
						Q[Nup+i] = t;
					}
				}
			}

			return R2 ;
		}


		{
			// Recursive call on NW
			size_t R2;
			{ // trans == FFLAS::FflasNoTrans
				typename Field::Element *Ac = A  + R*incCol;    // NE
				typename Field::Element *Ar = B;                // SW
				typename Field::Element *An = Ar + R*incCol;    // SE


				// Ar <- Ar.P
				FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
						Ndown, 0,(int) R, Ar, lda, P);
				// Ar <- Ar.U1^-1
				ftrsm( F, FFLAS::FflasRight, FFLAS::FflasUpper,
				       FFLAS::FflasNoTrans, Diag, Ndown, R,
				       F.one, A, lda, Ar, lda);
				// An <- An - Ar*Ac
				fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Ndown, colDim-R, R,
				       F.mOne, Ar, lda, Ac, lda, F.one, An, lda);

				// LU call on SE
				R2=LUdivine (F, Diag, trans, Ndown, N-R, An, lda,P+R, Q+Nup,
					     LuTag, cutoff);
				if (R2) {
					for (size_t i = R; i < R + R2; ++i)
						P[i] += R;
					// An <- An.P2
					FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
							Nup,(int) R, (int)(R+R2), A, lda, P);

				}
				else {
					if (LuTag == FFPACK::FfpackSingular)
						return 0;
				}

			}
			// Non zero row permutations
			for (size_t i = Nup; i < Nup + R2; ++i)
				Q[i] += Nup;
			if (R < Nup){ // Permutation of the 0 rows
				if (Diag == FFLAS::FflasNonUnit){
					for ( size_t i = 0, j = R ; i < R2; ++i, ++j){
						FFLAS::fcopy( F, colDim - j, B + i*incRow + j*incCol, incCol, A + j * (lda + 1),
							      incCol);
						typename Field::Element *Ai = B + i*incRow + j*incCol ;
						typename Field::Element *Aend = B + colDim*incCol ;
						for (; Ai != Aend + i*incRow ; Ai+=incCol)
							F.assign (*Ai, F.zero);
						///@todo std::swap ?
						size_t t = Q[j];
						Q[j]=Q[Nup+i];
						Q[Nup+i] = t;
					}
				}
				else { // Diag == FFLAS::FflasUnit
					for ( size_t i = 0, ii = R+1 ; i < R2; ++i, ++ii){
						if (ii < M)
							FFLAS::fcopy( F, colDim - ii,
								      B + i*incRow + ii*incCol, incCol,
								      A + (ii-1)*incRow + ii*incCol, incCol);
						else {
							// std::cout << "dangerous zone" << std::endl;
							FFLAS::fcopy( F, colDim - ii,
								      B + i*incRow + ii*incCol, incCol,
								      B + (ii-M-1)*incRow + ii*incCol, incCol);
						}

						typename Field::Element *Ai   = B + i*incRow + ii*incCol ;
						typename Field::Element *Aend = B + colDim*incCol ;
						for (; Ai != Aend + i*incRow ; Ai+=incCol)
							F.assign (*Ai, F.zero);
						size_t t = Q[ii-1];
						Q[ii-1]=Q[Nup+i];
						Q[Nup+i] = t;
					}
				}
			}
			return R + R2;
		}
	}
