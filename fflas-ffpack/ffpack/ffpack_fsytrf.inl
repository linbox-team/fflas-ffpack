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
            // TODO: maybe one day
        return 0;
	}

	template <class Field>
	inline size_t fsytrf_LOW_RPM_BC_Crout (const Field& F, const size_t N,
										  typename Field::Element_ptr A, const size_t lda,
										  typename Field::Element_ptr Dinv, const size_t incDinv,
										  size_t * P){
            //  TODO
        return 0;
    }
	template <class Field>
	inline size_t fsytrf_UP_RPM_BC_Crout (const Field& F, const size_t N,
										  typename Field::Element_ptr A, const size_t lda,
										  typename Field::Element_ptr Dinv, const size_t incDinv,
										  size_t * P){
        // std::cerr<<"*****************************************************"<<std::endl;
        // std::cerr<<"*****************************************************"<<std::endl;
        FFLAS::WriteMatrix(std::cerr<<"Base Case A = "<<std::endl,F,N,N,A,lda);	
		size_t * MathP = FFLAS::fflas_new<size_t>(N);
		int * twoBlocks = FFLAS::fflas_new<int>(N);
		for (size_t i=0; i<N; ++i){
			twoBlocks[i]=-1;
			MathP[i] = i;
		}
	
		size_t rank = 0;
        typename Field::Element two;
        F.init(two, 2UL);
        typename Field::Element_ptr CurrRow = A;
        for (size_t row = 0; row < N; row++, CurrRow += (lda+1)){

            // std::cerr<<"======================================================"<<std::endl;
            // FFLAS::WriteMatrix(std::cerr<<"A = "<<std::endl,F,N,N,A,lda);
				/* A =  [   U  | * | b | B ]  where U is rank x rank, b rank x 1
				 *      [------------------]
				 *      [      | 0 |   0   ]
				 *      [------------------]
				 *      [      | 0 |   c   ]
				 *      [      | 0 |   C   ]
				 */

				// tmp <-  D^-1 . b
			typename Field::Element_ptr tmp = FFLAS::fflas_new(F, rank);
            typename Field::Element_ptr bj = A+row, Dinvj=Dinv;

                //std::cerr<<"twoBlocks = ";
                //for (size_t h=0; h<rank; h++)
                    //std::cerr<<twoBlocks[h]<<" ";
                    //         std::cerr<<std::endl;
            for (size_t j=0; j<rank; ++j, bj+=lda, Dinvj+=incDinv){
                    //std::cerr<<"j = "<<j<<std::endl;
                if (twoBlocks[j]<0)
                    F.mul (tmp[j], *bj, *Dinvj);
                else{
                    F.mul (tmp[j+1], *bj, *Dinvj);
                    F.mul (tmp[j], *(bj+lda), *(Dinvj+1));
                    Dinvj+=incDinv;
                    bj+=lda;
                    j++;
                }
            }

           // c <- c - tmp^T . [ b | B ]
            // FFLAS::WriteMatrix(std::cerr<<"avant fgemv A = "<<std::endl,F,rank,N-row,A+row,lda);
            // FFLAS::WriteMatrix(std::cerr<<"avant fgemv tmp = "<<std::endl,F,1,rank,tmp,rank);
            // FFLAS::WriteMatrix(std::cerr<<"avant fgemv CurrRow = "<<std::endl,F,1,N-row,CurrRow,N-row);
            fgemv (F, FFLAS::FflasTrans, rank, N-row, F.mOne, A+row, lda, tmp, 1, F.one, CurrRow, 1);
            FFLAS::fflas_delete(tmp);
            // FFLAS::WriteMatrix(std::cerr<<"apres fgemv A = "<<std::endl,F,N,N,A,lda);
            size_t i = 0;
            while (F.isZero (*(CurrRow + i++)) && i<N-row);
			i--;
			if (!F.isZero (CurrRow[i])){
                // std::cerr<<"Found pivot "<< CurrRow[i]<<" at row = "<<row<<" col "<< i<<std::endl;
				    // found pivot
				if (!i){ // pivot is on the leading diagonal -> 1x1 diagonal block
                        //std::cerr<<"on the main diagonal"<<std::endl;
					F.inv (*Dinvj, CurrRow[i]);
					if (row>rank){ // some zero blocks are in the way
							// column rotations
							// On U
						cyclic_shift_col(F, A+rank, rank, row-rank+1, lda);
							// On P
						cyclic_shift_mathPerm(MathP+rank, row-rank+1);
                        // FFLAS::WritePermutation(std::cerr<<"MathP = ",MathP,N);
							// Moving pivot row
						F.assign (A[rank*(lda+1)], *CurrRow);
						FFLAS::fzero (F, row-rank, A+rank*(lda+1)+1, 1);
						FFLAS::fassign (F, N-row-1, CurrRow+1, 1, A+rank*lda+row+1, 1);
                        FFLAS::fzero (F, N-row, CurrRow,1); // CP : not sure we need it
					}
					rank++;
				} else { // off diagonal pivot -> forming a  2x2 diagonal block
					
                        //std::cerr<<"off diagonal pivot"<<std::endl;
                    
                        // Crout update of the row (row+i)
                    tmp = FFLAS::fflas_new(F, rank);
                    bj = A+row+i, Dinvj=Dinv;
                    for (size_t j=0; j<rank; ++j, bj+=lda, Dinvj+=incDinv){
                            //std::cerr<<"j = "<<j<<std::endl;
                        if (twoBlocks[j]<0)
                            F.mul (tmp[j], *bj, *Dinvj);
                        else{
                            F.mul (tmp[j+1], *bj, *Dinvj);
                            F.mul (tmp[j], *(bj+lda), *(Dinvj+1));
                            Dinvj+=incDinv;
                            bj+=lda;
                            j++;
                        }
                    } 
                        // c <- c - tmp^T . [ B ]
                    // FFLAS::WriteMatrix(std::cerr<<"avant fgemv A = "<<std::endl,F,rank,N-row-1,A+row+1,lda);
                    // FFLAS::WriteMatrix(std::cerr<<"avant fgemv tmp = "<<std::endl,F,1,rank,tmp,rank);
                    // FFLAS::WriteMatrix(std::cerr<<"avant fgemv CurrRow = "<<std::endl,F,1,N-row-1,CurrRow+i*lda+1,N-row );
                    FFLAS::fassign(F,i-1,CurrRow+lda+i,lda,CurrRow+i*lda+1,1);
                    fgemv (F, FFLAS::FflasTrans, rank, N-row-1, F.mOne, A+row+1, lda, tmp, 1, F.one, CurrRow+i*lda+1, 1);
                    FFLAS::fflas_delete(tmp);
                    // FFLAS::WriteMatrix(std::cerr<<"apres 2eme fgemv A = "<<std::endl,F,N,N,A,lda);
       
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
                         * A =  [   U  |        B'      ]  
						 *      [-----------------------]
						 *      [      | x y/2| 0  E' G']
						 *      [      | *  x | 0  0  C']
						 *      [-----------------------]
						 *      [             | 0  0  0 ]
						 *      [             | 0  D' F']
						 *      [             | 0  *  H']
						 */
                    size_t Ddim = i-1;
                    size_t Hdim = N-row-Ddim-2;
						// A[rank,rank] <- x
					F.assign (A[rank*(lda+1)], CurrRow[i]);
                    F.inv (*Dinvj, CurrRow[i]);
                        // A[rank,rank+1] <- y/2
					F.div (A[rank*(lda+1)+1], A[(row+i)*(lda+1)], two);

                         // 0 E' G' <- 0 E G
                    if (row==rank){ // Then we need to save C in a buffer
                            //std::cerr<<"Ddim = "<<Ddim<<" Hdim = "<<Hdim<<std::endl;
                        tmp = FFLAS::fflas_new(F, Hdim);
                        FFLAS::fassign (F, Hdim, CurrRow+i+1,1, tmp, 1);
                        // FFLAS::WriteMatrix(std::cerr<<"tmp = "<<std::endl,F,1,Hdim,tmp,Hdim);
                    } else{
                        tmp = CurrRow+i+1;
                        FFLAS::fzero(F, row-rank, A+rank*(lda+1)+2, 1);
                    }
                    FFLAS::fassign (F, Ddim, CurrRow+i*lda+1, 1, A+rank*lda+row+2,1);
                    FFLAS::fassign (F, Hdim, CurrRow+i*lda+2+Ddim, 1, A+rank*lda+N-Hdim,1);
                    
						// Moving D 1 row down and 1 col right
                    // FFLAS::WriteMatrix(std::cerr<<"avant moving D = "<<std::endl,F,N,N,A,lda);
                    for (size_t l=Ddim; l>0; l--)
                        FFLAS::fassign (F, 1, Ddim, CurrRow+l*lda+1,lda,CurrRow+(l+1)*lda+2,lda);
                        //FFLAS::WriteMatrix(std::cerr<<"apres moving D = "<<std::endl,F,N,N,A,lda);

						// Moving F 1 row down
                        //FFLAS::fassign (F, Ddim, Hdim, CurrRow+lda+i+1,lda,CurrRow+2*lda+i+1,lda);
					for (size_t l=Ddim; l>0; l--)
                        FFLAS::fassign (F, 1, Hdim, CurrRow+l*lda+i+1,lda,CurrRow+(l+1)*lda+i+1,lda);
                        //FFLAS::WriteMatrix(std::cerr<<"apres moving F = "<<std::endl,F,N,N,A,lda);

						// A[rank+1,rank+1] <- x
					F.assign (A[(rank+1)*(lda+1)], A[rank*(lda+1)]);
                    F.assign(*(Dinvj+1),*Dinvj);
                        // 0 0 C' <- 0 0 C
                        // TODO could be only 2 assigns to 0
					FFLAS::fzero (F, row-rank+Ddim, A+(rank+1)*(lda+1)+1, 1);
                        //std::cerr<<"row = "<<row<<" rank+1 = "<<rank+1<<std::endl;
                    if (row!=rank+1){
                        FFLAS::fassign (F, Hdim, tmp,1, A+(rank+1)*lda+row+i+1,1);
                    }
                        //FFLAS::WriteMatrix(std::cerr<<"apres moving C = "<<std::endl,F,N,N,A,lda);
                    if (row==rank)
                        FFLAS::fflas_delete(tmp);

                        // first row of [ D ET ] <-0
                        //std::cerr<<"row = "<<row<<" rank = "<<rank<<CurrRow[lda+1]<<std::endl;
                    if (row > rank)
                        FFLAS::fzero(F, Ddim+1, CurrRow+lda+1, 1);
                        // TODO: possibly need to zero out the row where x C used to be

                        // G'<- G' -y/2.x^-1.C'
                    FFLAS::WriteMatrix(std::cerr<<"avant update G' "<<std::endl,F,N,N,A,lda);
                    typename Field::Element x;
					F.init(x);
					F.mul (x, A[rank*(lda+1)+1], *Dinvj/*A[rank*(lda+1)]*/);
					F.negin(x);
					FFLAS::faxpy (F, Hdim, x, A+(rank+1)*lda+row+i+1, 1, A+rank*lda+N-Hdim, 1);

                    FFLAS::WriteMatrix(std::cerr<<"apres update G' "<<std::endl,F,N,N,A,lda);

                        // Rotate the columns of the upper part of U
                        //std::cerr<<"before cyclic shift row,i,rank="<<row<<" "<<i<<" "<<rank<<std::endl;
                    if (row==rank) // only need to rotate 2nd pivot
                        cyclic_shift_col(F, A+rank+1, rank, row+i-rank, lda);
                    else{ // need to rotate 2 pivots
                            //TODO: to both rotations at once
                        cyclic_shift_col(F, A+rank, rank, row-rank+1, lda);
                            //FFLAS::WriteMatrix(std::cerr<<"apres 1st cyclic' "<<std::endl,F,N,N,A,lda);
                        cyclic_shift_col(F, A+rank+1, rank, row+i-rank, lda);
                            //FFLAS::WriteMatrix(std::cerr<<"apres 2nd cyclic' "<<std::endl,F,N,N,A,lda);
                    }
                        // Update permutation
					size_t Prow = MathP[row];
					size_t Pi = MathP[row+i];
					std::copy(MathP+row+1,MathP+row+i,MathP+row+2);
					std::copy(MathP+rank,MathP+row,MathP+rank+2);
					MathP[rank]=Prow;
					MathP[rank+1]=Pi;

					twoBlocks[rank]=row;
					twoBlocks[rank+1]=row+i;
					row++;
					CurrRow+=(lda+1);

					rank+=2;
				}
			} // if no pivot found then keep going
		}
         FFLAS::WriteMatrix(std::cerr<<"done A = "<<std::endl,F,N,N,A,lda);
        // FFLAS::WritePermutation(std::cerr<<"MathPerm =",MathP,N);

		MathPerm2LAPACKPerm (P, MathP, N);
        // FFLAS::WritePermutation(std::cerr<<"LAPACKPerm =",MathP,N);
		for (size_t i=0; i<rank; i++){
			if (twoBlocks[i]>=0){
				P[i] = -P[i]-1;
				P[i+1] = -P[i+1]-1;
				// P[i] = -twoBlocks[i+1];
				// P[i+1] = -twoBlocks[i];
				i++;
			}
		}
        FFLAS::WritePermutation(std::cerr<<"LAPACKPerm =",P,N);
		FFLAS::fflas_delete (MathP, twoBlocks);
        // std::cerr<<"N = "<<N<<" rank = "<<rank<<std::endl;
		FFLAS::fzero (F, N-rank, N-rank, A+rank*(1+lda), lda);
		return (size_t) rank;
	}
    
	template <class Field>
    inline size_t fsytrf_UP_RPM (const Field& Fi, const size_t N,
                                 typename Field::Element_ptr A, const size_t lda,
                                 typename Field::Element_ptr Dinv, const size_t incDinv,
                                 size_t * P, size_t BCThreshold){
        if (N < BCThreshold){
            return fsytrf_UP_RPM_BC_Crout (Fi,N,A,lda,Dinv,incDinv,P);
		}
        std::vector<bool> twoBlocks(N,false);
        
        size_t N1 = N>>1;
        size_t N2 = N-N1;
		size_t * P1 = FFLAS::fflas_new<size_t >(N1);
		size_t R1,R2,R3;
        
            // A1 = P1^T [ U1^T ] D1 [ U1 V1 ] P1
		    //           [ V1^T ]
        R1 = fsytrf_UP_RPM (Fi, N1, A, lda, Dinv, incDinv, P1, BCThreshold);
        for (size_t i=0; i<R1; i++)
            if (int(P1[i])<0){
                P1[i]=-P1[i]-1;
                twoBlocks[i] = i? (!twoBlocks[i-1]) : true; // mark the first of the 2 positions
            }
        std::cerr<<"After R1 twoblocks =  "<<twoBlocks<<std::endl;
        FFLAS::WriteMatrix(std::cerr<<"A = "<<std::endl,Fi,N,N,A,lda);
        
        typename Field::Element_ptr A2 = A + N1;
        typename Field::Element_ptr A4 = A2 + N1*lda;
            // [ B1 ] <- P1^T A2
		    // [ B2 ]
		applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N2, 0, N1, A2, lda, P1);
        FFLAS::WriteMatrix(std::cerr<<"after P1^T A = "<<std::endl,Fi,N,N,A,lda);
  
        
        typename Field::Element_ptr B1 = A2;
        typename Field::Element_ptr B2 = B1 + R1*lda;
            /*     [ U1 V1 | B1 ]
             *     [    0  | B2 ]
             *     [ ------|--- ]
             *     [    0  | A4 ]
             *
             */
		    // C <- U1^-T B1
		ftrsm (Fi, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasTrans, FFLAS::FflasNonUnit, R1, N2, Fi.one, A, lda, B1, lda);
        FFLAS::WriteMatrix(std::cerr<<"after ftrsm A = "<<std::endl,Fi,N,N,A,lda);
            // F <- B2 - V1^T C
		fgemm (Fi, FFLAS::FflasTrans, FFLAS::FflasNoTrans, N1-R1, N2, R1, Fi.mOne, A + R1, lda, B1, lda, Fi.one, B2, lda);
        FFLAS::WriteMatrix(std::cerr<<"after fgemm A = "<<std::endl,Fi,N,N,A,lda);
            // G <- A4 - C^T D1^-1 C
            // E <- D1^-1 C (done simultaneously)
        fsyrk (Fi, FFLAS::FflasUpper, FFLAS::FflasTrans, N2, R1, Fi.mOne, B1, lda, A, lda+1, twoBlocks, Fi.one, A4, lda);
        FFLAS::WriteMatrix(std::cerr<<"after fsyrk A = "<<std::endl,Fi,N,N,A,lda);

            /*     [ U1 V1 | E ]
             *     [    0  | F ]
             *     [ ------|-- ]
             *     [    0  | G ]
             */
            // F = P2 [ L2 ] [ U2 V2 ] Q2
		    //        [ M2 ]
		size_t * P2 = FFLAS::fflas_new<size_t >(N1-R1);
		size_t * Q2 = FFLAS::fflas_new<size_t >(N2);
        typename Field::Element_ptr F=A2+R1*lda;
        FFLAS::WriteMatrix(std::cerr<<"F = "<<std::endl,Fi,N1-R1,N2,F,lda);
		R2 = _PLUQ (Fi, FFLAS::FflasUnit, N1-R1, N2, F, lda, P2, Q2, BCThreshold);
        FFLAS::WriteMatrix(std::cerr<<"After PLUQ R2 = "<<R2<<" A = "<<std::endl,Fi,N,N,A,lda);

        typename Field::Element_ptr H1 = A4, H2 = H1+R2, H3 = H2+R2*lda;

        if (R2){
                // [ G1   G2 ] <- Q2 G Q2^T
                // [ G2^T G3 ]
            applyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, N2, size_t(0), N2, A4, lda, Q2);
                //FFLAS::WriteMatrix(std::cerr<<"after G1 G2 A = "<<std::endl,Fi,N,N,A,lda);
            applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N2, size_t(0), N2, A4, lda, Q2);
                //FFLAS::WriteMatrix(std::cerr<<"after G2 G3 A = "<<std::endl,Fi,N,N,A,lda);
                // [ E1 E2 ] <- E Q2^T
            applyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, R1, size_t(0), N2, A2, lda, Q2);
                // [ V11 V12 ] <- V1 P2^T
            applyP (Fi, FFLAS::FflasRight, FFLAS::FflasNoTrans, R1, size_t(0), N1-R1, A+R1, lda, P2);
                // H1 <- upper tri such that U2^TxH1 + H1^T x U2= G1
            ftrssyr2k (Fi, FFLAS::FflasUpper, FFLAS::FflasUnit, R2, F, lda, A4, lda);
            FFLAS::WriteMatrix(std::cerr<<"after ftrssyr2k A = "<<std::endl,Fi,N,N,A,lda);
                // H2 <-  U2^-T (G2 - H1^T V2)
            FFLAS::ftrmm (Fi, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasTrans, FFLAS::FflasNonUnit, R2, N2-R2, Fi.mOne, A4, lda, F+R2, lda, Fi.one, H2, lda);
                //FFLAS::WriteMatrix(std::cerr<<"after ftrmm A = "<<std::endl,Fi,N,N,A,lda);
            FFLAS::ftrsm (Fi, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasTrans, FFLAS::FflasUnit, R2, N2-R2, Fi.one, F, lda, H2, lda);
                //FFLAS::WriteMatrix(std::cerr<<"after ftrsm A = "<<std::endl,Fi,N,N,A,lda);
                // H3 <- G3 - (V2^T H2 + H2^T V2)
            fsyr2k (Fi, FFLAS::FflasUpper, FFLAS::FflasTrans, N2-R2, R2, Fi.mOne, F+R2, lda, H2, lda, Fi.one, H3, lda);
            FFLAS::WriteMatrix(std::cerr<<"after fsyr2k A = "<<std::endl,Fi,N,N,A,lda);
                // U2' V2'  <-  D2 *  [U2 V2]
            typename Field::Element_ptr D2i = F, Dinvi = Dinv+R1*incDinv, H2i=H2;
            for (size_t i=0; i<R2; i++, Dinvi+=2, D2i+=lda+1, H2i+=lda){
                std::cerr<<"pivot = "<<D2i[0]<<std::endl;
                Fi.inv (*(Dinvi), *D2i);
                    //Fi.assign(*(Dinvi+1), *Dinvi);
                Fi.assign(*(Dinvi+1), *Dinvi);
                FFLAS::fscalin (Fi, N2-i-1, *D2i, D2i+1, 1);
                std::cerr<<"scaling H2 by "<<*D2i<<std::endl;
                FFLAS::WriteMatrix(std::cerr<<"H2 = "<<std::endl,Fi,1,N2-R2,H2i,lda);
                    //FFLAS::fscalin (Fi, N2-R2, *D2i, H2i, 1);
            }
        }
            /*     [ U1 V1 | E1      E2 ]
             *     [    0  | L2 \ U2 V2 ]
             *     [    0  | M2     0   ]
             *     [ ------|----------- ]
             *     [    0  | H1     H2  ]
             *     [    0  |        H3  ]
             */
        size_t * P3 = FFLAS::fflas_new<size_t >(N2-R2);
            // H3 = P3^T [ U3^T ] D3 [ U3 V3 ] P3
            //           [ V3^T ]
        FFLAS::WriteMatrix(std::cerr<<"Before R3 = "<<R3<<" A = "<<std::endl,Fi,N,N,A,lda);

        R3 = fsytrf_UP_RPM (Fi, N2-R2, H3, lda, Dinv+R1+2*R2, incDinv, P3, BCThreshold);

        FFLAS::WriteMatrix(std::cerr<<"After R3 = "<<R3<<" A = "<<std::endl,Fi,N,N,A,lda);
        for (size_t j=R1; j<R1+2*R2; j+=2)
            twoBlocks[j] = true;
        for (size_t i=0,j=R1+2*R2; i<R3; i++,j++)
            if (int(P3[i])<0){
                P3[i]=-P3[i]-1;
                twoBlocks[j] = j? (!twoBlocks[j-1]) : true; // mark the first of the 2 positions
            }
        
            // [ E21 E22 ]     [ E2 ]
            // [ V21 V22 ]  <- [ V2 ] P3^T
		    // [  0   0  ]     [  0 ]
		    // [ H21 H22 ]     [ H2 ]
        std::cerr<<"twoblocks =  "<<twoBlocks<<std::endl;
        FFLAS::WritePermutation(std::cerr<<"before applyP P3 =  ", P3, N2-R2);
		applyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, R1+R2, 0, N2-R2, A2+R2, lda, P3);
        FFLAS::WriteMatrix(std::cerr<<"After apply P3 = "<<R3<<" A = "<<std::endl,Fi,N,N,A,lda);
		applyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, R2,    0, N2-R2, H2, lda, P3);
        FFLAS::WriteMatrix(std::cerr<<"After apply P3 = "<<R3<<" A = "<<std::endl,Fi,N,N,A,lda);

		    // P <- Diag (P1 [ I_R1    ] , P3 [ I_R3    ])
		    //               [      P2 ]      [      P4 ]
		size_t* MathP = FFLAS::fflas_new<size_t>(N);
        for (size_t i=0; i<N; ++i)
            MathP[i] = i;
            /** MathP  <-[ [ I      ] x P1 |      ]   [ I_(N1+R2)      ]
             *             [   P2^T ]      |      ] x [           P3^T ]
             *             [ --------------|----- ]
             *             [               | Q2^T ]
             */
		composePermutationsLLM (MathP, P1, P2, R1, N1); // to be checked
        FFLAS::WritePermutation(std::cerr<<"composePermLLM MathP = ", MathP, N);
        LAPACKPerm2MathPerm (MathP+N1, Q2, N2);
        FFLAS::WritePermutation(std::cerr<<"Q2 = ", Q2, N2);
        FFLAS::WritePermutation(std::cerr<<"LAPACK2Math Q2 -> MathP = ", MathP, N);
        composePermutationsMLM (MathP, P3, N1+R2, N);
        FFLAS::WritePermutation(std::cerr<<"composePermMLM = ", MathP, N);

		FFLAS::fflas_delete( P1);
		FFLAS::fflas_delete( P2);
		FFLAS::fflas_delete( P3);
		for (size_t i=N1; i<N; ++i)
			MathP[i] += N1;

            /** Changing [ U1 V1 | E1      E21 E22 ] into [ U1 E11 E12  V1   E*   E*  ]
             *           [    0  | L2 \ U2 V21 V22 ]      [    U4  V41   0  V42  V43  ]
             *           [    0  | M2       0   0  ]      [         U3   0   0    V3  ]
             *           [ ------|---------------- ]      [              0   0     0  ]
             *           [    0  | H1      H21 H22 ]
             *           [    0  |          U3  V3 ]
             *           [    0  |               0 ]
             * where U4 is the 2R2 x 2R2 matrix formed by interleaving U2, L2^T and H1
             */
                // Rotating [V1 E1 E21]
        size_t dim = N1-R1+R2+R3;
        typename Field::Element_ptr tmpV1 = FFLAS::fflas_new(Fi, R1, dim);
        FFLAS::fassign(Fi, R1, dim, A+R1, lda, tmpV1, dim);
        FFLAS::WriteMatrix(std::cerr<<"tmpV1 =  "<<std::endl,Fi,R1,dim,tmpV1,dim);

        typename Field::Element_ptr Di=A+R1, tmpV1i=tmpV1;
        for (size_t i=0; i<R1; i++, tmpV1i+=dim, Di+=lda){
            FFLAS::fassign(Fi, R2, tmpV1i, 1, Di, 2);
            FFLAS::fassign(Fi, R2, tmpV1i+N1-R1, 1, Di+1, 2);
        }
        FFLAS::WriteMatrix(std::cerr<<"After for rotating [V1 E1 E21] A = "<<std::endl,Fi,N,N,A,lda);
        FFLAS::fassign(Fi, R1, R3, tmpV1+N1-R1+R2, dim, A+R1+2*R2, lda);
        FFLAS::fassign(Fi, R1, N1-R1-R2, tmpV1+R2, dim, A+R1+2*R2+R3, lda);
        FFLAS::WriteMatrix(std::cerr<<"After rotating [V1 E1 E21] A = "<<std::endl,Fi,N,N,A,lda);
        if (R2){

            
                // Interleaving L2^T, U2 and H1 into 
            typename Field::Element_ptr tmp = FFLAS::fflas_new(Fi, N1-R1, R2);
            FFLAS::fassign(Fi, N1-R1, R2, F, lda, tmp, R2);
            Di = A+R1*(lda+1);
            typename Field::Element_ptr tmpi = tmp;
            typename Field::Element_ptr H1i = H1;
            for (size_t i=0; i<R2; i++, Di+=2*lda+2, tmpi+=R2+1, H1i+=lda+1){
                   Fi.assign(*Di,*tmpi); // pivot
                   FFLAS::fassign (Fi, R2-i-1, tmpi+R2, R2, Di+2, 2); // copying L2
                   FFLAS::fassign (Fi, R2-i, H1i, 1, Di+1, 2); // copying H1
                   FFLAS::fassign (Fi, R2-i, tmpi, 1, Di+lda+1, 2); // copying U2
                   FFLAS::fzero (Fi, R2-i-1, Di+lda+2, 2); 
            }
// size_t hR2 =  (N1-R1)>>1; // Do a first pass with only N-R1 pivots (to avoid overwritting L2\U2)
            //     //if (hR2 %2) hR2--;
            // typename Field::Element_ptr Di = A+R1*(lda+1);
            // typename Field::Element_ptr Fptr = F;
            // typename Field::Element_ptr H1i = H1;
            // size_t i=0, dim=hR2;
            // std::cerr<<"dim = "<<dim<<std::endl;
            // FFLAS::WriteMatrix(std::cerr<<"-1/Before interleaving L2 U2 H1 A = "<<std::endl,Fi,N,N,A,lda);
            // for (; i<R2 && int(dim)> 0; i++, Di+=2*lda+2, Fptr+=lda+1, H1i+=lda, dim-=2){
            //     Fi.assign(*Di,*Fptr); // pivot
            //     FFLAS::fassign (Fi, dim-1, Fptr+lda, lda, Di+2, 2); // copying L2
            //     FFLAS::fassign (Fi, dim, H1i, 1, Di+1, 2); // copying H1
            //     FFLAS::fassign (Fi, dim, Fptr, 1, Di+lda+1, 2); // copying U2
            // }
            // FFLAS::WriteMatrix(std::cerr<<"0/ After interleaving L2 U2 H1 dim = "<<dim<<" A = "<<std::endl,Fi,N,N,A,lda);
            // typename Field::Element_ptr tmpM2 = FFLAS::fflas_new(Fi, R2, N1-R1-R2);
            //     // saving M2
            // for (size_t j=0; j<N1-R1-R2; j++)
            //     FFLAS::fassign(Fi, R2, F+R2*lda+j*lda, 1, tmpM2+j, N1-R1-R2);
            //     // finishing the interleave of L2^T U2 and H1:
            //     // 1/ copying the remaining of U2[1:i,i:R2]
            // std::cerr<<"i = "<<i<<" R2 = "<<R2<<std::endl;
            // FFLAS::fassign (Fi, i, R2-i, F+i, lda, A+(R1+1)*(lda+1)+2*i, 2*lda);
            // FFLAS::WriteMatrix(std::cerr<<"1/ After interleaving L2 U2 H1 A = "<<std::endl,Fi,N,N,A,lda);
            //     // 2/ copying the remaining of L2[i:R2,1:i]
            // for(size_t j=0; j<i; j++)
            //     FFLAS::fassign (Fi, R2-i, F+(i+1)*lda+j, lda, A+R1*(lda+1)+2*i+2*j*lda, 2);
            // FFLAS::WriteMatrix(std::cerr<<"2/ After interleaving L2 U2 H1 A = "<<std::endl,Fi,N,N,A,lda);
            //     // 3/ end of the iteration
            // dim = R2-i;
            // for (; i < R2; i++, Di+=2*lda+2, Fptr+=lda+1, H1i+=lda){
            //     Fi.assign(*Di,*Fptr); // pivot
            //     FFLAS::fassign (Fi, dim-1, Fptr+lda, lda, Di+2, 2); // copying L2
            //     FFLAS::fassign (Fi, dim, H1i, 1, Di+1, 2); // copying H1
            //     FFLAS::fassign (Fi, dim, Fptr, 1, Di+lda+1, 2); // copying U2
            // }
            FFLAS::WriteMatrix(std::cerr<<"After interleaving L2 U2 H1 A = "<<std::endl,Fi,N,N,A,lda);
                // Interleaving V21 and H21 into tmpV21
            typename Field::Element_ptr tmpV21 = FFLAS::fflas_new(Fi, 2*R2, R3);
            typename Field::Element_ptr V21i=F+R2, H2i = A4+R2;
            tmpi=tmpV21;
            for (size_t i=0; i<R2; i++, V21i+=lda, H2i+=lda,tmpi+=2*R3){
                FFLAS::fassign (Fi, R3, H2i, 1, tmpi, 1); // copying H21
                FFLAS::fassign (Fi, R3, V21i, 1, tmpi+R3, 1); // copying V21
            }
            // FFLAS::WriteMatrix(std::cerr<<" V21= "<<std::endl,Fi,R2,R3,F+R2,lda);
            // FFLAS::WriteMatrix(std::cerr<<" H21= "<<std::endl,Fi,R2,R3,A4+R2,lda);
            // FFLAS::WriteMatrix(std::cerr<<"After interleaving V21 H21 into tmpV21= "<<std::endl,
            //                    Fi,2*R2,R3,tmpV21,R3);
                //FFLAS::fassign(Fi, R2, R3, F+R2, lda, tmpV21, R3);

                // Interleaving M2^T and 0 into V42
            Di = A+R1*(lda+1)+2*R2+R3;
            tmpi=tmp+R2*R2;
            for (size_t i=0; i<R2; i++, Di+=2*lda, tmpi+=N1-R1-R2){
                FFLAS::fassign (Fi, N1-R2-R1, tmpi, 1, Di, 1); // copying M2^T
                FFLAS::fzero (Fi, N1-R2-R1, Di+lda, 1); // copying 0s
            }
            FFLAS::fflas_delete(tmp);
            FFLAS::WriteMatrix(std::cerr<<"After interleaving V42 A = "<<std::endl,Fi,N,N,A,lda);
                // Copying tmpV21 into V41
                //Di=A+R1*(lda+1)+2*R2;
            FFLAS::fassign(Fi, 2*R2, R3, tmpV21, R3, A+R1*(lda+1)+2*R2,lda);
            
            // for (size_t i=0; i<R2; i++, Di+=2*lda, H2i+=lda,tmpi+=R3){
            //     std::cerr<<" copying H2i = "<<*H2i<<" into Di = "<<*(Di)<<std::endl;
            //     FFLAS::fassign (Fi, R3, H2i, 1, Di, 1); // copying H21
            //     std::cerr<<" copying tmpi = "<<*tmpi<<" into Di+lda = "<<*(Di+lda)<<std::endl;
            //     FFLAS::fassign (Fi, R3, tmpi, 1, Di+lda, 1); // copying V21
            // }
            FFLAS::fflas_delete(tmpV21);

                //FFLAS::fflas_delete(tmpM2);

                // Interleaving V22 and H22 into V43
            dim = N2-R2-R3;
            tmp = FFLAS::fflas_new(Fi, R2, dim);
            FFLAS::fassign(Fi, R2, dim, F+R2+R3, lda, tmp, dim);
            Di=A+R1*lda+N1+R2+R3;
            H2i=A4+R2+R3;
            tmpi=tmp;
            for (size_t i=0; i<R2; i++, Di+=2*lda, H2i+=lda,tmpi+=dim){
                FFLAS::fassign (Fi, dim, H2i, 1, Di, 1); // copying H22
                FFLAS::fassign (Fi, dim, tmp+i*dim, 1, Di+lda, 1); // copying V22
            }
            FFLAS::fflas_delete(tmp);
        }
        
        if (R1+R2<N1){
            typename Field::Element_ptr DU3 = A+(R1+2*R2)*(lda+1);
                // Moving U3
            FFLAS::fassign(Fi, R3, R3, H3, lda, DU3, lda);
                // Zeros right of U3
                //FFLAS::fzero(Fi, R3, N1-R1, DU3+R3,lda);
            FFLAS::fzero(Fi, R3, N1-R1, DU3+R3,lda);
                // Moving V3
            FFLAS::fassign(Fi, R3, N2-R2-R3, H3+R3, lda, A2+R2+R3+(R1+2*R2)*lda, lda);
            //     // Rotating [V1 E1 E21] -> [E1 E21 V1]
            // typename Field::Element_ptr tmp = FFLAS::fflas_new(Fi, R1, N1-R1-R2);
            // FFLAS::fassign(Fi, R1, N1-R1-R2, A+R1+R2, lda, tmp, N1-R1-R2);
            // FFLAS::fassign(Fi, R1, R2+R3, A+N1, lda, A+R1+R2, lda);
            // FFLAS::fassign(Fi, R1, N1-R1-R2, tmp, N1-R1-R2, A+R1+2*R2+R3, lda);
//-------------------------
        }
            // Computing the permutation matrix
        size_t *tmpP=FFLAS::fflas_new<size_t>(N1-R1);
        FFLAS::WritePermutation(std::cerr<<"MathP = ", MathP, N);
        std::copy(MathP+R1, MathP+N1, tmpP);
        // FFLAS::WritePermutation(std::cerr<<"tmpP = ", tmpP, N1-R1);
        // FFLAS::WritePermutation(std::cerr<<"after copy MathP = ", MathP, N);
        size_t * MP1=MathP+R1, *tmpi = tmpP, * MP2=MathP+N1;
            // Applying Q = [e_1 e_{R2} e_2 e_{R2+1} ...] to MathP[R1:R1+2*R2]
//        std::cerr<<"tmpi= "<<*tmpi<<std::endl;
        std::cerr<<"R2 = "<<R2<<std::endl;
        for (; tmpi != tmpP+R2; tmpi++, MP2++, MP1++){
//            std::cerr<<"->tmpi = "<<*tmpi<<std::endl;
                //           std::cerr<<"->MP2 = "<<*MP2<<std::endl;
            *MP1 = *tmpi;
            *(++MP1) = *MP2;
        }
            // Moving P3
        for (size_t i=0; i<R3; i++, MP1++){
            *MP1=MathP[N1+R2+i];
        }
            //      std::cerr<<"tmpi= "<<*tmpi<<std::endl;
        FFLAS::WritePermutation(std::cerr<<"after for MathP = ", MathP, N);
            // end of the permutation
        // std::cerr<<"N1= "<<N1<<" R1 = "<<R1<<" R2  = "<<R2<<std::endl;
        // std::cerr<<"Avant tmpi= "<<*tmpi<<" tmpi+N1-R1 = "<<tmpi[N1-R1]<<" MP1  = "<<MP1[0]<<std::endl;

        std::copy(tmpi, tmpi+N1-R1-R2, MP1);
        // std::cerr<<"Apres tmpi= "<<*tmpi<<" tmpi+N1-R1 = "<<tmpi[N1-R1]<<" MP1  = "<<MP1[0]<<std::endl;
        
        // FFLAS::WritePermutation(std::cerr<<"after std::copy MAthP = ", MathP, N);

        MathPerm2LAPACKPerm (P, MathP, N);
        // FFLAS::WritePermutation(std::cerr<<"MathP2LAPACK P = ", P, N);
        FFLAS::fflas_delete( MathP);
		for (size_t i=0; i<R1+2*R2+R3; i++){
            std::cerr<<"twoBlocks["<<i<<"] ) "<<twoBlocks[i]<<std::endl;
			if (twoBlocks[i]){
				P[i] = -P[i]-1;
				P[i+1] = -P[i+1]-1;
				i++;
			}
		}
        FFLAS::WriteMatrix(std::cerr<<"done rec A = "<<std::endl,Fi,N,N,A,lda);
//        FFLAS::WriteMatrix(std::cerr<<"invD A = "<<std::endl,Fi,1,N,invD,N);
        return R1+2*R2+R3;
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

	template <class Field>
	inline size_t fsytrf_RPM (const Field& F, const FFLAS::FFLAS_UPLO UpLo, const size_t N,
                            typename Field::Element_ptr A, const size_t lda,
                            size_t * P, size_t threshold){
		typename Field::Element_ptr Dinv = FFLAS::fflas_new(F,N);
        size_t rank;
        if (UpLo==FFLAS::FflasUpper)
//            rank = fsytrf_UP_RPM_BC_Crout (F, N, A, lda, Dinv, 1, P, threshold);
            rank = fsytrf_UP_RPM (F, N, A, lda, Dinv, 1, P, threshold);
        else 
            rank = fsytrf_LOW_RPM_BC_Crout (F, N, A, lda, Dinv, 1, P);
                //rank = fsytrf_LOW_RPM (F, N, A, lda, Dinv, 1, P, P, threshold);
        
        FFLAS::WriteMatrix(std::cerr<<"After fsytrf_nonunit A = "<<std::endl,F,N,N,A, lda);
        FFLAS::WriteMatrix(std::cerr<<"After fsytrf_nonunit Dinv = "<<std::endl,F,1,N,Dinv, 1);
		size_t incA = (UpLo==FFLAS::FflasUpper) ? 1 : lda;
		for (size_t i=0; i<N; i++)
			FFLAS::fscalin (F, N-i-1, Dinv[i], A+i*(lda+1)+incA, incA);
		FFLAS::fflas_delete(Dinv);
		return rank;
	}

    template <class Field>
    inline void
    getTridiagonal (const Field& F, const size_t N, const size_t R,
                    typename Field::ConstElement_ptr A, const size_t lda, size_t * P,
                    typename Field::Element_ptr T, const size_t ldt){

        FFLAS::fzero(F,N,N,T,ldt);
        for (size_t i=0; i<R; i++)
            if ((int)P[i] <0){
                F.assign (T[i*(ldt+1)+1], A[i*(lda+1)]);
                F.assign (T[i*(ldt+1)+ldt], A[i*(lda+1)]);
                i++;
            } else {
                F.assign (T[i*(ldt+1)], A[i*(lda+1)]);
            }
    }
    
} //FFPACK

#endif // __FFLASFFPACK_ffpack_fsytrf_INL
