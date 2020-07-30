/* ffpack/ffpack_bruhatgen.inl
 * Copyright (C) 2020 Quentin Houssier
 *
 * Written by Quentin Houssier <quentin.houssier@ecl18.ec-lyon.fr>
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


#ifndef __FFLASFFPACK_ffpack_bruhatgen_inl
#define __FFLASFFPACK_ffpack_bruhatgen_inl

namespace FFPACK{
template<class Field>
inline size_t LTBruhatGen (const Field& Fi, const FFLAS::FFLAS_DIAG diag,
            const size_t N,
           typename Field::Element_ptr A, const size_t lda,
           size_t * P, size_t * Q)
    {
        if (N == 1) {Fi.assign(*(A+0), Fi.zero);
        return(0);}
        
        size_t N2 = N >> 1; // N = colonnes
        FFLAS::FFLAS_DIAG OppDiag =(diag==FFLAS::FflasUnit)?FFLAS::FflasNonUnit : FFLAS::FflasUnit;
        size_t * P1 = FFLAS::fflas_new<size_t>(N-N2);
        size_t * Q1 = FFLAS::fflas_new<size_t>(N2);
        // A1 = P1 [ L1 ] [ U1 V1 ] Q1
        //         [ M1 ]
        size_t r1 = PLUQ (Fi, diag, N-N2, N2, A, lda, P1, Q1);
        LAPACKPerm2MathPerm (P, P1, N-N2);
        LAPACKPerm2MathPerm (Q, Q1, N2);
        typename Field::Element_ptr A2 = A + N2;
        typename Field::Element_ptr A3 = A + (N-N2)*lda;
        typename Field::Element_ptr F = A2 + r1*lda;
        typename Field::Element_ptr G = A3 + r1;
        // [ B1 ] <- P1^T A2
        // [ B2 ]
        applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N-N2, size_t(0), N-N2, A2, lda, P1);
        // [ C1 C2 ] <- A3 Q1^T
        applyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, N2, size_t(0), N2, A3, lda, Q1);

        // D <- L1^-1 B1
        ftrsm (Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, r1, N-N2, Fi.one, A, lda, A2, lda);
        // E <- C1 U1^-1
        ftrsm (Fi, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, diag, N2, r1, Fi.one, A, lda, A3, lda);
        // F <- B2 - M1 D
        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, N-N2-r1, N-N2, r1, Fi.mOne, A + r1*lda, lda, A2, lda, Fi.one, A2+r1*lda, lda);
        // G <- C2 - E V1
        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, N2, N2-r1, r1, Fi.mOne, A3, lda, A+r1, lda, Fi.one, A3+r1, lda);
        // Expand L1\U1 into Bruhat generator
        applyP(Fi,FFLAS::FflasLeft,FFLAS::FflasTrans, N2, size_t(0), N-N2, A, lda, P1);
        applyP(Fi, FFLAS::FflasRight,FFLAS::FflasNoTrans, N-N2,size_t(0), N2, A, lda, Q1);
        //On stock D et E
        typename Field::Element_ptr D= FFLAS::fflas_new(Fi, r1, N-N2); 
        FFLAS::fassign(Fi, r1,N-N2, A2, lda, D, N-N2);
        size_t ldd = N-N2;
        typename Field::Element_ptr E= FFLAS::fflas_new(Fi, N2, r1); 
        FFLAS::fassign(Fi, N2,r1, A3, lda, E, r1);
        size_t lde = r1;
        FFLAS::fzero(Fi, r1, N-N2, A2,lda);
        FFLAS::fzero(Fi, N2, r1, A3, lda);
        //H <- P1  [ 0 ]
        //         [ F ]
        applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasTrans, N-N2, size_t(0),N-N2, A2, lda, P1);
        //I <- [ 0 G ] Q1
        applyP (Fi, FFLAS::FflasRight, FFLAS::FflasNoTrans, N2, size_t(0),N2, A3, lda, Q1);
        FFLAS::fflas_delete(P1,Q1);
        //A2 <- LT-Bruhat(H)
        
        size_t r2 = LTBruhatGen(Fi,diag,N-N2, A2,lda, P+r1, Q+r1);
        for (size_t i=r1;i<r1+r2;i++) Q[i]+=N2;
        //Restore raws of D into their position in A2
        for (size_t i=0;i<r1;i++){
            size_t row = P[i];
            FFLAS::fassign (Fi, N-N2-1-row, D+i*ldd, 1, A2+row*lda, 1);
        }
        FFLAS::fflas_delete(D);
        //A3 <- LT-Bruhat(I)
        size_t r3 = LTBruhatGen (Fi,diag, N2, A3, lda, P+r1+r2, Q+r1+r2);
        for (size_t i=r1+r2;i<r1+r2+r3;i++) P[i]+=N-N2;
        for (size_t i=0;i<r1;i++){
            size_t col = Q[i];
            FFLAS::fassign(Fi, N2-1-col, E+i, lde, A3+col,lda);
        }
        FFLAS::fflas_delete(E);
    return(r1+r2+r3);		
    
    }

template<class Field>
inline void getLTBruhatGen(const Field& Fi, const size_t N, const size_t r,const size_t * P, const size_t * Q,
                           typename Field::Element_ptr R, const size_t ldr){
    FFLAS::fzero(Fi, N, N, R,ldr);
    for(size_t i=0;i<r;i++)
        Fi.assign(R[P[i]*ldr+Q[i]],Fi.one);
}

template<class Field>
inline void getLTBruhatGen(const Field& Fi, const FFLAS::FFLAS_UPLO Uplo,const FFLAS::FFLAS_DIAG diag,
                           const size_t N, const size_t r, const size_t *P, const size_t * Q,
                           typename Field::ConstElement_ptr A, const size_t lda,
                           typename Field::Element_ptr T, const size_t ldt)

{   FFLAS::fzero(Fi, N, N, T, N);
    if (Uplo==FFLAS::FflasUpper) {     //U
        if(diag==FFLAS::FflasNonUnit){
            for(size_t i=0; i<r;i++){
                size_t row = P[i];
                size_t col = Q[i];
                FFLAS::fassign(Fi, N-1-row-col, A+row*lda+col,1,T+row*ldt+col,1);
                for(size_t j=0;j<i;j++){
                    Fi.assign(T[row*ldt+Q[j]],Fi.zero);
                }
            }
        } else {
            for(size_t i=0; i<r;i++){
                size_t row = P[i];
                size_t col = Q[i];
                FFLAS::fassign(Fi, N-2-row-col, A+row*lda+col+1,1,T+row*ldt+col+1,1);
                Fi.assign(T[row*ldt+col],Fi.one);
                for(size_t j=0;j<i;j++){
                    Fi.assign(T[row*ldt+Q[j]],Fi.zero);
                }
            }
        }
    } else{ //L
        if(diag==FFLAS::FflasNonUnit){
            for(size_t i=0; i<r;i++){
                size_t row = P[i];
                size_t col = Q[i];

                FFLAS::fassign(Fi, N-1-row-col, A+row*lda+col,lda,T+row*ldt+col,ldt);
                for(size_t j=0;j<i;j++){
                    Fi.assign(T[P[j]*ldt+col],Fi.zero);
                }
            }
        } else {
            for(size_t i=0; i<r;i++){
                size_t row = P[i];
                size_t col = Q[i];
                FFLAS::fassign(Fi, N-2-row-col, A+(row+1)*lda+col,lda,T+(row+1)*ldt+col,ldt);
                Fi.assign(T[row*ldt+col],Fi.one);
                for(size_t j=0;j<i;j++){
                    Fi.assign(T[P[j]*ldt+col],Fi.zero);
                }
            }
        }
    }
}
    

inline size_t LTQSorder(const size_t N, const size_t r,const size_t * P, const size_t * Q){
    std::vector<bool> rows(N,false);
    std::vector<bool> cols(N,false);
    
    for (size_t i=0;i<r;i++)
    {
        rows[P[i]]=true;
        cols[Q[i]]=true;
    }
    size_t s=0;
    size_t t=0;
    if (rows[0])
    {
        s=1;
        t=1;
    }

    for (size_t i=1;i<N;i++)
    {
        if (rows[i])   t+=1;
        if (cols[N-1-i]) t-=1;
        s = std::max(s,t);
    }
    return(s);
}

template<class Field>
inline size_t CompressToBlockBiDiagonal(const Field&Fi, const FFLAS::FFLAS_UPLO Uplo, size_t N, size_t s, size_t r, const size_t *P, const size_t *Q,  typename Field::Element_ptr A, size_t lda, typename Field::Element_ptr X, size_t ldx, size_t *K, size_t *M, size_t *T){
 
  typename Field::Element_ptr C=FFLAS::fflas_new(Fi, N, N);
  size_t ldc=N;
  FFLAS::fassign(Fi, N, N, A, lda, C, ldc);
  typename Field::Element_ptr D = X + 0;
  typename Field::Element_ptr S;
  FFLAS::FFLAS_SIDE Side;
  FFLAS::FFLAS_TRANSPOSE Trans;
  FFLAS::FFLAS_SIDE OppSide;
  FFLAS::FFLAS_TRANSPOSE OppTrans;
  const size_t * outer;
  const size_t * inner;
  size_t * Inv = FFLAS::fflas_new<size_t>(N);
  if(Uplo==FFLAS::FflasUpper){
    S = X + s*ldx;
    outer = Q;
    inner = P;
    Side = FFLAS::FflasLeft;
    Trans = FFLAS::FflasTrans;
    OppSide = FFLAS::FflasRight;
    OppTrans = FFLAS::FflasNoTrans;
      
  }
  else{
    S = X+s;
    outer = P;
    inner = Q;
    Side = FFLAS::FflasRight;
    Trans = FFLAS::FflasNoTrans;
    OppSide = FFLAS::FflasLeft;
    OppTrans = FFLAS::FflasTrans;
  }

   for (size_t i=0; i<r; i++){
         Inv [inner [i]] = i;
      }
   Bruhat2EchelonPermutation (N,r,outer,inner,M);
    size_t * MLap = FFLAS::fflas_new<size_t>(N);
    MathPerm2LAPACKPerm(MLap, M, N);
    applyP (Fi, Side, Trans, N, size_t(0), N, C, lda, MLap);
    size_t * pivot_outer_pos = FFLAS::fflas_new<size_t>(r);
  size_t * last_coeff = FFLAS::fflas_new<size_t>(r);
  for (size_t i=0;i<r;i++)
    {
      pivot_outer_pos[M[inner[i]]]= outer[i];
      last_coeff[M[inner[i]]] = N-2-inner[i];
    }
  if(Uplo==FFLAS::FflasUpper)
    FFLAS::fzero(Fi, 2*s,N,X,ldx);
  else
    FFLAS::fzero(Fi, N, 2*s, X,ldx);
  size_t CurrentBlockPos = 0;
  size_t NbBlocks = 0;
  size_t BlockPivot = 0;
  while (BlockPivot<r)
    { K[NbBlocks] = CurrentBlockPos;
      size_t BlockSize, NextBlockPos;
      if (BlockPivot+s >= r)
	{
	  BlockSize =r-BlockPivot;
	  NextBlockPos = N;
        } else{
	  BlockSize = s;
	  NextBlockPos = pivot_outer_pos[BlockPivot+BlockSize];}
	  
        if (Uplo==FFLAS::FflasUpper){//U
        FFLAS::fassign(Fi, BlockSize, NextBlockPos-CurrentBlockPos, C+CurrentBlockPos+BlockPivot*ldc,ldc,D+CurrentBlockPos,ldx);//On stock Di
	 
        }
        else{//L
        FFLAS::fassign(Fi, NextBlockPos-CurrentBlockPos, BlockSize, C+CurrentBlockPos*ldc+BlockPivot,ldc,D+CurrentBlockPos*ldx,ldx);//On stock Di
        }
	
      
      CurrentBlockPos=NextBlockPos;
      BlockPivot+=BlockSize; 
      NbBlocks++;
      
    }
   
     
  K[NbBlocks]=N;
  for (size_t i=0; i<r;i++)
    T[i]=i;
      //Construction de S
  for (size_t j=1;j<NbBlocks;j++)
    {if (Uplo==FFLAS::FflasUpper){//U
	FFLAS::fassign(Fi,s, K[j+1]-K[j], C+K[j]+(j-1)*s*ldc,ldc, S+K[j],ldx);
      }
      else{
	FFLAS::fassign(Fi, K[j+1]-K[j], s, C+K[j]*ldc+(j-1)*s,ldc, S+K[j]*ldx,ldx);
      }
      for(size_t t=0;t<s;t++)
        {
          if(last_coeff[(j-1)*s+t] >= K[j+1])//Si la colonne n'est pas nulle
            {
	      size_t l=0;
              while(last_coeff[j*s+l]>=K[j+1])
                {   
		  l++;
                }
	      if (Uplo==FFLAS::FflasUpper){//U
		FFLAS::fassign(Fi,last_coeff[(j-1)*s+t]-K[j+1]+1 , C+K[j+1]+((j-1)*s+t)*ldc,1, C+K[j+1]+(j*s+l)*ldc, 1);
	      }else{FFLAS::fassign(Fi,last_coeff[(j-1)*s+t]-K[j+1]+1 , C+K[j+1]*ldc+(j-1)*s+t,ldc, C+K[j+1]*ldc+j*s+l, ldc);}
	      T[(j-1)*s+t]= j*s+l;
	      last_coeff[j*s+l]= last_coeff[(j-1)*s+t];
	    }
        }
    }
 
  FFLAS::fflas_delete(C);
  FFLAS::fflas_delete(last_coeff);
  return(NbBlocks);

}
template<class Field>
inline void  ExpandBlockBiDiagonalToBruhat(const Field&Fi, const FFLAS::FFLAS_UPLO Uplo, size_t N, size_t s, size_t r, typename Field::Element_ptr A, size_t lda,typename Field::Element_ptr X, size_t ldx,size_t NbBlocks,size_t *K, size_t *M, size_t *T){

  FFLAS::fzero(Fi,N,N,A,lda);
  typename Field::Element_ptr D = X + 0;
  //We copy S
if (Uplo==FFLAS::FflasUpper)//U
  {typename Field::Element_ptr S = X + s*ldx;
    for (size_t j=1;j<NbBlocks;j++){
      FFLAS::fassign(Fi,s, K[j+1]-K[j],S+K[j],ldx, A+K[j]+(j-1)*s*lda,lda);}
  }
 else{
   typename Field::Element_ptr S = X+s;
   for (size_t j=1;j<NbBlocks;j++){
     FFLAS::fassign(Fi, K[j+1]-K[j], s,S+K[j]*ldx,ldx, A+K[j]*lda+(j-1)*s,lda);
    }
   
 }
 
//Extend S
 for (size_t i=1; i < NbBlocks-1; i++){
   for (size_t l=0; l<s; l++){
     if (T[(NbBlocks-i-2)*s+l] != (NbBlocks-i-2)*s+l){
	 if(Uplo==FFLAS::FflasUpper){
	   FFLAS::fassign(Fi, N-K[NbBlocks-i], A+K[NbBlocks-i]+T[(NbBlocks-i-2)*s+l]*lda, 1, A+K[NbBlocks-i]+((NbBlocks-i-2)*s+l)*lda,1);
	   FFLAS::fzero(Fi, N-K[NbBlocks-i], A+K[NbBlocks-i]+T[(NbBlocks-i-2)*s+l]*lda,1);
	 }else{
	   FFLAS::fassign(Fi, N-K[NbBlocks-i], A+K[NbBlocks-i]*lda+T[(NbBlocks-i-2)*s+l], lda, A+K[NbBlocks-i]*lda+(NbBlocks-i-2)*s+l,lda);
	   FFLAS::fzero(Fi, N-K[NbBlocks-i], A+K[NbBlocks-i]*lda+T[(NbBlocks-i-2)*s+l],lda);
	 }
       }
   }
   
 }
     //We copy D
     for (size_t i=0; i<NbBlocks-1;i++){
       if (Uplo==FFLAS::FflasUpper){//U
	 FFLAS::fassign(Fi, s, K[i+1]-K[i], D+K[i],ldx,A+K[i]+i*s*lda,lda);
	}
	else{//L
	  FFLAS::fassign(Fi, K[i+1]-K[i], s,D+K[i]*ldx,ldx, A+K[i]*lda+i*s,lda);
	}
     }
     if (Uplo==FFLAS::FflasUpper){//U
       FFLAS::fassign(Fi, r-(NbBlocks-1)*s, N-K[NbBlocks-1], D+K[NbBlocks-1],ldx,A+K[NbBlocks-1]+(NbBlocks-1)*s*lda,lda);
	}
	else{//L
	  FFLAS::fassign(Fi, N-K[NbBlocks-1], r-(NbBlocks-1)*s,D+K[NbBlocks-1]*ldx,ldx, A+K[NbBlocks-1]*lda+(NbBlocks-1)*s,lda);
	}
     //We Apply M-1
     size_t * MLap = FFLAS::fflas_new<size_t>(N);
    MathPerm2LAPACKPerm(MLap, M, N);
     if (Uplo==FFLAS::FflasUpper)//U
  {
    applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N, size_t(0), N, A, lda, MLap);
  }
 else{
  applyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, N, size_t(0), N, A, lda, MLap);
 }

}
 // Compute the permutation M (in MathPermutation format) such that LM^T is in column echelon form, where L is the
 // left factor of a Bruhat decomposition
      // M has to be by the caller of this function.
  void Bruhat2EchelonPermutation (size_t N,size_t R, const size_t* P,const size_t *Q, size_t* M){

      size_t * Pinv = FFLAS::fflas_new<size_t>(N);
      size_t * Ps = FFLAS::fflas_new<size_t>(R);

      for (size_t i=0; i<R; i++){
          Pinv [P [i]] = i;
          Ps[i] = P[i];
      }
      std::sort (Ps, Ps+R);
      std::vector<bool> ispivot(N,false);
      for (size_t i=0; i<R; i++){
          size_t piv = Q [Pinv [Ps[i]]];
          M[piv] = i;
          ispivot [piv]=true;
      }
      size_t curr=R;
      for (size_t i=0; i<N; ++i)
          if (!ispivot[i])
              M [i] = curr++;

      FFLAS::fflas_delete(Pinv,Ps);
   }

//Compute a table Tinv that gives us all the information we need to expand S. If T[i]=j means that line i of S has been fold into line j
// Then Tinv[j]=i means that line j was originally into line 1

size_t * TInverter (const size_t * T, size_t r)
    {
      size_t * Tinv = FFLAS::fflas_new<size_t>(r);
      std::vector<bool> IsFinal(r, true); //Stock the information if Line i is a line where no other line fold into
      for(size_t i=0;i<r;i++)
	{
	  if(IsFinal[i]) //First time we meet the line, means not other line fold into
	  {
	    Tinv[i]=i;
	    if (T[i]!=i)//if the line was fold into line T[i]
	      {
		Tinv[T[i]]=i;
		IsFinal[T[i]]=false;// Line T[i] is a line we have now met
	      }
	  }
	  else if(T[i]!=i)//A line was already fold into line i
	  {
	    Tinv[T[i]]=Tinv[i];
	    IsFinal[T[i]]=false;
	  }
      
	}
    
      return(Tinv);
    }

//Compute Rtranspose in the CRE decomposition, must be apply to the left of a matrix. R[col]=line
template <class Field>
inline void ComputeRPermutation (const Field&Fi, size_t N, size_t r, const size_t * P, const size_t * Q,
                                 size_t * R, const size_t * MU, const size_t * ML){
       for (size_t i=0;i<r;i++){
           R[MU[P[i]]] = ML[Q[i]];
       }
    }

/**
 * @brief Compute the product of a left-triangular quasi-separable matrix A, represented by a compact Bruhat generator, 
 * with a dense rectangular matrix B:  \f$ C \gets A \times B + beta C \f$
 *
 * @param F the base field
 * @param N the order of \p A
 * @param s the order of quasiseparability of \p A
 * @param r the number of pivots in the left-triangular par of the rank profile matrix of \p A
 * @param t the number of columns of \p B
 * @param P the row indices of the pivots of \p A
 * @param Q the column indices of the pivots of \p A
 * @param Xu the compact storage of U: Du blocks in the first s rows, Su blocks in the last s rows
 * @param ldxu the leading dimension of \p Xu
 * @param NbBlocksU the number of diagonal blocks in the compact storage of U
 * @param Ku the list of starting column positions for each block of the storage of U
 * @param Tu the folding matrix for the compact storage of U: \f$ Du + Tu  Su\f$ is in row echelon form
 * @param Mu a permutation matrix such that \f$ Mu (Du + Tu  Su)\f$ is the U factor of the Bruhat generator
 * @param Xl the compact storage of L: Dl blocks in the first s columns, Sl blocks in the last s columns
 * @param ldxl the leading dimension of \p Xl
 * @param NbBlocksL the number of diagonal blocks in the compact storage of L
 * @param Kl the list of starting row positions for each block of the storage of L
 * @param Tl the folding matrix for the compact storage of L: \f$ Dl + Sl  Tl \f$ is in column echelon form
 * @param Ml a permutation matrix such that \f$(Dl + Tl Sl) Ml \f$ is the L factor of the Bruhat generator
 * @param B an \f$ N \times t\f$ dense matrix
 * @param ldb leading dimension of \p B
 * @param beta scaling constant
 * @param [inout] C output matrix
 * @param ldc leading dimension of \p C
 */
template<class Field>
inline  void productBruhatxTS (const Field& Fi, size_t N, size_t s, size_t r, size_t t, const size_t *P, const size_t *Q,
                               typename Field::ConstElement_ptr Xu, size_t ldu, size_t NbBlocksU,
                               const size_t * Ku, const size_t *Tu , const size_t * MU,
                               typename Field::ConstElement_ptr Xl, size_t ldl, size_t NbBlocksL,
                               const size_t *Kl, const size_t *Tl, const size_t * ML,
                               typename  Field::Element_ptr B, size_t ldb,
                               const typename Field::Element beta,
                               typename Field::Element_ptr C, size_t ldc)
    {
        size_t * Tuinv = TInverter(Tu, r);
        size_t * Tlinv = TInverter(Tl, r);
        size_t k = N/s;  // Nb of slices of dimension s 
        size_t rs = N%s;
        if (rs) k++;
            /* A is split on a k x k grid of block-size s: 
             * columns slices have dimension (s, s, ..., rs) and row slices have dimension: (rs, s, s, ..., s)
             * so that the all blocks on the anti-diagonal have dimension sxs except the top-right one which is rs x rs
             * 
             * A = Left (C R E) where 
             *  - R is a permutation matrix (deduced from P, Q, Mu, Ml)
             *  - C = Dl + Sl x Tl is in column echelon form, where
             *    * Dl is block diagonal with slices of col. dim. (s,s,..,r%s) and starting at row position Kl[0],Kl[1]...Kl[NbBlocksL-1]
             *    * Sl is block sub-diagonal with slices of col. dim. (s,s,..,r%s) and starting at row position Kl[1]...Kl[NbBlocksL-1]
             *      Storage Xl = [ Dl[0] | Dl[1] | ... | Dl[NbBlocksL-2] | Dl[NbBlocksL-1] ] ^ T
             *                   [ Sl[0] | Sl[1] | ... | Sl[NbBlocksL-2] |                 ]
             *  - E = Du + Tu x Su  is in row echelon form, where
             *    * Du is block diagonal with slices of row dim. (s,s,..,r%s) and starting at column position Ku[0],Ku[1]...Ku[NbBlocksL-1]
             *    * Su is block sub-diagonal with slices of row. dim. (s,s,..,r%s) and starting at column position Ku[1]...Ku[NbBlocksL-1]
             *      Storage Xu = [ Du[0] | Du[1] | ... | Du[NbBlocksL-2] | Du[NbBlocksL-1] ]
             *                   [ Su[0] | Su[1] | ... | Su[NbBlocksL-2] |                 ]
             */
            
        size_t * R = FFLAS::fflas_new<size_t>(r);
        ComputeRPermutation(Fi, N, r, P, Q, R, MU, ML);
      
            // std::cerr<<"Entering CompactBruhat x TS"<<std::endl;
            // std::cerr<<"  Pivots: ";
            // for (size_t i=0; i<r; i++) std::cerr<<" ("<<P[i]<<", "<<Q[i]<<"),  ";
            // std::cerr<<std::endl;
            // FFLAS::WritePermutation(std::cerr<<"  Block structure of U: KU="<<std::endl,Ku, NbBlocksU+1)<<std::endl;
            // FFLAS::WritePermutation(std::cerr<<"  Block structure of L: KL="<<std::endl,Kl, NbBlocksL+1)<<std::endl;
            // FFLAS::WriteMatrix(std::cerr<<"  Xu = "<<std::endl,Fi, 2*s, N, Xu, ldu)<<std::endl;
            // FFLAS::WriteMatrix(std::cerr<<"  Xl = "<<std::endl,Fi, N, 2*s, Xl, ldl)<<std::endl;
        
            // std::cerr<<"N,s,rs,k = "<<N<<" "<<s<<" "<<rs<<" "<<k<<std::endl;

            //Gives the information about our position in XU (line = blocksu*s) and (Xl column= blocksl*s)
        size_t currbkU = 0;
        size_t currbkL = NbBlocksL;
      // The parial sums do not involve the last rs rows of C: skipping the last row slice of L when in the last s rows.
      if (Kl [NbBlocksL-1] >=  N-s) currbkL--;

      typename Field::Element_ptr SX= FFLAS::fflas_new(Fi, 2*s,t);
      typename Field::Element_ptr DX = FFLAS::fflas_new(Fi,2*s,t);
      typename Field::Element_ptr Z = FFLAS::fflas_new(Fi, r, t);
      typename Field::Element_ptr TlZ = FFLAS::fflas_new(Fi, r, t);
      FFLAS::fzero(Fi, r,t, Z,t);

      for (size_t i = 0; i < k-1; i++) { // Loop over the slices Bi of s rows of B except the last one

              // 1. Compute  X = Ei x Bi
              //    by storing Dx and Sx such that X = Dx + Tu x Sx
          
          FFLAS::fzero(Fi, 2*s, t, DX, t);
	  FFLAS::fzero(Fi, 2*s, t, SX, t);

	  if (Ku[currbkU+1]-i*s<s) { // current slice of E intersects two blocks of Du

              size_t fst_ncols = Ku[currbkU+1] - i*s;
              size_t snd_ncols = s-fst_ncols;
                  //We should not consider a line superior to r
              size_t bkU_nrows = (currbkU+1<NbBlocksU-1) ? s :  r-(currbkU+1)*s;
		
		//DX<-DjU*Bj
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,s, t, fst_ncols, Fi.one, Xu+i*s, ldu, B+i*s*ldb, ldb, Fi.zero, DX, t);
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,bkU_nrows, t, snd_ncols, Fi.one, Xu+Ku[currbkU+1], ldu, B+Ku[currbkU+1]*ldb, ldb, Fi.zero, DX+s*t, t);

                // FFLAS::WriteMatrix(std::cerr<<"  DX = "<<std::endl,Fi, s+bkU_nrows, t, DX, t)<<std::endl;
                // FFLAS::WriteMatrix(std::cout<<" avant accumulation Z="<<std::endl,Fi,r,t,Z,t)<<std::endl;

                    //Apply R and add to Z
		for (size_t l=0;l<s;l++){
                    FFLAS::faddin(Fi, t, DX+l*t, 1,Z+R[currbkU*s+l]*t ,1);
                }
		for (size_t l=0; l<bkU_nrows; l++) {
                    FFLAS::faddin(Fi, t, DX+(s+l)*t, 1, Z+R[(currbkU+1)*s+l]*t,1);
                }

                    //FFLAS::WriteMatrix(std::cout<<" apres accumulation Z="<<std::endl,Fi,r,t,Z,t)<<std::endl;
		//SX<-SjU*Bj 
		if(currbkU>0)
                    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,s, t, fst_ncols, Fi.one, Xu+i*s+s*ldu, ldu, B+i*s*ldb, ldb, Fi.zero, SX, t);

                fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,s, t, snd_ncols, Fi.one, Xu+Ku[currbkU+1]+s*ldu, ldu, B+Ku[currbkU+1]*ldb, ldb, Fi.zero, SX+s*t, t);
                
                    //FFLAS::WriteMatrix(std::cerr<<"  SX = "<<std::endl,Fi, 2*s, t, SX, t)<<std::endl;

		//Apply R*Tu and add to Z
		if(currbkU>0){
		    for(size_t l=0; l<s; l++)
		      FFLAS::faddin(Fi,t,SX+l*t,1,Z+R[Tuinv[(currbkU-1)*s+l]]*t,1);
		  }
		for (size_t l=0; l<s; l++){
                        //std::cerr<<"currbkU = "<<currbkU<<" currbkU*s+l = "<<currbkU*s+l<<" Tuinv[...] = "<<Tuinv[currbkU*s+l]<<std::endl;
                    FFLAS::faddin(Fi,t,SX+(l+s)*t,1,Z+R[Tuinv[currbkU*s+l]]*t,1);
                }
          } else { // s-dimensional column slice doesn't split a block in U
              size_t bkU_nrows = (currbkU<NbBlocksU-1) ? s : r-(currbkU)*s;
              
                  // DX<-DjU*Bj
              fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,bkU_nrows, t, s, Fi.one, Xu+i*s, ldu, B+i*s*ldb, ldb, Fi.zero, DX, t);
		//Apply R and add to Z
	
		for (size_t l=0;l<bkU_nrows;l++)
		  FFLAS::faddin(Fi, t, DX+l*t ,1,Z+R[currbkU*s+l]*t ,1);
		// SX<-SjU*Bj
		if(currbkU>0)
		  {fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,s, t, s, Fi.one, Xu+i*s+s*ldu, ldu, B+i*s*ldb, ldb, Fi.zero, SX, t);
		    //Apply R*Tu and add to Z
		    for (size_t l=0; l<s; l++)
		      FFLAS::faddin(Fi,t, SX+l*t, 1, Z+R[Tuinv[(currbkU-1)*s+l]]*t,1);
		   }
	      }

              //Compute Tl x Zj
	    for(size_t l=0;l<r;l++){ 
		FFLAS::fassign(Fi,t,Z+Tlinv[l]*t,1,TlZ+l*t,1);
            }
	    
            size_t bkL_ncols =  (currbkL<NbBlocksL) ? s : r-(currbkL-1)*s;
            size_t bkC_nrows = (i==k-2 && rs) ? rs : s; // last slice of C has rs rows
            size_t currbkrow = N - (i+1)*s - bkC_nrows;
                //std::cerr<<"i = "<<i<<" currbkrow = "<<currbkrow<<" bkL_ncols = "<<bkL_ncols<<std::endl;
	    if (currbkrow < Kl[currbkL-1]){ // If the current row slice intersects 2  blocks Dj
                    //DjL*Zj
                fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,Kl[currbkL-1]-currbkrow, t, s, Fi.one, Xl+currbkrow*ldl, ldl, Z+(currbkL-2)*s*t, t, beta, C+currbkrow*ldc, ldc);
	        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, N-(i+1)*s-Kl[currbkL-1],t, bkL_ncols, Fi.one, Xl+Kl[currbkL-1]*ldl, ldl, Z+(currbkL-1)*s*t, t, beta, C+Kl[currbkL-1]*ldc, ldc);
		
		//SjL*Zj
		if(currbkL>2)
                    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,Kl[currbkL-1]-currbkrow, t, s, Fi.one, Xl+currbkrow*ldl+s, ldl, TlZ+(currbkL-3)*s*t, t, Fi.one, C+currbkrow*ldc, ldc);
	        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,N-(i+1)*s-Kl[currbkL-1], t, s, Fi.one, Xl+Kl[currbkL-1]*ldl+s, ldl, TlZ+(currbkL-2)*s*t, t, Fi.one, C+Kl[currbkL-1]*ldc, ldc);	
            } else {
		 //DjL*Zj
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, bkC_nrows, t, bkL_ncols, Fi.one, Xl+currbkrow*ldl, ldl, Z+(currbkL-1)*s*t, t, beta, C+currbkrow*ldc, ldc);
		
		//SjL*Zj
		if (currbkL > 1)
                    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, bkC_nrows, t, s, Fi.one, Xl+currbkrow*ldl + s, ldl, TlZ+(currbkL-2)*s*t, t, Fi.one, C+currbkrow*ldc, ldc);
	      }
                //FFLAS::faddin(Fi,s,t,SX,t,C+(S-(i+1)*s)*ldc,ldc);
	    // FFLAS::WriteMatrix(std::cout<<"C="<<std::endl,Fi, N, t, C, t)<<std::endl;
	    //Next Block
	    if (Ku[currbkU+1] < (i+1)*s) // one block Dj has been completed in U
	      currbkU++;
	    if (currbkrow < Kl[currbkL-1]) // one block Dj has been completed in L
	      currbkL--;
	}
	FFLAS::fflas_delete(SX,DX,Z, TlZ);

            // TODO : move the following code into a dedicated routine: expandLCRE
        currbkU = 0;
	currbkL = NbBlocksL;  
	size_t row_pos = N;
	typename Field::Element_ptr trailing_term = FFLAS::fflas_new(Fi, s, t);
	typename Field::Element_ptr Er = FFLAS::fflas_new(Fi, r, s);
	typename Field::Element_ptr TlEr = FFLAS::fflas_new(Fi, r,s);
	typename Field::Element_ptr CRE = FFLAS::fflas_new(Fi, s,s);
//	typename Field::Element_ptr SCRE = FFLAS::fflas_new(Fi, s, s);
	for(size_t i=0; i<k; i++) { 
            size_t grid_dim = (i == k-1 && rs) ? rs : s;
            	    
	    //Compute (Left(CN-i+1REi)Bi) the trailing term
	       FFLAS::fzero(Fi, r, s, Er, s);
	       FFLAS::fzero(Fi, r,s, TlEr, s);
	       //Expand Ei and apply R
	    if(Ku[currbkU+1]-i*s<grid_dim) {
                size_t bkU_nrows =  (currbkU+1 == NbBlocksU-1) ? r - (currbkU+1)*s : s;

		if(currbkU>0)
		  {for (size_t l=0; l<s; l++)
		      { //S
			FFLAS::fassign(Fi,Ku[currbkU+1]-i*s,Xu+i*s+(s+l)*ldu,1, Er+R[Tuinv[(currbkU-1)*s+l]]*s,1);
		    
		  }}
		for (size_t l=0;l<s;l++)
		  { //S
		    FFLAS::fassign(Fi, i*s+grid_dim - Ku[currbkU+1], Xu+Ku[currbkU+1]+(s+l)*ldu,1, Er+R[Tuinv[currbkU*s+l]]*s+Ku[currbkU+1]-i*s,1);
		    //D
		    FFLAS::fassign(Fi,Ku[currbkU+1]-i*s, Xu+i*s+l*ldu,1 ,Er+R[currbkU*s+l]*s,1);}
		for (size_t l=0; l<bkU_nrows;l++)
		    FFLAS::fassign(Fi,  i*s+grid_dim -Ku[currbkU+1], Xu+Ku[currbkU+1]+l*ldu,1, Er+R[(currbkU+1)*s+l]*s+Ku[currbkU+1]-i*s,1);
		  
	      }
	    else {
                size_t bkU_nrows = (currbkU == NbBlocksU-1) ? r - currbkU*s : s;

                if(currbkU>0) {
                    for (size_t l=0; l<s; l++) { 
                        // FFLAS::WriteMatrix(std::cout<<"copying row "<<l<<"+s of Xu (D) = "<<std::endl,Fi, 1,grid_dim, Xu+i*s+(s+l)*ldu+l, ldu)<<std::endl;
                        FFLAS::fassign (Fi, grid_dim, Xu + i*s + (s+l)*ldu , 1, Er+R[Tuinv[(currbkU-1)*s+l]]*s,1);
                    }
                }
		for (size_t l=0; l<bkU_nrows; l++) {
                         // FFLAS::WriteMatrix(std::cout<<"copying row "<<l<<" of Xu = "<<std::endl,Fi, 1,grid_dim, Xu+i*s+l*ldu, ldu)<<std::endl;

                    FFLAS::fassign (Fi, grid_dim, Xu + i*s + l*ldu, 1, Er+R[currbkU*s+l]*s, 1);
                }
	      }
                //FFLAS::WriteMatrix(std::cout<<"Er="<<std::endl,Fi, r,s, Er, s)<<std::endl;

	    //Compute Left(CRE)
	    for (size_t j=0; j<r;j++)
	      {
		FFLAS::faddin(Fi,grid_dim,Er+Tlinv[j]*s,1,TlEr+j*s,1);
	      }
                // std::cerr<<"row_pos = "<<row_pos<<" currbkL-1 = "<<currbkL-1<<" Kl[currbkL-1] = "<<Kl[currbkL-1]<<" grid_dim = "<<grid_dim<<std::endl;

            size_t inner_dim = (currbkL<NbBlocksL) ? s : r-(NbBlocksL-1)*s;

            if (row_pos-Kl[currbkL-1]<grid_dim) { // the current row slice intersects 2  blocks Dj
                    //DN-j
                // std::cerr<<"ici"<<std::endl;
                if(currbkL>1){
                    // std::cerr<<"currbkL>1"<<std::endl;
                    fgemm (Fi,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans, Kl[currbkL-1]-(row_pos-grid_dim),grid_dim,s, Fi.one, Xl+(row_pos-grid_dim)*ldl,ldl,Er+(currbkL-2)*s *s,s,Fi.zero, CRE,s);
                }

                fgemm (Fi,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans, row_pos-Kl[currbkL-1],grid_dim, inner_dim, Fi.one, Xl+Kl[currbkL-1]*ldl,ldl,Er+(currbkL-1)*s *s,s, Fi.zero,CRE+(Kl[currbkL-1]-(row_pos-grid_dim))*s,s);

                    //SN-j
		if(currbkL>2)
		  fgemm (Fi,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans, Kl[currbkL-1]-(row_pos-grid_dim),grid_dim,s, Fi.one, Xl+(row_pos-grid_dim)*ldl+s,ldl,TlEr+(currbkL-3)*s *s,s, Fi.one,CRE,s);
		if(currbkL>1)
		  fgemm (Fi,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans, row_pos-Kl[currbkL-1],grid_dim,s, Fi.one, Xl+Kl[currbkL-1]*ldl+s,ldl,TlEr+(currbkL-2)*s *s,s, Fi.one,CRE+(Kl[currbkL-1]-(row_pos-grid_dim))*s,s);
	      }
	    else
	      { //DN-j
		  fgemm (Fi,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans, grid_dim,grid_dim, inner_dim, Fi.one, Xl+(row_pos-grid_dim)*ldl,ldl,Er+(currbkL-1)*s *s,s, Fi.zero,CRE,s);
		//SN-j
		if(currbkL>1)
		  fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,grid_dim,grid_dim,s, Fi.one, Xl+(row_pos-grid_dim)*ldl+s,ldl,TlEr+(currbkL-2)*s*s,s, Fi.one,CRE,s);
	      }
                //FFLAS::faddin(Fi, grid_dim, grid_dim, SCRE, s, CRE, s);
                // Left( )
            for (size_t i=0; i < grid_dim; i++)
                FFLAS::fzero(Fi, i+1, CRE + i*s + grid_dim-i-1, 1);

            // FFLAS::WriteMatrix(std::cout<<"CRE ="<<std::endl,Fi, grid_dim,grid_dim, CRE, s)<<std::endl;

            // FFLAS::WriteMatrix(std::cout<<"Befoe fgemm C="<<std::endl,Fi, N,t, C, ldc)<<std::endl;

	    fgemm(Fi,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans, grid_dim, t, grid_dim, Fi.one, CRE, s, B+i*s*t, t, (row_pos == N) ? beta : Fi.one, C+(row_pos-grid_dim)*ldc, ldc);
                //FFLAS::faddin(Fi,grid_dim,t,trailing_term,t,C+(row_pos-grid_dim)*ldc,ldc);
	    
	    //Next Block
	    row_pos -= grid_dim;
	    if (Ku[currbkU+1]-i*s < grid_dim)
	      currbkU++;
//	    if (row_pos-Kl[currbkL-1] < grid_dim)
	    if (row_pos < Kl[currbkL-1])
	      currbkL--;
        }
	FFLAS::fflas_delete (CRE, trailing_term, Er, TlEr);
	
    }
  
} //namespace FFPACK
#endif //_FFPACK_ffpack_bruhatgen_inl
