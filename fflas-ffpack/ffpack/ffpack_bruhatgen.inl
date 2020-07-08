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
inline void getLTBruhatGen(const Field& Fi, const size_t N, const size_t r,const size_t * P, const size_t * Q, typename Field::Element_ptr R, const size_t ldr)
{
  FFLAS::fzero(Fi, N, N, R,ldr);
    for(size_t i=0;i<r;i++){
            size_t row = P[i];
            size_t col = Q[i];
            Fi.assign(R[row*ldr+col],Fi.one);
        }
    
}
template<class Field>
inline void getLTBruhatGen(const Field& Fi, const FFLAS::FFLAS_UPLO Uplo,const FFLAS::FFLAS_DIAG diag,
                           const size_t N, const size_t r, const size_t *P, const size_t * Q,
                           typename Field::ConstElement_ptr A, const size_t lda,
                           typename Field::Element_ptr T, const size_t ldt)

{   FFLAS::fzero(Fi, N, N, T, N);
    //U
    if (Uplo==FFLAS::FflasUpper) {
      if(diag==FFLAS::FflasNonUnit){
          for(size_t i=0; i<r;i++){
            size_t row = P[i];
            size_t col = Q[i];
            FFLAS::fassign(Fi, N-1-row-col, A+row*lda+col,1,T+row*ldt+col,1);
            for(size_t j=0;j<i;j++){
                Fi.assign(T[row*ldt+Q[j]],Fi.zero);
            }
          }
      }
      else
        {
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
    }
    //L
    else
      {
        if(diag==FFLAS::FflasNonUnit){

          for(size_t i=0; i<r;i++){
            size_t row = P[i];
            size_t col = Q[i];

            FFLAS::fassign(Fi, N-1-row-col, A+row*lda+col,lda,T+row*ldt+col,ldt);
            for(size_t j=0;j<i;j++){
                Fi.assign(T[P[j]*ldt+col],Fi.zero);
            }
          
          }
        }

        else
          {
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
    Trans = FFLAS::FflasNoTrans;
    OppSide = FFLAS::FflasRight;
    OppTrans = FFLAS::FflasTrans;
      
  }
  else{
    S = X+s;
    outer = P;
    inner = Q;
    Side = FFLAS::FflasRight;
    Trans = FFLAS::FflasTrans;
    OppSide = FFLAS::FflasLeft;
    OppTrans = FFLAS::FflasNoTrans;
  }

   for (size_t i=0; i<r; i++){
         Inv [inner [i]] = i;
      }
   Bruhat2EchelonPermutation (N,r,outer,inner,M);
    size_t * MLap = FFLAS::fflas_new<size_t>(N);
    MathPerm2LAPACKPerm(MLap, M, N);
    applyP (Fi, Side, Trans, N, size_t(0), N, C, lda, MLap);
 
  size_t * last_coeff = FFLAS::fflas_new<size_t>(r);
  for (size_t i=0;i<r;i++)
    {
      last_coeff[i] = N-2-M[i];
    }
  if(Uplo==FFLAS::FflasUpper)
    FFLAS::fzero(Fi, 2*s,N,X,ldx);
  else
    FFLAS::fzero(Fi, N, 2*s, X,ldx);
  size_t CurrentBlockPos = 0;
  size_t NbBlocks = 0;
  size_t BlockPivot=0;
  while (BlockPivot<r)
    { K[NbBlocks] = CurrentBlockPos;
      size_t BlockSize, NextBlockPos;
      if (BlockPivot+s >= r)
	{
	  BlockSize =r-BlockPivot;
	  NextBlockPos = N;
        } else{
	  BlockSize = s;
	  NextBlockPos = outer[Inv[M[BlockPivot+BlockSize]]];}
	  
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
    applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasTrans, N, size_t(0), N, A, lda, MLap);
  }
 else{
  applyP (Fi, FFLAS::FflasRight, FFLAS::FflasNoTrans, N, size_t(0), N, A, lda, MLap);
 }

}
 // Compute M such that LM is in column echelon form, where L is the
      // left factor of a Bruhat decomposition
      // M is allocated in this function and should be deleted after using it.
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
          M[i] = piv;
          ispivot [piv]=true;
      }
      size_t curr=R;
      for (size_t i=0; i<N; ++i)
          if (!ispivot[i])
              M [curr++] = i;

      FFLAS::fflas_delete(Pinv,Ps);
   }
template<class Field>
 inline  void productBruhatxTS (const Field&Fi, size_t N, size_t s, size_t r, const size_t *P, const size_t *Q,  const typename Field::Element_ptr Xu,size_t ldu, size_t NbBlocksU, size_t * Ku, size_t *Tu ,const typename Field::Element_ptr Xl, size_t ldl, size_t NbBlocksL,size_t *Kl, size_t *Tl,typename  Field::Element_ptr B,size_t t, size_t ldb,typename Field::Element_ptr C, size_t ldc)
    {
      
    size_t * Tuinv = FFLAS::fflas_new<size_t>(r);
    std::vector<bool> IsFinalU(r,true);
    for(size_t i=0;i<r;i++)
      {
	if(Tu[i]==i && IsFinalU[i])
	  {
	    Tuinv[i]=i;
	  }
	else if(Tu[i]!=i &&IsFinalU[i])
	  {
	    Tuinv[Tu[i]]=i;
	    IsFinalU[Tu[i]]=false;
	  }
	else
	  {
	    Tuinv[Tu[i]]=Tuinv[i];
	    IsFinalU[Tu[i]]=false;
	  }
      }
    size_t * Tlinv = FFLAS::fflas_new<size_t>(r);
    std::vector<bool> IsFinalL(r,true);
    for(size_t i=0;i<r;i++)
      {
	if(Tl[i]==i && IsFinalL[i])
	  {
	    Tlinv[i]=i;
	  }
	else if(Tl[i]!=i &&IsFinalL[i])
	  {
	    Tlinv[Tl[i]]=i;
	    IsFinalL[Tl[i]]=false;
	  }
	else
	  {
	    Tuinv[Tl[i]]=Tlinv[i];
	    IsFinalL[Tl[i]]=false;
	  }
      }
  
    size_t S = 0;
    size_t k = 0;
    while( S+s<N)
      {
	S+=s;
	k++;
      }
    size_t rs = N-S;
    k++;
    size_t blocksu = 0;
    size_t blocksl=NbBlocksL;
    if(Kl[NbBlocksL-1]>S)
      blocksl-=1;
    size_t grid_sizeU;
    size_t grid_sizeL;
    typename Field::Element_ptr Sj= FFLAS::fflas_new(Fi, 2*s,t);
    typename Field::Element_ptr Xj = FFLAS::fflas_new(Fi,r,t);
    typename Field::Element_ptr Yj = FFLAS::fflas_new(Fi, r, t);
    typename Field::Element_ptr Z = FFLAS::fflas_new(Fi, r, t);
    typename Field::Element_ptr TlZ = FFLAS::fflas_new(Fi, r, t);
    typename Field::Element_ptr trailing_term = FFLAS::fflas_new(Fi, s, t);
    typename Field::Element_ptr E = FFLAS::fflas_new(Fi, r, s);
    typename Field::Element_ptr Er = FFLAS::fflas_new(Fi, r, s);
    typename Field::Element_ptr TlEr = FFLAS::fflas_new(Fi, r,s);
    typename Field::Element_ptr CRE = FFLAS::fflas_new(Fi, s,s);
    typename Field::Element_ptr SCRE = FFLAS::fflas_new(Fi, s, s);
    FFLAS::fzero(Fi, r,t, Z,t);
    for (size_t i=0;i<k-1;i++)
      {     FFLAS::fzero(Fi, 2*s,t,Sj,t);
	    FFLAS::fzero(Fi,r,t,Xj,t);
	    FFLAS::fzero(Fi,r,t,Yj,t);
	    FFLAS::fzero(Fi,r,t,TlZ,t);
	    if (Ku[blocksu+1]-i*s<s)
	      { //Dj*Bj
		if (blocksu+1<NbBlocksU)
		  grid_sizeU =s;
		else
		  grid_sizeU = r-(blocksu)*s;
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,s, t, Ku[blocksu+1]-i*s, Fi.one, Xu+i*s, ldu, B+i*s*ldb, ldb, Fi.zero, Xj+blocksu*s*t, t);
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,grid_sizeU, t, i*s+s-Ku[blocksu+1], Fi.one, Xu+Ku[blocksu+1], ldu, B+Ku[blocksu+1]*ldb, ldb, Fi.zero, Xj+(blocksu+1)*s*t, t);
		//Sj*Bj
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,s, t, Ku[blocksu+1]-i*s, Fi.one, Xu+i*s+s*ldu, ldu, B+i*s*ldb, ldb, Fi.zero, Sj, t);
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,s, t, i*s+s-Ku[blocksu+1], Fi.one, Xu+Ku[blocksu+1]+s*ldu, ldu, B+Ku[blocksu+1]*ldb, ldb, Fi.zero, Sj+s*t, t);
		//Apply Tu
		if(blocksu>0)
		  {
		    for(size_t l=0; l<s; l++)
		      FFLAS::faddin(Fi,t,Sj+l*t,1,Xj+Tuinv[(blocksu-1)*s+l]*t,1);
		  }
		for (size_t l=0; l<s; l++)
		  if(blocksu*s+l<r)
		    FFLAS::faddin(Fi,t,Sj+l*t,1,Xj+Tuinv[blocksu*s+l]*t,1);
	      }
	    else
	      { // Dj*Bj
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,s, t, s, Fi.one, Xu+i*s, ldu, B+i*s*ldb, ldb, Fi.zero, Xj+blocksu*s*t, t);
		// Sj*Bj
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,s, t, s, Fi.one, Xu+i*s+s*ldu, ldu, B+i*s*ldb, ldb, Fi.zero, Sj, t);
		//Apply Tu
		if(blocksu>0)
		  {
		    for(size_t l=0; l<s; l++)
		      FFLAS::faddin(Fi,t,Sj+l*t,1,Xj+Tuinv[(blocksu-1)*s+l]*t,1);
		  }
	      }
	    // Apply R to Xj
	    for(size_t j=0; j<r;j++)
	      {
		FFLAS::fassign(Fi,t, Xj+Q[j]*t,1, Yj+P[j]*t, 1);
	      }
	    //Zj
	    FFLAS::faddin(Fi,r,t,Yj,t,Z,t);
	    //Compute Tl*Zj
	    for(size_t l=0;l<r;l++)
	      {
		FFLAS::faddin(Fi,t,Z+l*t,1,TlZ+Tlinv[l]*t,1);
	      }
	    if (blocksl<NbBlocksL)
	      grid_sizeL= s;
	    else
	      grid_sizeL = r-(blocksl-1)*s;
	      
	    if (S-i*s-Kl[blocksl-1]<s)
	      { //Dj*Zj
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,Kl[blocksl-1]-(S-(i+1)*s), t, s, Fi.one, Xl+(S-(i+1)*s)*ldl, ldl, Z+(blocksl-2)*s*t, t, Fi.zero, C+(S-(i+1)*s)*ldc, ldc);
	        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, S-i*s-Kl[blocksl-1],t, grid_sizeL, Fi.one, Xl+Kl[blocksl-1]*ldl, ldl, Z+(blocksl-1)*s*t, t, Fi.zero, C+Kl[blocksl-1]*ldc, ldc);
		
		//Sj*Zj
	        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,Kl[blocksl-1]-(S-(i+1)*s), t, s, Fi.one, Xl+(S-(i+1)*s)*ldl+s, ldl, TlZ+(blocksl-1)*s*t, t, Fi.zero, Sj, t);
	        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,S-i*s-Kl[blocksl-1], t, s, Fi.one, Xl+Kl[blocksl-1]*ldl+s, ldl, TlZ+blocksl*s*t, t, Fi.zero, Sj+Kl[blocksl-1]-(S-(i+1)*s), t);	
	      }
	    else
	      {
		 //Dj*Bj
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, grid_sizeL, Fi.one, Xl+(S-(i+1)*s)*ldl, ldl, Z+(blocksl-1)*s*t, t, Fi.zero, C+(S-(i+1)*s)*ldc, ldc);
		
		//Sj*Bj
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, s, Fi.one, Xl+S*ldl+s, ldl, TlZ+blocksl*s*t, t, Fi.zero, Sj, t);
	      }
	    FFLAS::faddin(Fi,s,t,Sj,t,C+(N-S-s)*ldc,ldc);
	    //Next Block
	    if (Ku[blocksu+1]-S<s)
	      blocksu++;
	    if (N-S-Kl[blocksl-1]<s)
	      blocksl--;
	}
	blocksu = 0;
	blocksl = NbBlocksL;  
	size_t row_pos = N;
	size_t grid_sizerow;
	size_t grid_sizecol;
	for(size_t i=0; i<k;i++)
	  { if (i<k-1)
	       {
		 grid_sizecol = s;
	       }
	    else{
	      grid_sizecol = N-S;
	    }
	    
	    //Compute (Left(CiREN-i+1)BN-i+1) the trailing term
	       FFLAS::fzero(Fi, s, t, trailing_term,t);
	       FFLAS::fzero(Fi, r, s,E,s);
	       FFLAS::fzero(Fi, r, s, Er, s);
	       FFLAS::fzero(Fi, r,s, TlEr, s);
	       FFLAS::fzero(Fi, s,s, CRE,s);
	       FFLAS::fzero(Fi, s, s, SCRE,s);
	    if(Ku[blocksu+1]-i*s<grid_sizecol)
	      {
		if(blocksu>0)
		  {for (size_t l=0; l<s; l++)
		  {
		    FFLAS::fassign(Fi,s,Ku[blocksu+1]-i*s,Xu+i*s+s*ldu,ldu, E+Tuinv[(blocksu-1)*s+l]*s,s);
		    FFLAS::fassign(Fi, s, i*s+grid_sizecol - Ku[blocksu+1], Xu+Ku[blocksu+1]+s*ldu,ldu, E+Tuinv[blocksu*s+l]*s+Ku[blocksu+1]-i*s,s);
		  }}
		FFLAS::fassign(Fi,s,Ku[blocksu+1]-i*s, Xu+i*s,ldu, E+i*s*s,s);
		FFLAS::fassign(Fi, s, i*s+grid_sizecol -Ku[blocksu+1], Xu+Ku[blocksu+1],ldu, E+i*s*s+Ku[blocksu+1]-i*s,s);
	      }
	    else
	      { if(blocksu>0)
		  {for (size_t l=0; l<s; l++)
		  { 
		    FFLAS::fassign(Fi,s,grid_sizecol,Xu+i*s+s*ldu,ldu, E+Tuinv[(blocksu-1)*s+l]*s,s);
		 
		  }}
		FFLAS::fassign(Fi,s,grid_sizecol,Xu+i*s,ldu, E+i*s*s,s);
	      }
	    // Apply R to E
	    for(size_t j=0; j<r;j++)
	      {
		FFLAS::fassign(Fi,grid_sizecol, E+Q[j]*s,1, Er+P[j]*s, 1);
	      }
	    //Compute Left(CRE)
	    for (size_t j=0; j<r;j++)
	      {
		FFLAS::faddin(Fi,grid_sizecol,Er+j*s,1,TlEr+Tlinv[j]*s,1);
	      }
	    if (i<k)
	       {
		 grid_sizerow = s;
	       }
	    else{
	         grid_sizerow = N-S;
	    }
	    if (row_pos-Kl[blocksl-1]<grid_sizerow)
	      {if(blocksl>1)
		//DN-j
		  fgemm (Fi,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans, Kl[blocksl-1]-(row_pos-grid_sizerow),s,s, Fi.one, Xl+(row_pos-grid_sizerow)*ldl,ldl,Er+(blocksl-2)*s *s,s,Fi.zero, CRE,s);
		fgemm (Fi,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans, row_pos-Kl[blocksl-1],s,s, Fi.one, Xl+Kl[blocksl-1]*ldl,ldl,Er+(blocksl-1)*s *s,s, Fi.zero,CRE+(Kl[blocksl-1]-(row_pos-grid_sizerow))*s,s);
		//SN-j
		if(blocksl>2)
		  fgemm (Fi,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans, Kl[blocksl-1]-(row_pos-grid_sizerow),s,s, Fi.one, Xl+(row_pos-grid_sizerow)*ldl+s,ldl,TlEr+(blocksl-3)*s *s,s, Fi.zero,SCRE,s);
		if(blocksl>1)
		  fgemm (Fi,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans, row_pos-Kl[blocksl-1],s,s, Fi.one, Xl+Kl[blocksl-1]*ldl+s,ldl,TlEr+(blocksl-2)*s *s,s, Fi.zero,SCRE+(Kl[blocksl-1]-(row_pos-grid_sizerow))*s,s);
	      }
	    else
	      { //DN-j
		fgemm (Fi,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans, grid_sizerow,s,s, Fi.one, Xl+(row_pos-grid_sizerow)*ldl,ldl,Er+(blocksl-1)*s *s,s, Fi.zero,CRE,s);
		//SN-j
		fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,grid_sizerow,s,s, Fi.one, Xl+(N-S-s)*ldl+s,ldl,TlEr+(blocksl-2)*s*s,s, Fi.zero,SCRE,s);
	      }
	    FFLAS::faddin(Fi, grid_sizerow, grid_sizecol, SCRE, s, CRE, s);
	    for (size_t i=0; i<s; i++)
	      FFLAS::fzero(Fi, i+1, CRE + i*s + s-i-1, 1);
	    fgemm(Fi,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans, grid_sizerow, t, grid_sizecol, Fi.one, CRE, s, B+i*s*t, t, Fi.zero,trailing_term, t);
	    FFLAS::faddin(Fi,grid_sizerow,t,trailing_term,t,C+(row_pos-grid_sizerow)*ldc,ldc);
	    
	    //Next Block
	    row_pos -= grid_sizerow;
	    if (Ku[blocksu+1]-S<s)
	      blocksu++;
	    if (N-S-Kl[blocksl+1]<s)
	      blocksl--;
        }
  
	
    }
  
} //namespace FFPACK
#endif //_FFPACK_ffpack_bruhatgen_inl
