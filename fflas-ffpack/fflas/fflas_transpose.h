/*
 * Copyright (C) 2020 the FFLAS-FFPACK group
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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
/** @file fflas/fflas_transpose.h
 * @brief transpose the storage of the matrix (switch between row and col major mode)
 */


#ifndef __FFLASFFPACK_transpose_H
#define __FFLASFFPACK_transpose_H
#include "fflas-ffpack/utils/debug.h"
#include "fflas-ffpack/fflas/fflas.h"

#ifndef FFLAS_TRANSPOSE_BLOCKSIZE 
#define FFLAS_TRANSPOSE_BLOCKSIZE 32
#endif

#include "fflas-ffpack/fflas/fflas_simd.h"

namespace FFLAS {
  

  template <typename Field, typename Simd= Simd<typename Field::Element>, typename Enable=void>
  class BlockTransposeSIMD;

  template <typename Field>
  class BlockTransposeSIMD<Field,NoSimd<typename Field::Element>, void> {
  public:
    using Element          = typename Field::Element;
    using Element_ptr      = typename Field::Element_ptr;
    using ConstElement_ptr = typename Field::ConstElement_ptr;
    using Simd = NoSimd<typename Field::Element>;
       
    void transposein(const Field&F, Element_ptr A, size_t lda){transpose(F, A,lda,A,lda);}
    
    void transpose(const Field&F, ConstElement_ptr A, size_t lda, Element_ptr B, size_t ldb){*B=*A;}

    const constexpr size_t size() const { return 1;}
        
    void info() const {std::cerr<<"\n IN DEVELOPMENT: transpose with block SIMD: "<<typename NoSimd<typename Field::element>::type_string()<<" with vect_size="<<size()<<std::endl;}    
};

  template <typename Field, typename Simd, size_t>
  struct transpose_simd {
    void operator()(const Field&F, typename Field::ConstElement_ptr A, size_t lda, typename Field::Element_ptr B, size_t ldb);
  };
  

#define LD(i) R##i=Simd::loadu(A+lda*i)
#define ST(i) Simd::storeu(B+ldb*i,R##i)
#define PCK(i,j) Simd::unpacklohi(R##i,R##j,R##i,R##j);
  
  template <typename Field, typename Simd>
  struct transpose_simd<Field, Simd, 4> {
    inline void operator()(const Field&F, typename Field::ConstElement_ptr A, size_t lda, typename Field::Element_ptr B, size_t ldb){
      typename Simd::vect_t R0,R1,R2,R3;
      LD(0);LD(1);LD(2);LD(3);
      PCK(0,2); PCK(1,3);
      PCK(0,1); PCK(2,3);
      ST(0);ST(1);ST(2);ST(3);
    }
  };
  template <typename Field, typename Simd>
  struct transpose_simd<Field, Simd, 8> {
    inline void operator() (const Field&F, typename Field::ConstElement_ptr A, size_t lda, typename Field::Element_ptr B, size_t ldb){
      typename Simd::vect_t R0,R1,R2,R3,R4,R5,R6,R7; 
      LD(0);LD(1);LD(2);LD(3);LD(4);LD(5);LD(6);LD(7);
      PCK(0,4); PCK(1,5); PCK(2,6); PCK(3,7);
      PCK(0,2); PCK(1,3); PCK(4,6); PCK(5,7);
      PCK(0,1); PCK(2,3); PCK(4,5); PCK(6,7);
      ST(0);ST(1);ST(2);ST(3);ST(4);ST(5);ST(6);ST(7);
    }
  };
  template <typename Field, typename Simd>
  struct transpose_simd<Field, Simd, 16> {
    inline void operator()(const Field&F, typename Field::ConstElement_ptr A, size_t lda, typename Field::Element_ptr B, size_t ldb){
      typename Simd::vect_t R0,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15;
      LD(0);LD(1);LD(2);LD(3);LD(4);LD(5);LD(6);LD(7);LD(8);LD(9);LD(10);LD(11);LD(12);LD(13);LD(14);LD(15);
      PCK(0,8); PCK(1,9); PCK(2,10); PCK(3,11); PCK(4,12); PCK(5,13);  PCK(6,14);  PCK(7,15);
      PCK(0,4); PCK(1,5); PCK(2,6);  PCK(3,7);  PCK(8,12); PCK(9,13);  PCK(10,14); PCK(11,15);
      PCK(0,2); PCK(1,3); PCK(4,6);  PCK(5,7);  PCK(8,10); PCK(9,11);  PCK(12,14); PCK(13,15);
      PCK(0,1); PCK(2,3); PCK(4,5);  PCK(6,7);  PCK(8,9);  PCK(10,11); PCK(12,13); PCK(14,15);  
      ST(0);ST(1);ST(2);ST(3);ST(4);ST(5);ST(6);ST(7);ST(8);ST(9);ST(10);ST(11);ST(12);ST(13);ST(14);ST(15);
    }
  };
  
    template <typename Field, typename Simd>
  struct transpose_simd<Field, Simd, 32> {
    inline void operator()(const Field&F, typename Field::ConstElement_ptr A, size_t lda, typename Field::Element_ptr B, size_t ldb){
      typename Simd::vect_t R0,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15,R16,R17,R18,R19,R20,R21,R22,R23,R24,R25,R26,R27,R28,R29,R30,R31;
      LD(0);LD(1);LD(2);LD(3);LD(4);LD(5);LD(6);LD(7);LD(8);LD(9);LD(10);LD(11);LD(12);LD(13);LD(14);LD(15);
      LD(16);LD(17);LD(18);LD(19);LD(20);LD(21);LD(22);LD(23);LD(24);LD(25);LD(26);LD(27);LD(28);LD(29);LD(30);LD(31);
      PCK(0,16); PCK(1,17); PCK(2,18); PCK(3,19); PCK(4,20); PCK(5,21);  PCK(6,22);  PCK(7,23);  PCK(8,24);  PCK(9,25);  PCK(10,26); PCK(11,27); PCK(12,28); PCK(13,29); PCK(14,30); PCK(15,31);
      PCK(0,8);  PCK(1,9);  PCK(2,10); PCK(3,11); PCK(4,12); PCK(5,13);  PCK(6,14);  PCK(7,15);  PCK(16,24); PCK(17,25); PCK(18,26); PCK(19,27); PCK(20,28); PCK(21,29); PCK(22,30); PCK(23,31);
      PCK(0,4);  PCK(1,5);  PCK(2,6);  PCK(3,7);  PCK(8,12); PCK(9,13);  PCK(10,14); PCK(11,15); PCK(16,20); PCK(17,21); PCK(18,22); PCK(19,23); PCK(24,28); PCK(25,29); PCK(26,30); PCK(27,31);
      PCK(0,2);  PCK(1,3);  PCK(4,6);  PCK(5,7);  PCK(8,10); PCK(9,11);  PCK(12,14); PCK(13,15); PCK(16,18); PCK(17,19); PCK(20,22); PCK(21,23); PCK(24,26); PCK(25,27); PCK(28,30); PCK(29,31);
      PCK(0,1);  PCK(2,3);  PCK(4,5);  PCK(6,7);  PCK(8,9);  PCK(10,11); PCK(12,13); PCK(14,15); PCK(16,17); PCK(18,19); PCK(20,21); PCK(22,23); PCK(24,25); PCK(26,27); PCK(28,29); PCK(30,31);
      ST(0);ST(1);ST(2);ST(3);ST(4);ST(5);ST(6);ST(7);ST(8);ST(9);ST(10);ST(11);ST(12);ST(13);ST(14);ST(15);
      ST(16);ST(17);ST(18);ST(19);ST(20);ST(21);ST(22);ST(23);ST(24);ST(25);ST(26);ST(27);ST(28);ST(29);ST(30);ST(31);
    }
    };
  
#undef LD
#undef ST
#undef PCK  

  template <typename Field, typename Simd>
  class BlockTransposeSIMD<Field,Simd,
			   typename std::enable_if<FFLAS::support_simd<typename Field::Element>::value>::type>  {
  public:
    using Element          = typename Field::Element;
    using Element_ptr      = typename Field::Element_ptr;
    using ConstElement_ptr = typename Field::ConstElement_ptr;
    using vect_t           = typename Simd::vect_t;
    
    void transposein(const Field&F, Element_ptr A, size_t lda){transpose(F, A,lda,A,lda);}

    void transpose(const Field&F, ConstElement_ptr A, size_t lda, Element_ptr B, size_t ldb)  {transpose_simd<Field, Simd, Simd::vect_size>() (F,A,lda,B,ldb);}   
    
    void transposeLoop(const Field&F, ConstElement_ptr A, size_t lda, Element_ptr B, size_t ldb){
      vect_t R[Simd::vect_size];
      for (size_t i=0;i<Simd::vect_size;i++)
    	R[i] = Simd::loadu(A+i*lda);
      size_t w=Simd::vect_size>>1;
      size_t f=1;
      size_t idx0,idx1;
      for (;w>0; w>>=1, f<<=1){
    	idx0=0;idx1=w;
    	for (size_t i = 0; i < f; i++, idx0+=w,idx1+=w)
    	  for (size_t j = 0; j < w; j++, idx0++, idx1++)
    	    {
    	      Simd::unpacklohi(R[idx0],R[idx1],R[idx0],R[idx1]);
    	    }
      }      
      for (size_t i=0;i<Simd::vect_size;i++)
    	Simd::storeu(B+i*ldb, R[i]);
    }

     const constexpr size_t size() const { return Simd::vect_size;}
    
    void info() const {std::cerr<<"\n IN DEVELOPMENT: transpose with block SIMD: "<<Simd::type_string()<<" with vect_size="<<size()<<std::endl;}
  };

 

 


  template<class Field, size_t BLOCK>
  inline  typename Field::Element_ptr
  ftranspose_impl (const Field& F, const size_t m, const size_t n,
		   typename Field::ConstElement_ptr A, const size_t lda,
		   typename Field::Element_ptr B, const size_t ldb)
    {
      // rk: m<=ldb
      const size_t ls = BLOCK;
      typename Field::Element tmp; F.init(tmp);
      for (size_t i = 0; i < m; i+=ls){      
	for (size_t j =0; j < n; j+=ls)
	  for (size_t _i = i; _i < std::min(m, i+ls); _i++)
	    for (size_t _j = j; _j < std::min(n, j+ls); _j++){
	      *(B+_j*ldb+_i)=*(A+_i*lda+_j);
	    }
      }
      return B;
    }


  template<class Field, size_t BLOCK>
  inline  typename Field::Element_ptr
  ftranspose_impl_simd (const Field& F, const size_t m, const size_t n,
			typename Field::ConstElement_ptr A, const size_t lda,
			typename Field::Element_ptr B, const size_t ldb)
    {
      BlockTransposeSIMD<Field> BTS; //BTS.info();

      if (m%BTS.size() ||  n %BTS.size() || BLOCK % BTS.size())
	return ftranspose_impl<Field,BLOCK>(F,m,n,A,lda,B,ldb);
      
      // rk: m<=ldb
      const size_t ls = BLOCK;
      typename Field::Element tmp; F.init(tmp);
      for (size_t i = 0; i < m; i+=ls){      
	for (size_t j =0; j < n; j+=ls)
	  for (size_t _i = i; _i < std::min(m, i+ls); _i+=BTS.size())
	    for (size_t _j = j; _j < std::min(n, j+ls); _j+=BTS.size()){
	      BTS.transpose(F, A+_i*lda+_j, lda, B+_j*ldb+_i, ldb);
	    }
      }
      return B;
    }



  template<class Field, size_t BLOCK>
  inline  typename Field::Element_ptr
  ftransposein_impl (const Field& F, const size_t m, const size_t n,
		     typename Field::Element_ptr A, const size_t lda)
  {    
    // rk: m<=lda
    const size_t ls = BLOCK;
    typename Field::Element tmp; F.init(tmp);
    for (size_t i = 0; i < m; i+=ls){
      // these two loops are for diagonal blocks [i..i+ls,i..i+ls]
      for (size_t _i = i; _i < std::min(m, i+ls); _i++)
	for (size_t _j = _i+1; _j < std::min(n, i+ls); _j++){	
	  tmp= *(A+_i*lda+_j);
	  *(A+_i*lda+_j)=*(A+_j*lda+_i);
	  *(A+_j*lda+_i)=tmp;
	}
      // this loops is for off diagonal blocks
      for (size_t j =i+ls; j < n; j+=ls)
	// these two loops are for off diagonal blocks [i..i+ls,j..i+ls] and [j..j+ls,i..i+ls]
	// it might be usefull to copy these two ls x ls blocks into contiugous memory
	// -> this depends upon cache policy and mapping
	for (size_t _i = i; _i < std::min(m, i+ls); _i++)
	  for (size_t _j = j; _j < std::min(n, j+ls); _j++){
	    tmp= *(A+_i*lda+_j);
	    *(A+_i*lda+_j)=*(A+_j*lda+_i);
	    *(A+_j*lda+_i)=tmp;
	  }
    }
    return A;
  }
  template<class Field, size_t BLOCK>
  inline  typename Field::Element_ptr
  ftransposein_impl_simd (const Field& F, const size_t m, const size_t n,
			  typename Field::Element_ptr A, const size_t lda)
  {
    // rk: m<=lda
    BlockTransposeSIMD<Field> BTS; //BTS.info();

    if (m%BTS.size() ||  n %BTS.size() || BLOCK % BTS.size())
      return ftransposein_impl<Field,BLOCK>(F,m,n,A,lda);
      
    const size_t ls = BLOCK;
    typename Field::Element TMP[ls*ls];
    finit(F,ls,ls,TMP,ls);    
    for (size_t i = 0; i < m; i+=ls){
      // these two loops are for diagonal blocks [i..i+ls,i..i+ls]
      fassign(F, ls, ls, A+i*lda+i, lda, TMP, ls);       
      for (size_t _i = i; _i < std::min(m, i+ls); _i+=BTS.size()){
	BTS.transpose(F, A+_i*lda+_i, lda, A+_i*lda+_i, lda);
	for (size_t _j = _i+BTS.size(); _j < std::min(n, i+ls); _j+=BTS.size()){
	  BTS.transpose(F, A+_j*lda+_i, lda, A+_i*lda+_j, lda);
	  BTS.transpose(F, TMP+(_i-i)*ls+(_j-i), ls, A+_j*lda+_i, lda);
	}
      }      
      for (size_t j =i+ls; j < n; j+=ls){
	// these two loops are for off diagonal blocks [i..i+ls,j..i+ls] and [j..j+ls,i..i+ls]	
	// copy only the first block
	fassign(F, ls, ls, A+i*lda+j, lda, TMP, ls);       
	for (size_t _i = i; _i < std::min(m, i+ls); _i+=BTS.size())
	  for (size_t _j = j; _j < std::min(n, j+ls); _j+=BTS.size()){
	    BTS.transpose(F, A+_j*lda+_i, lda, A+_i*lda+_j, lda);
	    BTS.transpose(F, TMP+(_i-i)*ls+(_j-j), ls, A+_j*lda+_i, lda);
	  }
      }
    }
    return A;
  }
  
  template<class Field, size_t BLOCK>
  inline  typename Field::Element_ptr
  ftransposein_impl_simd_v2 (const Field& F, const size_t m, const size_t n,
			  typename Field::Element_ptr A, const size_t lda)
  {
    // rk: m<=lda
    BlockTransposeSIMD<Field> BTS; //BTS.info();

    if (m%BTS.size() ||  n %BTS.size() || BLOCK % BTS.size())
      return ftransposein_impl<Field,BLOCK>(F,m,n,A,lda);
      
    const size_t ls = BLOCK;
    typename Field::Element TMP1[ls*ls], TMP2[ls*ls];
    finit(F,ls,ls,TMP1,ls);
    finit(F,ls,ls,TMP2,ls);

    // This variant does not separate diagonal block ffrom the other ones.
    // -> since each block are copied in TMP1 and TMP2, we put back their transposed at the right position in the result
    // rk : diagonal blocks are transposed and written twice
    for (size_t i = 0; i < m; i+=ls){          
      for (size_t j =i; j < n; j+=ls){
	// these two loops are for off diagonal blocks [i..i+ls,j..i+ls] and [j..j+ls,i..i+ls]	
	// copy the two blocks
	fassign(F, ls, ls, A+i*lda+j, lda, TMP1, ls);
	fassign(F, ls, ls, A+j*lda+i, lda, TMP2, ls);       
	for (size_t _i = i; _i < std::min(m, i+ls); _i+=BTS.size())
	  for (size_t _j = j; _j < std::min(n, j+ls); _j+=BTS.size()){
	    BTS.transpose(F, TMP2+(_j-j)*ls+(_i-i), ls, A+_i*lda+_j, lda);
	    BTS.transpose(F, TMP1+(_i-i)*ls+(_j-j), ls, A+_j*lda+_i, lda);
	  }
      }
    }
    return A;
  }

  template<class Field>
  inline  typename Field::Element_ptr
  ftransposein (const Field& F, const size_t m, const size_t n,
		typename Field::Element_ptr A, const size_t lda)
  {
    FFLASFFPACK_check(m<=lda);
    //return ftransposein_impl_simd_v2<Field,4> (F,m,n,A,lda);
    return ftransposein_impl_simd<Field,FFLAS_TRANSPOSE_BLOCKSIZE> (F,m,n,A,lda);
  }

  template<class Field>
  inline  typename Field::Element_ptr
  ftranspose (const Field& F, const size_t m, const size_t n,
	      typename Field::ConstElement_ptr A, const size_t lda,
	      typename Field::Element_ptr B, const size_t ldb)
  {
    FFLASFFPACK_check(m<=ldb);
    return ftranspose_impl_simd<Field,FFLAS_TRANSPOSE_BLOCKSIZE> (F,m,n,A,lda,B,ldb);
  }


} // end of FFLAS namespace




#endif // __FFLASFFPACK_transpoose_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
