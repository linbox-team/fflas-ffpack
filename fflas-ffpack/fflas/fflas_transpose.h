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
  
  template <typename Field, typename Simd>
  class BlockTransposeSIMD<Field,Simd,
			   typename std::enable_if<FFLAS::support_simd<typename Field::Element>::value>::type>  {
  public:
    using Element          = typename Field::Element;
    using Element_ptr      = typename Field::Element_ptr;
    using ConstElement_ptr = typename Field::ConstElement_ptr;
    using vect_t           = typename Simd::vect_t;
    using unpacklohi       = typename Simd::unpacklohi;
    using loadu            = typename Simd::loadu;
    using storeu           = typename Simd::storeu;

    
    void transposein(Element_ptr A, size_t lda){
      transpose(A,lda,A,lda);
    }

    void transpose(ConstElement_ptr A, size_t lda, Element_ptr B, size_t ldb){
      vect_t R[Simd::vect_size];
      for (size_t i=0;i<Simd::vect_size;i++)
	loadu(R[i], A+i*lda);
      
      size_t w=Simd::vect_size>>1;
      size_t f=1;
      for (;w>0; w>>=1, f<<1)
	for (size_t i = 0; i < f; i++)
	  for (size_t j = 0; j < w; j)
	    {
	      unpacklohi(R[j+i*w],R[j+(i+1)*w],R[j+i*w],R[j+(i+1)*w]);
	    }
      
      for (size_t i=0;i<Simd::vect_size;i++)
	storeu(R[i], B+i*ldb);
    }

     const constexpr size_t size() const { return Simd::vect_size;}
    
    
};


  

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
      std::cerr<<"DEVELOPMENT: transpose with block SIMD\n";
      BlockTransposeSIMD<Field> BTS;
      if (m%BTS.size() ||  n %BTS.size() || BLOCK % BTS.size())
	return ftranspose_impl<Field,BLOCK>(F,m,n,A,lda,B,ldb);
      
      // rk: m<=ldb
      const size_t ls = BLOCK;
      typename Field::Element tmp; F.init(tmp);
      for (size_t i = 0; i < m; i+=ls){      
	for (size_t j =0; j < n; j+=ls)
	  for (size_t _i = i; _i < i+ls; _i+=BTS.size())
	    for (size_t _j = j; _j < j+ls; _j+=BTS.size()){
	      BTS.transpose(A+_i*lda+_j, lda, B+_j*ldb+_i, ldb);
	    }
      }
      return B;
    }


  
  
  template<class Field>
  inline  typename Field::Element_ptr
  ftransposein (const Field& F, const size_t m, const size_t n,
		typename Field::Element_ptr A, const size_t lda)
  {
    FFLASFFPACK_check(m<=lda);
    return ftransposein_impl<Field,FFLAS_TRANSPOSE_BLOCKSIZE> (F,m,n,A,lda);
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
