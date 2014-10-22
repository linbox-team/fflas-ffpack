/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_pfgemm.inl
 * Copyright (C) 2013 Jean Guillaume Dumas Clement Pernet Ziad Sultan
 *
 * Written by Jean Guillaume Dumas Clement Pernet Ziad Sultan
 * Time-stamp: <01 Oct 14 09:56:36 Jean-Guillaume.Dumas@imag.fr>
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

#ifndef __FFLASFFPACK_fflas_parfgemm_INL
#define __FFLASFFPACK_fflas_parfgemm_INL

#define __FFLASFFPACK_SEQPARTHRESHOLD 220
#define __FFLASFFPACK_DIMKPENALTY 1

#ifdef __FFLASFFPACK_USE_KAAPI
#include <kaapi++>
#include "fflas-ffpack/fflas/kaapi_routines.inl"
#endif
#ifdef __FFLASFFPACK_USE_OPENMP
#include <omp.h>
#endif

#include "fflas_blockcuts.inl"
#include "parallel.h"
#include "fflas-ffpack/utils/timer.h"

#ifdef __FFLASFFPACK_USE_OPENMP4
namespace FFLAS {

	template<class Field, class AlgoT, class FieldTrait>
	 typename Field::Element*
	 pfgemm_3D_rec_adapt( const Field& F,
			const FFLAS_TRANSPOSE ta,
			const FFLAS_TRANSPOSE tb,
			const size_t m,
			const size_t n,
			const size_t k,
			const typename Field::Element alpha,
			const typename Field::Element_ptr A, const size_t lda,
			const typename Field::Element_ptr B, const size_t ldb,
			const typename Field::Element beta,
			typename Field::Element * C, const size_t ldc, 
			MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel> & H){

	 typename Field::Element a = alpha;
	 typename Field::Element b = beta;

	 if (!m || !n) {return C;}
	 if (!k || F.isZero (alpha)){
		 fscalin(F, m, n, beta, C, ldc);
		 return C;
	 }
	 // threshold                                                                                                  

	 //      if(m==1 && n==1 && k==1)                                                                              
	 if (H.parseq.numthreads<=1 || std::min(m*n,std::min(m*k,k*n))<=__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD){
		 MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Sequential> SeqH(H);
		 return fgemm(F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, SeqH);
	 }
	  MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H1(H);
	  MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H2(H);
	  if(__FFLASFFPACK_DIMKPENALTY*m > k && m >= n) {
		 size_t M2= m>>1;
		 H1.parseq.numthreads /=2;
		 H2.parseq.numthreads = H.parseq.numthreads - H1.parseq.numthreads;
		 
		 typename Field::Element_ptr A1= A;
		 typename Field::Element_ptr A2= A+M2*lda;
		 typename Field::Element_ptr C1= C;
		 typename Field::Element_ptr C2= C+M2*ldc;
		 
		     // 2 multiply (1 split on dimension m)
		
#pragma omp task shared(F, A1, C1) depend(in:A1,B) depend(inout:C1)
		 pfgemm_3D_rec_adapt(F, ta, tb, M2, n, k, alpha, A1, lda, B, ldb, beta, C1, ldc, H1);
#pragma omp task shared(F, A2, C2) depend(inout:C2) depend(in:A2,B)
		 pfgemm_3D_rec_adapt(F, ta, tb, m-M2, n, k, alpha, A2, lda, B, ldb, beta, C2, ldc, H2);
#pragma omp taskwait
	 } else if (__FFLASFFPACK_DIMKPENALTY*n > k) {
		 size_t N2 = n>>1;
		 H1.parseq.numthreads /=2;
		 H2.parseq.numthreads = H.parseq.numthreads - H1.parseq.numthreads;
		 typename Field::Element_ptr B1= B;
		 typename Field::Element_ptr B2= B+N2;
		 
		 typename Field::Element_ptr C1= C;
		 typename Field::Element_ptr C2= C+N2;
		 
#pragma omp task shared(F, B1, C1) depend(in:A,B1) depend(inout:C1)
		 pfgemm_3D_rec_adapt(F, ta, tb, m, N2, k, a, A, lda, B1, ldb, b, C1, ldc, H1);
#pragma omp task shared(F, B2, C2) depend(in:A,B2) depend(inout:C2)
		 pfgemm_3D_rec_adapt(F, ta, tb, m, n-N2, k, a, A, lda, B2, ldb, b,C2, ldc, H2);
#pragma omp taskwait
	 } else {
		 size_t K2 = k>>1;
		 
		 typename Field::Element_ptr B1= B;
		 typename Field::Element_ptr B2= B+K2*ldb;
		 typename Field::Element_ptr A1= A;
		 typename Field::Element_ptr A2= A+K2;
		 typename Field::Element_ptr C2 = fflas_new (F, m, n,Alignment::PAGESIZE);
//#pragma omp task shared(F, A1, B1)                                                                  
		 H1.parseq.numthreads /= 2;
		 H2.parseq.numthreads = H.parseq.numthreads-H1.parseq.numthreads;
#pragma omp task shared(F, A1, C) depend(in:A1,B1) depend(inout:C)
		 pfgemm_3D_rec_adapt(F, ta, tb, m, n, K2, a, A1, lda, B1, ldb, b, C, ldc, H1);
#pragma omp task shared(F, A2, C2) depend(in:A2,B2) depend(inout:C2)
		 pfgemm_3D_rec_adapt(F, ta, tb, m, n, k-K2, a, A2, lda, B2, ldb, F.zero, C2, n, H2);
#pragma omp task shared(F) depend(inout:C) depend(in:C2)
		faddin(F, n, m, C2, n, C, ldc);
#pragma omp taskwait
		 fflas_delete(C2);
	  }
//#pragma omp taskwait                                                                                
	  return C;
 }
	
	template<class Field, class AlgoT, class FieldTrait>
	typename Field::Element*
	pfgemm_2D_rec_adapt( const Field& F,
			     const FFLAS_TRANSPOSE ta,
			     const FFLAS_TRANSPOSE tb,
			     const size_t m,
			     const size_t n,
			     const size_t k,
			     const typename Field::Element alpha,
			     const typename Field::Element_ptr A, const size_t lda,
			     const typename Field::Element_ptr B, const size_t ldb,
			     const typename Field::Element beta,
			     typename Field::Element * C, const size_t ldc, 
			     MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel> & H){

		typename Field::Element a = alpha;
		typename Field::Element b = beta;
		
		if (!m || !n) {return C;}
		if (!k || F.isZero (alpha)){
			fscalin(F, m, n, beta, C, ldc);
			return C;
		}
		if (H.parseq.numthreads<=1 || m*n<=__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD){
			MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Sequential> SeqH(H);
			return fgemm(F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, SeqH);
		}
		MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H1(H);
		MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H2(H);
		H1.parseq.numthreads /=2;
		H2.parseq.numthreads = H.parseq.numthreads - H1.parseq.numthreads;
		if(m >= n) {
			size_t M2= m>>1;
			typename Field::Element_ptr A1= A;
			typename Field::Element_ptr A2= A+M2*lda;
			typename Field::Element_ptr C1= C;
			typename Field::Element_ptr C2= C+M2*ldc;
#pragma omp task shared(F, A1, C1) depend(in:A1,B) depend(inout:C1)
			pfgemm_2D_rec_adapt(F, ta, tb, M2, n, k, alpha, A1, lda, B, ldb, beta, C1, ldc, H1);
#pragma omp task shared(F, A2, C2) depend(inout:C2) depend(in:A2,B)
			pfgemm_2D_rec_adapt(F, ta, tb, m-M2, n, k, alpha, A2, lda, B, ldb, beta, C2, ldc, H2);
#pragma omp taskwait
		} else {
			size_t N2 = n>>1;
			typename Field::Element_ptr B1= B;
			typename Field::Element_ptr B2= B+N2;
			typename Field::Element_ptr C1= C;
			typename Field::Element_ptr C2= C+N2;
#pragma omp task shared(F, B1, C1) depend(in:A,B1) depend(inout:C1)
			pfgemm_2D_rec_adapt(F, ta, tb, m, N2, k, a, A, lda, B1, ldb, b, C1, ldc, H1);
#pragma omp task shared(F, B2, C2) depend(in:A,B2) depend(inout:C2)
			pfgemm_2D_rec_adapt(F, ta, tb, m, n-N2, k, a, A, lda, B2, ldb, b,C2, ldc, H2);
#pragma omp taskwait
		}
	
//#pragma omp taskwait                                                                                
	return C;
	}

	template<class Field, class AlgoT, class FieldTrait>
	typename Field::Element*
	pfgemm_2D_rec( const Field& F,
			const FFLAS_TRANSPOSE ta,
			const FFLAS_TRANSPOSE tb,
			const size_t m,
			const size_t n,
			const size_t k,
			const typename Field::Element alpha,
			const typename Field::Element_ptr A, const size_t lda,
			const typename Field::Element_ptr B, const size_t ldb,
			const typename Field::Element beta,
			typename Field::Element * C, const size_t ldc, 
			MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel> & H){

	 size_t a = alpha;
	 size_t b = beta;

	 if (!m || !n) {return C;}
	 if (!k || F.isZero (alpha)){
		 fscalin(F, m, n, beta, C, ldc);
		 return C;
	 }
	 // threshold                                                                                                  

	 //      if(m==1 && n==1 && k==1)                                                                              
	 if(H.parseq.numthreads<=1|| m*n<=__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD){
		 MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Sequential> SeqH(H);
		 return fgemm(F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, SeqH);
	 } else 
	 {
		 size_t M2= m>>1;
		 size_t N2= n>>1;
		 
		 typename Field::Element_ptr A1= A;
		 typename Field::Element_ptr A2= A+M2*lda;
		 typename Field::Element_ptr B1= B;
		 typename Field::Element_ptr B2= B+N2;

		 typename Field::Element_ptr C11= C;
		 typename Field::Element_ptr C21= C+M2*ldc;
		 typename Field::Element_ptr C12= C+N2;
		 typename Field::Element_ptr C22= C+N2+M2*ldc;

		 MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H1(H);
		 MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H2(H);
		 MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H3(H);
		 MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H4(H);
		 MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H5(H);
		 MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H6(H);
		 MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H7(H);
		 MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H8(H);
		 int nt = H.parseq.numthreads;
		 int nt_rec = nt/4;
		 int nt_mod = nt%4;
		 H1.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0)); 
		 H2.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0)); 
		 H3.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0)); 
		 H4.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0)); 
		 // H1.parseq.numthreads = std::max(nt/4,1);
		 // H2.parseq.numthreads = std::max(nt/4,1);
		 // H3.parseq.numthreads = std::max(nt/4,1);
		 // H4.parseq.numthreads = std::max(1,nt - H1.parseq.numthreads -H2.parseq.numthreads -H3.parseq.numthreads);
#pragma omp task shared(F, A1, B1, C11) depend(inout:C11) depend(in:A1,B1)
		 pfgemm_2D_rec(F, ta, tb, M2, N2, k, alpha, A1, lda, B1, ldb, beta, C11, ldc, H1);

#pragma omp task shared(F, A1, B2, C12) depend(inout:C12) depend(in:A1,B2)
		 pfgemm_2D_rec(F, ta, tb, M2, n-N2, k, alpha, A1, lda, B2, ldb, beta, C12, ldc, H2);

#pragma omp task shared(F, A2, B1, C21) depend(inout:C21) depend(in:A2,B1)
		 pfgemm_2D_rec(F, ta, tb, m-M2, N2, k, a, A2, lda, B1, ldb, b, C21, ldc, H3);

#pragma omp task shared(F, A2, B2, C22) depend(inout:C22) depend(in:A2,B2)
		 pfgemm_2D_rec(F, ta, tb, m-M2, n-N2, k, a, A2, lda, B2, ldb, b,C22, ldc, H4);

	 }
#pragma omp taskwait
	 return C;
 }



	template<class Field, class AlgoT, class FieldTrait>
	typename Field::Element_ptr
	pfgemm_3D_rec2_V2 (const Field& F,
			   const FFLAS_TRANSPOSE ta,
			   const FFLAS_TRANSPOSE tb,
			   const size_t m,
			   const size_t n,
			   const size_t k,
			   const typename Field::Element alpha,
			   const typename Field::Element_ptr A, const size_t lda,
			   const typename Field::Element_ptr B, const size_t ldb,
			   const typename Field::Element beta,
			   typename Field::Element_ptr C, const size_t ldc, 
			   MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel> & H){
		/*
	size_t a = 1.0;
	size_t b = 0.0;
		*/
	//	int wino = WWINO;
	if (!m || !n) {return C;}
	if (!k || F.isZero (alpha)){
		fscalin(F, m, n, beta, C, ldc);
		return C;
	}
	if(H.parseq.numthreads <= 1|| std::min(m*n,std::min(m*k,k*n))<=__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD){	// threshold
		FFLAS::MMHelper<Field, AlgoT, FieldTrait,FFLAS::ParSeqHelper::Sequential> WH (H);
		return fgemm(F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, WH);
	}
	else
	{
		typename Field::Element a = alpha;
		typename Field::Element b = 0;

		size_t M2= m>>1;
		size_t N2= n>>1;
		size_t K2= k>>1;
		typename Field::Element_ptr A11= A;
		typename Field::Element_ptr A12= A+K2;
		typename Field::Element_ptr A21= A+M2*lda;
		typename Field::Element_ptr A22= A+K2+M2*lda;

		typename Field::Element_ptr B11= B;
		typename Field::Element_ptr B12= B+N2;
		typename Field::Element_ptr B21= B+K2*ldb;
		typename Field::Element_ptr B22= B+N2+K2*ldb;

		typename Field::Element_ptr C11= C;
		typename Field::Element_ptr C_11 = fflas_new (F, M2, N2,Alignment::PAGESIZE);
		
		typename Field::Element_ptr C12= C+N2;
		typename Field::Element_ptr C_12 = fflas_new (F, M2, n-N2,Alignment::PAGESIZE);
		
		typename Field::Element_ptr C21= C+M2*ldc;
		typename Field::Element_ptr C_21 = fflas_new (F, m-M2, N2,Alignment::PAGESIZE);		

		typename Field::Element_ptr C22= C+N2+M2*ldc;
		typename Field::Element_ptr C_22 = fflas_new (F, m-M2, n-N2,Alignment::PAGESIZE);

		// 1/ 8 multiply in parallel
		    //omp_set_task_affinity(omp_get_locality_domain_num_for( C11));
		
		MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H1(H);
		MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H2(H);
		MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H3(H);
		MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H4(H);
		MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H5(H);
		MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H6(H);
		MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H7(H);
		MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H8(H);
		int nt = H.parseq.numthreads;
		int nt_rec = nt/8;
		int nt_mod = nt % 8 ;
		H1.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0));
		H2.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0)); 
		H3.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0)); 
		H4.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0)); 
		H5.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0)); 
		H6.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0)); 
		H7.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0)); 
		H8.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0)); 

                #pragma omp task shared(F, A11, B11, C11) depend(out:C11) depend(in:A11,B11)
		pfgemm_3D_rec2_V2(F, ta, tb, M2, N2, K2, alpha, A11, lda, B11, ldb, beta, C11, ldc, H1);
		    //omp_set_task_affinity(omp_get_locality_domain_num_for( C_11));
                #pragma omp task shared(F, A12, B21, C_11) depend(out:C_11) depend(in:A12,B21)
		pfgemm_3D_rec2_V2(F, ta, tb, M2, N2, k-K2, a, A12, lda, B21, ldb, b,C_11, N2, H2);
		    //omp_set_task_affinity(omp_get_locality_domain_num_for( C12));
                #pragma omp task shared(F, A12, B22, C12) depend(out:C12) depend(in:A12,B22)
		pfgemm_3D_rec2_V2(F, ta, tb, M2, n-N2, k-K2, alpha, A12, lda, B22, ldb, beta, C12, ldc, H3);
		    //omp_set_task_affinity(omp_get_locality_domain_num_for( C_12));
                #pragma omp task shared(F, A11, B12, C_12) depend(out:C_12) depend(in:A11,B12)
		pfgemm_3D_rec2_V2(F, ta, tb, M2, n-N2, K2, a, A11, lda, B12, ldb, b, C_12, n-N2, H4);
		    //omp_set_task_affinity(omp_get_locality_domain_num_for( C21));
                #pragma omp task shared(F,A22, B21, C21) depend(out:C21) depend(in:A22,B21)
		pfgemm_3D_rec2_V2(F, ta, tb, m-M2, N2, k-K2, alpha, A22, lda, B21, ldb, beta, C21, ldc, H5);
		    //omp_set_task_affinity(omp_get_locality_domain_num_for( C_21));
                #pragma omp task shared(F, A21, B11, C_21) depend(out:C_21) depend(in:A21, B11)
		pfgemm_3D_rec2_V2(F, ta, tb, m-M2, N2, K2, a, A21, lda, B11, ldb, b,C_21, N2, H6);
		    //omp_set_task_affinity(omp_get_locality_domain_num_for( C22));
                #pragma omp task shared(F, A21, B12, C22) depend(out:C22) depend(in:A21,B12)
		pfgemm_3D_rec2_V2(F, ta, tb, m-M2, n-N2, K2, alpha, A21, lda, B12, ldb, beta, C22, ldc, H7);
		    //omp_set_task_affinity(omp_get_locality_domain_num_for( C_22));
                #pragma omp task shared(F, A22, B22, C_22) depend(out:C_22) depend(in:A22,B22)
		pfgemm_3D_rec2_V2(F, ta, tb, m-M2, n-N2, k-K2, a, A22, lda, B22, ldb, b,C_22, n-N2, H8);

		// 2/ final add
		    //omp_set_task_affinity(omp_get_locality_domain_num_for( C11));
                #pragma omp task shared(F, C11, C_11) depend(inout:C11) depend(in:C_11)
		faddin(F, M2, N2, C_11, N2, C11, ldc);
		    //omp_set_task_affinity(omp_get_locality_domain_num_for( C12));
                #pragma omp task shared(F, C12, C_12) depend(inout:C12) depend(in:C_12)
		faddin(F, M2, n-N2, C_12, n-N2, C12, ldc);
		    //omp_set_task_affinity(omp_get_locality_domain_num_for( C21));
                #pragma omp task shared(F, C21, C_21) depend(inout:C21) depend(in:C_21)
		faddin(F, m-M2, N2, C_21, N2, C21, ldc);
		    //omp_set_task_affinity(omp_get_locality_domain_num_for( C22));
                #pragma omp task shared(F, C22, C_22) depend(inout:C22) depend(in:C_22)
		faddin(F, m-M2, n-N2, C_22, n-N2, C22, ldc);

                #pragma omp taskwait	
		
        	FFLAS::fflas_delete (C_11);
		FFLAS::fflas_delete (C_12);
		FFLAS::fflas_delete (C_21);
		FFLAS::fflas_delete (C_22);
	}
	return C;
}

	template<class Field, class AlgoT, class FieldTrait>
	typename Field::Element*
	pfgemm_3D_rec_V2( const Field& F,
			const FFLAS_TRANSPOSE ta,
			const FFLAS_TRANSPOSE tb,
			const size_t m,
			const size_t n,
			const size_t k,
			const typename Field::Element alpha,
			const typename Field::Element_ptr A, const size_t lda,
			const typename Field::Element_ptr B, const size_t ldb,
		  const typename Field::Element beta,
		  typename Field::Element_ptr C, const size_t ldc, 
		  MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel> & H){
	
	
	if (!m || !n) {return C;}
	if (!k || F.isZero (alpha)){
		fscalin(F, m, n, beta, C, ldc);
		return C;
	}

	if(H.parseq.numthreads <= 1|| std::min(m*n,std::min(m*k,k*n))<=__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD){	// threshold
		FFLAS::MMHelper<Field, AlgoT, FieldTrait,FFLAS::ParSeqHelper::Sequential> WH (H);
		return fgemm(F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, WH);
	}else{
		size_t M2= m>>1;
		size_t N2= n>>1;
		size_t K2= k>>1;
		typename Field::Element_ptr A11= A;
		typename Field::Element_ptr A12= A+K2;
		typename Field::Element_ptr A21= A+M2*lda;
		typename Field::Element_ptr A22= A+K2+M2*lda;

		typename Field::Element_ptr B11= B;
		typename Field::Element_ptr B12= B+N2;
		typename Field::Element_ptr B21= B+K2*ldb;
		typename Field::Element_ptr B22= B+N2+K2*ldb;

		typename Field::Element_ptr C11= C;
		typename Field::Element_ptr C12= C+N2;
		typename Field::Element_ptr C21= C+M2*ldc;
		typename Field::Element_ptr C22= C+N2+M2*ldc;

		MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H1(H);
		MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H2(H);
		MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H3(H);
		MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel> H4(H);
		int nt = H.parseq.numthreads;
		int nt_rec = nt/4;
		int nt_mod = nt%4;
		H1.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0)); 
		H2.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0)); 
		H3.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0)); 
		H4.parseq.numthreads = std::max(1,nt_rec + ((nt_mod-- > 0)?1:0)); 

                // 1/ 4 multiply
                #pragma omp task shared(F, A11, B11, C11) depend(inout:C11) depend(in:A11, B11)
		pfgemm_3D_rec_V2(F, ta, tb, M2, N2, K2, alpha, A11, lda, B11, ldb, beta, C11, ldc, H1);
                #pragma omp task shared(F, A12, B22, C12) depend(inout:C12 ) depend(in:A12, B22)
		pfgemm_3D_rec_V2(F, ta, tb, M2, n-N2, k-K2, alpha, A12, lda, B22, ldb, beta, C12, ldc, H2);
                #pragma omp task shared(F, A22, B21, C21) depend(inout:C21) depend(in:A22, B21)
		pfgemm_3D_rec_V2(F, ta, tb, m-M2, N2, k-K2, alpha, A22, lda, B21, ldb, beta, C21, ldc, H3);
                #pragma omp task shared(F, A21, B12, C22) depend(inout:C22) depend(in:A21, B12)
		pfgemm_3D_rec_V2(F, ta, tb, m-M2, n-N2, K2, alpha, A21, lda, B12, ldb, beta, C22, ldc, H4);

                // 2/ 4 add+multiply
                #pragma omp task shared(F, A12, B21, C11) depend(inout:C11) depend(in:A12, B21)
		pfgemm_3D_rec_V2(F, ta, tb, M2, N2, k-K2, alpha, A12, lda, B21, ldb, F.one, C11, ldc, H1);
                #pragma omp task shared(F, A11, B12, C12)  depend(inout:C12) depend(in:A11, B12)
		pfgemm_3D_rec_V2(F, ta, tb, M2, n-N2, K2, alpha, A11, lda, B12, ldb, F.one, C12, ldc, H2);
                #pragma omp task shared(F, A21, B11, C21)  depend(inout:C21) depend(in:B11, A21)
		pfgemm_3D_rec_V2(F, ta, tb, m-M2, N2, K2, alpha, A21, lda, B11, ldb, F.one, C21, ldc, H3);
                #pragma omp task shared(F, A22, B22, C22)  depend(inout:C22) depend(in:A22, B22)
		pfgemm_3D_rec_V2(F, ta, tb, m-M2, n-N2, k-K2, alpha, A22, lda, B22, ldb, F.one, C22, ldc, H4);

                #pragma omp taskwait
	}
	return C;
}


	template<class Field, class AlgoT, class FieldTrait>
	inline typename Field::Element_ptr
	fgemm( const Field& F,
		const FFLAS::FFLAS_TRANSPOSE ta,
		const FFLAS::FFLAS_TRANSPOSE tb,
		const size_t m,
		const size_t n,
		const size_t k,
		const typename Field::Element alpha,
		typename Field::Element_ptr A, const size_t lda,
		typename Field::Element_ptr B, const size_t ldb,
		const typename Field::Element beta,
		typename Field::Element_ptr C, const size_t ldc,
		MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel> & H) 
	{

            if ((ta != FFLAS::FflasNoTrans) || (tb != FFLAS::FflasNoTrans)) {
                std::cerr << "*** ERROR ***: pfgemm ^T NOT YET IMPLEMENTED" << std::endl;
                return C;
            }

	    switch (H.parseq.method){
		case THREE_D_ADAPT: // Splitting 1 dimension at a time recursively: the largest one
			return pfgemm_3D_rec_adapt (F, ta, tb, m, n, k ,alpha, A, lda, B, ldb, beta, C, ldc, H);
		case TWO_D_ADAPT: // Splitting 1 dimension at a time recursively: the largest one
			return pfgemm_2D_rec_adapt (F, ta, tb, m, n, k ,alpha, A, lda, B, ldb, beta, C, ldc, H);
		case TWO_D: // Splitting the outer dimensions m and n recursively
			return pfgemm_2D_rec (F, ta, tb, m, n, k ,alpha, A, lda, B, ldb, beta, C, ldc, H);
		case THREE_D_INPLACE: // Splitting the three dimensions recursively, without temp and with synchro
			return pfgemm_3D_rec_V2(F, ta, tb, m, n, k ,alpha, A, lda, B, ldb, beta, C, ldc, H);
		case THREE_D: // Splitting the three dimensions recursively, with temp alloc and fewer synchro
			return pfgemm_3D_rec2_V2(F, ta, tb, m, n, k ,alpha, A, lda, B, ldb, beta, C, ldc, H);
		default: // 2D iterative: splitting the outer dimensions m and n iteratively 
			H.parseq.numthreads = std::min(H.parseq.numthreads, std::max(1,int(m*n/(__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD))));
			ForStrategy2D iter(m,n,H.parseq.method,H.parseq.numthreads);
			// if (H.recLevel < 0) 
			// 	H.recLevel = Protected::WinogradSteps (F, min3(iter.rowBlockSize,k,iter.colBlockSize));
			for (iter.begin(); ! iter.end(); ++iter){
				MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Sequential> SeqH (H);
				TASK(READ(A,B,F), NOWRITE(), READWRITE(C), fgemm, F, ta, tb, iter.iend-iter.ibeg, iter.jend-iter.jbeg, k, alpha, A + iter.ibeg*lda, lda, B +iter.jbeg, ldb, beta, C+ iter.ibeg*ldc+iter.jbeg, ldc, SeqH);
			}
			WAIT;
	    }
	    return C;
	}

} // FFLAS
#else
namespace FFLAS{
	template<class Field, class AlgoT, class FieldTrait>
	inline typename Field::Element_ptr
	fgemm( const Field& F,
		const FFLAS::FFLAS_TRANSPOSE ta,
		const FFLAS::FFLAS_TRANSPOSE tb,
		const size_t m,
		const size_t n,
		const size_t k,
		const typename Field::Element alpha,
		typename Field::Element_ptr A, const size_t lda,
		typename Field::Element_ptr B, const size_t ldb,
		const typename Field::Element beta,
		typename Field::Element_ptr C, const size_t ldc,
		MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel> & H) 
	{

            if ((ta != FFLAS::FflasNoTrans) || (tb != FFLAS::FflasNoTrans)) {
                std::cerr << "*** ERROR ***: pfgemm ^T NOT YET IMPLEMENTED" << std::endl;
                return C;
            }

	    /*	    std::cout<<H.parseq.numthreads<<" "<<std::endl;
	    std::cout<<(m*n/(__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD))<<std::endl;
	    std::cout<<"computation"<<std::endl;
	    */
		// Threshold: no need to slice to blocks smaller than __FFLASFFPACK_SEQPARTHRESHOLD
	    H.parseq.numthreads = std::min(H.parseq.numthreads, std::max(1,(int)(m*n/(__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD))));
	    ForStrategy2D iter(m,n,H.parseq.method,H.parseq.numthreads);
	    if (H.recLevel < 0) 
		    H.recLevel = Protected::WinogradSteps (F, min3(iter.rowBlockSize,k,iter.colBlockSize));
	    for (iter.begin(); ! iter.end(); ++iter){
		MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Sequential> SeqH (H);
		TASK(READ(A,B,F), NOWRITE(), READWRITE(C), fgemm, F, ta, tb, iter.iend-iter.ibeg, iter.jend-iter.jbeg, k, alpha, A + iter.ibeg*lda, lda, B +iter.jbeg, ldb, beta, C+ iter.ibeg*ldc+iter.jbeg, ldc, SeqH);
	}
	    
	WAIT;
	return C;
	}
}
#endif // __FFLASFFPACK_USE_OPENMP4
#endif // __FFLASFFPACK_fflas_parfgemm_INL

