/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_pfgemm.inl
 * Copyright (C) 2013 Jean Guillaume Dumas Clement Pernet Ziad Sultan
 *
 * Written by Jean Guillaume Dumas Clement Pernet Ziad Sultan
 * Time-stamp: <04 Sep 14 14:08:06 Jean-Guillaume.Dumas@imag.fr>
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


#ifdef __FFLASFFPACK_USE_KAAPI
#include <kaapi++>
#include "kaapi_routines.inl"
#endif
#ifdef __FFLASFFPACK_USE_OPENMP
#include <omp.h>
#endif

#include "fflas_blockcuts.inl"
#include "parallel.h"


namespace FFLAS {


	template<class Field>
	typename Field::Element*
	pfgemm_1D_rec( const Field& F,
		       const FFLAS_TRANSPOSE ta,
		       const FFLAS_TRANSPOSE tb,
		       const size_t m,
		       const size_t n,
		       const size_t k,
		       const typename Field::Element alpha,
		       const typename Field::Element_ptr A, const size_t lda,
		       const typename Field::Element_ptr B, const size_t ldb,
		       const typename Field::Element beta,
		       typename Field::Element * C, const size_t ldc, size_t seuil){

		size_t a = alpha;
		size_t b = beta;

		if (!m || !n) {return C;}
		if (!k || F.isZero (alpha)){
			fscalin(F, m, n, beta, C, ldc);
			return C;
		}
		// threshold                                                                                                                                                                        

		//      if(m==1 && n==1 && k==1)                                                                                                                                                    
		if(m<seuil || n<seuil || k<seuil)
			return fgemm(F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
		else if(m==std::max(std::max(m,n),k))
			{
				size_t M2= m>>1;

				typename Field::Element_ptr A1= A;
				typename Field::Element_ptr A2= A+M2*lda;
				/*                                                                                                                                                                  
                typename Field::Element_ptr B11= B;                                                                                                                                          
                typename Field::Element_ptr B12= B+N2;                                                                                                                                       
                typename Field::Element_ptr B21= B+K2*ldb;                                                                                                                                   
                typename Field::Element_ptr B22= B+N2+K2*ldb;                                                                                                                                
				*/

				typename Field::Element_ptr C1= C;
				typename Field::Element_ptr C2= C+M2*ldc;

				// 2 multiply (1 split on dimension m)                                                                                                                              
#pragma omp task shared(F, A1, C1) depend(in:A1,B) depend(inout:C1)
				pfgemm_1D_rec(F, ta, tb, M2, n, k, alpha, A1, lda, B, ldb, beta, C1, ldc, seuil);

#pragma omp task shared(F, A2, C2) depend(inout:C2) depend(in:A2,B)
				pfgemm_1D_rec(F, ta, tb, m-M2, n, k, alpha, A2, lda, B, ldb, beta, C2, ldc,seuil);

				//                #pragma omp taskwait                                                                                                                              
			}
		else if(n==std::max(n,k))
			{

				size_t N2 = n>>1;
				typename Field::Element_ptr B1= B;
				typename Field::Element_ptr B2= B+N2;

				typename Field::Element_ptr C1= C;
				typename Field::Element_ptr C2= C+N2;

#pragma omp task shared(F, B1, C1) depend(in:A,B1) depend(inout:C1)
				pfgemm_1D_rec(F, ta, tb, m, N2, k, a, A, lda, B1, ldb, b, C1, ldc,seuil);

#pragma omp task shared(F, B2, C2) depend(in:A,B2) depend(inout:C2)
				pfgemm_1D_rec(F, ta, tb, m, n-N2, k, a, A, lda, B2, ldb, b,C2, ldc,seuil);

				//                #pragma omp taskwait                                                                                                                              
			}
		else if(k==std::max(n,k))
			{
				size_t K2 = k>>1;


				typename Field::Element_ptr B1= B;
				typename Field::Element_ptr B2= B+K2*ldb;

				typename Field::Element_ptr A1= A;
				typename Field::Element_ptr A2= A+K2;
				// 2 variantes possibles avec tmp et avec synchro(depend)

				//#pragma omp task shared(F, A1, B1)                                                                                                                                
				pfgemm_1D_rec(F, ta, tb, m, n, K2, a, A1, lda, B1, ldb, b, C, ldc,seuil);
				//              #ragma omp taskwait                                                                                                                                 
				//              #pragma omp task shared(F, A2, B2, C)                                                                                                               
				b=1;
				pfgemm_1D_rec(F, ta, tb, m, n, k-K2, a, A2, lda, B2, ldb, b, C, ldc,seuil);


			}
                         #pragma omp taskwait
		return C;
	}



 template<class Field>
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
			typename Field::Element * C, const size_t ldc, size_t seuil){

	 size_t a = alpha;
	 size_t b = beta;

	 if (!m || !n) {return C;}
	 if (!k || F.isZero (alpha)){
		 fscalin(F, m, n, beta, C, ldc);
		 return C;
	 }
	 // threshold                                                                                                  

	 //      if(m==1 && n==1 && k==1)                                                                              
	 if(m<seuil || n<seuil || k<seuil)
		 return fgemm(F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
	 else 
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


#pragma omp task shared(F, A1, B1, C11) depend(inout:C11) depend(in:A1,B1)
		 pfgemm_1D_rec(F, ta, tb, M2, N2, k, alpha, A1, lda, B1, ldb, beta, C11, ldc, seuil);

#pragma omp task shared(F, A1, B2, C12) depend(inout:C12) depend(in:A1,B2)
		 pfgemm_1D_rec(F, ta, tb, M2, n-N2, k, alpha, A1, lda, B2, ldb, beta, C12, ldc,seuil);

#pragma omp task shared(F, A2, B1, C21) depend(inout:C21) depend(in:A2,B1)
		 pfgemm_1D_rec(F, ta, tb, m-M2, N2, k, a, A2, lda, B1, ldb, b, C21, ldc,seuil);

#pragma omp task shared(F, A2, B2, C22) depend(inout:C22) depend(in:A2,B2)
		 pfgemm_1D_rec(F, ta, tb, m-M2, n-N2, k, a, A2, lda, B2, ldb, b,C22, ldc,seuil);

	 }
#pragma omp taskwait
	 return C;
 }



template<class Field>
typename Field::Element*
pfgemm_3D_rec( const Field& F,
		const FFLAS_TRANSPOSE ta,
		const FFLAS_TRANSPOSE tb,
		const size_t m,
		const size_t n,
		const size_t k,
		const typename Field::Element alpha,
		const typename Field::Element_ptr A, const size_t lda,
		const typename Field::Element_ptr B, const size_t ldb,
		const typename Field::Element beta,
	       typename Field::Element_ptr C, const size_t ldc, size_t seuil, size_t * x){

	size_t a = 1;
	size_t b = 1;
	
	if (!m || !n) {return C;}
	if (!k || F.isZero (alpha)){
		fscalin(F, m, n, beta, C, ldc);
		return C;
	}
	// threshold
	
	if(m<seuil || n<seuil || k<seuil)
		return fgemm(F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
	else
	{
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

		//handling dependencies
		size_t * x1;
		size_t * x2;
		size_t * x3;
		size_t * x4;

// 1/ 4 multiply
                #pragma omp task shared(F, A11, B11, C11) depend(inout:x1)
		pfgemm_3D_rec(F, ta, tb, M2, N2, K2, alpha, A11, lda, B11, ldb, beta, C11, ldc, seuil, x1);

                #pragma omp task shared(F, A12, B22, C12) depend(out:x2)
		pfgemm_3D_rec(F, ta, tb, M2, n-N2, k-K2, alpha, A12, lda, B22, ldb, beta, C12, ldc,seuil, x2);

                #pragma omp task shared(F, A22, B21, C21) depend(out:x3)
		pfgemm_3D_rec(F, ta, tb, m-M2, N2, k-K2, alpha, A22, lda, B21, ldb, beta, C21, ldc,seuil, x3);

                #pragma omp task shared(F, A21, B12, C22) depend(out:x4)
		pfgemm_3D_rec(F, ta, tb, m-M2, n-N2, K2, alpha, A21, lda, B12, ldb, beta, C22, ldc,seuil, x4);
// Sync
//		#pragma omp taskwait

// 2/ 4 add+multiply
		
                #pragma omp task shared(F, A12, B21, C11)  depend(out:x1)
		pfgemm_3D_rec(F, ta, tb, M2, N2, k-K2, a, A12, lda, B21, ldb, b,C11, ldc,seuil, x1);

		#pragma omp task shared(F, A11, B12, C12)  depend(out:x2)
		pfgemm_3D_rec(F, ta, tb, M2, n-N2, K2, a, A11, lda, B12, ldb, b,C12, ldc,seuil, x2);

		#pragma omp task shared(F, A21, B11, C21)  depend(out:x3)
		pfgemm_3D_rec(F, ta, tb, m-M2, N2, K2, a, A21, lda, B11, ldb, b,C21, ldc,seuil, x3);

		#pragma omp task shared(F, A22, B22, C22)  depend(out:x4)
		pfgemm_3D_rec(F, ta, tb, m-M2, n-N2, k-K2, a, A22, lda, B22, ldb, b,C22, ldc,seuil, x4);

		#pragma omp taskwait

	}
	
		return C;
}

template<class Field>
typename Field::Element_ptr
pfgemm_3D_rec2( const Field& F,
		const FFLAS_TRANSPOSE ta,
		const FFLAS_TRANSPOSE tb,
		const size_t m,
		const size_t n,
		const size_t k,
		const typename Field::Element alpha,
		const typename Field::Element_ptr A, const size_t lda,
		const typename Field::Element_ptr B, const size_t ldb,
		const typename Field::Element beta,
		typename Field::Element_ptr C, const size_t ldc, size_t seuil, size_t *x){

	size_t a = 1.0;
	size_t b = 0.0;
	
	if (!m || !n) {return C;}
	if (!k || F.isZero (alpha)){
		fscalin(F, m, n, beta, C, ldc);
		return C;
	}
	if(m<seuil || n<seuil || k<seuil){	// threshold
		//                #pragma omp task shared(F, C)
		return fgemm(F, ta, tb, m, n, k, a, A, lda, B, ldb, b, C, ldc);
		//                #pragma omp taskwait
	}
	else
	{
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
		typename Field::Element_ptr C_11 = fflas_new (F, M2, N2);

		typename Field::Element_ptr C12= C+N2;
		typename Field::Element_ptr C_12 = fflas_new (F, M2, n-N2);
		
		typename Field::Element_ptr C21= C+M2*ldc;
		typename Field::Element_ptr C_21 = fflas_new (F, m-M2, N2);		

		typename Field::Element_ptr C22= C+N2+M2*ldc;
		typename Field::Element_ptr C_22 = fflas_new (F, m-M2, n-N2);

		//dependencies handling		
		size_t * x1;
		size_t * x2;
		size_t * x3;
		size_t * x4;

		// 1/ 8 multiply in parallel
                #pragma omp task shared(F, A11, B11, C11) depend(in:x1)
		pfgemm_3D_rec2(F, ta, tb, M2, N2, K2, a, A11, lda, B11, ldb, b,C11, ldc, seuil, x1);
                #pragma omp task shared(F, A12, B21, C_11) depend(in:x1)
		pfgemm_3D_rec2(F, ta, tb, M2, N2, k-K2, a, A12, lda, B21, ldb, b,C_11, N2,seuil,x1);

#pragma omp task shared(F, A12, B22, C12) depend(in:x2)
		pfgemm_3D_rec2(F, ta, tb, M2, n-N2, k-K2, a, A12, lda, B22, ldb, b,C12, ldc,seuil,x2);
		#pragma omp task shared(F, A11, B12, C_12) depend(in:x2)
		pfgemm_3D_rec2(F, ta, tb, M2, n-N2, K2, a, A11, lda, B12, ldb, b, C_12, n-N2,seuil,x2);

                #pragma omp task shared(F,A22, B21, C21) depend(in:x3)
		pfgemm_3D_rec2(F, ta, tb, m-M2, N2, k-K2, a, A22, lda, B21, ldb, b,C21, ldc,seuil,x3);
		#pragma omp task shared(F, A21, B11, C_21) depend(in:x3)
		pfgemm_3D_rec2(F, ta, tb, m-M2, N2, K2, a, A21, lda, B11, ldb, b,C_21, N2,seuil,x3);

                #pragma omp task shared(F, A21, B12, C22) depend(in:x4)
		pfgemm_3D_rec2(F, ta, tb, m-M2, n-N2, K2, a, A21, lda, B12, ldb, b,C22, ldc,seuil,x4);
		#pragma omp task shared(F, A22, B22, C_22) depend(in:x4)
		pfgemm_3D_rec2(F, ta, tb, m-M2, n-N2, k-K2, a, A22, lda, B22, ldb, b,C_22, n-N2,seuil,x4);

		//#pragma omp taskwait
		
		// 2/ final add
		// modifier les add, 

#pragma omp task shared(F, C11, C_11) depend(inout:x1)
		faddin(F, M2, N2, C_11, N2, C11, ldc);
#pragma omp task shared(F, C12, C_12) depend(inout:x2)
		faddin(F, M2, n-N2, C_12, n-N2, C12, ldc);
#pragma omp task shared(F, C21, C_21) depend(inout:x3)
		faddin(F, m-M2, N2, C_21, N2, C21, ldc);
#pragma omp task shared(F, C22, C_22) depend(inout:x4)
		faddin(F, m-M2, n-N2, C_22, n-N2, C22, ldc);
#pragma omp taskwait	
		/*
		for (size_t i=0; i<M2; i++){
#pragma omp task shared(F, C11, C_11) depend(inout:x1)
			{
			for(size_t j=0; j<N2; j++)
				F.add (C11[i*ldc+j], C11[i*ldc+j], C_11[i*N2+j]);
			}
#pragma omp task shared(F, C12, C_12) depend(inout:x2)
			{
			for(size_t j=0; j<n-N2; j++)
				F.add (C12[i*ldc+j], C12[i*ldc+j], C_12[i*(n-N2)+j]);
			}
		}
		for (size_t i=0; i<m-M2; i++){
#pragma omp task shared(F, C21, C_21) depend(inout:x3)
			{
			for (size_t j=0; j<N2; j++)
				F.add (C21[i*ldc+j], C21[i*ldc+j], C_21[i*N2+j]);
			}
#pragma omp task shared(F, C22, C_22) depend(inout:x4)
			{
			for (size_t j=0; j<n-N2; j++)
				F.add (C22[i*ldc+j], C22[i*ldc+j], C_22[i*(n-N2)+j]);
		        }
		}
		*/
		/*
		//                #pragma omp task shared(F, C11, C_11)
		fadd(F, M2, C11, ldc, C_11, N2, C11, ldc);
		//                #pragma omp task shared(F, C12, C_12)
		fadd(F, M2, C12, ldc, C_12, n-N2, C12, ldc);
		//                #pragma omp task shared(F, C21, C_21)
		fadd(F, m-M2, C21, ldc, C_21, N2, C21, ldc);
		//                #pragma omp task shared(F, C22, C_22)
		fadd(F, m-M2, C22, ldc, C_22, n-N2, C22, ldc);
		*/

		

        	FFLAS::fflas_delete (C_11);
		FFLAS::fflas_delete (C_12);
		FFLAS::fflas_delete (C_21);
		FFLAS::fflas_delete (C_22);

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
            
		ForStrategy2D iter(m,n,H.parseq.method,H.parseq.numthreads);
		if (H.recLevel < 0) {
			H.recLevel = Protected::WinogradSteps (F, min3(iter.rowBlockSize,k,iter.colBlockSize));
		}
		for (iter.begin(); ! iter.end(); ++iter){
			MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Sequential> SeqH (H);
			TASK(READ(A,B,F), NOWRITE(), READWRITE(C), fgemm, F, ta, tb, iter.iend-iter.ibeg, iter.jend-iter.jbeg, k, alpha, A + iter.ibeg*lda, lda, B +iter.jbeg, ldb, beta, C+ iter.ibeg*ldc+iter.jbeg, ldc, SeqH);
		}

		WAIT;
		//		BARRIER;

		return C;
	}

} // FFLAS

#endif // __FFLASFFPACK_fflas_parfgemm_INL

