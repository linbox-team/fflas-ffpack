/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_pfgemm.inl
 * Copyright (C) 2013 Jean Guillaume Dumas Clement Pernet Ziad Sultan
 *
 * Written by Jean Guillaume Dumas Clement Pernet Ziad Sultan
 * Time-stamp: <09 Dec 14 10:02:11 Jean-Guillaume.Dumas@imag.fr>
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


//#endif

#include "fflas_blockcuts.inl"
#include "parallel.h"
#include "fflas-ffpack/utils/timer.h"

//#ifdef __FFLASFFPACK_USE_DATAFLOW
#include "fflas_pfgemm_variants.inl"

//#else

namespace FFLAS {

	template<class Field, class AlgoT, class FieldTrait>
	inline typename Field::Element_ptr
	fgemm_iter( const Field& F,
		    ForStrategy2D iter,
		    const FFLAS::FFLAS_TRANSPOSE ta,
		    const FFLAS::FFLAS_TRANSPOSE tb,
		    const size_t m,
		    const size_t n,
		    const size_t k,
		    const typename Field::Element alpha,
		    typename Field::ConstElement_ptr A, const size_t lda,
		    typename Field::ConstElement_ptr B, const size_t ldb,
		    const typename Field::Element beta,
		    typename Field::Element_ptr C, const size_t ldc,
		    MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Sequential> & H){

		A = A+iter.ibeg*lda;
		B = B+iter.jbeg;
		C = C+iter.ibeg*ldc+iter.jbeg;

		fgemm( F, ta, tb, iter.iend-iter.ibeg, iter.jend-iter.jbeg, k, alpha, A, lda, B, ldb, beta, C, ldc, H);
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
		typename Field::ConstElement_ptr A, const size_t lda,
		typename Field::ConstElement_ptr B, const size_t ldb,
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
			H.parseq.numthreads = std::min(H.parseq.numthreads, std::max((size_t)1,(size_t)(m*n/(__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD))));
			/*
			typename Field::Element_ptr CC = 0;
			typename Field::ConstElement_ptr AA = 0;
			typename Field::ConstElement_ptr BB = 0;
			*/
			//			ForStrategy2D iter(m,n,H.parseq);

			//			for (iter.begin(); ! iter.end(); ++iter){
				//				CC = C+iter.ibeg*ldc+iter.jbeg;
			TASKFOR2D(iter,m,n,H.parseq){
				
				MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Sequential> SeqH (H);

				TASK(READ(A,B,F), NOWRITE(), READWRITE(C), fgemm, F, ta, tb, iter.iend-iter.ibeg, iter.jend-iter.jbeg, k, alpha, A+iter.ibeg*lda, lda, B+iter.jbeg, ldb, beta, C+iter.ibeg*ldc+iter.jbeg, ldc, SeqH);
				//				BEGINTASK(TASK(READ(A,B,F), NOWRITE(), READWRITE(CC))
				//				typename Field::Element_ptr CCC=0;
				/*				
#pragma omp task shared(AA, BB, CC, F) depend(in:AA, BB) depend(out:CC)
				{
					  AA = A+iter.ibeg*lda;
					  BB = B+iter.jbeg;

					  CC = C+iter.ibeg*ldc+iter.jbeg;
					  cout<<" iter.ibeg "<<iter.ibeg<<"iter.jbeg "<<iter.jbeg<<endl;
					  fgemm(F, ta, tb, iter.iend-iter.ibeg, iter.jend-iter.jbeg, k, alpha, AA, lda, BB, ldb, beta, CC, ldc, SeqH);
				}
					  //		ENDTASK
					  */
			}
			WAIT;
	    }
		return C;
}



}
//#endif //DATAFLOW
#endif // __FFLASFFPACK_fflas_parfgemm_INL

