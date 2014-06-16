/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_pfgemm.inl
 * Copyright (C) 2013 Jean Guillaume Dumas Clement Pernet Ziad Sultan
 *
 * Written by Jean Guillaume Dumas Clement Pernet Ziad Sultan
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

#ifndef __FFLASFFPACK_fflas_pfgmm_INL
#define __FFLASFFPACK_fflas_pfgmm_INL


#ifdef __FFLASFFPACK_USE_KAAPI
#include <kaapi++>
#endif
#ifdef __FFLASFFPACK_USE_OPENMP
#include <omp.h>
#endif

#include "fflas-ffpack/ffpack/parallel.h"


namespace FFLAS {



#ifdef __FFLASFFPACK_USE_KAAPI
	template<class Field>			
	struct Taskfgemm : public ka::Task<15>::Signature<
		Field,
		FFLAS_TRANSPOSE,
		FFLAS_TRANSPOSE,
		size_t ,
		size_t ,
		size_t ,
		typename Field::Element,
		ka::R<typename Field::Element>,
		size_t ,
		ka::R<typename Field::Element>,
		size_t ,
		typename Field::Element,
		ka::RW<typename Field::Element>,
		size_t ,
		size_t
		>{};
	template<class Field>			
	struct Taskfgemmw : public ka::Task<14>::Signature<
		Field,
		FFLAS_TRANSPOSE,
		FFLAS_TRANSPOSE,
		size_t ,
		size_t ,
		size_t ,
		typename Field::Element,
		ka::R<typename Field::Element>,
		size_t ,
		ka::R<typename Field::Element>,
		size_t ,
		typename Field::Element,
		ka::RW<typename Field::Element>,
		size_t 
		>{};
}

template<class Field>
struct TaskBodyCPU<FFLAS::Taskfgemm<Field> >{
	void operator()(const Field& F,
			const FFLAS::FFLAS_TRANSPOSE ta,
			const FFLAS::FFLAS_TRANSPOSE tb,
			size_t BlockRowDim,
			size_t BlockColDim,
			size_t k,
			const typename Field::Element alpha,
			ka::pointer_r<typename Field::Element> A,
			const size_t lda,
			ka::pointer_r<typename Field::Element> B,
			const size_t ldb,
			const typename Field::Element beta,
			ka::pointer_rw<typename Field::Element> C, const size_t ldc,
			const size_t w)
	{
		FFLAS::fgemm( F, ta, tb, BlockRowDim, BlockColDim, k, alpha, A.ptr(), lda, B.ptr() , ldb,
			      beta, C.ptr(), ldc, w);
	}
};

template<class Field>
struct TaskBodyCPU<FFLAS::Taskfgemmw<Field> >{
	void operator()(const Field& F,
			const FFLAS::FFLAS_TRANSPOSE ta,
			const FFLAS::FFLAS_TRANSPOSE tb,
			size_t BlockRowDim,
			size_t BlockColDim,
			size_t k,
			const typename Field::Element alpha,
			ka::pointer_r<typename Field::Element> A,
			const size_t lda,
			ka::pointer_r<typename Field::Element> B,
			const size_t ldb,
			const typename Field::Element beta,
			ka::pointer_rw<typename Field::Element> C, const size_t ldc)
	{
		FFLAS::fgemm( F, ta, tb, BlockRowDim, BlockColDim, k, alpha, A.ptr(), lda, B.ptr() , ldb,
			      beta, C.ptr(), ldc);
	}
};

namespace FFLAS {
	
#endif
	
	template<class Field>
	inline typename Field::Element*
	pfgemm( const Field& F,
		const FFLAS::FFLAS_TRANSPOSE ta,
		const FFLAS::FFLAS_TRANSPOSE tb,
		const size_t m,
		const size_t n,
		const size_t k,
		const typename Field::Element alpha,
		const typename Field::Element* A, const size_t lda,
		const typename Field::Element* B, const size_t ldb,
		const typename Field::Element beta,
		typename Field::Element* C, const size_t ldc,
		const size_t w,
		const FFLAS::CuttingStrategy method,
		const int maxThreads
		){
        ForStrategy2D iter(m,n,method,maxThreads);
        for (iter.begin(); ! iter.end(); ++iter){
		TASK(READ(A,B,F), NOWRITE(), READWRITE(C), fgemm, F, ta, tb, iter.iend-iter.ibeg, iter.jend-iter.jbeg, k, alpha, A + iter.ibeg*lda, lda, B +iter.jbeg, ldb, beta, C+ iter.ibeg*ldc+iter.jbeg, ldc, w);
		//			TASK(read(A,B,F, NcolBlocks, NrowBlocks), nowrite(), readwrite(C), fgemm, F, ta, tb, BlockRowDim, BlockColDim, k, alpha, A + RBLOCKSIZE * i*lda, lda, B + CBLOCKSIZE * j, ldb, beta, C+ RBLOCKSIZE*i*ldc+j*CBLOCKSIZE, ldc, w );

        }

		WAIT;
//		BARRIER;

		return C;
	}
	
	template<class Field>
	inline typename Field::Element*
	pfgemm( const Field& F,
		const FFLAS_TRANSPOSE ta,
		const FFLAS_TRANSPOSE tb,
		const size_t m,
		const size_t n,
		const size_t k,
		const typename Field::Element alpha,
		const typename Field::Element* A, const size_t lda,
		const typename Field::Element* B, const size_t ldb,
		const typename Field::Element beta,
		typename Field::Element* C, const size_t ldc,
        const CuttingStrategy method,
        const int maxThreads
            ){
        ForStrategy2D iter(m,n,method,maxThreads);
        for (iter.begin(); ! iter.end(); ++iter){
		TASK(READ(A,B,F), NOWRITE(), READWRITE(C), fgemm, F, ta, tb, iter.iend-iter.ibeg, iter.jend-iter.jbeg, k, alpha, A + iter.ibeg*lda, lda, B +iter.jbeg, ldb, beta, C+ iter.ibeg*ldc+iter.jbeg, ldc);
		/*
#ifdef __FFLASFFPACK_USE_OPENMP
#pragma omp task shared (A, B, C, F)
            fgemm( F, ta, tb, iter.iend-iter.ibeg, iter.jend-iter.jbeg, k, alpha, A + iter.ibeg*lda, lda, B +iter.jbeg, ldb, beta, C+ iter.ibeg*ldc+iter.jbeg, ldc);
#endif
#ifdef __FFLASFFPACK_USE_KAAPI
			ka::Spawn<Taskfgemmw<Field> >()( F, ta, tb, iter.iend-iter.ibeg, iter.jend-iter.jbeg, k, alpha, A + iter.ibeg*lda, lda, B +iter.jbeg, ldb, beta, C+ iter.ibeg*ldc+iter.jbeg, ldc);
#endif
		}
#ifdef __FFLASFFPACK_USE_OPENMP
       #pragma omp taskwait
#endif
#ifdef __FFLASFFPACK_USE_KAAPI
		ka::Sync();
#endif
		*/
	}
	WAIT;
	return C;
	}
	
} // FFLAS                                                                                                                   

#endif // __FFLASFFPACK_fflas_pfgmm_INL  

