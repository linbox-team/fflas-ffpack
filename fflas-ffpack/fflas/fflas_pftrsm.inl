/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_pftrsm.inl
 * Copyright (C) 2013 Ziad Sultan
 *
 * Written by Ziad Sultan  < Ziad.Sultan@imag.fr >
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


#ifndef __FFLASFFPACK_fflas_pftrsm_INL
#define __FFLASFFPACK_fflas_pftrsm_INL


#ifdef __FFLASFFPACK_USE_OPENMP
#include <omp.h>
#endif
#ifdef __FFLASFFPACK_USE_KAAPI
#include <kaapi++>
#endif

#include "fflas-ffpack/fflas/parallel.h"

namespace FFLAS {


#ifdef __FFLASFFPACK_USE_KAAPI

	template<class Field>
	struct Taskftrsm: public ka::Task<12>::Signature<
		Field , /* Field F */
		FFLAS::FFLAS_SIDE ,
		FFLAS::FFLAS_UPLO ,
		FFLAS::FFLAS_TRANSPOSE ,
		FFLAS::FFLAS_DIAG ,
		size_t ,   /* size : M */
		size_t ,   /* size : N */
		typename Field::Element ,
		ka::R<typename Field::Element >, /* Matrix A */
		size_t , /* lda */
		ka::RW<typename Field::Element >, /* Matrix B */
		size_t  /* ldb */
		>{};
}

template<class Field>
struct TaskBodyCPU<FFLAS::Taskftrsm<Field> > {
	void operator()(const Field & F, const FFLAS::FFLAS_SIDE Side,
			const FFLAS::FFLAS_UPLO Uplo,
			const FFLAS::FFLAS_TRANSPOSE TransA,
			const FFLAS::FFLAS_DIAG Diag,
			const size_t M, const size_t N,
			const typename Field::Element alpha,
			ka::pointer_r<typename Field::Element > A, const size_t lda,
			ka::pointer_rw<typename Field::Element > B, const size_t ldb )
	{

		ftrsm(F, Side, Uplo, TransA, Diag, M, N, alpha, A.ptr(), lda, B.ptr(), ldb);
	}
};

	namespace FFLAS {

#endif
	//#ifdef __FFLASFFPACK_USE_OPENMP
	template<class Field>
	inline typename Field::Element*
	pftrsm( const Field& F,
		const FFLAS::FFLAS_SIDE Side,
		const FFLAS::FFLAS_UPLO UpLo,
		const FFLAS::FFLAS_TRANSPOSE TA,
		const FFLAS::FFLAS_DIAG Diag,
		const size_t m,
		const size_t n,
		const typename Field::Element alpha,
#ifdef __FFLAS__TRSM_READONLY
		const
#endif
		typename Field::Element* A, const size_t lda,
		typename Field::Element* B, const size_t ldb,
		const FFLAS::CuttingStrategy method,
                const size_t numThreads)
	{
		if(Side == FflasRight){
			/*
			size_t BLOCKSIZE=std::max(m/numThreads,(size_t)1); // There is always 2 TRSM taking place in parallel
			size_t NBlocks = m/BLOCKSIZE;
			size_t LastBlockSize = m % BLOCKSIZE;
			if (LastBlockSize)
				NBlocks++;
			else
				LastBlockSize=BLOCKSIZE;
			//#pragma omp parallel for default (none) shared(A,B,F,NBlocks, LastBlockSize, BLOCKSIZE)
			for (size_t t = 0; t < NBlocks; ++t){
				size_t i = t % NBlocks;
				size_t BlockDim = BLOCKSIZE;
				if (i == NBlocks-1)
					BlockDim = LastBlockSize;

				std::cout<<" BlockDim : "<<BlockDim<<std::endl;
                                #pragma omp task shared (A, B, F)
				ftrsm(F, Side, UpLo, TA, Diag, BlockDim, n, alpha, A , lda, B + BLOCKSIZE*i*ldb, ldb);
				}
							
			*/
		ForStrategy1D iter(m, method, numThreads);
		for (iter.begin(); ! iter.end(); ++iter)
			{
				TASK(READ(F, A), NOWRITE(), READWRITE(B), ftrsm, F, Side, UpLo, TA, Diag, iter.iend-iter.ibeg, n, alpha, A, lda, B + iter.ibeg*ldb, ldb);
				//	std::cout<<" iter.iend-iter.ibeg : "<<iter.iend-iter.ibeg<<std::endl;
			//	ftrsm( F, Side, UpLo, TA, Diag, iter.iend-iter.ibeg, n, alpha, A, lda,  B + iter.ibeg*ldb, ldb);
			}

		}
		else {
			/*
			size_t BLOCKSIZE=std::max(2*n/numthreads,(size_t)1); // There is always 2 TRSM taking place in parallel
			size_t NBlocks = n/BLOCKSIZE;
			size_t LastBlockSize = n % BLOCKSIZE;
			if (LastBlockSize)
				NBlocks++;
			else
				LastBlockSize=BLOCKSIZE;
			//#pragma omp parallel for default (none) shared(A,B,F,NBlocks, LastBlockSize, BLOCKSIZE)
			for (size_t t = 0; t < NBlocks; ++t)
				{
					size_t j = t % NBlocks;
					size_t BlockDim = BLOCKSIZE;
					if (j == NBlocks-1)
						BlockDim = LastBlockSize;
                                        #pragma omp task shared (A, B, F)
					ftrsm(F, Side, UpLo, TA, Diag, m, BlockDim, alpha, A , lda, B + BLOCKSIZE*j, ldb);
					}
			*/
			ForStrategy1D iter(n, method, numThreads);
			for (iter.begin(); ! iter.end(); ++iter)
				TASK(READ(F, A), NOWRITE(), READWRITE(B), ftrsm, F, Side, UpLo, TA, Diag, iter.iend-iter.ibeg, n, alpha, A , lda, B + iter.ibeg, ldb);
			//			WAIT;
			//#pragma omp taskwait
		}
		WAIT;		      
		return B;
	}

} // FFLAS

// #else
// namespace FFLAS {

// 	template<class Field>
// 	typename Field::Element*
// 	pftrsm( const Field& F,
// 		const FFLAS::FFLAS_SIDE Side,
// 		const FFLAS::FFLAS_UPLO UpLo,
// 		const FFLAS::FFLAS_TRANSPOSE TA,
// 		const FFLAS::FFLAS_DIAG Diag,
// 		const size_t m,
// 		const size_t n,
// 		const typename Field::Element alpha,
// 		const typename Field::Element* A, const size_t lda,
// 		typename Field::Element* B, const size_t ldb,
// 		const FFLAS::CuttingStrategy method,
//                 const int maxThreads)
// 	{
// 		return ftrsm(F,Side,UpLo,TA,Diag,m,n,alpha,A,lda,B,ldb);
// 	}
// } //FFLAS

// #endif // __FFLASFFPACK_USE_OPENMP

#endif // __FFLASFFPACK_fflas_pftrsm_INL 
