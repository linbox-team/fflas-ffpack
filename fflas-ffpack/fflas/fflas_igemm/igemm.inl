/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2013,2014  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 * the code is inspired and adapted from the Eigen library
 * modified by BB <brice.boyer@lip6.fr>
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

#ifndef __FFLASFFPACK_fflas_igemm_igemm_INL
#define __FFLASFFPACK_fflas_igemm_igemm_INL

#include "fflas-ffpack/utils/fflas_memory.h"

namespace FFLAS { namespace Protected {

	// Assume matrices A,B,C are stored in column major order
	void igemm_colmajor(size_t rows, size_t cols, size_t depth, int64_t* C, size_t ldc, int64_t* A, size_t lda, int64_t* B, size_t ldb)
	{

		using simd = Simd<int64_t> ;
		size_t mc,kc,nc;
		mc=rows;
		nc=cols;
		kc=depth;
		FFLAS::details::BlockingFactor(mc,nc,kc);
		size_t sizeA = mc*kc;
		size_t sizeB = kc*cols;
		size_t sizeW = simd::vect_size*kc*_nr; // store data duplicated by the number of elements fitting in vector register

		// these data must be (simd::alignment byte aligned
		int64_t *blockA, *blockB, *blockW;


		blockA = fflas_new<int64_t>(sizeA, (Alignment)simd::alignment);
		blockB = fflas_new<int64_t>(sizeB, (Alignment)simd::alignment);
		blockW = fflas_new<int64_t>(sizeW, (Alignment)simd::alignment);


		// For each horizontal panel of B, and corresponding vertical panel of A
		for(size_t k2=0; k2<depth; k2+=kc){

			const size_t actual_kc = std::min(k2+kc,depth)-k2;
			FFLASFFPACK_check(kc <= depth);

			// pack horizontal panel of B into sequential memory (L2 cache)
			FFLAS::details::pack_rhs<_nr>(blockB, B+k2, ldb, actual_kc, cols);

			// For each mc x kc block of the lhs's vertical panel...
			for(size_t i2=0; i2<rows; i2+=mc){

				const size_t actual_mc = std::min(i2+mc,rows)-i2;


				FFLASFFPACK_check(mc <= rows);
				// pack a chunk of the vertical panel of A into a sequential memory (L1 cache)
				FFLAS::details::pack_lhs<_mr>(blockA, A+i2+k2*lda, lda, actual_mc, actual_kc);
				// call block*panel kernel
				FFLAS::details::igebp(actual_mc, cols, actual_kc, C+i2, ldc, blockA, actual_kc, blockB, actual_kc, blockW);
			}
		}

		fflas_delete(blockA);
		fflas_delete(blockB);
		fflas_delete(blockW);

	}

	void igemm(size_t rows, size_t cols, size_t depth, int64_t* C, size_t ldc, int64_t* A, size_t lda, int64_t* B, size_t ldb)
	{
		igemm_colmajor(rows, cols, depth, C, ldc, A, lda, B, ldb);
	}

} // Protected
} // FFLAS


// igemm

namespace FFLAS {
	inline void igemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
			   const int M, const int N, const int K,
			   const int64_t alpha, const int64_t *A, const int lda, const int64_t *B, const int ldb,
			   const int64_t beta, int64_t *C, const int ldc)
	{
		if (beta != 0 || alpha != 1) {
			std::cout << __func__ << " alpha, beta are 1 and 0" << std::endl;
			exit(-1);
		}
		if (TransA !=  CblasNoTrans || TransB !=  CblasNoTrans) {
			std::cout << __func__ << " transpose not surpported ";
			exit(-1);
		}
		if (Order != CblasColMajor) {
			std::cout << __func__ << " transpose not surpported ";
			exit(-1);
		}
		Protected::igemm(M,N,K,C,ldc,const_cast<int64_t*>(A),lda,const_cast<int64_t*>(B),ldb);
	}


} // FFLAS

#endif // __FFLASFFPACK_fflas_igemm_igemm_INL

