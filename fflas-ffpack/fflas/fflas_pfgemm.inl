/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* fflas/fflas_pfgemm.inl
 * Copyright (C) 2013 ??
 *
 * Written by ??
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

#ifdef __FFLAS_USE_OMP
#include <omp.h>
#endif
#ifdef __FFLAS_USE_KAAPI
#include <kaapi++.h>
#endif 

namespace FFLAS {

#ifdef __FFLAS_USE_KAAPI
	class KaapifgemmPar{
		operator (){
			pfgemm();
		}
	class KaapifgemmSeq{
		operator (){
			fgemm();
		}
	};
#endif

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
        const size_t w,
        const CuttingStrategy method = BLOCK_THREADS
        ){

    size_t RBLOCKSIZE, CBLOCKSIZE;
    BlockCuts(RBLOCKSIZE, CBLOCKSIZE, m, n, method);

    size_t NrowBlocks = m/RBLOCKSIZE;
    size_t LastrowBlockSize = m % RBLOCKSIZE;
    if (LastrowBlockSize)
        NrowBlocks++;
    else
        LastrowBlockSize = RBLOCKSIZE;
    size_t NcolBlocks = n/CBLOCKSIZE;
    size_t LastcolBlockSize = n % CBLOCKSIZE;
    if (LastcolBlockSize)
        NcolBlocks++;
    else
        LastcolBlockSize = CBLOCKSIZE;
    
    const size_t BLOCKS = NrowBlocks*NcolBlocks;
    

#pragma omp parallel for default (none) shared (A, B, C, F, RBLOCKSIZE, CBLOCKSIZE, NcolBlocks, NrowBlocks, LastcolBlockSize, LastrowBlockSize)
    for (size_t t = 0; t < BLOCKS; ++t){
        size_t i = t / NcolBlocks;
        size_t j = t % NcolBlocks;
        size_t BlockRowDim = RBLOCKSIZE;
        if (i == NrowBlocks-1)
            BlockRowDim = LastrowBlockSize;
        size_t BlockColDim = CBLOCKSIZE;
        if (j == NcolBlocks-1)
            BlockColDim = LastcolBlockSize;

#ifdef __FFLAS_USE_OMP
        fgemm( F, ta, tb, BlockRowDim, BlockColDim, k, alpha, A + RBLOCKSIZE * i*lda, lda, B + CBLOCKSIZE * j, ldb, beta, C+ RBLOCKSIZE*i*ldc+j*CBLOCKSIZE, ldc, w);
#endif
#ifdef __FFLAS_USE_KAAPI
	spawn<fgemm>
#endif
    }
    return C;
}

} // FFLAS

#endif // __FFLASFFPACK_fflas_pfgmm_INL
