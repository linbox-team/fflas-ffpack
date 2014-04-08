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

#ifdef __FFLASFFPACK_USE_OPENMP
#include <omp.h>

namespace FFLAS {

	template<class Field>
	typename Field::Element*
	pftrsm( const Field& F,
		     const FFLAS_SIDE Side,
		     const FFLAS_UPLO UpLo,
		     const FFLAS_TRANSPOSE TA,
		     const FFLAS_DIAG Diag,
		     const size_t m,
		     const size_t n,
		     const typename Field::Element alpha,
		     const typename Field::Element* A, const size_t lda,
		     typename Field::Element* B, const size_t ldb)
	{
		int numthreads = omp_get_max_threads();
		if(Side == FflasRight){
			size_t BLOCKSIZE=std::max(2*m/numthreads,(size_t)1); // There is always 2 TRSM taking place in parallel
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
                                #pragma omp task shared (A, B, F)
				ftrsm(F, Side, UpLo, TA, Diag, BlockDim, n, alpha, A , lda, B + BLOCKSIZE*i*ldb, ldb);
			}
		}
		else {
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
		}
                #pragma omp taskwait

		return B;
	}


} // FFLAS
#else
namespace FFLAS {

	template<class Field>
	typename Field::Element*
	pftrsm( const Field& F,
		     const FFLAS_SIDE Side,
		     const FFLAS_UPLO UpLo,
		     const FFLAS_TRANSPOSE TA,
		     const FFLAS_DIAG Diag,
		     const size_t m,
		     const size_t n,
		     const typename Field::Element alpha,
		     const typename Field::Element* A, const size_t lda,
		     typename Field::Element* B, const size_t ldb)
	{
		return ftrsm(F,Side,UpLo,TA,Diag,m,n,alpha,A,lda,B,ldb);
	}
} //FFLAS

#endif // __FFLASFFPACK_USE_OPENMP


