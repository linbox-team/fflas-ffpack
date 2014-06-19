/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_pftrsm.inl
 * Copyright (C) 2013 Ziad Sultan
 *
 * Written by Ziad Sultan  < Ziad.Sultan@imag.fr >
 * Time-stamp: <19 Jun 14 13:39:29 Jean-Guillaume.Dumas@imag.fr>
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
#include "fflas-ffpack/fflas/kaapi_routines.inl"
#endif

#include "fflas-ffpack/fflas/parallel.h"

namespace FFLAS {

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
            ForStrategy1D iter(m, method, numThreads);
            for (iter.begin(); ! iter.end(); ++iter)
				TASK(READ(F, A), NOWRITE(), READWRITE(B), ftrsm, F, Side, UpLo, TA, Diag, iter.iend-iter.ibeg, n, alpha, A, lda, B + iter.ibeg*ldb, ldb);
		} else {
			ForStrategy1D iter(n, method, numThreads);
			for (iter.begin(); ! iter.end(); ++iter) {
				TASK(READ(F, A), NOWRITE(), READWRITE(B), ftrsm, F, Side, UpLo, TA, Diag, m, iter.iend-iter.ibeg, alpha, A , lda, B + iter.ibeg, ldb);
            }
		}
		WAIT;		      
		return B;
	}

} // FFLAS


#endif // __FFLASFFPACK_fflas_pftrsm_INL 
