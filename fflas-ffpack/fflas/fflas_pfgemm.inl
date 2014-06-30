/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_pfgemm.inl
 * Copyright (C) 2013 Jean Guillaume Dumas Clement Pernet Ziad Sultan
 *
 * Written by Jean Guillaume Dumas Clement Pernet Ziad Sultan
 * Time-stamp: <30 Jun 14 15:30:43 Jean-Guillaume.Dumas@imag.fr>
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


	template<class Field, class AlgoT, class FieldTrait>
	inline typename Field::Element*
	fgemm( const Field& F,
		const FFLAS::FFLAS_TRANSPOSE ta,
		const FFLAS::FFLAS_TRANSPOSE tb,
		const size_t m,
		const size_t n,
		const size_t k,
		const typename Field::Element alpha,
		typename Field::Element* A, const size_t lda,
		typename Field::Element* B, const size_t ldb,
		const typename Field::Element beta,
		typename Field::Element* C, const size_t ldc,
		MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel> & H) 
	{
		ForStrategy2D iter(m,n,H.parseq.method,H.parseq.numthreads);
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

