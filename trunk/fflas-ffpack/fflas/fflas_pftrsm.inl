/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_pftrsm.inl
 * Copyright (C) 2013 Ziad Sultan
 *
 * Written by Ziad Sultan  < Ziad.Sultan@imag.fr >
 * Time-stamp: <09 Dec 14 10:01:06 Jean-Guillaume.Dumas@imag.fr>
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

#define PTRSM_HYBRID_THRESHOLD 256
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
	inline typename Field::Element_ptr
	ftrsm( const Field& F,
		const FFLAS::FFLAS_SIDE Side,
		const FFLAS::FFLAS_UPLO UpLo,
		const FFLAS::FFLAS_TRANSPOSE TA,
		const FFLAS::FFLAS_DIAG Diag,
		const size_t m,
		const size_t n,
		const typename Field::Element alpha,
#ifdef __FFLAS__TRSM_READONLY
		typename Field::ConstElement_ptr
#else
		typename Field::Element_ptr
#endif
		 A, const size_t lda,
		typename Field::Element_ptr B, const size_t ldb,
		TRSMHelper <StructureHelper::Iterative, ParSeqHelper::Parallel> & H)
		// const FFLAS::CuttingStrategy method,
                // const size_t numThreads)
	{
		if(Side == FflasRight){
			ForStrategy1D iter(m, H.parseq);
			for (iter.begin(); ! iter.end(); ++iter) {
				TRSMHelper<StructureHelper::Recursive, ParSeqHelper::Sequential> SeqH (H);
				TASK(READ(F, A), NOWRITE(), READWRITE(B), ftrsm, F, Side, UpLo, TA, Diag, iter.iend-iter.ibeg, n, alpha, A, lda, B + iter.ibeg*ldb, ldb, SeqH);
			}
		} else {
			ForStrategy1D iter(n, H.parseq);
			for (iter.begin(); ! iter.end(); ++iter) {
				TRSMHelper<StructureHelper::Recursive, ParSeqHelper::Sequential> SeqH (H);
				TASK(READ(F, A), NOWRITE(), READWRITE(B), ftrsm, F, Side, UpLo, TA, Diag, m, iter.iend-iter.ibeg, alpha, A , lda, B + iter.ibeg, ldb, SeqH);
			}
		}
		WAIT;
		return B;
	}
	template<class Field>
	inline typename Field::Element_ptr
	ftrsm( const Field& F,
		const FFLAS::FFLAS_SIDE Side,
		const FFLAS::FFLAS_UPLO UpLo,
		const FFLAS::FFLAS_TRANSPOSE TA,
		const FFLAS::FFLAS_DIAG Diag,
		const size_t m,
		const size_t n,
		const typename Field::Element alpha,
#ifdef __FFLAS__TRSM_READONLY
		typename Field::ConstElement_ptr
#else
		typename Field::Element_ptr
#endif
		 A, const size_t lda,
		typename Field::Element_ptr B, const size_t ldb,
		TRSMHelper <StructureHelper::Hybrid, ParSeqHelper::Parallel> & H)
		// const FFLAS::CuttingStrategy method,
                // const size_t numThreads)
	{
		if(Side == FflasRight){
			int nt = H.parseq.numthreads;
			int nt_it,nt_rec;
			if ((int)m/PTRSM_HYBRID_THRESHOLD < nt){
				nt_it = (int)ceil(double(m)/PTRSM_HYBRID_THRESHOLD);
				nt_rec = (int)ceil(double(nt)/nt_it);
			} else { nt_it = nt; nt_rec = 1;}
			ForStrategy1D iter(m, ParSeqHelper::Parallel((size_t)nt_it,H.parseq.method));
			for (iter.begin(); ! iter.end(); ++iter) {
				ParSeqHelper::Parallel psh(nt_rec,CuttingStrategy::TWO_D_ADAPT);
				TRSMHelper<StructureHelper::Recursive, ParSeqHelper::Parallel> SeqH (psh);
				std::cerr<<"trsm_rec nt = "<<nt_rec<<std::endl;
				TASK(READ(F, A), NOWRITE(), READWRITE(B), ftrsm, F, Side, UpLo, TA, Diag, iter.iend-iter.ibeg, n, alpha, A, lda, B + iter.ibeg*ldb, ldb, SeqH);
			}
		} else {
			int nt = H.parseq.numthreads;
			int nt_it=nt;
			int nt_rec=1;
			while(nt_it*PTRSM_HYBRID_THRESHOLD >= (int)n){
				nt_it>>=1;
				nt_rec<<=1;
			}
			nt_it<<=1;
			nt_rec>>=1;
			// if ((int)n/PTRSM_HYBRID_THRESHOLD < nt){
			// 	nt_it = std::min(nt,(int)ceil(double(n)/PTRSM_HYBRID_THRESHOLD));
			// 	nt_rec = ceil(double(nt)/nt_it);
			// } else { nt_it = nt; nt_rec = 1;}

			ForStrategy1D iter(n, ParSeqHelper::Parallel((size_t)nt_it,H.parseq.method));
			for (iter.begin(); ! iter.end(); ++iter) {
				ParSeqHelper::Parallel psh(nt_rec, CuttingStrategy::TWO_D_ADAPT);
				TRSMHelper<StructureHelper::Recursive, ParSeqHelper::Parallel> SeqH (psh);
				    //std::cerr<<"trsm_rec nt = "<<nt_rec<<std::endl;
				TASK(READ(F, A), NOWRITE(), READWRITE(B), ftrsm, F, Side, UpLo, TA, Diag, m, iter.iend-iter.ibeg, alpha, A , lda, B + iter.ibeg, ldb, SeqH);
			}
		}
		WAIT;
		return B;
	}

} // FFLAS


#endif // __FFLASFFPACK_fflas_pftrsm_INL
