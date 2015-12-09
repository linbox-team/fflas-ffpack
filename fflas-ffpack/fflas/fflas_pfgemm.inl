/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_pfgemm.inl
 * Copyright (C) 2013 Jean Guillaume Dumas Clement Pernet Ziad Sultan
 *
 * Written by Jean Guillaume Dumas Clement Pernet Ziad Sultan
 * Time-stamp: <27 Nov 15 14:07:46 Jean-Guillaume.Dumas@imag.fr>
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

#ifndef __FFLASFFPACK_fflas_pfgemm_INL
#define __FFLASFFPACK_fflas_pgemm_INL

#define __FFLASFFPACK_SEQPARTHRESHOLD 220
#define __FFLASFFPACK_DIMKPENALTY 1

#ifdef __FFLASFFPACK_USE_KAAPI
#include <kaapi++>
#include "fflas-ffpack/fflas/kaapi_routines.inl"
#endif
#ifdef __FFLASFFPACK_USE_OPENMP
#include <omp.h>
#endif


#include "fflas-ffpack/paladin/blockcuts.inl"
#include "fflas-ffpack/paladin/parallel.h"
#include "fflas-ffpack/utils/timer.h"


#include "fflas-ffpack/paladin/pfgemm_variants.inl"

namespace FFLAS {

	template<class Field, class AlgoT, class ModeTrait, class Strat, class Param>
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
	       MMHelper<Field, AlgoT, ModeTrait, ParSeqHelper::Parallel<Strat,Param> > & H) 
	{

        if ((ta != FFLAS::FflasNoTrans) || (tb != FFLAS::FflasNoTrans)) {
            std::cerr << "*** ERROR ***: pfgemm ^T NOT YET IMPLEMENTED" << std::endl;
            return C;
        }
        
	pfgemm (F, ta, tb, m, n, k ,alpha, A, lda, B, ldb, beta, C, ldc, H);
	
        // if (H.parseq.method() == RECURSIVE) {
        //     switch (H.parseq.strategy()){
        //         case StrategyParameter::THREE_D_ADAPT: // Splitting 1 dimension at a time recursively: the largest one
        //             return pfgemm_3D_rec_adapt (F, ta, tb, m, n, k ,alpha, A, lda, B, ldb, beta, C, ldc, H);
        //         case StrategyParameter::TWO_D_ADAPT: // Splitting 1 dimension at a time recursively: the largest one
        //             return pfgemm_2D_rec_adapt (F, ta, tb, m, n, k ,alpha, A, lda, B, ldb, beta, C, ldc, H);
        //         case StrategyParameter::TWO_D: // Splitting the outer dimensions m and n recursively
        //             return pfgemm_2D_rec (F, ta, tb, m, n, k ,alpha, A, lda, B, ldb, beta, C, ldc, H);
        //         case StrategyParameter::THREE_D_INPLACE: // Splitting the three dimensions recursively, without temp and with synchro
        //             return pfgemm_3D_rec_V2(F, ta, tb, m, n, k ,alpha, A, lda, B, ldb, beta, C, ldc, H);
        //         case StrategyParameter::THREE_D: // Splitting the three dimensions recursively, with temp alloc and fewer synchro
        //             return pfgemm_3D_rec2_V2(F, ta, tb, m, n, k ,alpha, A, lda, B, ldb, beta, C, ldc, H);
        //         default:
        //             return pfgemm_2D_rec_adapt (F, ta, tb, m, n, k ,alpha, A, lda, B, ldb, beta, C, ldc, H);
        //     }
        // } else {
        //     H.parseq.set_numthreads( std::min(H.parseq.numthreads(), std::max((size_t)1,(size_t)(m*n/(__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD)))) );
                
//             MMHelper<Field, AlgoT, ModeTrait, ParSeqHelper::Sequential> SeqH (H);
//             SYNCH_GROUP( 
//             {FORBLOCK2D(iter,m,n,H.parseq,
//                    TASK( MODE(
//                        READ(A[iter.ibegin()*lda],B[iter.jbegin()]) 
//                        CONSTREFERENCE(F, SeqH) 
//                        READWRITE(C[iter.ibegin()*ldc+iter.jbegin()])), 
//                          fgemm( F, ta, tb, iter.iend()-iter.ibegin(), iter.jend()-iter.jbegin(), k, alpha, A+iter.ibegin()*lda, lda, B+iter.jbegin(), ldb, beta, C+iter.ibegin()*ldc+iter.jbegin(), ldc, SeqH);
// 			     //	 std::cout<<" "<<iter.iend()-iter.ibegin()<<std::endl;
// );
//                    );
//             });
	return C;
	}
}

#endif // __FFLASFFPACK_fflas_pfgemm_INL

