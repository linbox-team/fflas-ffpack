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
#define __FFLASFFPACK_fflas_pfgemm_INL

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

    template<class Field, class ModeTrait, class Strat, class Param>
    inline typename  std::enable_if<!std::is_same<ModeTrait,ModeCategories::ConvertTo<ElementCategories::RNSElementTag> >::value,typename Field::Element_ptr>::type
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
           MMHelper<Field, MMHelperAlgo::Winograd, ModeTrait, ParSeqHelper::Parallel<Strat,Param> > & H)
    {
        return pfgemm (F, ta, tb, m, n, k ,alpha, A, lda, B, ldb, beta, C, ldc, H);
    }

    // template<class Field, class ModeTrait, class Strat, class Param>
    // inline typename  std::enable_if<!std::is_same<ModeTrait,ModeCategories::ConvertTo<ElementCategories::RNSElementTag> >::value,typename Field::Element_ptr>::type
    // fgemm( const Field& F,
    //        const FFLAS::FFLAS_TRANSPOSE ta,
    //        const FFLAS::FFLAS_TRANSPOSE tb,
    //        const size_t m,
    //        const size_t n,
    //        const size_t k,
    //        const typename Field::Element alpha,
    //        typename Field::ConstElement_ptr A, const size_t lda,
    //        typename Field::ConstElement_ptr B, const size_t ldb,
    //        const typename Field::Element beta,
    //        typename Field::Element_ptr C, const size_t ldc,
    //        MMHelper<Field, MMHelperAlgo::WinogradPar, ModeTrait, ParSeqHelper::Parallel<Strat,Param> > & H)
    // {
    // 	std::cerr<<"coucou"<<std::endl;
    // 	return BLAS3::WinoPar(F, ta, tb, m, n, k ,alpha, A, lda, B, ldb, beta, C, ldc, H);
    // }
}

#endif // __FFLASFFPACK_fflas_pfgemm_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
