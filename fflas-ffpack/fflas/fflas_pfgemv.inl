/* fflas/fflas_pfgemv.inl
 * Copyright (C) 2019 
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

#ifndef __FFLASFFPACK_fflas_pfgemv_INL
#define __FFLASFFPACK_fflas_pfgemv_INL

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


#include "fflas-ffpack/paladin/pfgemv.inl"

namespace FFLAS {
/*
    template<class Field, class ModeTrait, class Strat, class Param>
    inline typename  std::enable_if<!std::is_same<ModeTrait,ModeCategories::ConvertTo<ElementCategories::RNSElementTag> >::value,typename Field::Element_ptr>::type
    fgemv( const Field& F,
           const FFLAS_TRANSPOSE ta,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr X, const size_t incX,
           const typename Field::Element beta,
           typename Field::Element_ptr Y, const size_t incY,
           MMHelper<Field, MMHelperAlgo::Winograd, ModeTrait, ParSeqHelper::Parallel<Strat,Param> > & H)
    {
        return pfgemv (F, ta, M, N, alpha, A, lda, X, incX, beta, Y, incY, H);
    }
*/
    template<class Field, class Cut, class Param>
    typename Field::Element_ptr
    fgemv(const Field& F,
           const FFLAS_TRANSPOSE ta,
           const size_t m,
           const size_t n,
           const typename Field::Element alpha,
           const typename Field::ConstElement_ptr A, const size_t lda,
           const typename Field::ConstElement_ptr X, const size_t incX,
           const typename Field::Element beta,
           typename Field::Element_ptr Y, const size_t incY,
           ParSeqHelper::Parallel<Cut,Param> parH){

        MMHelper<Field, MMHelperAlgo::Auto, typename FFLAS::ModeTraits<Field>::value, ParSeqHelper::Parallel<Cut,Param> > pH (F,m,n,1,parH);
        return fgemv(F, ta, m, n, alpha, A, lda, X, incX, beta, Y, incY, pH);
    }

    template<class Field>
    typename Field::Element_ptr
    fgemv(const Field& F,
           const FFLAS_TRANSPOSE ta,
           const size_t m,
           const size_t n,
           const typename Field::Element alpha,
           const typename Field::ConstElement_ptr A, const size_t lda,
           const typename Field::ConstElement_ptr X, const size_t incX,
           const typename Field::Element beta,
           typename Field::Element_ptr Y, const size_t incY,
           ParSeqHelper::Sequential seqH ){

        return fgemv(F, ta, m, n, alpha, A, lda, X, incX, beta, Y, incY);
    }

}

#endif // __FFLASFFPACK_fflas_pfgemv_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
