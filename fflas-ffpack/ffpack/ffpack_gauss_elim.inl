/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* ffpack.inl
 * Copyright (C) 2014 FFLAS-FFACK group
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

#ifndef __FFLASFFPACK_ffpack_INL
#define __FFLASFFPACK_ffpack_INL

#include "fflas-ffpack/fflas/fflas.h"
using namespace FFLAS;
//#include  "mpi.h"
//#include "fflas-ffpack/paladin/pfgemv.inl"


namespace FFPACK {


template <class Field>
typename Field::Element_ptr
gauss_elim( const Field& F, const size_t M, const size_t k,
	   typename Field::Element_ptr A, const size_t lda,
	   typename Field::Element_ptr x, const int incx,
	   typename Field::ConstElement_ptr b, const int incb )
{

	size_t *P = FFLAS::fflas_new<size_t>(M);
	size_t *rowP = FFLAS::fflas_new<size_t>(M);

	if (LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, M, A, lda, P, rowP) < M){
		std::cerr<<"SINGULAR MATRIX"<<std::endl;
	}
	else{

		typename Field::Element_ptr r, rr;
		r = FFLAS::fflas_new(F,M,incx);
		FFLAS::fassign( F, M, r, incx, b, incb );
		for(size_t i=1; i<k; i++){
			FFLAS::fgemv(F, FFLAS::FflasNoTrans, M, k, F.one, A, lda, x, incx, F.zero, x, incx);
			fsub (F, M, incx, r, incx, x, incx, r, incx);
			F.assign(x, fdot (F, M, b, incb, r, incx));			
		}

	}
		FFLAS::fflas_delete( P);
		FFLAS::fflas_delete( rowP);
		return x;
}



} // FFPACK

#endif // __FFLASFFPACK_ffpack_INL
