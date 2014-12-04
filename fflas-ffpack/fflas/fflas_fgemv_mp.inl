/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2014 FFLAS-FFPACK group
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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

#ifndef __FFLASFFPACK_fgemv_mp_INL
#define __FFLASFFPACK_fgemv_mp_INL
// activate only if FFLAS-FFPACK haves multiprecision integer
#ifdef __FFLASFFPACK_HAVE_INTEGER
#include "fflas-ffpack/field/rns-integer-mod.h"

namespace FFLAS {

	// specialization of the fgemv function for the field RNSIntegerMod<rns_double>
	//template<>
	inline FFPACK::rns_double::Element_ptr
	fgemv (const FFPACK::RNSIntegerMod<FFPACK::rns_double>& F, const FFLAS_TRANSPOSE ta,
	       const size_t M, const size_t N,
	       const FFPACK::rns_double::Element alpha,
	       FFPACK::rns_double::ConstElement_ptr A, const size_t lda,
	       FFPACK::rns_double::ConstElement_ptr X, const size_t incX,
	        const FFPACK::rns_double::Element beta,
	       FFPACK::rns_double::Element_ptr Y, const size_t incY)
	{
		if (M!=0 && N !=0){
			for (size_t i=0;i<F.size();i++)
				fgemv(F.rns()._field_rns[i], ta,
				      M, N,
				      alpha._ptr[i*alpha._stride],
				      A._ptr+i*A._stride, lda,
				      X._ptr+i*X._stride, incX,
				      beta._ptr[i*beta._stride],
				      Y._ptr+i*Y._stride, incY
				      );
			size_t Ydim = (ta == FflasNoTrans)?M:N;
			freduce (F, Ydim, Y, incY);
		}
		return Y;
	} 
} // end namespace FFLAS 

#endif // __FFLASFFPACK_HAVE_INTEGER

#endif
