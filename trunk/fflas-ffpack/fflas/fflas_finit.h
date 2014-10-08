/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_finit.inl
 * Copyright (C) 2014 FFLAS FFPACK group
 *
 * Written by  BB<bboyer@ncsu.edu>
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

#ifndef __FFLASFFPACK_fflas_init_H
#define __FFLASFFPACK_fflas_init_H

#include "fflas-ffpack/fflas/fflas_simd_functions.h"
#include "fflas-ffpack/field/field-traits.h"
#include "fflas_finit.inl"

namespace FFLAS {

	/***************************/
	/*         LEVEL 1         */
	/***************************/

	template<class Field>
	void
	finit (const Field & F, const size_t m,
	       typename Field::ConstElement_ptr  B, const size_t incY,
	       typename Field::Element_ptr A, const size_t incX
	      )
	{
		return details::finit(F,m,B,incY,A,incX,typename FieldTraits<Field>::category());
	}

	template<class Field>
	void
	finit (const Field & F, const size_t m,
	       typename Field::Element_ptr A, const size_t incX
	      )
	{
		return details::finit(F,m,A,incX,typename FieldTraits<Field>::category());
	}

	// OOOPS
	template<class Field, class ConstOtherElement_ptr>
	void
	finit (const Field& F, const size_t n,
	       ConstOtherElement_ptr Y, const size_t incY,
	       typename Field::Element_ptr X, const size_t incX)
	{
		typename Field::Element_ptr Xi = X ;
		ConstOtherElement_ptr Yi = Y ;

		if (incX == 1 && incY == 1)
			for (; Yi < Y + n ; ++Xi, ++Yi)
				F.init( *Xi , *Yi);
		else
			for (; Yi < Y+n*incY; Xi+=incX, Yi += incY )
				F.init( *Xi , *Yi);
	}


	/***************************/
	/*         LEVEL 2         */
	/***************************/


	template<class Field>
	void
	finit (const Field& F, const size_t m , const size_t n,
	       typename Field::Element_ptr A, const size_t lda)
	{
		if (n == lda)
			finit(F,n*m,A,1);
		else
			for (size_t i = 0 ; i < m ; ++i)
				finit(F,n,A+i*lda,1);
		return;
	}

	template<class Field, class OtherElement_ptr>
	void
	finit (const Field& F, const size_t m , const size_t n,
	       const OtherElement_ptr B, const size_t ldb,
	       typename Field::Element_ptr A, const size_t lda)
	{
		if (n == lda && n == ldb)
			finit(F,n*m,B,1,A,1);
		else
			for (size_t i = 0 ; i < m ; ++i)
				finit(F,n,B+i*ldb,1,A+i*lda,1);
		return;
	}

} // end of namespace FFLAS

#endif // __FFLASFFPACK_fflas_init_H
