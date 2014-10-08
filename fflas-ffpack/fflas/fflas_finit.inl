/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* fflas/fflas_finit.inl
 * Copyright (C) 2014 Pascal Giorgi
 *
 * Written by Pascal Giorgi <Pascal.Giorgi@lirmm.fr>
 * BB<bboyer@ncsu.edu>
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
#ifndef __FFLASFFPACK_fflas_init_INL
#define __FFLASFFPACK_fflas_init_INL


#include "fflas-ffpack/fflas/fflas_fcopy.h"

#define FFLASFFPACK_COPY_INIT 100

namespace FFLAS { namespace details {


	// specialised
	template<class Field>
	typename std::enable_if<FFLAS::support_simd_mod<typename Field::Element>::value, void>::type
	finit (const Field & F, const size_t m,
	       typename Field::Element_ptr A, const size_t incX
	       , FieldCategories::ModularTag
	      )
	{
		if(incX == 1) {
			typename Field::Element p;
			p=(typename Field::Element)F.cardinality();
			typename Field::Element invp=1/(typename Field::Element) p;
			vectorised::modp<!FieldTraits<Field>::balanced,false>(A,A,m,p,invp,F.minElement(),F.maxElement());
		}
		else { /*  faster with copy, use incX=1, copy back ? */
			if (m < FFLASFFPACK_COPY_INIT) {
				typename Field::Element_ptr Xi = A ;
				for (; Xi < A+m*incX; Xi+=incX )
					F.init( *Xi , *Xi);
			}
			else {
				typename Field::Element_ptr Ac = fflas_new (F,m,1) ;
				fcopy(F,m,A,incX,Ac,1);
				finit(F,m,Ac,1,FieldCategories::ModularTag());
				fcopy(F,m,Ac,1,A,incX);
				fflas_delete (Ac);
			}
		}
	}

	template<class Field>
	typename std::enable_if<!FFLAS::support_simd_mod<typename Field::Element>::value, void>::type
	finit (const Field & F, const size_t m,
	       typename Field::Element_ptr A, const size_t incX
	       , FieldCategories::ModularTag
	      )
	{ /* ??? ( faster with copy, use incX=1, copy back ? */
		typename Field::Element_ptr  Xi = A ;
		for (; Xi < A+m*incX; Xi+=incX )
			F.init( *Xi , *Xi);
	}

	template<class Field>
	void
	finit (const Field & F, const size_t m,
	       typename Field::Element_ptr A, const size_t incX
	       , FieldCategories::GenericTag
	      )
	{
		typename Field::Element_ptr Xi = A ;
		for (; Xi < A+m*incX; Xi+=incX )
			F.init( *Xi , *Xi);
	}

	template<class Field>
	void
	finit (const Field & F, const size_t m,
	       typename Field::Element_ptr A, const size_t incX
	       , FieldCategories::UnparametricTag
	      )
	{
		return;
	}

	template<class Field>
	typename std::enable_if<FFLAS::support_simd_mod<typename Field::Element>::value, void>::type
	finit (const Field & F, const size_t m,
	       typename Field::ConstElement_ptr  B, const size_t incY,
	       typename Field::Element_ptr A, const size_t incX
	       , FieldCategories::ModularTag
	       )
	{
		if(incX == 1 && incY == 1) {
			typename Field::Element p;
			p=(typename Field::Element)F.cardinality();
			double invp=1./(double)p;
			vectorised::modp<!FieldTraits<Field>::balanced,false>(A,B,m,p,invp,F.minElement(),F.maxElement());
		}
		else {
			typename Field::Element_ptr Xi = A ;
			typename Field::ConstElement_ptr Yi = B ;
			for (; Xi < A+m*incX; Xi+=incX, Yi += incY )
				F.init( *Xi , *Yi);
		}
	}

	template<class Field>
	typename std::enable_if<!FFLAS::support_simd_mod<typename Field::Element>::value, void>::type
	finit (const Field & F, const size_t m,
	       typename Field::ConstElement_ptr  B, const size_t incY,
	       typename Field::Element_ptr A, const size_t incX
	       , FieldCategories::ModularTag
	       )
	{

		typename Field::Element_ptr Xi = A ;
		typename Field::ConstElement_ptr Yi = B ;
		for (; Xi < A+m*incX; Xi+=incX, Yi += incY )
			F.init( *Xi , *Yi);
	}

	template<class Field>
	void
	finit (const Field & F, const size_t m,
	       typename Field::ConstElement_ptr  B, const size_t incY,
	       typename Field::Element_ptr A, const size_t incX
	       , FieldCategories::GenericTag
	       )
	{

		typename Field::Element_ptr Xi = A ;
		typename Field::ConstElement_ptr Yi = B ;
		for (; Xi < A+m*incX; Xi+=incX, Yi += incY )
			F.init( *Xi , *Yi);
	}

	template<class Field>
	void
	finit (const Field & F, const size_t m,
	       typename Field::ConstElement_ptr  B, const size_t incY,
	       typename Field::Element_ptr A, const size_t incX
	       , FieldCategories::UnparametricTag
	      )
	{
		return;
	}


} // details
} // FFLAS

#endif // __FFLASFFPACK_fflas_init_INL

