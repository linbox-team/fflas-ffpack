/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2014 FFLAS-FFPACK group
 *
 * Written by BB <bbboyer@ncsu.edu>
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

#ifndef __FFLASFFPACK_fadd_INL
#define __FFLASFFPACK_fadd_INL


namespace FFLAS { namespace details {

	/**** Specialised ****/

	template <class Field, bool ADD>
	typename std::enable_if<FFLAS::support_simd<typename Field::Element>::value, void>::type
	fadd (const Field & F,  const size_t N,
	      typename Field::ConstElement_ptr A, const size_t inca,
	      typename Field::ConstElement_ptr B, const size_t incb,
	      typename Field::Element_ptr C, const size_t incc
	      , FieldCategories::ModularTag
	     )
	{
		if (inca == 1 && incb == 1 && incc == 1) {
			typename Field::Element p = (typename Field::Element) F.characteristic();
			if (ADD)
				FFLAS::vectorised::addp<!FieldTraits<Field>::balanced>(C,A,B,N,p,F.minElement(),F.maxElement());
			else
				FFLAS::vectorised::subp<!FieldTraits<Field>::balanced>(C,A,B,N,p,F.minElement(),F.maxElement());
		}
		else {
			for (size_t i=0; i<N; i++)
				if (ADD)
					F.add (C[i*incc], A[i*inca], B[i*incb]);
				else
					F.sub (C[i*incc], A[i*inca], B[i*incb]);
		}
	}

	template <class Field, bool ADD>
	typename std::enable_if<!FFLAS::support_simd<typename Field::Element>::value, void>::type
	fadd (const Field & F,  const size_t N,
	      typename Field::ConstElement_ptr A, const size_t inca,
	      typename Field::ConstElement_ptr B, const size_t incb,
	      typename Field::Element_ptr C, const size_t incc
	      , FieldCategories::ModularTag
	     )
	{
		if (inca == 1 && incb == 1 && incc == 1) {
				for (size_t i=0; i<N; i++)
				if (ADD)
					F.add (C[i], A[i], B[i]);
				else
					F.sub (C[i], A[i], B[i]);
		}
		else {
			for (size_t i=0; i<N; i++)
				if (ADD)
				F.add (C[i*incc], A[i*inca], B[i*incb]);
				else
				F.sub (C[i*incc], A[i*inca], B[i*incb]);
		}
	}



	template <class Field, bool ADD>
	void
	fadd (const Field & F,  const size_t N,
	      const double* A, const size_t inca,
	      const double* B, const size_t incb,
	      double* C, const size_t incc
	      , FieldCategories::GenericTag
	      )
	{
		if (inca == 1 && incb == 1 && incc == 1) {
			for (size_t i=0; i<N; i++) {
				if (ADD)
				F.add (C[i], A[i], B[i]);
				else
				F.sub (C[i], A[i], B[i]);
			}
		}
		else {
			for (size_t i=0; i<N; i++)
				if (ADD)
				F.add (C[i*incc], A[i*inca], B[i*incb]);
				else
				F.add (C[i*incc], A[i*inca], B[i*incb]);
		}
	}

	template <class Field, bool ADD>
	void
	fadd (const Field & F,  const size_t N,
	      typename Field::ConstElement_ptr A, const size_t inca,
	      typename Field::ConstElement_ptr B, const size_t incb,
	      typename Field::Element_ptr C, const size_t incc
	      , FieldCategories::UnparametricTag
	     )
	{
		for (size_t i=0; i<N; i++)
			if (ADD)
				C[i] = A[i] + B[i];
			else
				C[i] = A[i] - B[i];
	}



} // details
} // FFLAS


#endif // __FFLASFFPACK_fscal_INL
