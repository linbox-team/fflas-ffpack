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

#include "fflas-ffpack/fflas/fflas_simd_functions.h"

namespace FFLAS {

	/***************************/
	/*         LEVEL 1         */
	/***************************/


	/**** Specialised ****/


	template <>
	void
	fadd (const FFPACK:: Modular<double> & F,  const size_t N,
	      const double* A, const size_t inca,
	      const double* B, const size_t incb,
	      double* C, const size_t incc)
	{
		if (inca == 1 && incb == 1 && incc == 1) {
			double p = (double)F.characteristic();
			vectorised::addp<true>(C,A,B,N,p,0,p-1);

		}
		else {
			for (size_t i=0; i<N; i++)
				F.add (C[i*incc], A[i*inca], B[i*incb]);
		}
	}


	template <>
	void
	fadd (const FFPACK:: ModularBalanced<double> & F,  const size_t N,
	      const double* A, const size_t inca,
	      const double* B, const size_t incb,
	      double* C, const size_t incc)
	{
		if (inca == 1 && incb == 1 && incc == 1) {
			double p = (double)F.characteristic();
			vectorised::addp<false>(C,A,B,N,p,F.minElement(), F.maxElement());
		}
		else {
			for (size_t i=0; i<N; i++)
				F.add (C[i*incc], A[i*inca], B[i*incb]);
		}
	}

	template <>
	void
	fadd (const FFPACK:: Modular<float> & F,  const size_t N,
	      const float* A, const size_t inca,
	      const float* B, const size_t incb,
	      float* C, const size_t incc)
	{
		if (inca == 1 && incb == 1 && incc == 1) {
			float p = (float)F.characteristic();
			vectorised::addp<true>(C,A,B,N,p,0,p-1);
		}
		else {
			for (size_t i=0; i<N; i++)
				F.add (C[i*incc], A[i*inca], B[i*incb]);
		}
	}

	template <>
	void
	fsub (const FFPACK:: Modular<double> & F,  const size_t N,
	      const double* A, const size_t inca,
	      const double* B, const size_t incb,
	      double* C, const size_t incc)
	{
		if (inca == 1 && incb == 1 && incc == 1) {
			double p = (double)F.characteristic();
			vectorised::subp<true>(C,A,B,N,p,0,p-1);

		}
		else {
			for (size_t i=0; i<N; i++)
				F.sub (C[i*incc], A[i*inca], B[i*incb]);
		}
	}


	template <>
	void
	fsub (const FFPACK:: ModularBalanced<double> & F,  const size_t N,
	      const double* A, const size_t inca,
	      const double* B, const size_t incb,
	      double* C, const size_t incc)
	{
		if (inca == 1 && incb == 1 && incc == 1) {
			double p = (double)F.characteristic();
			vectorised::subp<false>(C,A,B,N,p,F.minElement(), F.maxElement());

		}
		else {
			for (size_t i=0; i<N; i++)
				F.sub (C[i*incc], A[i*inca], B[i*incb]);
		}
	}

	template <>
	void
	fsub (const FFPACK:: Modular<float> & F,  const size_t N,
	      const float* A, const size_t inca,
	      const float* B, const size_t incb,
	      float* C, const size_t incc)
	{
		if (inca == 1 && incb == 1 && incc == 1) {
			float p = (float)F.characteristic();
			vectorised::subp<true>(C,A,B,N,p,0,p-1);

		}
		else {
			for (size_t i=0; i<N; i++)
				F.sub (C[i*incc], A[i*inca], B[i*incb]);
		}
	}


	template <>
	void
	fsub (const FFPACK:: ModularBalanced<float> & F,  const size_t N,
	      const float* A, const size_t inca,
	      const float* B, const size_t incb,
	      float* C, const size_t incc)
	{
		if (inca == 1 && incb == 1 && incc == 1) {
			float p = (float)F.characteristic();
			vectorised::subp<false>(C,A,B,N,p,F.minElement(), F.maxElement());

		}
		else {
			for (size_t i=0; i<N; i++)
				F.sub (C[i*incc], A[i*inca], B[i*incb]);
		}
	}

#if 0
	template <>
	void
	fadd (const FFPACK:: UnparametricField<float> & F,  const size_t N,
	      const float* A, const size_t inca,
	      const float* B, const size_t incb,
	      float* C, const size_t incc)
	{
		if (inca == 1 && incb == 1 && incc == 1) {
			vectorised::addp(C,A,B,N);
		}
		else
			for (size_t i=0; i<N; i++)
				F.add (C[i*incc], A[i*inca], B[i*incb]);
	}
#endif

	/****   Generic   ****/

	template <class Field>
	void
	fadd (const Field& F,  const size_t N,
	      typename Field::ConstElement_ptr A, const size_t inca,
	      typename Field::ConstElement_ptr B, const size_t incb,
	      typename Field::Element_ptr C, const size_t incc)
	{
		if (inca == 1 && incb == 1 && incc == 1)
			for (size_t i=0; i<N; i++)
				F.add (C[i], A[i], B[i]);
		else
			for (size_t i=0; i<N; i++)
				F.add (C[i*incc], A[i*inca], B[i*incb]);
	}

	template <class Field>
	void
	faddin (const Field& F,  const size_t N,
		typename Field::ConstElement_ptr B, const size_t incb,
		typename Field::Element_ptr C, const size_t incc)
	{
		fadd(F,N,B,incb,C,incc,C,incc);
		return;

#if 0
		if (incb == 1 && incc == 1)
			for (size_t i=0; i<N; i++)
				F.addin (C[i], B[i]);
		else
			for (size_t i=0; i<N; i++)
				F.addin (C[i*incc], B[i*incb]);
#endif
	}

	template <class Field>
	void
	fsub (const Field& F,  const size_t N,
	      typename Field::ConstElement_ptr A, const size_t inca,
	      typename Field::ConstElement_ptr B, const size_t incb,
	      typename Field::Element_ptr C, const size_t incc)
	{
		if (inca == 1 && incb == 1 && incc == 1)
			for (size_t i=0; i<N; i++)
				F.sub (C[i], A[i], B[i]);
		else
			for (size_t i=0; i<N; i++)
				F.sub (C[i*incc], A[i*inca], B[i*incb]);

	}

	template <class Field>
	void
	fsubin (const Field& F,  const size_t N,
		typename Field::ConstElement_ptr B, const size_t incb,
		typename Field::Element_ptr C, const size_t incc)
	{

		fsub(F,N,C,incc,B,incb,C,incc);
	}


	// C = A + a B
	template <class Field>
	void
	fadd (const Field& F, const size_t N,
	      typename Field::ConstElement_ptr A, const size_t inca,
	      const typename Field::Element alpha,
	      typename Field::ConstElement_ptr B, const size_t incb,
	      typename Field::Element_ptr C, const size_t incc)
	{
		if (C == A && inca == incc)
			return faxpy(F,N,alpha,B,incb,C,incc);
		if (F.isOne(alpha))
			return fadd(F,N,A,inca,B,incb,C,incc);
		if (F.isMOne(alpha)){
			return fsub(F,N,A,inca,B,incb,C,incc);
		}
		if (F.isZero(alpha))
			return fcopy(F,N,A,inca,C,incc);

		if (inca == 1 && incb == 1 && incc == 1) {
			for (size_t i = 0 ; i < N ; ++i) {
				//!@todo optimise here
				F.mul(C[i],alpha,B[i]);
				F.addin(C[i],A[i]);
			}
			return;
		}

		typename Field::ConstElement_ptr Ai = A, Bi = B;
		typename Field::Element_ptr Ci = C;
		for (; Ai < A+N*inca; Ai+=inca, Bi+=incb, Ci+=incc) {
			F.mul(*Ci,alpha,*Bi);
			F.addin (*Ci, *Ai);
		}
	}


	/***************************/
	/*         LEVEL 2         */
	/***************************/



	template <class Field>
	void
	fadd (const Field& F, const size_t M, const size_t N,
	      typename Field::ConstElement_ptr A, const size_t lda,
	      typename Field::ConstElement_ptr B, const size_t ldb,
	      typename Field::Element_ptr C, const size_t ldc)
	{
		if (N == lda && N == ldb && N == ldc)
			return fadd(F,M*N,A,1,B,1,C,1);
		typename Field::ConstElement_ptr Ai = A, Bi = B;
		typename Field::Element_ptr Ci = C;
		for (; Ai < A+M*lda; Ai+=lda, Bi+=ldb, Ci+=ldc)
			fadd(F,N,Ai,1,Bi,1,Ci,1);
	}

	template <class Field>
	void
	fsub (const Field& F, const size_t M, const size_t N,
	      typename Field::ConstElement_ptr A, const size_t lda,
	      typename Field::ConstElement_ptr B, const size_t ldb,
	      typename Field::Element_ptr C, const size_t ldc)
	{
		if (N == lda && N == ldb && N == ldc)
			return fsub(F,M*N,A,1,B,1,C,1);
		typename Field::ConstElement_ptr Ai = A, Bi = B;
		typename Field::Element_ptr Ci = C;
		for (; Ai < A+M*lda; Ai+=lda, Bi+=ldb, Ci+=ldc)
			fsub(F,N,Ai,1,Bi,1,Ci,1);
	}

	template <class Field>
	void
	faddin (const Field& F, const size_t M, const size_t N,
		typename Field::ConstElement_ptr B, const size_t ldb,
		typename Field::Element_ptr C, const size_t ldc)
	{
		if (N == ldb && N == ldc)
			return faddin(F,M*N,B,1,C,1);
		const typename Field::Element  *Bi = B;
		typename Field::Element_ptr Ci = C;
		for (; Bi < B+M*ldb;  Bi+=ldb, Ci+=ldc)
			faddin(F,N,Bi,1,Ci,1);
	}

	template <class Field>
	void
	fsubin (const Field& F, const size_t M, const size_t N,
		typename Field::ConstElement_ptr B, const size_t ldb,
		typename Field::Element_ptr C, const size_t ldc)
	{
		if (N == ldb && N == ldc)
			return fsubin(F,M*N,B,1,C,1);
		typename Field::ConstElement_ptr Bi = B;
		typename Field::Element_ptr Ci = C;
		for (; Bi < B+M*ldb;  Bi+=ldb, Ci+=ldc)
			fsubin(F,N,Bi,1,Ci,1);
	}


	// C = A + a B
	template <class Field>
	void
	fadd (const Field& F, const size_t M, const size_t N,
	      typename Field::ConstElement_ptr A, const size_t lda,
	      const typename Field::Element alpha,
	      typename Field::ConstElement_ptr B, const size_t ldb,
	      typename Field::Element_ptr C, const size_t ldc)
	{
		if (C == A && lda == ldc)
			return faxpy(F,M,N,alpha,B,ldb,C,ldc);
		if (F.isOne(alpha))
			return fadd(F,M,N,A,lda,B,ldb,C,ldc);
		if (F.isMOne(alpha))
			return fsub(F,M,N,A,lda,B,ldb,C,ldc);
		if (F.isZero(alpha))
			return fcopy(F,M,N,A,lda,C,ldc);

		if (N == lda && N == ldb && N == ldc)
			return fadd(F,M*N,A,1,alpha,B,1,C,1);

		typename Field::ConstElement_ptr Ai = A, Bi = B;
		typename Field::Element_ptr Ci = C;
		for (; Ai < A+M*lda; Ai+=lda, Bi+=ldb, Ci+=ldc)
			for (size_t i=0; i<N; i++) {
				F.mul(Ci[i],alpha,Bi[i]);
				F.addin (Ci[i], Ai[i]);
			}
	}


} // FFLAS

#endif // __FFLASFFPACK_fscal_INL
