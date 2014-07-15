/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* ffpack_echelon.h
 * Copyright (C) 2009, 2010 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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

#ifndef __FFLASFFPACK_ffpack_echelon_forms_INL
#define __FFLASFFPACK_ffpack_echelon_forms_INL

template <class Field>
size_t FFPACK::ColumnEchelonForm (const Field& F, const size_t M, const size_t N,
				  typename Field::Element_ptr A, const size_t lda,
				  size_t* P, size_t* Qt, const bool transform)
{

	size_t r;
	r = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, N, A, lda, P, Qt);

	if (transform){
		ftrtri (F, FFLAS::FflasUpper, FFLAS::FflasNonUnit, r, A, lda);
		ftrmm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, r, N-r,
		       F.mOne, A, lda, A+r, lda);
	}

	return r;
}

template <class Field>
size_t FFPACK::RowEchelonForm (const Field& F, const size_t M, const size_t N,
			       typename Field::Element_ptr A, const size_t lda,
			       size_t* P, size_t* Qt, const bool transform)
{

	size_t r;
	r = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasTrans,  M, N, A, lda, P, Qt);

	if (transform){

		ftrtri (F, FFLAS::FflasLower, FFLAS::FflasNonUnit, r, A, lda);
		ftrmm (F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, M-r, r,
		       F.mOne, A, lda, A+r*lda, lda);
	}

	return r;
}

template <class Field>
size_t
FFPACK::ReducedColumnEchelonForm (const Field& F, const size_t M, const size_t N,
				  typename Field::Element_ptr A, const size_t lda,
				  size_t* P, size_t* Qt, const bool transform)
{

	size_t r;
	r = ColumnEchelonForm (F, M, N, A, lda, P, Qt, transform);
	// M = Q^T M
	for (size_t i=0; i<r; ++i){
		if ( Qt[i]> (size_t) i ){
			FFLAS::fswap( F, i,
			       A + Qt[i]*lda, 1,
			       A + i*lda, 1 );
		}
	}
	if (transform){
		ftrtri (F, FFLAS::FflasLower, FFLAS::FflasUnit, r, A, lda);
		ftrmm (F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, M-r, r,
		       F.one, A, lda, A+r*lda, lda);
		ftrtrm (F, FFLAS::FflasNonUnit, r, A, lda);
	} else {
		ftrsm (F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, M-r, r,
		       F.one, A, lda, A+r*lda, lda);
		for (size_t i=0; i<r; i++){
			for (size_t j=0; j<N; j++)
				F.assign (*(A+i*lda+j), F.zero);
			F.assign (*(A + i*(lda+1)), F.one);
		}
		applyP(F, FFLAS::FflasLeft, FFLAS::FflasTrans,
		       r, 0,(int) r, A, lda, Qt);
	}

	return r;
}

template <class Field>
size_t
FFPACK::ReducedRowEchelonForm (const Field& F, const size_t M, const size_t N,
			       typename Field::Element_ptr A, const size_t lda,
			       size_t* P, size_t* Qt, const bool transform)
{

	size_t r;
	r = RowEchelonForm (F, M, N, A, lda, P, Qt, transform);
	// M = M Q
	for (size_t i=0; i<r; ++i){
		if ( Qt[i]> i ){
			FFLAS::fswap( F, i,
			       A + Qt[i], lda,
			       A + i, lda );
		}
	}
	if (transform){
		ftrtri (F, FFLAS::FflasUpper, FFLAS::FflasUnit, r, A, lda);
		ftrmm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasUnit, r, N-r,
		       F.one, A, lda, A+r, lda);
		ftrtrm (F, FFLAS::FflasUnit, r, A, lda);
	} else {
		ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasUnit, r, N-r,
		       F.one, A, lda, A+r, lda);
		for (size_t i=0; i<r; i++){
			for (size_t j=0; j<M; j++)
				F.assign (*(A+j*lda+i), F.zero);
			F.assign (*(A + i*(lda+1)), F.one);
		}
		applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
		       r, 0, (int)r, A, lda, Qt);
	}
	return r;
}

/*
 * Warning, this implementation is currently broken:
 * the LAPACK permutation mechanism can not be used here as is
 * More work required on the construction of the permutation P...
 */
template <class Field>
size_t
FFPACK::REF (const Field& F, const size_t M, const size_t N,
	     typename Field::Element_ptr A, const size_t lda,
	     const size_t colbeg, const size_t rowbeg, const size_t colsize,
	     size_t* Qt, size_t* P)
{
	if (colsize == 1){
		for (size_t i=rowbeg; i<M; ++i){
			if (!F.isZero(*(A+i*lda+colbeg))){
				Qt[rowbeg] = i;
				if (i!= rowbeg){
					F.assign(*(A+rowbeg*lda+colbeg),*(A+i*lda+colbeg));
					F.assign(*(A+i*lda+colbeg), F.zero);
				}
				typename Field::Element invpiv;
				F.inv(invpiv, *(A+rowbeg*lda + colbeg));
				F.assign(*(A+rowbeg*lda+colbeg), invpiv);
				F.negin(invpiv);
				// for (size_t j=0; j<rowbeg; ++j)
					// F.mulin (*(A+j*lda+colbeg), invpiv);
				FFLAS::fscalin(F,rowbeg,invpiv,A+colbeg,lda);
				// for (size_t j=rowbeg+1; j<M; ++j)
					// F.mulin (*(A+j*lda+colbeg), invpiv);
					FFLAS::fscalin(F,M-rowbeg-1,invpiv,A+colbeg,lda);
				return 1;
			}
		}
		Qt[rowbeg]=colbeg;
		return 0;
	}
	size_t recsize = colsize / 2;

	// Recurive call on slice A*1
	size_t r1 = REF(F, M, N, A, lda, colbeg, rowbeg, recsize, Qt, P);

	typename Field::Element_ptr A11 = A+colbeg;
	typename Field::Element_ptr A12 = A11+recsize;
	typename Field::Element_ptr A22 = A12+rowbeg*lda;
	typename Field::Element_ptr A21 = A11+rowbeg*lda;
	typename Field::Element_ptr A31 = A21+r1*lda;
	typename Field::Element_ptr A32 = A22+r1*lda;

	/**
	 *  ---------------------
	 * | I  | A11 | A12 |    |
	 * |----|-----|-----|----|
	 * |    |I | *| A22 |    |
	 * |    |0 | 0| A22 |    |
	 * |----|-----|-----|----|
	 * |    | 0   | A32 |    |
	 * |----|-----|-----|----|
	 *
	 * where the transformation matrix is stored at the pivot column position
	 */
	// Apply row permutation on A*2
	applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, colsize - recsize, rowbeg, rowbeg+r1, A12, lda, Qt);

	// A12 <- A12 - A11 * A22
	fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, rowbeg, colsize - recsize, r1,
	       F.one, A11, lda, A22, lda, F.one, A12, lda);

	// A32 <- A32 - A31 * A22
	fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-rowbeg-r1, colsize - recsize, r1,
	       F.one, A31, lda, A22, lda, F.one, A32, lda);

	// A22 <- A21*A22
	typename Field::Element_ptr tmp = fflas_new (F, r1, colsize-recsize);
	for (size_t i = 0; i < r1; ++i)
		fcopy (F, colsize-recsize, A22+i*lda, 1, tmp+i*(colsize-recsize), 1);
	fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, r1, colsize-recsize, r1,
	       F.one, A21, lda, tmp, colsize-recsize, F.zero, A22, lda);
	fflas_delete (tmp);

	// Recurive call on slice A*2
	size_t r2 = REF(F, M, N, A, lda, colbeg + recsize, rowbeg + r1,
			colsize - recsize, Qt, P);

	// Apply permutation on A*1
	applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, r1, rowbeg+r1, rowbeg+r1+r2, A11, lda, Qt);

	typename Field::Element_ptr U11 = A11;
	typename Field::Element_ptr U12 = A12;
	typename Field::Element_ptr U21 = A31;
	typename Field::Element_ptr U22 = A32;
	typename Field::Element_ptr U31 = U21+r2*lda;
	typename Field::Element_ptr U32 = U31+recsize;

	// U11 <- U11 + U12 * U21
	fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, rowbeg+r1, r1, r2,
	       F.one, U12, lda, U21, lda, F.one, U11, lda);

	// U31 <- U31 + U32 * U21
	fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-rowbeg-r1-r2, r1, r2,
	       F.one, U32, lda, U21, lda, F.one, U31, lda);

	// U21 <- U22*U21
	tmp = fflas_new (F, r2, r1);
	for (size_t i = 0; i < r2; ++i)
		fcopy (F, r1, U21+i*lda, 1, tmp+i*r1, 1);

	fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, r2, r1, r2,
	       F.one, U22, lda, tmp, r1, F.zero, U21, lda);
	fflas_delete(tmp);

	//Permute the non pivot columns to the end
	if (r1 < recsize){
		size_t ncol = recsize -r1;
		size_t nrow = rowbeg + r1;
		typename Field::Element_ptr NZ1 = A11+r1;

		tmp = fflas_new (F, nrow, ncol);
		for (size_t i=0; i < nrow; ++i)
			fcopy (F, ncol, NZ1 + i*lda, 1, tmp+i*ncol, 1);
		for (size_t i=0; i < M; ++i)
			// Risky copy with overlap, but safe with the naive
			// implementation of fcopy
			//! @bug safe ???
			fcopy (F, r2, A12 + i*lda, 1, NZ1+i*lda, 1);
		NZ1 +=  r2;
		for (size_t i=0; i<nrow; ++i)
			fcopy (F, ncol, tmp+i*ncol,1, NZ1 + i*lda, 1);
		fflas_delete (tmp);

		for (size_t i=rowbeg+r1; i<M; ++i)
			for (size_t j=0; j<recsize-r1; ++j)
				F.assign(*(NZ1+i*lda+j), F.zero);
		// size_t * temp = new size_t[recsize-r1];
		// for (size_t i=0,j = colbeg+r1; j<colbeg+recsize; ++i,++j)
		//  	temp[i] = P[j];
		// for (size_t  i = colbeg+recsize, j = colbeg+r1; i<colbeg+recsize+r2; ++i,++j)
		// 	P[j] = P[i];
		// for (size_t i=0,j = colbeg+r1+r2; i<recsize-r1; ++i,++j)
		// 	P[j] = temp[i];
		// delete temp;
		for (size_t  i = colbeg+recsize, j = colbeg+r1; i<colbeg+recsize+r2; ++i,++j){
			size_t t = P[i];
			P[i] = P[j];
			P[j] = t;
			//P[j]=P[i];
		}

	}

	return r1+r2;
}

namespace FFPACK {
	template <class Field>
	size_t
	ReducedRowEchelonForm2 (const Field& F, const size_t M, const size_t N,
				typename Field::Element_ptr A, const size_t lda,
				size_t* P, size_t* Qt, const bool transform /*= true */)
	{
		for (size_t i=0; i<N; ++i)
			Qt[i] = i;
		return REF (F, M, N, A, lda, 0, 0, N, P, Qt);

	}
} // FFPACK

#endif  // __FFLASFFPACK_ffpack_echelon_forms_INL
