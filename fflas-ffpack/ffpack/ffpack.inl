/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* ffpack.inl
 * Copyright (C) 2014 FFLAS-FFACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 * BB <bbboyer@ncsu.edu>
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

namespace FFPACK {


template <class Field>
size_t
Rank (const Field& F, const size_t M, const size_t N,
	  typename Field::Element_ptr A, const size_t lda)
{
	if (M == 0 and  N  == 0)
		return 0 ;

	size_t *P = FFLAS::fflas_new<size_t>(N);
	size_t *Q = FFLAS::fflas_new<size_t>(M);
	size_t R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, N, A, lda, P, Q);
	FFLAS::fflas_delete( Q);
	FFLAS::fflas_delete( P);
	return R;
}

template <class Field>
bool
IsSingular (const Field& F, const size_t M, const size_t N,
			typename Field::Element_ptr A, const size_t lda)
{
	if ( (M==0) and (N==0) ) return  false;
	if ( (M==0) or (N==0) )	return  true;
	if ( M != N ) return  true ;


		size_t *P = FFLAS::fflas_new<size_t>(N);
		size_t *Q = FFLAS::fflas_new<size_t>(M);
		bool singular = !LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, N, A, lda, P, Q, FfpackSingular);

		FFLAS::fflas_delete( P);
		FFLAS::fflas_delete( Q);
		return singular;
	}

template <class Field>
typename Field::Element
Det( const Field& F, const size_t M, const size_t N,
	 typename Field::Element_ptr A, const size_t lda)
{
	if ( (M==0) and (N==0) )
		return  F.one ;
	if ( (M==0) or (N==0) )
		return  F.zero ;
	if ( M != N )
		return  F.zero ;

	typename Field::Element det; F.init(det);
	bool singular;
	size_t *P = FFLAS::fflas_new<size_t>(N);
	size_t *Q = FFLAS::fflas_new<size_t>(M);
	singular  = !LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans,  M, N,
						   A, lda, P, Q, FfpackSingular);
	if (singular){
		F.assign(det,F.zero);
		FFLAS::fflas_delete( P);
		FFLAS::fflas_delete( Q);
		return det;
	}
	else{
		F.assign(det,F.one);
		typename Field::Element_ptr Ai=A;
		for (; Ai < A+ M*lda+N; Ai+=lda+1 )
			F.mulin( det, *Ai );
		int count=0;
		for (size_t i=0;i<N;++i)
			if (P[i] != i) ++count;

		if ((count&1) == 1)
			F.negin(det);
	}
	FFLAS::fflas_delete( P);
	FFLAS::fflas_delete( Q);
	return det;
}

template <class Field>
typename Field::Element_ptr
Solve( const Field& F, const size_t M,
	   typename Field::Element_ptr A, const size_t lda,
	   typename Field::Element_ptr x, const int incx,
	   typename Field::ConstElement_ptr b, const int incb )
{

	size_t *P = FFLAS::fflas_new<size_t>(M);
	size_t *rowP = FFLAS::fflas_new<size_t>(M);

	if (LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, M, A, lda, P, rowP) < M){
		std::cerr<<"SINGULAR MATRIX"<<std::endl;
		FFLAS::fflas_delete( P);
		FFLAS::fflas_delete( rowP);
		return x;
	}
	else{
		FFLAS::fassign( F, M, b, incb, x, incx );

		ftrsv(F,  FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, M,
			  A, lda , x, incx);
		ftrsv(F,  FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, M,
			  A, lda , x, incx);
		applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
				M, 0,(int) M, x, incx, P );
		FFLAS::fflas_delete( rowP);
		FFLAS::fflas_delete( P);

		return x;

	}
}

template <class Field>
size_t NullSpaceBasis (const Field& F, const FFLAS::FFLAS_SIDE Side,
					   const size_t M, const size_t N,
					   typename Field::Element_ptr A, const size_t lda,
					   typename Field::Element_ptr& NS, size_t& ldn,
					   size_t& NSdim)
{
	if (Side == FFLAS::FflasRight) { // Right NullSpace
		size_t* P = FFLAS::fflas_new<size_t>(N);
		size_t* Qt = FFLAS::fflas_new<size_t>(M);

		size_t R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, N, A, lda, P, Qt);
		delete [] Qt;

		ldn = N-R;
		NSdim = ldn;

		if (NSdim == 0) {
			FFLAS::fflas_delete( P);
			NS = NULL ;
			return NSdim ;
		}

		NS = FFLAS::fflas_new (F, N, ldn);

		if (R == 0) {
			FFLAS::fflas_delete( P);
			FFLAS::fidentity(F,N,ldn,NS,ldn);
			return NSdim;
		}

		FFLAS::fassign (F, R, ldn,  A + R,  lda, NS , ldn );

		ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, R, ldn,
			   F.mOne, A, lda, NS, ldn);

		FFLAS::fidentity(F,NSdim,NSdim,NS+R*ldn,ldn);

		applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
				NSdim, 0,(int) R, NS, ldn, P);

		delete [] P;

		return NSdim;
	}
	else { // Left NullSpace
		size_t* P = FFLAS::fflas_new<size_t>(M);
		size_t* Qt = FFLAS::fflas_new<size_t>(N);

		size_t R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasTrans, M, N, A, lda, P, Qt);
		delete [] Qt;

		ldn = M;
		NSdim = M-R;

		if (NSdim == 0) {
			FFLAS::fflas_delete( P);
			NS = NULL;
			return NSdim;
		}

		NS = FFLAS::fflas_new (F, NSdim, ldn);

		if (R == 0) {
			FFLAS::fflas_delete( P);
			FFLAS::fidentity(F,NSdim,ldn,NS,ldn);
			return NSdim;
		}

		FFLAS::fassign (F, NSdim, R, A + R *lda, lda, NS, ldn);
		ftrsm (F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, NSdim, R, F.mOne, A, lda, NS, ldn);

		FFLAS::fidentity(F,NSdim,NSdim,NS+R,ldn);
		applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, NSdim, 0,(int) R, NS, ldn, P);

		delete [] P;
		return NSdim;
	}
}

template <class Field>
size_t RowRankProfile (const Field& F, const size_t M, const size_t N,
					   typename Field::Element_ptr A, const size_t lda,
					   size_t* &rkprofile)
{
	size_t *P = FFLAS::fflas_new<size_t>(N);
	size_t *Q = FFLAS::fflas_new<size_t>(M);
	size_t R;

	R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, N, A, lda, P, Q);
	rkprofile = FFLAS::fflas_new<size_t>(R);
	
	for (size_t i=0; i<R; ++i)
		rkprofile[i] = Q[i];
	FFLAS::fflas_delete( P);
	FFLAS::fflas_delete( Q);
	return R;
}

template <class Field>
size_t ColumnRankProfile (const Field& F, const size_t M, const size_t N,
						  typename Field::Element_ptr A, const size_t lda,
						  size_t* &rkprofile)
{
	size_t *P = FFLAS::fflas_new<size_t>(M);
	size_t *Q = FFLAS::fflas_new<size_t>(N);
	size_t R;
	
	R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasTrans, M, N, A, lda, P, Q);
	rkprofile = FFLAS::fflas_new<size_t>(R);
	
	for (size_t i=0; i<R; ++i)
		rkprofile[i] = Q[i];
	FFLAS::fflas_delete( P);
	FFLAS::fflas_delete( Q);
	return R;
}

template <class Field>
size_t RowRankProfileSubmatrixIndices (const Field& F,
									   const size_t M, const size_t N,
									   typename Field::Element_ptr A,
									   const size_t lda,
									   size_t*& rowindices,
									   size_t*& colindices,
									   size_t& R)
{
	size_t *P = FFLAS::fflas_new<size_t>(N);
	size_t *Q = FFLAS::fflas_new<size_t>(M);
	
	R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, M, N, A, lda, P, Q);
	rowindices = FFLAS::fflas_new<size_t>(M);
	colindices = FFLAS::fflas_new<size_t>(N);
	for (size_t i=0; i<R; ++i){
		rowindices [i] = Q [i];
	}
	for (size_t i=0; i<N; ++i)
		colindices [i] = i;
	size_t tmp;
	for (size_t i=0; i<R; ++i){
		if (i != P[i]){
			tmp = colindices[i];
			colindices[i] = colindices[P[i]];
			colindices[P[i]] = tmp;
		}
	}
	
	FFLAS::fflas_delete( P);
	FFLAS::fflas_delete( Q);
	
	return R;
}
	
template <class Field>
size_t ColRankProfileSubmatrixIndices (const Field& F,
									   const size_t M, const size_t N,
									   typename Field::Element_ptr A,
									   const size_t lda,
									   size_t*& rowindices,
									   size_t*& colindices,
									   size_t& R)
{
	size_t *P = FFLAS::fflas_new<size_t>(M);
	size_t *Q = FFLAS::fflas_new<size_t>(N);
	
	R = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasTrans, M, N, A, lda, P, Q);
	rowindices = FFLAS::fflas_new<size_t>(M);
	colindices = FFLAS::fflas_new<size_t>(N);
	for (size_t i=0; i<R; ++i)
		colindices [i] = Q [i];
	
	for (size_t i=0; i<N; ++i)
		rowindices [i] = i;
	
	size_t tmp;
	for (size_t i=0; i<R; ++i){
		if (i != P[i]){
			tmp = rowindices[i];
			rowindices[i] = rowindices[P[i]];
			rowindices[P[i]] = tmp;
		}
	}
	FFLAS::fflas_delete( P);
	FFLAS::fflas_delete( Q);
	
	return R;
}
	
template <class Field>
size_t RowRankProfileSubmatrix (const Field& F,
								const size_t M, const size_t N,
								typename Field::Element_ptr A,
								const size_t lda,
								typename Field::Element_ptr& X, size_t& R)
{
	
	size_t * rowindices, * colindices;
	
	typename Field::Element_ptr A2 = FFLAS::fflas_new (F, M, N) ;
	FFLAS::fassign(F,M,N,A,lda,A2,N);

	RowRankProfileSubmatrixIndices (F, M, N, A2, N, rowindices, colindices, R);

	X = FFLAS::fflas_new (F, R, R);
	for (size_t i=0; i<R; ++i)
		for (size_t j=0; j<R; ++j)
			F.assign (*(X + i*R + j), *(A + rowindices[i]*lda + colindices[j]));
	FFLAS::fflas_delete (A2);
	FFLAS::fflas_delete( rowindices);
	FFLAS::fflas_delete( colindices);
	return R;
}

template <class Field>
size_t ColRankProfileSubmatrix (const Field& F, const size_t M, const size_t N,
								typename Field::Element_ptr A, const size_t lda,
								typename Field::Element_ptr& X, size_t& R)
{
	
	size_t * rowindices, * colindices;

	typename Field::Element_ptr A2 = FFLAS::fflas_new (F, M, N);
	FFLAS::fassign(F,M,N,A,lda,A2,N);

	ColRankProfileSubmatrixIndices (F, M, N, A2, N, rowindices, colindices, R);

	X = FFLAS::fflas_new (F, R, R);
	for (size_t i=0; i<R; ++i)
		for (size_t j=0; j<R; ++j)
			F.assign (*(X + i*R + j), *(A + rowindices[i]*lda + colindices[j]));
	FFLAS::fflas_delete (A2);
	FFLAS::fflas_delete( colindices);
	FFLAS::fflas_delete( rowindices);
	return R;
}
	
template <class Field>
typename Field::Element_ptr
LQUPtoInverseOfFullRankMinor( const Field& F, const size_t rank,
							  typename Field::Element_ptr A_factors, const size_t lda,
							  const size_t* QtPointer,
							  typename Field::Element_ptr X, const size_t ldx)
{

		// upper entries are okay, just need to move up bottom ones
	const size_t* srcRow = QtPointer;
	for (size_t row=0; row<rank; row++, srcRow++)
		if (*srcRow != row) {
			typename Field::Element_ptr oldRow = A_factors + (*srcRow) * lda;
			typename Field::Element_ptr newRow = A_factors + row * lda;
			for (size_t col=0; col<row; col++, oldRow++, newRow++)
				F.assign(*newRow, *oldRow);
		}
	
		// X <- (Qt.L.Q)^(-1)
		//invL( F, rank, A_factors, lda, X, ldx);
	ftrtri (F, FFLAS::FflasLower, FFLAS::FflasUnit, rank, A_factors, lda);
	FFLAS::fassign(F,rank,rank,X,ldx,A_factors,lda);
	
		// X = U^-1.X
	ftrsm( F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans,
		   FFLAS::FflasNonUnit, rank, rank, F.one, A_factors, lda, X, ldx);
	
	return X;
	
}

template<class Field>
void
solveLB( const Field& F, const FFLAS::FFLAS_SIDE Side,
		 const size_t M, const size_t N, const size_t R,
		 typename Field::Element_ptr L, const size_t ldl,
		 const size_t * Q,
		 typename Field::Element_ptr B, const size_t ldb )
{

	size_t LM = (Side == FFLAS::FflasRight)?N:M;
	int i = (int)R ;
	for (; i--; ){ // much faster for
		if (  Q[i] > (size_t) i){
				//for (size_t j=0; j<=Q[i]; ++j)
				//F.init( *(L+Q[i]+j*ldl), 0 );
				//std::cerr<<"1 deplacement "<<i<<"<-->"<<Q[i]<<endl;
			FFLAS::fassign( F, LM-Q[i]-1, L+(Q[i]+1)*ldl+i, ldl , L+Q[i]*(ldl+1)+ldl,ldl);
			for ( size_t j=Q[i]*ldl; j<LM*ldl; j+=ldl)
				F.assign( *(L+i+j), F.zero );
		}
	}
	ftrsm( F, Side, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, M, N, F.one, L, ldl , B, ldb);
		//write_field(F,std::cerr<<"dans solveLB "<<endl,L,N,N,ldl);
		// Undo the permutation of L
	for (size_t ii=0; ii<R; ++ii){
		if ( Q[ii] > (size_t) ii){
				//for (size_t j=0; j<=Q[ii]; ++j)
				//F.init( *(L+Q[ii]+j*ldl), 0 );
			FFLAS::fassign( F, LM-Q[ii]-1, L+Q[ii]*(ldl+1)+ldl,ldl, L+(Q[ii]+1)*ldl+ii, ldl );
			for ( size_t j=Q[ii]*ldl; j<LM*ldl; j+=ldl)
				F.assign( *(L+Q[ii]+j), F.zero );
		}
	}
}

template<class Field>
void
solveLB2( const Field& F, const FFLAS::FFLAS_SIDE Side,
		  const size_t M, const size_t N, const size_t R,
		  typename Field::Element_ptr L, const size_t ldl,
		  const size_t * Q,
		  typename Field::Element_ptr B, const size_t ldb )
{
	typename Field::Element_ptr Lcurr, Rcurr, Bcurr;
	size_t ib,  Ldim;
	int k;
	if ( Side == FFLAS::FflasLeft ){
		size_t j = 0;
		while ( j<R ) {
			ib = Q[j];
			k = (int)ib ;
			while ((j<R) && ( (int) Q[j] == k)  ) {k++;j++;}
			Ldim = (size_t)k-ib;
			Lcurr = L + j-Ldim + ib*ldl;
			Bcurr = B + ib*ldb;
			Rcurr = Lcurr + Ldim*ldl;

			ftrsm( F, Side, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, Ldim, N, F.one,
				   Lcurr, ldl , Bcurr, ldb );

			fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-(size_t)k, N, Ldim, F.mOne,
				   Rcurr , ldl, Bcurr, ldb, F.one, Bcurr+Ldim*ldb, ldb);
		}
	}
	else{ // Side == FFLAS::FflasRight
		int j=(int)R-1;
		while ( j >= 0 ) {
			ib = Q[j];
			k = (int) ib;
			while ( (j >= 0) &&  ( (int)Q[j] == k)  ) {--k;--j;}
			Ldim = ib-(size_t)k;
			Lcurr = L + j+1 + (k+1)*(int)ldl;
			Bcurr = B + ib+1;
			Rcurr = Lcurr + Ldim*ldl;

			fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M,  Ldim, N-ib-1, F.mOne,
				   Bcurr, ldb, Rcurr, ldl,  F.one, Bcurr-Ldim, ldb);

			ftrsm (F, Side, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, M, Ldim, F.one,
				   Lcurr, ldl , Bcurr-Ldim, ldb );
		}
	}
}

} // FFPACK

#endif // __FFLASFFPACK_ffpack_INL