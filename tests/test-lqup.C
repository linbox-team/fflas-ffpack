/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Clément Pernet
 * This file is Free Software and part of FFLAS-FFPACK.
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


//--------------------------------------------------------------------------
//          Test for the lsp factorisation
//--------------------------------------------------------------------------
// usage: test-lsp p A n, for n lsp factorization
// of A over Z/pZ
//-------------------------------------------------------------------------

using namespace std;


//#define __LUDIVINE_CUTOFF 1
#include <iostream>
#include <iomanip>
#include "Matio.h"
#include "utils/timer.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "test-utils.h"

#include "fflas-ffpack/utils/args-parser.h"

using namespace FFPACK;


/*! Tests the LUdivine routine.
 * @tparam Field Field
 * @tparam Diag  Unit diagonal in L ?
 * @tparam Trans ?
 * @param F field
 * @param A Matrix (preallocated)
 * @param r rank of A
 * @param m rows
 * @param n cols
 * @param lda leading dim of A
 * @return 0 iff correct, 1 otherwise
 */
template<class Field, FFLAS::FFLAS_DIAG diag, FFLAS::FFLAS_TRANSPOSE trans>
bool test_lu(const Field & F,
	     const typename Field::Element * A,
	     size_t r,
	     size_t m, size_t n, size_t lda)
{
	bool fail = false;
	typedef typename Field::Element Element ;
	Element * B = new Element[m*lda] ;
	// memcpy(B,A,m*lda*sizeof(Element)); // probably faster than ::fcopy !
	FFLAS::fcopy(F,m,n,B,lda,A,lda);

	size_t maxP, maxQ ;

	if (trans == FFLAS::FflasTrans){
		maxP = m;
		maxQ = n;
	}
	else{ // trans == FFLAS::FflasNoTrans
		maxP = n;
		maxQ = m;
	}

	size_t * P = new size_t[maxP] ;
	size_t * Q = new size_t[maxQ] ;

	size_t R = FFPACK::LUdivine (F, diag, trans, m, n, B, lda, P, Q,
				     FFPACK::FfpackLQUP);

	if (R != r) {
		std::cout << "rank is wrong (expected " << R << " but got " << r << ")" << std::endl;
		delete[] B ;
		delete[] P ;
		delete[] Q ;
		return fail = true;
	}

	Element * C = new Element[m*n]; // compute C=LQUP and check C == A
	/*  Build L,U */
	Element * L, *U;
	if (trans == FFLAS::FflasNoTrans){
		L = new Element[m*m];
		U = new Element[m*n];

		Element zero,one;
		F.init(zero,0.0);
		F.init(one,1.0);
		/*  build U */
		for (size_t i=0; i<R; ++i){
			for (size_t j=0; j<i; ++j)
				F.assign ( *(U + i*n + j), zero);
			for (size_t j=i+1; j<n; ++j)
				F.assign (*(U + i*n + j), *(B+ i*lda+j));
		}
		for (size_t i=R;i<m; ++i) {
			for (size_t j=0; j<n; ++j)
				F.assign(*(U+i*n+j), zero);
		}
		/*  build L */
		for ( size_t i=0; i<m; ++i ){
			size_t j=0;
			for (; j< ((i<R)?i:R) ; ++j )
				F.assign( *(L + i*m+j), *(B+i*lda+j));
			for (; j<m; ++j )
				F.assign( *(L+i*m+j), zero);
		}

		// write_field(F,cerr<<"L = "<<endl,L,m,m,m);
		//write_field(F,cerr<<"U = "<<endl,U,m,n,n);
		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				m,0,(int)R, L, m, Q);
		for ( size_t i=0; i<m; ++i )
			F.assign(*(L+i*(m+1)), one);

		/*  reconstruct the diagonal */
		//write_field(F,cerr<<"L = "<<endl,L,m,m,m);
		//write_field(F,cerr<<"U = "<<endl,U,m,n,n);
		if (diag == FFLAS::FflasNonUnit) {
			for ( size_t i=0; i<R; ++i )
				F.assign (*(U+i*(n+1)), *(B+i*(lda+1)));
		}
		else{ // diag == FFLAS::FflasUnit
			for ( size_t i=0; i<R; ++i ){
				*(L+Q[i]*(m+1)) = *(B+Q[i]*lda+i);
				F.assign (*(U+i*(n+1)),one);
			}
		}
		//write_field(F,cerr<<"L = "<<endl,L,m,m,m);
		//write_field(F,cerr<<"U = "<<endl,U,m,n,n);

		/*  Compute LQUP */
		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				m,0,(int) R, U, n, P);
		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
				n,0,(int)R, U, n, Q);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
			      m,n,m, 1.0, L,m, U,n, 0.0, C,n);
		//delete[] A;
	}
	else { /*  trans == FFLAS::FflasTrans */

		L = new Element[m*n];
		U = new Element[n*n];


		Element zero,one;
		F.init(zero,0.0);
		F.init(one,1.0);
		/*  build L */
		for (size_t i=0; i<R; ++i){
			for (size_t j=0; j<i; ++j)
				F.assign ( *(L + i + j*n), zero);
			for (size_t j=i+1; j<m; ++j)
				F.assign (*(L + i + j*n), *(B+ i+j*lda));
		}
		for (size_t i=R;i<n; ++i) {
			for (size_t j=0; j<m; ++j)
				F.assign(*(L+i+j*n), zero);
		}
		/*  build U */
		for ( size_t i=0; i<n; ++i ){
			size_t j=0;
			for (;  j< ((i<R)?i:R) ; ++j )
				F.assign( *(U + i+j*n), *(B+i+j*lda));
			for (; j<n; ++j )
				F.assign( *(U+i+j*n), zero);
		}
		//write_field(F,cerr<<"L = "<<endl,L,m,n,n);
		//write_field(F,cerr<<"U = "<<endl,U,n,n,n);

		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans,
				n,0,(int)R, U, n, Q);

		for (size_t i=0; i<n; ++i)
			F.assign (*(U+i*(n+1)),one);

		/*  reconstruct the diagonal */
		if (diag == FFLAS::FflasNonUnit) {
			for ( size_t i=0; i<R; ++i )
				F.assign (*(L+i*(n+1)), *(B+i*(lda+1)));
		}
		else{ // diag == FFLAS::FflasUnit
			for ( size_t i=0; i<R; ++i ){
				*(U+Q[i]*(n+1)) = *(B+Q[i]+i*lda);
				F.assign (*(L+i*(n+1)),one);
			}
		}
		//write_field(F,cerr<<"L = "<<endl,L,m,n,n);
		//write_field(F,cerr<<"U = "<<endl,U,n,n,n);

		/*  Compute LQUP */
		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
				n,0,(int)R, L, n, P);
		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				m,0,(int)R, L, n, Q);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
			      m,n,n, 1.0, L,n, U,n, 0.0, C,n);
	}
	/*  check equality */
	for (size_t i=0; i<m; ++i) {
		for (size_t j=0; j<n; ++j)
			if (!F.areEqual (*(A+i*lda+j), *(C+i*n+j))){
				std::cerr << " A["<<i<<","<<j<<"]    = " << (*(A+i*lda+j))
				<< " PLUQ["<<i<<","<<j<<"] = " << (*(C+i*n+j))
				<< endl;
				fail|=true;
			}
	}

	delete[] P;
	delete[] L;
	delete[] U;
	delete[] Q;
	delete[] B;
	delete[] C;
	return fail;


}

/*! Tests the LUpdate routine.
 * @tparam Field Field
 * @tparam Diag  Unit diagonal in L ?
 * @tparam Trans ?
 * @param F field
 * @param A Matrix (preallocated)
 * @param r rank of A
 * @param B Matrix (preallocated)
 * @param m rows in A
 * @param n cols in A (and B)
 * @param k rows in B
 * @param lda leading dim of A (and B)
 * @return 0 iff correct, 1 otherwise
 */
template<class Field, FFLAS::FFLAS_DIAG diag, FFLAS::FFLAS_TRANSPOSE trans>
bool test_lu_append(const Field & F,
		    const typename Field::Element * A,
		    const typename Field::Element * B,
		    size_t m, size_t n, size_t k, size_t lda)
{
	FFLASFFPACK_check(n<=lda);

	bool fail = false;
	size_t M = m + k ;
	typedef typename Field::Element Element ;
	Element * Acop = new Element[m*lda] ;
	FFLAS::fcopy(F,m,n,Acop,lda,A,lda) ;

	Element * Bcop = new Element[k*lda] ;
	FFLAS::fcopy(F,k,n,Bcop,lda,B,lda) ;

	Element * Append = new Element[M*lda];
	FFLAS::fcopy(F,m,n,Append,lda,A,lda) ;
	FFLAS::fcopy(F,k,n,Append+m*lda,lda,B,lda) ;

#if 0 /*  paranoid check */
	for (size_t i = 0 ; i < m ; ++i) {
		for (size_t j = 0 ; j < n ; ++j) {
			FFLASFFPACK_check(Append[i*lda+j]==A[i*lda+j]);
		}
	}
	for (size_t i = 0 ; i < k ; ++i) {
		for (size_t j = 0 ; j < n ; ++j) {
			FFLASFFPACK_check(Append[(i+m)*lda+j]==B[i*lda+j]);
		}
	}
#endif

	Element * Afull = new Element[M*lda];
	FFLAS::fcopy(F,M,n,Afull,lda,Append,lda) ;
	// FFLAS::fcopy(F,m,n,Afull,lda,A,lda) ;
	// FFLAS::fcopy(F,k,n,Afull+m*lda,lda,B,lda) ;

#if 0
std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
	for (size_t i = 0 ; i < m ; ++i) {
		for (size_t j = 0 ; j < n ; ++j) {
			std::cout << Append[i*lda+j] << "(" << A[i*lda+j] << ") " ;
		} std::cout << std::endl;
	}
std::cout << "-----------------------------------" << std::endl;
	for (size_t i = 0 ; i < k ; ++i) {
		for (size_t j = 0 ; j < n ; ++j) {
			std::cout << Append[(i+m)*lda+j] ;
			std::cout << "(" << B[i*lda+j] << ") "  ;
		}std::cout << std::endl;
	}
std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::flush << std::endl;
#endif



#if 0
	for (size_t i = 0 ; i < m ; ++i)
		for (size_t j = 0 ; j < n ; ++j)
			FFLASFFPACK_check(Acop[i*lda+j]==A[i*lda+j]);
	for (size_t i = 0 ; i < k ; ++i)
		for (size_t j = 0 ; j < n ; ++j)
			FFLASFFPACK_check(Bcop[i*lda+j]==B[i*lda+j]);
	for (size_t i = 0 ; i < M ; ++i)
		for (size_t j = 0 ; j < n ; ++j)
			if (i < m)
				FFLASFFPACK_check(Afull[i*lda+j]==A[i*lda+j]);
			else
				FFLASFFPACK_check(Afull[i*lda+j]==B[(i-m)*lda+j]);
#endif




	size_t maxP, maxQ ;

	if (trans == FFLAS::FflasTrans){
		maxP = M;
		maxQ = n;
	}
	else{ // trans == FFLAS::FflasNoTrans
		maxP = n;
		maxQ = M;
	}

	size_t * P = new size_t[maxP] ;
	size_t * Q = new size_t[maxQ] ;

	size_t * PP = new size_t[maxP] ;
	size_t * QQ = new size_t[maxQ] ;

	/* valgrind says the following leaks. Just incroyable. */
	size_t R  = FFPACK::LUdivine (F, diag, trans, M, n, Append, lda, PP, QQ,
				      FFPACK::FfpackLQUP);

	size_t R1 = FFPACK::LUdivine (F, diag, trans, m, n, Acop,   lda, P, Q,
				      FFPACK::FfpackLQUP);

	size_t R2 = FFPACK::LUpdate  (F,diag,trans,m,n,Acop,lda,R1,k,Bcop,lda,P,Q,
				      FFPACK::FfpackLQUP);
#if 0
	std::cout << "P := [ " ;
	for (size_t i = 0 ; i < maxP ; ++i)
		std::cout << P[i] << " " ;
	std::cout << ']' << std::endl;
	std::cout << "Q := [ ";
	for (size_t i = 0 ; i < maxQ ; ++i)
		std::cout << Q[i] << " " ;
	std::cout << ']' << std::endl;
	std::cout << "PP := [ ";
	for (size_t i = 0 ; i < maxP ; ++i)
		std::cout << PP[i] << " " ;
	std::cout << ']' << std::endl;
	std::cout << "QQ := [ ";
	for (size_t i = 0 ; i < maxQ ; ++i)
		std::cout << QQ[i] << " " ;
	std::cout << ']' << std::endl;
#endif

	if (R2 != R) {
		std::cout << "error, bad rank " << R2 << " <> " << R << " (expected) " << std::endl;
		delete[] Bcop ;
		delete[] Acop ;
		delete[] Append ;
		delete[] PP;
		delete[] QQ;
		delete[] P ;
		delete[] Q ;
		return fail=true;

	}

	// compute C=LQUP and check C == A
	Element * C = new Element[M*lda];
	/*  Build L,U */
	Element * L, *U;
	if (trans == FFLAS::FflasNoTrans){
		L = new Element[M*M];
		U = new Element[M*n];

		Element zero,one;
		F.init(zero,0.0);
		F.init(one,1.0);
		/*  build U */
		for (size_t i=0; i<R; ++i){
			for (size_t j=0; j<i; ++j)
				F.assign ( *(U + i*n + j), zero);
			for (size_t j=i+1; j<n; ++j)
				if (i < m)
					F.assign (*(U + i*n + j), *(Acop+ i*lda+j));
				else
					F.assign (*(U + i*n + j), *(Bcop+ (i-m)*lda+j));
		}

		for (size_t i=R;i<M; ++i) {
			for (size_t j=0; j<n; ++j)
				F.assign(*(U+i*n+j), zero);
		}
		/*  build L */
		for ( size_t i=0; i<M; ++i ){
			size_t j=0;
			for (; j< ((i<R)?i:R) ; ++j ) {
				if (i<m)
					F.assign( *(L + i*M+j), *(Acop+i*lda+j));
				else
					F.assign( *(L + i*M+j), *(Bcop+(i-m)*lda+j));
			}
			for (; j<M; ++j )
				F.assign( *(L+i*M+j), zero);
		}

		// write_field(F,cerr<<"L = "<<endl,L,m,m,m);
		//write_field(F,cerr<<"U = "<<endl,U,m,n,n);
		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				M,0,(int)R, L, M, Q);
		for ( size_t i=0; i<M; ++i )
			F.assign(*(L+i*(M+1)), one);

		/*  reconstruct the diagonal */
		//write_field(F,cerr<<"L = "<<endl,L,m,m,m);
		//write_field(F,cerr<<"U = "<<endl,U,m,n,n);
		if (diag == FFLAS::FflasNonUnit) {
			for ( size_t i=0; i<R; ++i ) {
				if (i<m)
					F.assign (*(U+i*(n+1)), *(Acop+i*(lda+1)));
				else
					F.assign (*(U+i*(n+1)), *(Bcop+(i-m)*(lda+1)));
			}
		}
		else{ // diag == FFLAS::FflasUnit
			for ( size_t i=0; i<R; ++i ){
				if (Q[i] < m)
					*(L+Q[i]*(M+1)) = *(Acop+Q[i]*lda+i);
				else
					*(L+Q[i]*(M+1)) = *(Bcop+(Q[i]-m)*lda+i);
				F.assign (*(U+i*(n+1)),one);
			}
		}
		// write_field(F,cerr<<"L = "<<endl,L,(int)M,(int)M,(int)M);
		// write_field(F,cerr<<"U = "<<endl,U,(int)M,(int)n,(int)n);

		/*  Compute LQUP */
		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				M,0,(int) R, U, n, P);
		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
				n,0,(int)R, U, n, Q);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
			      M,n,M, 1.0, L,M, U,n, 0.0, C,lda);
		//delete[] A;
	}
#if 0 /*  not working */
	else { /*  trans == FFLAS::FflasTrans */

		L = new Element[M*n];
		U = new Element[n*n];


		Element zero,one;
		F.init(zero,0.0);
		F.init(one,1.0);
		/*  build L */
		for (size_t i=0; i<R; ++i){
			for (size_t j=0; j<i; ++j)
				F.assign ( *(L + i + j*n), zero);
			for (size_t j=i+1; j<M; ++j) {
				if (i < m)
					F.assign (*(L + i + j*n), *(Acop+ i+j*lda));
				else
					F.assign (*(L + i + j*n), *(Bcop+ (i-m)+j*lda));
			}
		}
		for (size_t i=R;i<n; ++i) {
			for (size_t j=0; j<M; ++j)
				F.assign(*(L+i+j*n), zero);
		}
		/*  build U */
		for ( size_t i=0; i<n; ++i ){
			size_t j=0;
			for (;  j< ((i<R)?i:R) ; ++j ) {
				if (i < m)
					F.assign( *(U + i+j*n), *(Acop+i+j*lda));
				else
					F.assign( *(U + i+j*n), *(Bcop+(i-m)+j*lda));
			}
			for (; j<n; ++j )
				F.assign( *(U+i+j*n), zero);
		}
		//write_field(F,cerr<<"L = "<<endl,L,m,n,n);
		//write_field(F,cerr<<"U = "<<endl,U,n,n,n);

		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans,
				n,0,(int)R, U, n, Q);

		for (size_t i=0; i<n; ++i)
			F.assign (*(U+i*(n+1)),one);

		/*  reconstruct the diagonal */
		if (diag == FFLAS::FflasNonUnit) {
			for ( size_t i=0; i<R; ++i ) {
				if (i < m)
					F.assign (*(L+i*(n+1)), *(Acop+i*(lda+1)));
				else
					F.assign (*(L+i*(n+1)), *(Bcop+(i-m)*(lda+1)));
			}
		}
		else{ // diag == FFLAS::FflasUnit
			for ( size_t i=0; i<R; ++i ){
				if (i<m)
					*(U+Q[i]*(n+1)) = *(Acop+Q[i]+i*lda);
				else
					*(U+Q[i]*(n+1)) = *(Bcop+Q[i]+(i-m)*lda);
				F.assign (*(L+i*(n+1)),one);
			}
		}
		// write_field(F,cerr<<"L = "<<endl,L,(int)M,(int)n,(int)n);
		// write_field(F,cerr<<"U = "<<endl,U,(int)n,(int)n,(int)n);

		/*  Compute LQUP */
		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
				n,0,(int)R, L, n, P);
		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				M,0,(int)R, L, n, Q);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
			      M,n,n, 1.0, L,n, U,n, 0.0, C,lda);
	}
#endif
#if 0 /*  check CC == LL UU */
	Element * LL, *UU;
	Element * CC = new Element[M*lda];
	if (trans == FFLAS::FflasNoTrans){
		LL = new Element[M*M];
		UU = new Element[M*n];

		Element zero,one;
		F.init(zero,0.0);
		F.init(one,1.0);
		/*  build U */
		for (size_t i=0; i<R; ++i){
			for (size_t j=0; j<i; ++j)
				F.assign ( *(UU + i*n + j), zero);
			for (size_t j=i+1; j<n; ++j)
				F.assign (*(UU + i*n + j), *(Append+ i*lda+j));
		}

		for (size_t i=R;i<M; ++i) {
			for (size_t j=0; j<n; ++j)
				F.assign(*(UU+i*n+j), zero);
		}
		/*  build L */
		for ( size_t i=0; i<M; ++i ){
			size_t j=0;
			for (; j< ((i<R)?i:R) ; ++j ) {
				F.assign( *(LL + i*M+j), *(Append+i*lda+j));
			}
			for (; j<M; ++j )
				F.assign( *(LL+i*M+j), zero);
		}

		// write_field(F,cerr<<"LL = "<<endl,LL,m,m,m);
		//write_field(F,cerr<<"UU = "<<endl,UU,m,n,n);
		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				M,0,(int)R, LL, M, Q);
		for ( size_t i=0; i<M; ++i )
			F.assign(*(LL+i*(M+1)), one);

		/*  reconstruct the diagonal */
		//write_field(F,cerr<<"LL = "<<endl,LL,m,m,m);
		//write_field(F,cerr<<"UU = "<<endl,UU,m,n,n);
		if (diag == FFLAS::FflasNonUnit) {
			for ( size_t i=0; i<R; ++i ) {
				F.assign (*(UU+i*(n+1)), *(Append+i*(lda+1)));
			}
		}
		else{ // diag == FFLAS::FflasUnit
			for ( size_t i=0; i<R; ++i ){
				*(LL+Q[i]*(M+1)) = *(Append+Q[i]*lda+i);
				F.assign (*(UU+i*(n+1)),one);
			}
		}
		write_field(F,cerr<<"L = "<<endl,LL,(int)M,(int)M,(int)M);
		write_field(F,cerr<<"U = "<<endl,UU,(int)M,(int)n,(int)n);

		/*  Compute LQUP */
		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				M,0,(int) R, UU, n, P);
		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
				n,0,(int)R, UU, n, Q);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
			      M,n,M, 1.0, LL,M, UU,n, 0.0, CC,lda);
		//delete[] A;
	}
	else { /*  trans == FFLAS::FflasTrans */

		LL = new Element[M*n];
		UU = new Element[n*n];


		Element zero,one;
		F.init(zero,0.0);
		F.init(one,1.0);
		/*  build L */
		for (size_t i=0; i<R; ++i){
			for (size_t j=0; j<i; ++j)
				F.assign ( *(LL + i + j*n), zero);
			for (size_t j=i+1; j<M; ++j) {
				F.assign (*(LL + i + j*n), *(Append+ i+j*lda));
			}
		}
		for (size_t i=R;i<n; ++i) {
			for (size_t j=0; j<M; ++j)
				F.assign(*(LL+i+j*n), zero);
		}
		/*  build UU */
		for ( size_t i=0; i<n; ++i ){
			size_t j=0;
			for (;  j< ((i<R)?i:R) ; ++j ) {
				F.assign( *(UU + i+j*n), *(Append+i+j*lda));
			}
			for (; j<n; ++j )
				F.assign( *(UU+i+j*n), zero);
		}
		// 		write_field(F,cerr<<"LL = "<<endl,LL,m,n,n);
		// 		write_field(F,cerr<<"UU = "<<endl,UU,n,n,n);

		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans,
				n,0,(int)R, UU, n, Q);

		for (size_t i=0; i<n; ++i)
			F.assign (*(UU+i*(n+1)),one);

		/*  reconstruct the diagonal */
		if (diag == FFLAS::FflasNonUnit) {
			for ( size_t i=0; i<R; ++i ) {
				F.assign (*(LL+i*(n+1)), *(Append+i*(lda+1)));
			}
		}
		else{ // diag == FFLAS::FflasUnit
			for ( size_t i=0; i<R; ++i ){
				*(UU+Q[i]*(n+1)) = *(Append+Q[i]+i*lda);
				F.assign (*(LL+i*(n+1)),one);
			}
		}
		write_field(F,cerr<<"LL = "<<endl,LL,(int)M,(int)n,(int)n);
		write_field(F,cerr<<"UU = "<<endl,UU,(int)n,(int)n,(int)n);

		/*  Compute LQUP */
		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
				n,0,(int)R, LL, n, P);
		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				M,0,(int)R, LL, n, Q);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
			      M,n,n, 1.0, LL,n, UU,n, 0.0, CC,lda);
	}
	for (size_t i=0; i<M; ++i) {
		for (size_t j=0; j<n; ++j)
			if (!F.areEqual (*(Afull+i*lda+j), *(CC+i*lda+j))){
				std::cerr << " A["<<i<<","<<j<<"]    = " << (*(Afull+i*lda+j))
				<< " PLUQ["<<i<<","<<j<<"] = " << (*(CC+i*lda+j))
				<< endl << "xxxx" << endl;
				fail|=true;
			}
	}

#endif

	/*  check equality */
	for (size_t i=0; i<M; ++i) {
		for (size_t j=0; j<n; ++j)
			if (!F.areEqual (*(Afull+i*lda+j), *(C+i*lda+j))){
				std::cerr << " A["<<i<<","<<j<<"]    = " << (*(Afull+i*lda+j))
				<< " PLUQ(append)["<<i<<","<<j<<"] = " << (*(C+i*lda+j))
				<< endl;
				fail|=true;
			}
	}

	delete[] PP;
	delete[] P;
	delete[] L;
	delete[] U;
	delete[] Q;
	delete[] QQ;
	delete[] Acop;
	delete[] Bcop;
	delete[] Append;
	delete[] Afull;
	delete[] C;

	return fail;


}



template<class Field, FFLAS::FFLAS_DIAG diag, FFLAS::FFLAS_TRANSPOSE trans>
bool launch_test(const Field & F,
		 size_t r,
		 size_t m, size_t n)
{
	typedef typename Field::Element Element ;
	bool fail = false ;
	{ /*  user given and lda bigger */
		size_t lda = n+10 ;
		Element * A = new Element[m*lda];
		RandomMatrixWithRank(F,A,r,m,n,lda);
		fail |= test_lu<Field,diag,trans>(F,A,r,m,n,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
	}
	{ /*  user given and lda bigger. Rank is max */
		size_t lda = n+10 ;
		size_t R = std::min(m,n);
		Element * A = new Element[m*lda];
		RandomMatrixWithRank(F,A,R,m,n,lda);
		fail |= test_lu<Field,diag,trans>(F,A,R,m,n,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
	}
	{ /*  user given and lda bigger. Rank is min */
		size_t lda = n+10 ;
		size_t R = 0;
		Element * A = new Element[m*lda];
		RandomMatrixWithRank(F,A,R,m,n,lda);
		fail |= test_lu<Field,diag,trans>(F,A,R,m,n,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
	}
	{ /*  square  */
		size_t M = std::max(m,n);
		size_t N = M ;
		size_t R = M/2 ;
		size_t lda = N+10 ;
		Element * A = new Element[M*lda];
		RandomMatrixWithRank(F,A,R,M,N,lda);
		fail |= test_lu<Field,diag,trans>(F,A,R,M,N,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
	}
	{ /*  wide  */
		size_t M = std::max(m,n);
		size_t N = 2*M ;
		size_t R = 3*M/4 ;
		size_t lda = N+5 ;
		Element * A = new Element[M*lda];
		RandomMatrixWithRank(F,A,R,M,N,lda);
		fail |= test_lu<Field,diag,trans>(F,A,R,M,N,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
	}
	{ /*  narrow  */
		size_t M = std::max(m,n);
		size_t N = M/2 ;
		size_t R = 3*M/8 ;
		size_t lda = N+5 ;
		Element * A = new Element[M*lda];
		RandomMatrixWithRank(F,A,R,M,N,lda);
		fail |= test_lu<Field,diag,trans>(F,A,R,M,N,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
	}

	return fail;
}

template<class Field, FFLAS::FFLAS_DIAG diag, FFLAS::FFLAS_TRANSPOSE trans>
bool launch_test_append(const Field & F,
			size_t r,
			size_t m, size_t n)
{
	typedef typename Field::Element Element ;
	bool fail = false ;
	{ /*  user given and lda bigger */
		size_t lda = n+10 ;
		size_t k = m/2+1 ;
		Element * A = new Element[m*lda];
		Element * B = new Element[k*lda];
		RandomMatrixWithRank(F,A,r,m,n,lda);
		RandomMatrixWithRank(F,B,k/2+1,k,n,lda);
		fail |= test_lu_append<Field,diag,trans>(F,A,B,m,n,k,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
		delete[] B ;
	}
	{ /*  user given and lda bigger. Rank is max */
		size_t lda = n+10 ;
		size_t R = std::min(m,n);
		size_t k = m/2+1 ;
		Element * A = new Element[m*lda];
		Element * B = new Element[k*lda];
		RandomMatrixWithRank(F,A,R,m,n,lda);
		RandomMatrixWithRank(F,B,k/2+1,k,n,lda);
		fail |= test_lu_append<Field,diag,trans>(F,A,B,m,n,k,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
		delete[] B ;
	}
	{ /*  user given and lda bigger. Appended Rank is min */
		size_t lda = n+10 ;
		size_t R = std::min(m,n);
		size_t k = m/2+1 ;
		Element * A = new Element[m*lda];
		Element * B = new Element[k*lda];
		RandomMatrixWithRank(F,A,R,m,n,lda);
		RandomMatrixWithRank(F,B,0,k,n,lda);
		fail |= test_lu_append<Field,diag,trans>(F,A,B,m,n,k,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
		delete[] B ;
	}
	{ /*  user given and lda bigger. Rank is min */
		size_t lda = n+10 ;
		size_t R = 0;
		size_t k = m/2+1 ;
		Element * A = new Element[m*lda];
		Element * B = new Element[k*lda];
		RandomMatrixWithRank(F,A,R,m,n,lda);
		RandomMatrixWithRank(F,B,k/2+1,k,n,lda);
		fail |= test_lu_append<Field,diag,trans>(F,A,B,m,n,k,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
		delete[] B ;
	}
	{ /*  square  */
		size_t M = std::max(m,n);
		size_t N = M ;
		size_t R = M/2 ;
		size_t lda = N+10 ;
		size_t k = R ;
		Element * A = new Element[M*lda];
		Element * B = new Element[k*lda];
		RandomMatrixWithRank(F,A,R,M,N,lda);
		RandomMatrixWithRank(F,B,R/2,k,N,lda);
		fail |= test_lu_append<Field,diag,trans>(F,A,B,M,N,k,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
		delete[] B ;
	}
	{ /*  wide  */
		size_t M = std::max(m,n);
		size_t N = 2*M ;
		size_t R = M/2 ;
		size_t k = R ;
		size_t lda = N+10 ;
		Element * A = new Element[M*lda];
		Element * B = new Element[k*lda];
		RandomMatrixWithRank(F,A,R,M,N,lda);
		RandomMatrixWithRank(F,B,k/2,k,N,lda);
		fail |= test_lu_append<Field,diag,trans>(F,A,B,M,N,k,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
		delete[] B ;
	}
#if 0 /*  leak here */
	{ /*  narrow  */
		size_t M = std::max(m,n);
		size_t N = M/2 ;
		size_t R = M/3 ;
		size_t k = N ;
		size_t lda = N+10 ;
		Element * A = new Element[M*lda];
		Element * B = new Element[k*lda];
		RandomMatrixWithRank(F,A,R,M,N,lda);
		RandomMatrixWithRank(F,A,std::min(k/2,M/2),k,N,lda);
		fail |= test_lu_append<Field,diag,trans>(F,A,B,M,N,k,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
		delete[] B ;
	}
#endif

	return fail;
}

int main(int argc, char** argv)
{
	cerr<<setprecision(20);
	int p = 101;
	size_t m = 50;
	size_t n = 50;
	size_t r = 20;
	int iter = 2 ;
	bool fail = false;

	static Argument as[] = {
		{ 'p', "-p P", "Set the field characteristic.",         TYPE_INT , &p },
		{ 'n', "-n N", "Set the number of cols in the matrix.", TYPE_INT , &n },
		{ 'm', "-m N", "Set the number of rows in the matrix.", TYPE_INT , &m },
		{ 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iter },
		// { 'f', "-f file", "Set input file", TYPE_STR, &file },
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(argc,argv,as);

	{
		typedef ModularBalanced<double> Field;
		typedef Field::Element Element;
		Field F(p);

		for (int i = 0 ; i < iter ; ++i) {
			fail |= launch_test<Field,FFLAS::FflasUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			fail |= launch_test<Field,FFLAS::FflasUnit,FFLAS::FflasTrans>(F,r,m,n);
			fail |= launch_test<Field,FFLAS::FflasNonUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			fail |= launch_test<Field,FFLAS::FflasNonUnit,FFLAS::FflasTrans>(F,r,m,n);

#if 1 /*  may be bogus */
			fail |= launch_test_append<Field,FFLAS::FflasUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			fail |= launch_test_append<Field,FFLAS::FflasNonUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			// fail |= launch_test_append<Field,FFLAS::FflasUnit,FFLAS::FflasTrans>(F,r,m,n);
			// fail |= launch_test_append<Field,FFLAS::FflasNonUnit,FFLAS::FflasTrans>(F,r,m,n);
#endif
		}
	}

	{
		typedef Modular<double> Field;
		typedef Field::Element Element;
		Field F(p);

		for (int i = 0 ; i < iter ; ++i) {
			fail |= launch_test<Field,FFLAS::FflasUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			fail |= launch_test<Field,FFLAS::FflasUnit,FFLAS::FflasTrans>(F,r,m,n);
			fail |= launch_test<Field,FFLAS::FflasNonUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			fail |= launch_test<Field,FFLAS::FflasNonUnit,FFLAS::FflasTrans>(F,r,m,n);

#if 1 /*  may be bogus */
			fail |= launch_test_append<Field,FFLAS::FflasUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			fail |= launch_test_append<Field,FFLAS::FflasNonUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			// fail |= launch_test_append<Field,FFLAS::FflasUnit,FFLAS::FflasTrans>(F,r,m,n);
			// fail |= launch_test_append<Field,FFLAS::FflasNonUnit,FFLAS::FflasTrans>(F,r,m,n);
#endif

		}
	}

	{
		typedef ModularBalanced<float> Field;
		typedef Field::Element Element;
		Field F(p);

		for (int i = 0 ; i < iter ; ++i) {
			fail |= launch_test<Field,FFLAS::FflasUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			fail |= launch_test<Field,FFLAS::FflasUnit,FFLAS::FflasTrans>(F,r,m,n);
			fail |= launch_test<Field,FFLAS::FflasNonUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			fail |= launch_test<Field,FFLAS::FflasNonUnit,FFLAS::FflasTrans>(F,r,m,n);

#if 1 /*  may be bogus */
			fail |= launch_test_append<Field,FFLAS::FflasUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			fail |= launch_test_append<Field,FFLAS::FflasNonUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			// fail |= launch_test_append<Field,FFLAS::FflasUnit,FFLAS::FflasTrans>(F,r,m,n);
			// fail |= launch_test_append<Field,FFLAS::FflasNonUnit,FFLAS::FflasTrans>(F,r,m,n);
#endif

		}
	}

	{
		typedef Modular<float> Field;
		typedef Field::Element Element;
		Field F(p);

		for (int i = 0 ; i < iter ; ++i) {
			fail |= launch_test<Field,FFLAS::FflasUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			fail |= launch_test<Field,FFLAS::FflasUnit,FFLAS::FflasTrans>(F,r,m,n);
			fail |= launch_test<Field,FFLAS::FflasNonUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			fail |= launch_test<Field,FFLAS::FflasNonUnit,FFLAS::FflasTrans>(F,r,m,n);

#if 1 /*  may be bogus */
			fail |= launch_test_append<Field,FFLAS::FflasUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			fail |= launch_test_append<Field,FFLAS::FflasNonUnit,FFLAS::FflasNoTrans>(F,r,m,n);
			// fail |= launch_test_append<Field,FFLAS::FflasUnit,FFLAS::FflasTrans>(F,r,m,n);
			// fail |= launch_test_append<Field,FFLAS::FflasNonUnit,FFLAS::FflasTrans>(F,r,m,n);
#endif


		}
	}

	return fail ;

}
