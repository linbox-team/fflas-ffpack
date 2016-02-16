/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

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


//-------------------------------------------------------------------------
//      Test suite for the Gaussian elimination routines: LUdivine and PLUQ
//-------------------------------------------------------------------------

#define  __FFLASFFPACK_SEQUENTIAL
#define __LUDIVINE_CUTOFF 1
#include "fflas-ffpack/fflas-ffpack-config.h"
#include <givaro/modular-balanced.h>
#include <iostream>
#include <iomanip>

#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "test-utils.h"

#include "fflas-ffpack/utils/args-parser.h"

using namespace std;
using namespace FFPACK;


/*! Tests the LUdivine routine.
 * @tparam Field Field
 * @tparam Diag  Unit diagonal in U 
 * @tparam Trans 
 * @param F field
 * @param A Matrix (preallocated)
 * @param r rank of A
 * @param m rows
 * @param n cols
 * @param lda leading dim of A
 * @return 0 iff correct, 1 otherwise
 */
template<class Field, FFLAS::FFLAS_DIAG diag, FFLAS::FFLAS_TRANSPOSE trans>
bool test_LUdivine(const Field & F,
				   typename Field::ConstElement_ptr A, size_t lda,
				   size_t r, size_t m, size_t n)
{
	bool fail = false;
	typedef typename Field::Element_ptr Element_ptr ;
	typedef typename Field::Element Element ;
	Element_ptr B = FFLAS::fflas_new(F,m,lda) ;
	FFLAS::fassign(F,m,n,A,lda,B,lda);

	size_t maxP, maxQ ;

	if (trans == FFLAS::FflasTrans){
		maxP = m;
		maxQ = n;
	}
	else{ // trans == FFLAS::FflasNoTrans
		maxP = n;
		maxQ = m;
	}

	size_t * P = FFLAS::fflas_new<size_t>(maxP) ;
	size_t * Q = FFLAS::fflas_new<size_t>(maxQ) ;
	
	size_t R = FFPACK::LUdivine (F, diag, trans, m, n, B, lda, P, Q);

	if (R != r) {
		std::cout << "rank is wrong (expecting " << r << " but got " << R << ")" << std::endl;
		FFLAS::fflas_delete( B );
		FFLAS::fflas_delete( P );
		FFLAS::fflas_delete( Q );
		return fail = true;
	}

	Element_ptr X = FFLAS::fflas_new(F, m, n); // compute X=CUP and check X == A
		/*  Build L,U */
	Element_ptr L, U;
	if (trans == FFLAS::FflasNoTrans){
		L = FFLAS::fflas_new(F, m, m);
		U = FFLAS::fflas_new(F, m, n);

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

			/*  reconstruct the diagonal */
		if (diag == FFLAS::FflasNonUnit) {
			for ( size_t i=0; i<R; ++i ){
				F.assign (*(U+i*(n+1)), *(B+i*(lda+1)));
				F.assign (*(L+Q[i]*m+i), F.one);
			}
		}
		else{ // diag == FFLAS::FflasUnit
			for ( size_t i=0; i<R; ++i ){
				F.assign (*(L+Q[i]*m+i), *(B+Q[i]*lda+i));
				F.assign (*(U+i*(n+1)),one);
			}
		}
		
			/*  Compute CUP */
		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
						m,0,(int) R, U, n, P);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
					  m,n,R, 1.0, L,m, U,n, 0.0, X,n);
	}
	else { /*  trans == FFLAS::FflasTrans */

		L = FFLAS::fflas_new(F, m, n);
		U = FFLAS::fflas_new(F, n, n);


		typename Field::Element zero,one;
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

			/*  Compute LQUP */
		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
						n,0,(int)R, L, n, P);
		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
						m,0,(int)R, L, n, Q);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
					  m,n,n, 1.0, L,n, U,n, 0.0, X,n);
	}
		/*  check equality */
	for (size_t i=0; i<m; ++i) {
		for (size_t j=0; j<n; ++j)
			if (!F.areEqual (*(A+i*lda+j), *(X+i*n+j))){
				std::cerr << std::endl<<" A["<<i<<","<<j<<"] = " << (*(A+i*lda+j))
						  << " LQUP["<<i<<","<<j<<"] = " << (*(X+i*n+j));
				fail|=true;
			}
	}
		// if (fail){
		// 	write_field(F,cerr<<"A = "<<endl,A,m,n,lda);
		// 	write_field(F,cerr<<"LU = "<<endl,B,m,n,lda);
		// 	write_field(F,cerr<<"L = "<<endl,L,m,m,m);
		// 	write_field(F,cerr<<"U = "<<endl,U,m,n,n);
		// }

	FFLAS::fflas_delete( P);
	FFLAS::fflas_delete( L);
	FFLAS::fflas_delete( U);
	FFLAS::fflas_delete( Q);
	FFLAS::fflas_delete( B);
	FFLAS::fflas_delete( X);
	return fail;


}


/*! Verifies that B = PLUQ where A stores [L\U]
 * @tparam Field Field
 * @tparam Diag  Unit diagonal in U
 * @param F field
 * @param A Matrix (preallocated)
 * @param r rank of A
 * @param m rows
 * @param n cols
 * @param lda leading dim of A
 * @return 0 iff correct, 1 otherwise
 */
template<class Field, FFLAS::FFLAS_DIAG diag>
bool verifPLUQ (const Field & F, typename Field::ConstElement_ptr A, size_t lda, 
				typename Field::Element_ptr PLUQ, size_t ldpluq,
				size_t * P, size_t * Q, size_t m, size_t n, size_t R)
{


	typename Field::Element_ptr X = FFLAS::fflas_new (F, m, n);
	typename Field::Element_ptr L = FFLAS::fflas_new (F, m, R);
	typename Field::Element_ptr	U = FFLAS::fflas_new (F, R, n);
	FFLAS::fzero(F, m, R, L, R);
	FFLAS::fzero(F, R, n, U, n);
	
	typename Field::Element zero,one;
	F.init(zero,0.0);
	F.init(one,1.0);
	FFPACK::getTriangular(F, FFLAS::FflasUpper, diag, m,n,R, PLUQ, ldpluq, U, n, true);
	FFPACK::getTriangular(F, FFLAS::FflasLower, (diag==FFLAS::FflasNonUnit)?FFLAS::FflasUnit:FFLAS::FflasNonUnit, 
						  m,n,R, PLUQ, ldpluq, L, R, true);
	FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, R,0,m, L, R, P);
	FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, R,0,n, U, n, Q);
	FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,R, F.one, L,R, U,n, F.zero, X,n);

	bool fail = false;
	for(size_t i=0; i<m; ++i)
		for (size_t j=0; j<n; ++j)
			if (!F.areEqual (*(A+i*lda+j), *(X+i*n+j))){
				std::cerr << std::endl<<" A ["<<i<<","<<j<<"] = " << (*(A+i*lda+j))
						  << " PLUQ ["<<i<<","<<j<<"] = " << (*(X+i*n+j))
						  << std::endl;
				fail=true;
			}
		//write_field(F, std::cerr<<"X = "<<std::endl,X, m,n,n);
	FFLAS::fflas_delete( U);
	FFLAS::fflas_delete( L);
	FFLAS::fflas_delete( X);
	return fail;
}
/*! Tests the LUdivine routine.
 * @tparam Field Field
 * @tparam Diag  Unit diagonal in U 
 * @tparam Trans 
 * @param F field
 * @param A Matrix (preallocated)
 * @param r rank of A
 * @param m rows
 * @param n cols
 * @param lda leading dim of A
 * @return 0 iff correct, 1 otherwise
 */
template<class Field, FFLAS::FFLAS_DIAG diag>
bool test_pluq (const Field & F,
				typename Field::ConstElement_ptr A,
				size_t r, size_t m, size_t n, size_t lda)
{
	bool fail = false;
	typedef typename Field::Element_ptr Element_ptr ;
	Element_ptr B = FFLAS::fflas_new(F,m,lda) ;
	FFLAS::fassign(F,m,n,A,lda,B,lda);

	size_t * P = FFLAS::fflas_new<size_t> (m);
	size_t * Q = FFLAS::fflas_new<size_t> (n);
	
		//write_field(F,std::cerr<<"\n B = \n",B,m,n,lda);
	size_t R = FFPACK::PLUQ (F, diag, m, n, B, lda, P, Q);
		//write_field(F,std::cerr<<"\n PLUQ = \n",B,m,n,lda);

	if (R != r) {
		std::cout << "rank is wrong (expected " << r << " but got " << R << ")" << std::endl;
		FFLAS::fflas_delete (B);
		FFLAS::fflas_delete (P);
		FFLAS::fflas_delete (Q);
		return fail = true;
	}
	fail |=  verifPLUQ<Field,diag> (F,A, lda, B, lda, P, Q, m, n, r);
	FFLAS::fflas_delete (B);
	FFLAS::fflas_delete(P);
	FFLAS::fflas_delete(Q);
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
// template<class Field, FFLAS::FFLAS_DIAG diag, FFLAS::FFLAS_TRANSPOSE trans>
// bool test_lu_append(const Field & F,
// 		    const typename Field::Element_ptr A,
// 		    const typename Field::Element_ptr B,
// 		    size_t m, size_t n, size_t k, size_t lda)
// {
// 	FFLASFFPACK_check(n<=lda);

// 	bool fail = false;
// 	size_t M = m + k ;
// 	typedef typename Field::Element Element ;
// 	Element_ptr Acop = FFLAS::fflas_new(F, m, lda) ;
// 	FFLAS::fassign(F,m,n,A,lda,Acop,lda) ;

// 	Element_ptr Bcop = FFLAS::fflas_new(F, k, lda) ;
// 	FFLAS::fassign(F,k,n,B,lda,Bcop,lda) ;

// 	Element_ptr Append = FFLAS::fflas_new (F, M, lda);
// 	FFLAS::fassign(F,m,n,A,lda,Append,lda) ;
// 	FFLAS::fassign(F,k,n,B,lda,Append+m*lda,lda) ;

// #if 0 /*  paranoid check */
// 	for (size_t i = 0 ; i < m ; ++i) {
// 		for (size_t j = 0 ; j < n ; ++j) {
// 			FFLASFFPACK_check(Append[i*lda+j]==A[i*lda+j]);
// 		}
// 	}
// 	for (size_t i = 0 ; i < k ; ++i) {
// 		for (size_t j = 0 ; j < n ; ++j) {
// 			FFLASFFPACK_check(Append[(i+m)*lda+j]==B[i*lda+j]);
// 		}
// 	}
// #endif

// 	Element_ptr Afull = FFLAS::fflas_new(F, M, lda);
// 	FFLAS::fassign(F,M,n,Append,lda,Afull,lda) ;
// 	// FFLAS::fassign(F,m,n,A,lda,Afull,lda) ;
// 	// FFLAS::fassign(F,k,n,B,lda,Afull+m*lda,lda) ;

// #if 0
// std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
// 	for (size_t i = 0 ; i < m ; ++i) {
// 		for (size_t j = 0 ; j < n ; ++j) {
// 			std::cout << Append[i*lda+j] << "(" << A[i*lda+j] << ") " ;
// 		} std::cout << std::endl;
// 	}
// std::cout << "-----------------------------------" << std::endl;
// 	for (size_t i = 0 ; i < k ; ++i) {
// 		for (size_t j = 0 ; j < n ; ++j) {
// 			std::cout << Append[(i+m)*lda+j] ;
// 			std::cout << "(" << B[i*lda+j] << ") "  ;
// 		}std::cout << std::endl;
// 	}
// std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::flush << std::endl;
// #endif



// #if 0
// 	for (size_t i = 0 ; i < m ; ++i)
// 		for (size_t j = 0 ; j < n ; ++j)
// 			FFLASFFPACK_check(Acop[i*lda+j]==A[i*lda+j]);
// 	for (size_t i = 0 ; i < k ; ++i)
// 		for (size_t j = 0 ; j < n ; ++j)
// 			FFLASFFPACK_check(Bcop[i*lda+j]==B[i*lda+j]);
// 	for (size_t i = 0 ; i < M ; ++i)
// 		for (size_t j = 0 ; j < n ; ++j)
// 			if (i < m)
// 				FFLASFFPACK_check(Afull[i*lda+j]==A[i*lda+j]);
// 			else
// 				FFLASFFPACK_check(Afull[i*lda+j]==B[(i-m)*lda+j]);
// #endif




// 	size_t maxP, maxQ ;

// 	if (trans == FFLAS::FflasTrans){
// 		maxP = M;
// 		maxQ = n;
// 	}
// 	else{ // trans == FFLAS::FflasNoTrans
// 		maxP = n;
// 		maxQ = M;
// 	}

// 	size_t * P = FFLAS::fflas_new<size_t>(maxP) ;
// 	size_t * Q = FFLAS::fflas_new<size_t>(maxQ) ;

// 	size_t * PP = FFLAS::fflas_new<size_t>(maxP) ;
// 	size_t * QQ = FFLAS::fflas_new<size_t>(maxQ) ;

// 	/* valgrind says the following leaks. Just incroyable. */
// 	size_t R  = FFPACK::LUdivine (F, diag, trans, M, n, Append, lda, PP, QQ);

// 	size_t R1 = FFPACK::LUdivine (F, diag, trans, m, n, Acop,   lda, P, Q);

// 	size_t R2 = FFPACK::LUpdate  (F,diag,trans,m,n,Acop,lda,R1,k,Bcop,lda,P,Q,
// 				      FFPACK::FfpackLQUP);
// #if 0
// 	std::cout << "P := [ " ;
// 	for (size_t i = 0 ; i < maxP ; ++i)
// 		std::cout << P[i] << " " ;
// 	std::cout << ']' << std::endl;
// 	std::cout << "Q := [ ";
// 	for (size_t i = 0 ; i < maxQ ; ++i)
// 		std::cout << Q[i] << " " ;
// 	std::cout << ']' << std::endl;
// 	std::cout << "PP := [ ";
// 	for (size_t i = 0 ; i < maxP ; ++i)
// 		std::cout << PP[i] << " " ;
// 	std::cout << ']' << std::endl;
// 	std::cout << "QQ := [ ";
// 	for (size_t i = 0 ; i < maxQ ; ++i)
// 		std::cout << QQ[i] << " " ;
// 	std::cout << ']' << std::endl;
// #endif

// 	if (R2 != R) {
// 		std::cout << "error, bad rank " << R2 << " <> " << R << " (expected) " << std::endl;
// 		FFLAS::fflas_delete( Bcop );
// 		FFLAS::fflas_delete( Acop );
// 		FFLAS::fflas_delete( Append );
// 		FFLAS::fflas_delete( PP);
// 		FFLAS::fflas_delete( QQ);
// 		FFLAS::fflas_delete( P );
// 		FFLAS::fflas_delete( Q );
// 		return fail=true;

// 	}

// 	// compute C=LQUP and check C == A
// 	Element_ptr C = FFLAS::fflas_new (F, M, lda);
// 	/*  Build L,U */
// 	Element_ptr L, U;
// 	if (trans == FFLAS::FflasNoTrans){
// 		L = FFLAS::fflas_new(F, M, M);
// 		U = FFLAS::fflas_new(F, M, n);

// 		typename Field::Element zero,one;
// 		F.init(zero,0.0);
// 		F.init(one,1.0);
// 		/*  build U */
// 		for (size_t i=0; i<R; ++i){
// 			for (size_t j=0; j<i; ++j)
// 				F.assign ( *(U + i*n + j), zero);
// 			for (size_t j=i+1; j<n; ++j)
// 				if (i < m)
// 					F.assign (*(U + i*n + j), *(Acop+ i*lda+j));
// 				else
// 					F.assign (*(U + i*n + j), *(Bcop+ (i-m)*lda+j));
// 		}

// 		for (size_t i=R;i<M; ++i) {
// 			for (size_t j=0; j<n; ++j)
// 				F.assign(*(U+i*n+j), zero);
// 		}
// 		/*  build L */
// 		for ( size_t i=0; i<M; ++i ){
// 			size_t j=0;
// 			for (; j< ((i<R)?i:R) ; ++j ) {
// 				if (i<m)
// 					F.assign( *(L + i*M+j), *(Acop+i*lda+j));
// 				else
// 					F.assign( *(L + i*M+j), *(Bcop+(i-m)*lda+j));
// 			}
// 			for (; j<M; ++j )
// 				F.assign( *(L+i*M+j), zero);
// 		}

// 		// write_field(F,cerr<<"L = "<<endl,L,m,m,m);
// 		//write_field(F,cerr<<"U = "<<endl,U,m,n,n);
// 		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
// 				M,0,(int)R, L, M, Q);
// 		for ( size_t i=0; i<M; ++i )
// 			F.assign(*(L+i*(M+1)), one);

// 		/*  reconstruct the diagonal */
// 		//write_field(F,cerr<<"L = "<<endl,L,m,m,m);
// 		//write_field(F,cerr<<"U = "<<endl,U,m,n,n);
// 		if (diag == FFLAS::FflasNonUnit) {
// 			for ( size_t i=0; i<R; ++i ) {
// 				if (i<m)
// 					F.assign (*(U+i*(n+1)), *(Acop+i*(lda+1)));
// 				else
// 					F.assign (*(U+i*(n+1)), *(Bcop+(i-m)*(lda+1)));
// 			}
// 		}
// 		else{ // diag == FFLAS::FflasUnit
// 			for ( size_t i=0; i<R; ++i ){
// 				if (Q[i] < m)
// 					*(L+Q[i]*(M+1)) = *(Acop+Q[i]*lda+i);
// 				else
// 					*(L+Q[i]*(M+1)) = *(Bcop+(Q[i]-m)*lda+i);
// 				F.assign (*(U+i*(n+1)),one);
// 			}
// 		}
// 		// write_field(F,cerr<<"L = "<<endl,L,(int)M,(int)M,(int)M);
// 		// write_field(F,cerr<<"U = "<<endl,U,(int)M,(int)n,(int)n);

// 		/*  Compute LQUP */
// 		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
// 				M,0,(int) R, U, n, P);
// 		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
// 				n,0,(int)R, U, n, Q);
// 		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
// 			      M,n,M, 1.0, L,M, U,n, 0.0, C,lda);
// 		//FFLAS::fflas_delete( A);
// 	}
// #if 0 /*  not working */
// 	else { /*  trans == FFLAS::FflasTrans */

// 		L = FFLAS::fflas_new(F, M, n);
// 		U = FFLAS::fflas_new(F, n, n);


// 		typename Field::Element zero,one;
// 		F.init(zero,0.0);
// 		F.init(one,1.0);
// 		/*  build L */
// 		for (size_t i=0; i<R; ++i){
// 			for (size_t j=0; j<i; ++j)
// 				F.assign ( *(L + i + j*n), zero);
// 			for (size_t j=i+1; j<M; ++j) {
// 				if (i < m)
// 					F.assign (*(L + i + j*n), *(Acop+ i+j*lda));
// 				else
// 					F.assign (*(L + i + j*n), *(Bcop+ (i-m)+j*lda));
// 			}
// 		}
// 		for (size_t i=R;i<n; ++i) {
// 			for (size_t j=0; j<M; ++j)
// 				F.assign(*(L+i+j*n), zero);
// 		}
// 		/*  build U */
// 		for ( size_t i=0; i<n; ++i ){
// 			size_t j=0;
// 			for (;  j< ((i<R)?i:R) ; ++j ) {
// 				if (i < m)
// 					F.assign( *(U + i+j*n), *(Acop+i+j*lda));
// 				else
// 					F.assign( *(U + i+j*n), *(Bcop+(i-m)+j*lda));
// 			}
// 			for (; j<n; ++j )
// 				F.assign( *(U+i+j*n), zero);
// 		}
// 		//write_field(F,cerr<<"L = "<<endl,L,m,n,n);
// 		//write_field(F,cerr<<"U = "<<endl,U,n,n,n);

// 		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans,
// 				n,0,(int)R, U, n, Q);

// 		for (size_t i=0; i<n; ++i)
// 			F.assign (*(U+i*(n+1)),one);

// 		/*  reconstruct the diagonal */
// 		if (diag == FFLAS::FflasNonUnit) {
// 			for ( size_t i=0; i<R; ++i ) {
// 				if (i < m)
// 					F.assign (*(L+i*(n+1)), *(Acop+i*(lda+1)));
// 				else
// 					F.assign (*(L+i*(n+1)), *(Bcop+(i-m)*(lda+1)));
// 			}
// 		}
// 		else{ // diag == FFLAS::FflasUnit
// 			for ( size_t i=0; i<R; ++i ){
// 				if (i<m)
// 					*(U+Q[i]*(n+1)) = *(Acop+Q[i]+i*lda);
// 				else
// 					*(U+Q[i]*(n+1)) = *(Bcop+Q[i]+(i-m)*lda);
// 				F.assign (*(L+i*(n+1)),one);
// 			}
// 		}
// 		// write_field(F,cerr<<"L = "<<endl,L,(int)M,(int)n,(int)n);
// 		// write_field(F,cerr<<"U = "<<endl,U,(int)n,(int)n,(int)n);

// 		/*  Compute LQUP */
// 		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
// 				n,0,(int)R, L, n, P);
// 		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
// 				M,0,(int)R, L, n, Q);
// 		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
// 			      M,n,n, 1.0, L,n, U,n, 0.0, C,lda);
// 	}
// #endif
// #if 0 /*  check CC == LL UU */
// 	Element_ptr LL, UU;
// 	Element_ptr CC = FFLAS::fflas_new (F, M, lda);
// 	if (trans == FFLAS::FflasNoTrans){
// 		LL = FFLAS::fflas_new (F, M, M);
// 		UU = FFLAS::fflas_new (F, M, n);

// 		Element zero,one;
// 		F.init(zero,0.0);
// 		F.init(one,1.0);
// 		/*  build U */
// 		for (size_t i=0; i<R; ++i){
// 			for (size_t j=0; j<i; ++j)
// 				F.assign ( *(UU + i*n + j), zero);
// 			for (size_t j=i+1; j<n; ++j)
// 				F.assign (*(UU + i*n + j), *(Append+ i*lda+j));
// 		}

// 		for (size_t i=R;i<M; ++i) {
// 			for (size_t j=0; j<n; ++j)
// 				F.assign(*(UU+i*n+j), zero);
// 		}
// 		/*  build L */
// 		for ( size_t i=0; i<M; ++i ){
// 			size_t j=0;
// 			for (; j< ((i<R)?i:R) ; ++j ) {
// 				F.assign( *(LL + i*M+j), *(Append+i*lda+j));
// 			}
// 			for (; j<M; ++j )
// 				F.assign( *(LL+i*M+j), zero);
// 		}

// 		// write_field(F,cerr<<"LL = "<<endl,LL,m,m,m);
// 		//write_field(F,cerr<<"UU = "<<endl,UU,m,n,n);
// 		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
// 				M,0,(int)R, LL, M, Q);
// 		for ( size_t i=0; i<M; ++i )
// 			F.assign(*(LL+i*(M+1)), one);

// 		/*  reconstruct the diagonal */
// 		//write_field(F,cerr<<"LL = "<<endl,LL,m,m,m);
// 		//write_field(F,cerr<<"UU = "<<endl,UU,m,n,n);
// 		if (diag == FFLAS::FflasNonUnit) {
// 			for ( size_t i=0; i<R; ++i ) {
// 				F.assign (*(UU+i*(n+1)), *(Append+i*(lda+1)));
// 			}
// 		}
// 		else{ // diag == FFLAS::FflasUnit
// 			for ( size_t i=0; i<R; ++i ){
// 				*(LL+Q[i]*(M+1)) = *(Append+Q[i]*lda+i);
// 				F.assign (*(UU+i*(n+1)),one);
// 			}
// 		}
// 		write_field(F,cerr<<"L = "<<endl,LL,(int)M,(int)M,(int)M);
// 		write_field(F,cerr<<"U = "<<endl,UU,(int)M,(int)n,(int)n);

// 		/*  Compute LQUP */
// 		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
// 				M,0,(int) R, UU, n, P);
// 		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
// 				n,0,(int)R, UU, n, Q);
// 		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
// 			      M,n,M, 1.0, LL,M, UU,n, 0.0, CC,lda);
// 		//FFLAS::fflas_delete( A);
// 	}
// 	else { /*  trans == FFLAS::FflasTrans */

// 		LL = FFLAS::fflas_new(F, M, n);
// 		UU = FFLAS::fflas_new(F, n, n);


// 		typename Field::Element zero,one;
// 		F.init(zero,0.0);
// 		F.init(one,1.0);
// 		/*  build L */
// 		for (size_t i=0; i<R; ++i){
// 			for (size_t j=0; j<i; ++j)
// 				F.assign ( *(LL + i + j*n), zero);
// 			for (size_t j=i+1; j<M; ++j) {
// 				F.assign (*(LL + i + j*n), *(Append+ i+j*lda));
// 			}
// 		}
// 		for (size_t i=R;i<n; ++i) {
// 			for (size_t j=0; j<M; ++j)
// 				F.assign(*(LL+i+j*n), zero);
// 		}
// 		/*  build UU */
// 		for ( size_t i=0; i<n; ++i ){
// 			size_t j=0;
// 			for (;  j< ((i<R)?i:R) ; ++j ) {
// 				F.assign( *(UU + i+j*n), *(Append+i+j*lda));
// 			}
// 			for (; j<n; ++j )
// 				F.assign( *(UU+i+j*n), zero);
// 		}
// 		// 		write_field(F,cerr<<"LL = "<<endl,LL,m,n,n);
// 		// 		write_field(F,cerr<<"UU = "<<endl,UU,n,n,n);

// 		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans,
// 				n,0,(int)R, UU, n, Q);

// 		for (size_t i=0; i<n; ++i)
// 			F.assign (*(UU+i*(n+1)),one);

// 		/*  reconstruct the diagonal */
// 		if (diag == FFLAS::FflasNonUnit) {
// 			for ( size_t i=0; i<R; ++i ) {
// 				F.assign (*(LL+i*(n+1)), *(Append+i*(lda+1)));
// 			}
// 		}
// 		else{ // diag == FFLAS::FflasUnit
// 			for ( size_t i=0; i<R; ++i ){
// 				*(UU+Q[i]*(n+1)) = *(Append+Q[i]+i*lda);
// 				F.assign (*(LL+i*(n+1)),one);
// 			}
// 		}
// 		write_field(F,cerr<<"LL = "<<endl,LL,(int)M,(int)n,(int)n);
// 		write_field(F,cerr<<"UU = "<<endl,UU,(int)n,(int)n,(int)n);

// 		/*  Compute LQUP */
// 		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
// 				n,0,(int)R, LL, n, P);
// 		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
// 				M,0,(int)R, LL, n, Q);
// 		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
// 			      M,n,n, 1.0, LL,n, UU,n, 0.0, CC,lda);
// 	}
// 	for (size_t i=0; i<M; ++i) {
// 		for (size_t j=0; j<n; ++j)
// 			if (!F.areEqual (*(Afull+i*lda+j), *(CC+i*lda+j))){
// 				std::cerr << " A["<<i<<","<<j<<"]    = " << (*(Afull+i*lda+j))
// 				<< " LQUP["<<i<<","<<j<<"] = " << (*(CC+i*lda+j))
// 				<< endl << "xxxx" << endl;
// 				fail|=true;
// 			}
// 	}

// #endif

// 	/*  check equality */
// 	for (size_t i=0; i<M; ++i) {
// 		for (size_t j=0; j<n; ++j)
// 			if (!F.areEqual (*(Afull+i*lda+j), *(C+i*lda+j))){
// 				std::cerr << " A["<<i<<","<<j<<"]    = " << (*(Afull+i*lda+j))
// 				<< " LQUP(append)["<<i<<","<<j<<"] = " << (*(C+i*lda+j))
// 				<< endl;
// 				fail|=true;
// 			}
// 	}

// 	FFLAS::fflas_delete( PP);
// 	FFLAS::fflas_delete( P);
// 	FFLAS::fflas_delete( L);
// 	FFLAS::fflas_delete( U);
// 	FFLAS::fflas_delete( Q);
// 	FFLAS::fflas_delete( QQ);
// 	FFLAS::fflas_delete( Acop);
// 	FFLAS::fflas_delete( Bcop);
// 	FFLAS::fflas_delete( Append);
// 	FFLAS::fflas_delete( Afull);
// 	FFLAS::fflas_delete( C);

// 	return fail;


// }



template<class Field, FFLAS::FFLAS_DIAG diag, FFLAS::FFLAS_TRANSPOSE trans>
bool launch_test(const Field & F,
				 size_t r,
				 size_t m, size_t n)
{
		//typedef typename Field::Element Element ;
	typedef typename Field::Element_ptr Element_ptr ;
	bool fail = false ;
	{ /*  user given and lda bigger */
		size_t lda = n+10 ;
		Element_ptr A = FFLAS::fflas_new (F, m, lda);
		RandomMatrixWithRankandRandomRPM(F,A,lda,r,m,n);
		fail |= test_LUdivine<Field,diag,trans>(F,A,lda,r,m,n);
		fail |= test_pluq<Field,diag>(F,A,r,m,n,lda);
		if (fail) std::cout << "failed at big lda" << std::endl;
		FFLAS::fflas_delete( A );
	}
	{ /*  user given and lda bigger. Rank is max */
		size_t lda = n+10 ;
		size_t R = std::min(m,n);
		Element_ptr A = FFLAS::fflas_new (F, m, lda);
		RandomMatrixWithRankandRandomRPM(F,A,lda,R,m,n);
		fail |= test_LUdivine<Field,diag,trans>(F,A,lda,R,m,n);
		fail |= test_pluq<Field,diag>(F,A,R,m,n,lda);
		if (fail) std::cout << "failed at big lda max rank" << std::endl;
		FFLAS::fflas_delete( A );
	}
	{ /*  user given and lda bigger. Rank is min */
		size_t lda = n+10 ;
		size_t R = 0;
		Element_ptr A = FFLAS::fflas_new (F, m, lda);
		RandomMatrixWithRankandRandomRPM(F,A,lda,R,m,n);
		fail |= test_LUdivine<Field,diag,trans>(F,A,lda,R,m,n);
		fail |= test_pluq<Field,diag>(F,A,R,m,n,lda);
		if (fail) std::cout << "failed at big lda, rank 0" << std::endl;
		FFLAS::fflas_delete( A );
	}
	{ /*  square  */
		size_t M = std::max(m,n);
		size_t N = M ;
		size_t R = M/2 ;
		size_t lda = N+10 ;
		Element_ptr A = FFLAS::fflas_new (F, M, lda);
		RandomMatrixWithRankandRandomRPM(F,A,lda,R,M,N);
		fail |= test_LUdivine<Field,diag,trans>(F,A,lda,R,M,N);
		fail |= test_pluq<Field,diag>(F,A,R,M,N,lda);
		if (fail) std::cout << "failed at square" << std::endl;
		FFLAS::fflas_delete( A );
	}
	{ /*  wide  */
		size_t M = std::max(m,n);
		size_t N = 2*M ;
		size_t R = 3*M/4 ;
		size_t lda = N+5 ;
		Element_ptr A = FFLAS::fflas_new (F, M, lda);
		RandomMatrixWithRankandRandomRPM(F,A,lda,R,M,N);
		fail |= test_LUdivine<Field,diag,trans>(F,A,lda,R,M,N);
		fail |= test_pluq<Field,diag>(F,A,R,M,N,lda);
		if (fail) std::cout << "failed at wide" << std::endl;
		FFLAS::fflas_delete( A );
	}
	{ /*  narrow  */
		size_t M = std::max(m,n);
		size_t N = M/2 ;
		size_t R = 3*M/8 ;
		size_t lda = N+5 ;
		Element_ptr A = FFLAS::fflas_new (F, M, lda);
		RandomMatrixWithRankandRandomRPM(F,A,lda,R,M,N);
		fail |= test_LUdivine<Field,diag,trans>(F,A,lda,R,M,N);
		fail |= test_pluq<Field,diag>(F,A,R,M,N,lda);
		if (fail) std::cout << "failed at narrow" << std::endl;
		FFLAS::fflas_delete( A );
	}
	return !fail;
}

// template<class Field, FFLAS::FFLAS_DIAG diag, FFLAS::FFLAS_TRANSPOSE trans>
// bool launch_test_append(const Field & F,
// 			size_t r,
// 			size_t m, size_t n)
// {
// 	typedef typename Field::Element Element ;
// 	bool fail = false ;
// 	{ /*  user given and lda bigger */
// 		size_t lda = n+10 ;
// 		size_t k = m/2+1 ;
// 		Element_ptr A = FFLAS::fflas_new (F, m, lda);
// 		Element_ptr B = FFLAS::fflas_new (F, k, lda);
// 		RandomMatrixWithRank(F,A,lda,r,m,n);
// 		RandomMatrixWithRank(F,B,lda,k/2+1,k,n);
// 		fail |= test_lu_append<Field,diag,trans>(F,A,B,m,n,k,lda);
// 		if (fail) std::cout << "failed" << std::endl;
// 		FFLAS::fflas_delete( A );
// 		FFLAS::fflas_delete( B );
// 	}
// 	{ /*  user given and lda bigger. Rank is max */
// 		size_t lda = n+10 ;
// 		size_t R = std::min(m,n);
// 		size_t k = m/2+1 ;
// 		Element_ptr A = FFLAS::fflas_new (F, m, lda);
// 		Element_ptr B = FFLAS::fflas_new (F, k, lda);
// 		RandomMatrixWithRank(F,A,lda,R,m,n);
// 		RandomMatrixWithRank(F,B,lda,k/2+1,k,n);
// 		fail |= test_lu_append<Field,diag,trans>(F,A,B,m,n,k,lda);
// 		if (fail) std::cout << "failed" << std::endl;
// 		FFLAS::fflas_delete( A );
// 		FFLAS::fflas_delete( B );
// 	}
// 	{ /*  user given and lda bigger. Appended Rank is min */
// 		size_t lda = n+10 ;
// 		size_t R = std::min(m,n);
// 		size_t k = m/2+1 ;
// 		Element_ptr A = FFLAS::fflas_new (F, m, lda);
// 		Element_ptr B = FFLAS::fflas_new (F, k, lda);
// 		RandomMatrixWithRank(F,A,lda,R,m,n);
// 		RandomMatrixWithRank(F,B,lda,0,k,n);
// 		fail |= test_lu_append<Field,diag,trans>(F,A,B,m,n,k,lda);
// 		if (fail) std::cout << "failed" << std::endl;
// 		FFLAS::fflas_delete( A );
// 		FFLAS::fflas_delete( B );
// 	}
// 	{ /*  user given and lda bigger. Rank is min */
// 		size_t lda = n+10 ;
// 		size_t R = 0;
// 		size_t k = m/2+1 ;
// 		Element_ptr A = FFLAS::fflas_new (F, m, lda);
// 		Element_ptr B = FFLAS::fflas_new (F, k, lda);
// 		RandomMatrixWithRank(F,A,lda,R,m,n);
// 		RandomMatrixWithRank(F,B,lda,k/2+1,k,n);
// 		fail |= test_lu_append<Field,diag,trans>(F,A,B,m,n,k,lda);
// 		if (fail) std::cout << "failed" << std::endl;
// 		FFLAS::fflas_delete( A );
// 		FFLAS::fflas_delete( B );
// 	}
// 	{ /*  square  */
// 		size_t M = std::max(m,n);
// 		size_t N = M ;
// 		size_t R = M/2 ;
// 		size_t lda = N+10 ;
// 		size_t k = R ;
// 		Element_ptr A = FFLAS::fflas_new (F, M, lda);
// 		Element_ptr B = FFLAS::fflas_new (F, k, lda);
// 		RandomMatrixWithRank(F,A,lda,R,M,N);
// 		RandomMatrixWithRank(F,B,lda,R/2,k,N);
// 		fail |= test_lu_append<Field,diag,trans>(F,A,B,M,N,k,lda);
// 		if (fail) std::cout << "failed" << std::endl;
// 		FFLAS::fflas_delete( A );
// 		FFLAS::fflas_delete( B );
// 	}
// 	{ /*  wide  */
// 		size_t M = std::max(m,n);
// 		size_t N = 2*M ;
// 		size_t R = M/2 ;
// 		size_t k = R ;
// 		size_t lda = N+10 ;
// 		Element_ptr A = FFLAS::fflas_new (F, M, lda);
// 		Element_ptr B = FFLAS::fflas_new (F, k, lda);
// 		RandomMatrixWithRank(F,A,lda,R,M,N);
// 		RandomMatrixWithRank(F,B,lda,k/2,k,N);
// 		fail |= test_lu_append<Field,diag,trans>(F,A,B,M,N,k,lda);
// 		if (fail) std::cout << "failed" << std::endl;
// 		FFLAS::fflas_delete( A );
// 		FFLAS::fflas_delete( B );
// 	}
// 	//! @bug leaks :
// #if 0 /*  leak here */
// 	{ /*  narrow  */
// 		size_t M = std::max(m,n);
// 		size_t N = M/2 ;
// 		size_t R = M/3 ;
// 		size_t k = N ;
// 		size_t lda = N+10 ;
// 		Element_ptr A = FFLAS::fflas_new (F, M, lda);
// 		Element_ptr B = FFLAS::fflas_new (F, k, lda);
// 		RandomMatrixWithRank(F,A,lda,R,M,N);
// 		RandomMatrixWithRank(F,A,lda,std::min(k/2,M/2),k,N);
// 		fail |= test_lu_append<Field,diag,trans>(F,A,B,M,N,k,lda);
// 		if (fail) std::cout << "failed" << std::endl;
// 		FFLAS::fflas_delete( A );
// 		FFLAS::fflas_delete( B );
// 	}
// #endif

// 	return fail;
// }


template<class Field>
bool run_with_field(Givaro::Integer q, uint64_t b, size_t m, size_t n, size_t r, size_t iters){
	bool ok = true ;
	int nbit=(int)iters;
	
	while (ok &&  nbit){
			// choose Field 
		Field* F= chooseField<Field>(q,b);
		if (F==nullptr)
			return true;
		std::ostringstream oss;
		F->write(oss);
		
		std::cout.fill('.');
		std::cout<<"Checking ";
		std::cout.width(40);
		std::cout<<oss.str();
		std::cout<<" ... ";


		ok&= launch_test<Field,FFLAS::FflasUnit,FFLAS::FflasNoTrans>    (*F,r,m,n);
		ok&= launch_test<Field,FFLAS::FflasUnit,FFLAS::FflasTrans>      (*F,r,m,n);
		ok&= launch_test<Field,FFLAS::FflasNonUnit,FFLAS::FflasNoTrans> (*F,r,m,n);
		ok&= launch_test<Field,FFLAS::FflasNonUnit,FFLAS::FflasTrans>   (*F,r,m,n);		

#if 0 /*  may be bogus */
		ok&= launch_test_append<Field,FFLAS::FflasUnit,FFLAS::FflasNoTrans>   (*F,r,m,n);
		ok&= launch_test_append<Field,FFLAS::FflasNonUnit,FFLAS::FflasNoTrans>(*F,r,m,n);
		ok&= launch_test_append<Field,FFLAS::FflasUnit,FFLAS::FflasTrans>     (*F,r,m,n);
		ok&= launch_test_append<Field,FFLAS::FflasNonUnit,FFLAS::FflasTrans>  (*F,r,m,n);
#endif
		nbit--;
		if ( !ok )
				//std::cout << "\033[1;31mFAILED\033[0m "<<std::endl;		
			std::cout << "FAILED "<<std::endl;		
		else
				//std::cout << "\033[1;32mPASSED\033[0m "<<std::endl;
			std::cout << "PASSED "<<std::endl;
		delete F;
	}
	return ok;
}

int main(int argc, char** argv)
{
	cerr<<setprecision(20);
	static Givaro::Integer q=-1;
	static size_t b=0;
	static size_t m=120;
	static size_t n=120;
	static size_t r=80;
	static size_t iters=2;
	static bool loop=false;
	static Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
		{ 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
		{ 'm', "-m M", "Set the row dimension of the matrix.",      TYPE_INT , &m },
		{ 'n', "-n N", "Set the column dimension of the matrix.", TYPE_INT , &n },
		{ 'r', "-r R", "Set the rank.", TYPE_INT , &r },
		{ 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
		{ 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(argc,argv,as);

	if (r > std::min (m,n)) 
		r = std::min (m, n);

	bool ok=true;
	do{
		ok&=run_with_field<Givaro::Modular<float> >           (q,b,m,n,r,iters);
		ok&=run_with_field<Givaro::Modular<double> >          (q,b,m,n,r,iters);
		ok&=run_with_field<Givaro::ModularBalanced<float> >   (q,b,m,n,r,iters);
		ok&=run_with_field<Givaro::ModularBalanced<double> >  (q,b,m,n,r,iters);
		ok&=run_with_field<Givaro::Modular<int32_t> >         (q,b,m,n,r,iters);
		ok&=run_with_field<Givaro::ModularBalanced<int32_t> > (q,b,m,n,r,iters);
		ok&=run_with_field<Givaro::Modular<int64_t> >         (q,b,m,n,r,iters);
		ok&=run_with_field<Givaro::ModularBalanced<int64_t> > (q,b,m,n,r,iters);
		ok&=run_with_field<Givaro::Modular<Givaro::Integer> > (q,(b?b:512),m/6,n/6,r/6,iters);		
	} while (loop && ok);

	return !ok;
}
