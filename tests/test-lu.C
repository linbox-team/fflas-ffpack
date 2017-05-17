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

// #define MONOTONIC_CYLCES
// #define MONOTONIC_MOREPIVOTS
// #define MONOTONIC_FEWPIVOTS
#ifdef MONOTONIC_CYLCES
  #define MONOTONIC_APPLYP
#endif
#ifdef MONOTONIC_MOREPIVOTS
  #define MONOTONIC_APPLYP
#endif
#ifdef MONOTONIC_FEWPIVOTS
  #define MONOTONIC_APPLYP
#endif

#define BASECASE_K 37 // Forcing a lower base case to be able to test a few recursive steps with smallish dimensions


#define  __FFLASFFPACK_SEQUENTIAL
#define __LUDIVINE_CUTOFF 1
#include "fflas-ffpack/fflas-ffpack-config.h"
#include <givaro/modular-balanced.h>
#include <iostream>
#include <iomanip>
Givaro::Timer tperm, tgemm, tBC, ttrsm,trest,timtot;
size_t mvcnt = 0;
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "test-utils.h"

#include "fflas-ffpack/utils/args-parser.h"

using namespace std;
using namespace FFPACK;
using namespace FFLAS;

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
template<class Field, FFLAS::FFLAS_DIAG diag, FFLAS_TRANSPOSE trans>
bool test_LUdivine(const Field & F,
				   typename Field::ConstElement_ptr A, size_t lda,
				   size_t r, size_t m, size_t n)
{
	bool fail = false;
	typename Field::Element_ptr B = fflas_new(F,m,lda) ;
	fassign(F,m,n,A,lda,B,lda);

	size_t maxP, maxQ ;

	if (trans == FflasTrans){
		maxP = m;
		maxQ = n;
	}
	else{ // trans == FflasNoTrans
		maxP = n;
		maxQ = m;
	}

	size_t * P = fflas_new<size_t>(maxP) ;
	size_t * Q = fflas_new<size_t>(maxQ) ;
	
	size_t R = LUdivine (F, diag, trans, m, n, B, lda, P, Q);

	if (R != r) {
		std::cout << "rank is wrong (expecting " << r << " but got " << R << ")" << std::endl;
		fflas_delete( B );
		fflas_delete( P );
		fflas_delete( Q );
		return fail = true;
	}

		/*  Build L,U */
	typename Field::Element_ptr L = fflas_new(F, m, R);
	typename Field::Element_ptr U = fflas_new(F, R, n);

	if (trans == FflasNoTrans){
		getTriangular (F, FflasUpper, diag, m, n, R, B, lda, U, n, true);
		getEchelonForm (F, FflasLower, (diag==FflasNonUnit)?FflasUnit:FflasNonUnit, m, n, R, Q, B, lda, L, R, true);
		applyP (F, FflasRight, FflasNoTrans, R, 0, R, U, n, P);
	} else {
		getTriangular (F, FflasLower, diag, m, n, R, B, lda, L, R, true);
		getEchelonForm (F, FflasUpper, (diag==FflasNonUnit)?FflasUnit:FflasNonUnit, m, n, R, Q, B, lda, U, n, true);
		applyP (F, FflasLeft, FflasTrans, R, 0, R, L, R, P);
	}
	fgemm (F, FflasNoTrans, FflasNoTrans, m, n, R, 1.0, L, R, U, n, 0.0, B, lda);

	fail |= !fequal(F, m, n, A, lda, B, lda);
	
	if (fail){
		write_field(F,cerr<<"A = "<<endl,A,m,n,lda);
		write_field(F,cerr<<"LU = "<<endl,B,m,n,lda);
		write_field(F,cerr<<"L = "<<endl,L,m,R,R);
		write_field(F,cerr<<"U = "<<endl,U,R,n,n);
	}

	fflas_delete( P);
	fflas_delete( L);
	fflas_delete( U);
	fflas_delete( Q);
	fflas_delete( B);
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
template<class Field, FFLAS_DIAG diag>
bool verifPLUQ (const Field & F, typename Field::ConstElement_ptr A, size_t lda, 
				typename Field::Element_ptr PLUQ, size_t ldpluq,
				size_t * P, size_t * Q, size_t m, size_t n, size_t R)
{


	typename Field::Element_ptr X = fflas_new (F, m, n);
	typename Field::Element_ptr L = fflas_new (F, m, R);
	typename Field::Element_ptr	U = fflas_new (F, R, n);
	fzero(F, m, R, L, R);
	fzero(F, R, n, U, n);
	
	typename Field::Element zero,one;
	F.init(zero,0.0);
	F.init(one,1.0);
	// write_field(F,std::cerr<<"PLUQ = "<<std::endl,PLUQ,m,n,ldpluq);
	getTriangular(F, FflasUpper, diag, m,n,R, PLUQ, ldpluq, U, n, true);
	getTriangular(F, FflasLower, (diag==FflasNonUnit)?FflasUnit:FflasNonUnit, 
						  m,n,R, PLUQ, ldpluq, L, R, true);
	applyP( F, FflasLeft, FflasTrans, R,0,m, L, R, P);
	applyP (F, FflasRight, FflasNoTrans, R,0,n, U, n, Q);
	fgemm (F, FflasNoTrans, FflasNoTrans, m,n,R, F.one, L,R, U,n, F.zero, X,n);

	// write_perm(std::cerr<<"P = ",P,m);
	// write_perm(std::cerr<<"Q = ",Q,n);
	// write_field(F,std::cerr<<"L = "<<std::endl,L,m,R,R);
	// write_field(F,std::cerr<<"U = "<<std::endl,U,R,n,n);
	

	bool fail = false;
	for(size_t i=0; i<m; ++i)
		for (size_t j=0; j<n; ++j)
			if (!F.areEqual (*(A+i*lda+j), *(X+i*n+j))){
				std::cerr << std::endl<<" A ["<<i<<","<<j<<"] = " << (*(A+i*lda+j))
						  << " PLUQ ["<<i<<","<<j<<"] = " << (*(X+i*n+j));
				fail=true;
			}
		//write_field(F, std::cerr<<"X = "<<std::endl,X, m,n,n);
	if (fail)
		std::cerr << std::endl;

	fflas_delete( U);
	fflas_delete( L);
	fflas_delete( X);
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
template<class Field, FFLAS_DIAG diag, class RandIter>
bool test_pluq (const Field & F,
				typename Field::ConstElement_ptr A,
				size_t r, size_t m, size_t n, size_t lda, RandIter& G)
{
	bool fail = false;
	typedef typename Field::Element_ptr Element_ptr ;
	Element_ptr B = fflas_new(F,m,lda) ;
	fassign(F,m,n,A,lda,B,lda);

	size_t * P = fflas_new<size_t> (m);
	size_t * Q = fflas_new<size_t> (n);
	
    ForceCheck_PLUQ<Field> checker (G,m,n,B,lda);

	size_t R = PLUQ (F, diag, m, n, B, lda, P, Q);

    try {
        checker.check(B,lda,diag,R,P,Q);
    } catch(FailurePLUQCheck &e) {
        std::cout << m << 'x' << n << " pluq verification failed!\n";
    }

	if (R != r) {
		std::cout << "rank is wrong (expected " << r << " but got " << R << ")" << std::endl;
		fflas_delete (B);
		fflas_delete (P);
		fflas_delete (Q);
		return fail = true;
	}
	fail |=  verifPLUQ<Field,diag> (F,A, lda, B, lda, P, Q, m, n, r);
	fflas_delete (B);
	fflas_delete(P);
	fflas_delete(Q);
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
// template<class Field, FFLAS_DIAG diag, FFLAS_TRANSPOSE trans>
// bool test_lu_append(const Field & F,
// 		    const typename Field::Element_ptr A,
// 		    const typename Field::Element_ptr B,
// 		    size_t m, size_t n, size_t k, size_t lda)
// {
// 	FFLASFFPACK_check(n<=lda);

// 	bool fail = false;
// 	size_t M = m + k ;
// 	typedef typename Field::Element Element ;
// 	Element_ptr Acop = fflas_new(F, m, lda) ;
// 	fassign(F,m,n,A,lda,Acop,lda) ;

// 	Element_ptr Bcop = fflas_new(F, k, lda) ;
// 	fassign(F,k,n,B,lda,Bcop,lda) ;

// 	Element_ptr Append = fflas_new (F, M, lda);
// 	fassign(F,m,n,A,lda,Append,lda) ;
// 	fassign(F,k,n,B,lda,Append+m*lda,lda) ;

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

// 	Element_ptr Afull = fflas_new(F, M, lda);
// 	fassign(F,M,n,Append,lda,Afull,lda) ;
// 	// fassign(F,m,n,A,lda,Afull,lda) ;
// 	// fassign(F,k,n,B,lda,Afull+m*lda,lda) ;

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

// 	if (trans == FflasTrans){
// 		maxP = M;
// 		maxQ = n;
// 	}
// 	else{ // trans == FflasNoTrans
// 		maxP = n;
// 		maxQ = M;
// 	}

// 	size_t * P = fflas_new<size_t>(maxP) ;
// 	size_t * Q = fflas_new<size_t>(maxQ) ;

// 	size_t * PP = fflas_new<size_t>(maxP) ;
// 	size_t * QQ = fflas_new<size_t>(maxQ) ;

// 	/* valgrind says the following leaks. Just incroyable. */
// 	size_t R  = LUdivine (F, diag, trans, M, n, Append, lda, PP, QQ);

// 	size_t R1 = LUdivine (F, diag, trans, m, n, Acop,   lda, P, Q);

// 	size_t R2 = LUpdate  (F,diag,trans,m,n,Acop,lda,R1,k,Bcop,lda,P,Q,
// 				      FfpackLQUP);
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
// 		fflas_delete( Bcop );
// 		fflas_delete( Acop );
// 		fflas_delete( Append );
// 		fflas_delete( PP);
// 		fflas_delete( QQ);
// 		fflas_delete( P );
// 		fflas_delete( Q );
// 		return fail=true;

// 	}

// 	// compute C=LQUP and check C == A
// 	Element_ptr C = fflas_new (F, M, lda);
// 	/*  Build L,U */
// 	Element_ptr L, U;
// 	if (trans == FflasNoTrans){
// 		L = fflas_new(F, M, M);
// 		U = fflas_new(F, M, n);

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
// 		applyP( F, FflasRight, FflasNoTrans,
// 				M,0,(int)R, L, M, Q);
// 		for ( size_t i=0; i<M; ++i )
// 			F.assign(*(L+i*(M+1)), one);

// 		/*  reconstruct the diagonal */
// 		//write_field(F,cerr<<"L = "<<endl,L,m,m,m);
// 		//write_field(F,cerr<<"U = "<<endl,U,m,n,n);
// 		if (diag == FflasNonUnit) {
// 			for ( size_t i=0; i<R; ++i ) {
// 				if (i<m)
// 					F.assign (*(U+i*(n+1)), *(Acop+i*(lda+1)));
// 				else
// 					F.assign (*(U+i*(n+1)), *(Bcop+(i-m)*(lda+1)));
// 			}
// 		}
// 		else{ // diag == FflasUnit
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
// 		applyP (F, FflasRight, FflasNoTrans,
// 				M,0,(int) R, U, n, P);
// 		applyP (F, FflasLeft, FflasTrans,
// 				n,0,(int)R, U, n, Q);
// 		fgemm (F, FflasNoTrans, FflasNoTrans,
// 			      M,n,M, 1.0, L,M, U,n, 0.0, C,lda);
// 		//fflas_delete( A);
// 	}
// #if 0 /*  not working */
// 	else { /*  trans == FflasTrans */

// 		L = fflas_new(F, M, n);
// 		U = fflas_new(F, n, n);


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

// 		applyP( F, FflasLeft, FflasTrans,
// 				n,0,(int)R, U, n, Q);

// 		for (size_t i=0; i<n; ++i)
// 			F.assign (*(U+i*(n+1)),one);

// 		/*  reconstruct the diagonal */
// 		if (diag == FflasNonUnit) {
// 			for ( size_t i=0; i<R; ++i ) {
// 				if (i < m)
// 					F.assign (*(L+i*(n+1)), *(Acop+i*(lda+1)));
// 				else
// 					F.assign (*(L+i*(n+1)), *(Bcop+(i-m)*(lda+1)));
// 			}
// 		}
// 		else{ // diag == FflasUnit
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
// 		applyP (F, FflasLeft, FflasTrans,
// 				n,0,(int)R, L, n, P);
// 		applyP (F, FflasRight, FflasNoTrans,
// 				M,0,(int)R, L, n, Q);
// 		fgemm (F, FflasNoTrans, FflasNoTrans,
// 			      M,n,n, 1.0, L,n, U,n, 0.0, C,lda);
// 	}
// #endif
// #if 0 /*  check CC == LL UU */
// 	Element_ptr LL, UU;
// 	Element_ptr CC = fflas_new (F, M, lda);
// 	if (trans == FflasNoTrans){
// 		LL = fflas_new (F, M, M);
// 		UU = fflas_new (F, M, n);

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
// 		applyP( F, FflasRight, FflasNoTrans,
// 				M,0,(int)R, LL, M, Q);
// 		for ( size_t i=0; i<M; ++i )
// 			F.assign(*(LL+i*(M+1)), one);

// 		/*  reconstruct the diagonal */
// 		//write_field(F,cerr<<"LL = "<<endl,LL,m,m,m);
// 		//write_field(F,cerr<<"UU = "<<endl,UU,m,n,n);
// 		if (diag == FflasNonUnit) {
// 			for ( size_t i=0; i<R; ++i ) {
// 				F.assign (*(UU+i*(n+1)), *(Append+i*(lda+1)));
// 			}
// 		}
// 		else{ // diag == FflasUnit
// 			for ( size_t i=0; i<R; ++i ){
// 				*(LL+Q[i]*(M+1)) = *(Append+Q[i]*lda+i);
// 				F.assign (*(UU+i*(n+1)),one);
// 			}
// 		}
// 		write_field(F,cerr<<"L = "<<endl,LL,(int)M,(int)M,(int)M);
// 		write_field(F,cerr<<"U = "<<endl,UU,(int)M,(int)n,(int)n);

// 		/*  Compute LQUP */
// 		applyP (F, FflasRight, FflasNoTrans,
// 				M,0,(int) R, UU, n, P);
// 		applyP (F, FflasLeft, FflasTrans,
// 				n,0,(int)R, UU, n, Q);
// 		fgemm (F, FflasNoTrans, FflasNoTrans,
// 			      M,n,M, 1.0, LL,M, UU,n, 0.0, CC,lda);
// 		//fflas_delete( A);
// 	}
// 	else { /*  trans == FflasTrans */

// 		LL = fflas_new(F, M, n);
// 		UU = fflas_new(F, n, n);


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

// 		applyP( F, FflasLeft, FflasTrans,
// 				n,0,(int)R, UU, n, Q);

// 		for (size_t i=0; i<n; ++i)
// 			F.assign (*(UU+i*(n+1)),one);

// 		/*  reconstruct the diagonal */
// 		if (diag == FflasNonUnit) {
// 			for ( size_t i=0; i<R; ++i ) {
// 				F.assign (*(LL+i*(n+1)), *(Append+i*(lda+1)));
// 			}
// 		}
// 		else{ // diag == FflasUnit
// 			for ( size_t i=0; i<R; ++i ){
// 				*(UU+Q[i]*(n+1)) = *(Append+Q[i]+i*lda);
// 				F.assign (*(LL+i*(n+1)),one);
// 			}
// 		}
// 		write_field(F,cerr<<"LL = "<<endl,LL,(int)M,(int)n,(int)n);
// 		write_field(F,cerr<<"UU = "<<endl,UU,(int)n,(int)n,(int)n);

// 		/*  Compute LQUP */
// 		applyP (F, FflasLeft, FflasTrans,
// 				n,0,(int)R, LL, n, P);
// 		applyP (F, FflasRight, FflasNoTrans,
// 				M,0,(int)R, LL, n, Q);
// 		fgemm (F, FflasNoTrans, FflasNoTrans,
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

// 	fflas_delete( PP);
// 	fflas_delete( P);
// 	fflas_delete( L);
// 	fflas_delete( U);
// 	fflas_delete( Q);
// 	fflas_delete( QQ);
// 	fflas_delete( Acop);
// 	fflas_delete( Bcop);
// 	fflas_delete( Append);
// 	fflas_delete( Afull);
// 	fflas_delete( C);

// 	return fail;


// }

template<class Field, FFLAS_DIAG diag, FFLAS_TRANSPOSE trans, class RandIter>
bool launch_test(const Field & F,
				 size_t r,
				 size_t m, size_t n, RandIter& G)
{
		//typedef typename Field::Element Element ;
	typedef typename Field::Element_ptr Element_ptr ;
	bool fail = false ;
	{ /*  user given and lda bigger */
		size_t lda = n+10 ;
		Element_ptr A = fflas_new (F, m, lda);
		RandomMatrixWithRankandRandomRPM(F,m,n,r,A,lda,G);
		fail |= test_LUdivine<Field,diag,trans>(F,A,lda,r,m,n);
		RandomMatrixWithRankandRandomRPM(F,m,n,r,A,lda,G);
		fail |= test_pluq<Field,diag>(F,A,r,m,n,lda,G);
		if (fail) std::cout << "failed at big lda" << std::endl;
		fflas_delete( A );
	}
	{ /*  user given and lda bigger. Rank is max */
		size_t lda = n+10 ;
		size_t R = std::min(m,n);
		Element_ptr A = fflas_new (F, m, lda);
		RandomMatrixWithRankandRandomRPM(F,m,n,R,A,lda,G);
		fail |= test_LUdivine<Field,diag,trans>(F,A,lda,R,m,n);
		RandomMatrixWithRankandRandomRPM(F,m,n,R,A,lda,G);
		fail |= test_pluq<Field,diag>(F,A,R,m,n,lda,G);
		if (fail) std::cout << "failed at big lda max rank" << std::endl;
		fflas_delete( A );
	}
	{ /*  user given and lda bigger. Rank is min */
		size_t lda = n+10 ;
		size_t R = 0;
		Element_ptr A = fflas_new (F, m, lda);
		RandomMatrixWithRankandRandomRPM(F,m,n,R,A,lda,G);
		fail |= test_LUdivine<Field,diag,trans>(F,A,lda,R,m,n);
		RandomMatrixWithRankandRandomRPM(F,m,n,R,A,lda,G);
		fail |= test_pluq<Field,diag>(F,A,R,m,n,lda,G);
		if (fail) std::cout << "failed at big lda, rank 0" << std::endl;
		fflas_delete( A );
	}
	{ /*  square  */
		size_t M = std::max(m,n);
		size_t N = M ;
		size_t R = M/2 ;
		size_t lda = N+10 ;
		Element_ptr A = fflas_new (F, M, lda);
		RandomMatrixWithRankandRandomRPM(F,M,N,R,A,lda,G);
		fail |= test_LUdivine<Field,diag,trans>(F,A,lda,R,M,N);
		RandomMatrixWithRankandRandomRPM(F,M,N,R,A,lda,G);
		fail |= test_pluq<Field,diag>(F,A,R,M,N,lda,G);
		if (fail) std::cout << "failed at square" << std::endl;
		fflas_delete( A );
	}
	{ /*  wide  */
		size_t M = std::max(m,n);
		size_t N = 2*M ;
		size_t R = 3*M/4 ;
		size_t lda = N+5 ;
		Element_ptr A = fflas_new (F, M, lda);
		RandomMatrixWithRankandRandomRPM(F,M,N,R,A,lda,G);
		fail |= test_LUdivine<Field,diag,trans>(F,A,lda,R,M,N);
		RandomMatrixWithRankandRandomRPM(F,M,N,R,A,lda,G);
		fail |= test_pluq<Field,diag>(F,A,R,M,N,lda,G);
		if (fail) std::cout << "failed at wide" << std::endl;
		fflas_delete( A );
	}
	{ /*  narrow  */
		size_t M = std::max(m,n);
		size_t N = M/2 ;
		size_t R = 3*M/8 ;
		size_t lda = N+5 ;
		Element_ptr A = fflas_new (F, M, lda);
		RandomMatrixWithRankandRandomRPM(F,M,N,R,A,lda,G);
		fail |= test_LUdivine<Field,diag,trans>(F,A,lda,R,M,N);
		RandomMatrixWithRankandRandomRPM(F,M,N,R,A,lda,G);
		fail |= test_pluq<Field,diag>(F,A,R,M,N,lda,G);
		if (fail) std::cout << "failed at narrow" << std::endl;
		fflas_delete( A );
	}
 	return !fail;
 }

// template<class Field, FFLAS_DIAG diag, FFLAS_TRANSPOSE trans>
// bool launch_test_append(const Field & F,
// 			size_t r,
// 			size_t m, size_t n)
// {
// 	typedef typename Field::Element Element ;
// 	bool fail = false ;
// 	{ /*  user given and lda bigger */
// 		size_t lda = n+10 ;
// 		size_t k = m/2+1 ;
// 		Element_ptr A = fflas_new (F, m, lda);
// 		Element_ptr B = fflas_new (F, k, lda);
// 		RandomMatrixWithRank(F,A,lda,r,m,n,G);
// 		RandomMatrixWithRank(F,B,lda,k/2+1,k,n,G);
// 		fail |= test_lu_append<Field,diag,trans>(F,A,B,m,n,k,lda);
// 		if (fail) std::cout << "failed" << std::endl;
// 		fflas_delete( A );
// 		fflas_delete( B );
// 	}
// 	{ /*  user given and lda bigger. Rank is max */
// 		size_t lda = n+10 ;
// 		size_t R = std::min(m,n);
// 		size_t k = m/2+1 ;
// 		Element_ptr A = fflas_new (F, m, lda);
// 		Element_ptr B = fflas_new (F, k, lda);
// 		RandomMatrixWithRank(F,m,n,R,A,lda,G);
// 		RandomMatrixWithRank(F,B,lda,k/2+1,k,n,G);
// 		fail |= test_lu_append<Field,diag,trans>(F,A,B,m,n,k,lda);
// 		if (fail) std::cout << "failed" << std::endl;
// 		fflas_delete( A );
// 		fflas_delete( B );
// 	}
// 	{ /*  user given and lda bigger. Appended Rank is min */
// 		size_t lda = n+10 ;
// 		size_t R = std::min(m,n);
// 		size_t k = m/2+1 ;
// 		Element_ptr A = fflas_new (F, m, lda);
// 		Element_ptr B = fflas_new (F, k, lda);
// 		RandomMatrixWithRank(F,m,n,R,A,lda,G);
// 		RandomMatrixWithRank(F,B,lda,0,k,n,G);
// 		fail |= test_lu_append<Field,diag,trans>(F,A,B,m,n,k,lda);
// 		if (fail) std::cout << "failed" << std::endl;
// 		fflas_delete( A );
// 		fflas_delete( B );
// 	}
// 	{ /*  user given and lda bigger. Rank is min */
// 		size_t lda = n+10 ;
// 		size_t R = 0;
// 		size_t k = m/2+1 ;
// 		Element_ptr A = fflas_new (F, m, lda);
// 		Element_ptr B = fflas_new (F, k, lda);
// 		RandomMatrixWithRank(F,m,n,R,A,lda,G);
// 		RandomMatrixWithRank(F,B,lda,k/2+1,k,n,G);
// 		fail |= test_lu_append<Field,diag,trans>(F,A,B,m,n,k,lda);
// 		if (fail) std::cout << "failed" << std::endl;
// 		fflas_delete( A );
// 		fflas_delete( B );
// 	}
// 	{ /*  square  */
// 		size_t M = std::max(m,n);
// 		size_t N = M ;
// 		size_t R = M/2 ;
// 		size_t lda = N+10 ;
// 		size_t k = R ;
// 		Element_ptr A = fflas_new (F, M, lda);
// 		Element_ptr B = fflas_new (F, k, lda);
// 		RandomMatrixWithRank(F,A,lda,R,M,N,G);
// 		RandomMatrixWithRank(F,B,lda,R/2,k,N,G);
// 		fail |= test_lu_append<Field,diag,trans>(F,A,B,M,N,k,lda);
// 		if (fail) std::cout << "failed" << std::endl;
// 		fflas_delete( A );
// 		fflas_delete( B );
// 	}
// 	{ /*  wide  */
// 		size_t M = std::max(m,n);
// 		size_t N = 2*M ;
// 		size_t R = M/2 ;
// 		size_t k = R ;
// 		size_t lda = N+10 ;
// 		Element_ptr A = fflas_new (F, M, lda);
// 		Element_ptr B = fflas_new (F, k, lda);
// 		RandomMatrixWithRank(F,A,lda,R,M,N,G);
// 		RandomMatrixWithRank(F,B,lda,k/2,k,N,G);
// 		fail |= test_lu_append<Field,diag,trans>(F,A,B,M,N,k,lda);
// 		if (fail) std::cout << "failed" << std::endl;
// 		fflas_delete( A );
// 		fflas_delete( B );
// 	}
// 	//! @bug leaks :
// #if 0 /*  leak here */
// 	{ /*  narrow  */
// 		size_t M = std::max(m,n);
// 		size_t N = M/2 ;
// 		size_t R = M/3 ;
// 		size_t k = N ;
// 		size_t lda = N+10 ;
// 		Element_ptr A = fflas_new (F, M, lda);
// 		Element_ptr B = fflas_new (F, k, lda);
// 		RandomMatrixWithRank(F,A,lda,R,M,N,G);
// 		RandomMatrixWithRank(F,A,lda,std::min(k/2,M/2),k,N,G);
// 		fail |= test_lu_append<Field,diag,trans>(F,A,B,M,N,k,lda);
// 		if (fail) std::cout << "failed" << std::endl;
// 		fflas_delete( A );
// 		fflas_delete( B );
// 	}
// #endif

// 	return fail;
// }


template<class Field>
bool run_with_field(Givaro::Integer q, uint64_t b, size_t m, size_t n, size_t r, size_t iters, uint64_t seed){
	bool ok = true ;
	int nbit=(int)iters;
	
	while (ok &&  nbit){
			// choose Field 
		Field* F= chooseField<Field>(q,b);
		if (F==nullptr)
			return true;
		typename Field::RandIter G(*F,0,seed);
		std::ostringstream oss;
		F->write(oss);
		
		std::cout.fill('.');
		std::cout<<"Checking ";
		std::cout.width(40);
		std::cout<<oss.str();
		std::cout<<" ... ";


		ok&= launch_test<Field,FflasUnit,FflasNoTrans>    (*F,r,m,n,G);
		ok&= launch_test<Field,FflasUnit,FflasTrans>      (*F,r,m,n,G);
		ok&= launch_test<Field,FflasNonUnit,FflasNoTrans> (*F,r,m,n,G);
		ok&= launch_test<Field,FflasNonUnit,FflasTrans>   (*F,r,m,n,G);

#if 0 /*  may be bogus */
		ok&= launch_test_append<Field,FflasUnit,FflasNoTrans>   (*F,r,m,n,G);
		ok&= launch_test_append<Field,FflasNonUnit,FflasNoTrans>(*F,r,m,n,G);
		ok&= launch_test_append<Field,FflasUnit,FflasTrans>     (*F,r,m,n,G);
		ok&= launch_test_append<Field,FflasNonUnit,FflasTrans>  (*F,r,m,n,G);
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
	Givaro::Integer q=-1;
	size_t b=0;
	size_t m=120;
	size_t n=120;
	size_t r=70;
	size_t iters=3;
	bool loop=false;
	size_t seed=time(NULL);
	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
		{ 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
		{ 'm', "-m M", "Set the row dimension of the matrix.",      TYPE_INT , &m },
		{ 'n', "-n N", "Set the column dimension of the matrix.", TYPE_INT , &n },
		{ 'r', "-r R", "Set the rank.", TYPE_INT , &r },
		{ 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
		{ 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
		{ 's', "-s seed", "Set seed for the random generator", TYPE_INT, &seed },
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,as);

	if (r > std::min (m,n)) 
		r = std::min (m, n);

	srand(seed);

	bool ok=true;
	do{
		ok&=run_with_field<Givaro::Modular<float> >           (q,b,m,n,r,iters,seed);
		ok&=run_with_field<Givaro::Modular<double> >          (q,b,m,n,r,iters,seed);
		ok&=run_with_field<Givaro::ModularBalanced<float> >   (q,b,m,n,r,iters,seed);
		ok&=run_with_field<Givaro::ModularBalanced<double> >  (q,b,m,n,r,iters,seed);
		ok&=run_with_field<Givaro::Modular<int32_t> >         (q,b,m,n,r,iters,seed);
		ok&=run_with_field<Givaro::ModularBalanced<int32_t> > (q,b,m,n,r,iters,seed);
		ok&=run_with_field<Givaro::Modular<int64_t> >         (q,b,m,n,r,iters,seed);
		ok&=run_with_field<Givaro::ModularBalanced<int64_t> > (q,b,m,n,r,iters,seed);
		ok&=run_with_field<Givaro::Modular<Givaro::Integer> > (q,5,m/6,n/6,r/6,iters,seed);
		ok&=run_with_field<Givaro::Modular<Givaro::Integer> > (q,(b?b:512),m/6,n/6,r/6,iters,seed);
	} while (loop && ok);

	return !ok;
}
