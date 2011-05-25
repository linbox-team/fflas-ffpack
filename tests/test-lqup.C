/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) Fflas-Ffpack
 * Written by Clément Pernet
 * This file is Free Software and part of Fflas-Ffpack.
 * See COPYING for license information.
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
#include "timer.h"
#include "fflas-ffpack/field/modular-balanced.h"
// #include "fflas-ffpack/field/nonzero-randiter.h"
#include "fflas-ffpack/ffpack/ffpack.h"

#include "fflas-ffpack/utils/args-parser.h"

using namespace FFPACK;

/*! @brief  Random Matrix.
 * Creates a \c m x \c n matrix with random entries.
 * @param F field
 * @param A pointer to the matrix (preallocated to at least \c m x \c lda field elements)
 * @param m number of rows in \p A
 * @param n number of cols in \p A
 * @param lda leading dimension of \p A
 * @return pointer to \c A.
 */
template<class Field>
typename Field::Element * RandomMatrix(const Field & F,
				      typename Field::Element * A,
				      size_t m, size_t n, size_t lda)
{
	typedef typename Field::RandIter Randiter ;
	Randiter R(F);
	for (size_t i=0 ; i<m ; ++i)
		for (size_t j= 0; j<n ;++j)
			R.random( A[i*lda+j] );
	return A;

}

/*! Random integer in range.
 * @param a min bound
 * @param b max bound
 * @return a random integer in [a,b[  */
int RandInt(int a, int b)
{
	int x = a ;
	x += rand()%(b-a+1);
	fflaflas_check(x<b && x>=a);
	return x ;
}
/*! Random integer in range.
 * @param a min bound
 * @param b max bound
 * @return a random integer in [a,b[  */
size_t RandInt(size_t a, size_t b)
{
	size_t x = a ;
	x += rand()%(b-a);
	fflaflas_check(x<b && x>=a);
	return x ;
}

/*! @brief  Random Matrix with prescribed rank.
 * Creates a \c m x \c n matrix with random entries and rank \c r.
 * @param F field
 * @param A pointer to the matrix (preallocated to at least \c m x \c lda field elements)
 * @param r rank of the matrix to build
 * @param m number of rows in \p A
 * @param n number of cols in \p A
 * @param lda leading dimension of \p A
 * @return pointer to \c A.
 */
template<class Field>
typename Field::Element * RandomMatrixWithRank(const Field & F,
				      typename Field::Element * A,
				      size_t r,
				      size_t m, size_t n, size_t lda)
{
	typedef typename Field::RandIter Randiter ;
	typedef typename Field::Element  Element ;
	Randiter R(F);
	NonzeroRandIter<Field,Randiter> nzR(F,R);

	size_t * P = new size_t[n];
	size_t * Q = new size_t[m];
	for (size_t i = 0 ; i < m ; ++i ) Q[i] = 0;
	for (size_t i = 0 ; i < n ; ++i ) P[i] = 0;

	Element * U = new Element[m*n];
	Element * L = new Element[m*m];


	/*  Create L, lower invertible */
	for (size_t i=0 ; i<m ; ++i)
		for (size_t j= 0; j<i ;++j)
			R.random( L[i*m+j] );

	for (size_t i=0 ; i<m ; ++i)
		nzR.random( L[i*m+i] );

	for (size_t i=0 ; i<m ; ++i)
		for (size_t j= i+1; j<m ;++j)
			 F.init(L[i*m+j],0);


	/*  Create U, upper or rank r */
	for (size_t i=0 ; i<r ; ++i)
		for (size_t j= i+1; j<r ;++j)
			R.random( U[i*n+j] );
	for (size_t i=0 ; i<r ; ++i)
			nzR.random( U[i*n+i] );
	for (size_t i=0 ; i<r ; ++i)
		for (size_t j= 0 ; j<i ;++j)
			 F.init(U[i*n+j],0);

	for (size_t i=r ; i<m ; ++i)
		for (size_t j= 0 ; j<n ;++j)
			 F.init(U[i*n+j],0);

	for (size_t i=0 ; i<r ; ++i)
		for (size_t j= r ; j<n ;++j)
			R.random( U[i*n+j] );

	/*  Create a random P,Q */

	for (size_t i = 0 ; i < n ; ++i)
		P[i] = i + RandInt(0UL,n-i);
	for (size_t i = 0 ; i < m ; ++i)
		Q[i] = i + RandInt(0UL,m-i);

	/*  compute product */

	FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
			m,0,(int)n, U, n, P);
	FFPACK::applyP (F, FFLAS::FflasLeft,  FFLAS::FflasNoTrans,
		       	m,0,(int)m, L, m, Q);
	FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
		      m,n,m, 1.0, L,m, U,n, 0.0, A,lda);
	//! @todo compute LU with ftrtr

	delete[] P;
	delete[] L;
	delete[] U;
	delete[] Q;

	return A;

}

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
 * @return true iff correct
 */
template<class Field, FFLAS::FFLAS_DIAG diag, FFLAS::FFLAS_TRANSPOSE trans>
bool test_lu(const Field & F,
	     const typename Field::Element * A,
	     size_t r,
	     size_t m, size_t n, size_t lda)
{
	typedef typename Field::Element Element ;
	Element * B = new Element[m*lda] ;
	memcpy(B,A,m*lda*sizeof(Element)); // probably faster than ::fcopy !

	size_t maxP, maxQ ;

	if (trans == FFLAS::FflasTrans){
		maxP = m;
		maxQ = n;
	} else{
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
		return true;
	}

	Element * C = new Element[m*n]; // compute C=LQUP and check C == A
	Element * L, *U;
	if (trans == FFLAS::FflasNoTrans){
		L = new Element[m*m];
		U = new Element[m*n];

		Element zero,one;
		F.init(zero,0.0);
		F.init(one,1.0);
		for (size_t i=0; i<R; ++i){
			for (size_t j=0; j<i; ++j)
				F.assign ( *(U + i*n + j), zero);
			for (size_t j=i+1; j<n; ++j)
				F.assign (*(U + i*n + j), *(B+ i*lda+j));
		}
		for (size_t i=R;i<m; ++i)
			for (size_t j=0; j<n; ++j)
				F.assign(*(U+i*n+j), zero);
		for ( size_t i=0; i<m; ++i ){
			size_t j=0;
			for (; j< ((i<R)?i:R) ; ++j )
				F.assign( *(L + i*m+j), *(B+i*lda+j));
			for (; j<m; ++j )
				F.assign( *(L+i*m+j), zero);
		}

		// write_field(F,cerr<<"L = "<<endl,L,m,m,m);
// 		write_field(F,cerr<<"U = "<<endl,U,m,n,n);
		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
			       	m,0,(int)R, L, m, Q);
		for ( size_t i=0; i<m; ++i )
			F.assign(*(L+i*(m+1)), one);

 		//write_field(F,cerr<<"L = "<<endl,L,m,m,m);
 		//write_field(F,cerr<<"U = "<<endl,U,m,n,n);
		if (diag == FFLAS::FflasNonUnit)
			for ( size_t i=0; i<R; ++i )
				F.assign (*(U+i*(n+1)), *(B+i*(lda+1)));

		else{
			for ( size_t i=0; i<R; ++i ){
				*(L+Q[i]*(m+1)) = *(B+Q[i]*lda+i);
				F.assign (*(U+i*(n+1)),one);
			}
		}
// 		write_field(F,cerr<<"L = "<<endl,L,m,m,m);
// 		write_field(F,cerr<<"U = "<<endl,U,m,n,n);

		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
			       	m,0,(int) R, U, n, P);
		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
			       	n,0,(int)R, U, n, Q);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
			      m,n,m, 1.0, L,m, U,n, 0.0, C,n);
		//delete[] A;
	} else {

		L = new Element[m*n];
		U = new Element[n*n];


		Element zero,one;
		F.init(zero,0.0);
		F.init(one,1.0);
		for (size_t i=0; i<R; ++i){
			for (size_t j=0; j<i; ++j)
				F.assign ( *(L + i + j*n), zero);
			for (size_t j=i+1; j<m; ++j)
				F.assign (*(L + i + j*n), *(B+ i+j*lda));
		}

		for (size_t i=R;i<n; ++i)
			for (size_t j=0; j<m; ++j)
				F.assign(*(L+i+j*n), zero);
		for ( size_t i=0; i<n; ++i ){
			size_t j=0;
			for (;  j< ((i<R)?i:R) ; ++j )
				F.assign( *(U + i+j*n), *(B+i+j*lda));
			for (; j<n; ++j )
				F.assign( *(U+i+j*n), zero);
		}
// 		write_field(F,cerr<<"L = "<<endl,L,m,n,n);
// 		write_field(F,cerr<<"U = "<<endl,U,n,n,n);

		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans,
			       	n,0,(int)R, U, n, Q);


		for (size_t i=0; i<n; ++i)
			F.assign (*(U+i*(n+1)),one);

		if (diag == FFLAS::FflasNonUnit)
			for ( size_t i=0; i<R; ++i )
				F.assign (*(L+i*(n+1)), *(B+i*(lda+1)));
		else{
			for ( size_t i=0; i<R; ++i ){
				*(U+Q[i]*(n+1)) = *(B+Q[i]+i*lda);
				F.assign (*(L+i*(n+1)),one);
			}
		}
// 		write_field(F,cerr<<"L = "<<endl,L,m,n,n);
// 		write_field(F,cerr<<"U = "<<endl,U,n,n,n);

		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
			       	n,0,(int)R, L, n, P);
		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
			       	m,0,(int)R, L, n, Q);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
			      m,n,n, 1.0, L,n, U,n, 0.0, C,n);
	}
	bool fail = false;
	for (size_t i=0; i<m; ++i)
		for (size_t j=0; j<n; ++j)
			if (!F.areEqual (*(A+i*lda+j), *(C+i*n+j))){
				std::cerr << " A["<<i<<","<<j<<"]    = " << (*(A+i*lda+j))
					  << " PLUQ["<<i<<","<<j<<"] = " << (*(C+i*n+j))
					  << endl;
				fail=true;
			}

	delete[] P;
	delete[] L;
	delete[] U;
	delete[] Q;
	delete[] B;
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
		fail &= test_lu<Field,diag,trans>(F,A,r,m,n,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
	}
	{ /*  user given and lda bigger. Rank is max */
		size_t lda = n+10 ;
		size_t R = std::min(m,n);
		Element * A = new Element[m*lda];
		RandomMatrixWithRank(F,A,R,m,n,lda);
		fail &= test_lu<Field,diag,trans>(F,A,R,m,n,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
	}
	{ /*  user given and lda bigger. Rank is min */
		size_t lda = n+10 ;
		size_t R = 0;
		Element * A = new Element[m*lda];
		RandomMatrixWithRank(F,A,R,m,n,lda);
		fail &= test_lu<Field,diag,trans>(F,A,R,m,n,lda);
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
		fail &= test_lu<Field,diag,trans>(F,A,R,M,N,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
	}
	{ /*  wide  */
		size_t M = std::max(m,n);
		size_t N = 2*M ;
		size_t R = M/2 ;
		size_t lda = N+10 ;
		Element * A = new Element[M*lda];
		RandomMatrixWithRank(F,A,R,M,N,lda);
		fail &= test_lu<Field,diag,trans>(F,A,R,M,N,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
	}
	{ /*  narrow  */
		size_t M = std::max(m,n);
		size_t N = M/2 ;
		size_t R = M/3 ;
		size_t lda = N+10 ;
		Element * A = new Element[M*lda];
		RandomMatrixWithRank(F,A,R,M,N,lda);
		fail &= test_lu<Field,diag,trans>(F,A,R,M,N,lda);
		if (fail) std::cout << "failed" << std::endl;
		delete[] A ;
	}

	return fail;
}

int main(int argc, char** argv){
	cerr<<setprecision(20);
	int p = 101;
	int m = 50;
	int n = 50;
	int r = 20;
	int iter = 5 ;
	bool fail = false;

	static Argument as[] = {
		{ 'p', "-p P", "Set the field characteristic.",         TYPE_INT , &p },
		{ 'n', "-n N", "Set the number of cols in the matrix.", TYPE_INT , &n },
		{ 'm', "-m N", "Set the number of rows in the matrix.", TYPE_INT , &m },
		{ 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iter },
		// { 'f', "-f file", "Set input file", TYPE_STR, &file },
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,as);

	typedef ModularBalanced<double> Field;

	typedef Field::Element Element;

	Field F(p);

	for (int i = 0 ; i < iter ; ++i) {
		fail &= launch_test<Field,FFLAS::FflasUnit,FFLAS::FflasNoTrans>(F,r,m,n);
		fail &= launch_test<Field,FFLAS::FflasUnit,FFLAS::FflasTrans>(F,r,m,n);
		fail &= launch_test<Field,FFLAS::FflasNonUnit,FFLAS::FflasNoTrans>(F,r,m,n);
		fail &= launch_test<Field,FFLAS::FflasNonUnit,FFLAS::FflasTrans>(F,r,m,n);
	}

	return fail ;

}
