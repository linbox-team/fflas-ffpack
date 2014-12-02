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



//#define __LUDIVINE_CUTOFF 1
#include <iostream>
#include <iomanip>
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/field/modular-integer.h" 
#include "fflas-ffpack/ffpack/ffpack.h" 
#include "test-utils.h"

#include "fflas-ffpack/utils/args-parser.h"

using namespace std;
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
	     const typename Field::Element_ptr A,
	     size_t r,
	     size_t m, size_t n, size_t lda)
{

	bool fail = false;
	typedef typename Field::Element Element ;
	typedef typename Field::Element_ptr Element_ptr ;
	Element_ptr  B = FFLAS::fflas_new(F,m,lda) ;
	// memcpy(B,A,m*lda*sizeof(Element)); // probably faster than ::fassign !
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


	FFLAS::Timer chrono;
	chrono.start();
	size_t R = FFPACK::LUdivine (F, diag, trans, m, n, B, lda, P, Q,
				     FFPACK::FfpackLQUP);
	chrono.stop();
	std::cout<<"Ludivine done"<<std::endl;

	if (R != r) { 
		std::cout << "rank is wrong (expected " << r << " but got " << R << ")" << std::endl;
		FFLAS::fflas_delete(B) ;
		FFLAS::fflas_delete(P) ;
		FFLAS::fflas_delete(Q) ;
		return fail = true; 
	}

	Element_ptr C = FFLAS::fflas_new(F,m,n); // compute C=LQUP and check C == A
	/*  Build L,U */
	Element_ptr L, U;

	if (trans == FFLAS::FflasNoTrans){
		L = FFLAS::fflas_new(F,m,m);
		U = FFLAS::fflas_new(F,m,n);

		Element zero,one;
		F.init(zero,0.0);
		F.init(one,1.0);
		/*  build U */
		for (size_t i=0; i<R; ++i){
			for (size_t j=0; j<i; ++j)
				F.assign (*(U + i*n + j), zero);
			for (size_t j=i+1; j<n; ++j)
				F.assign (*(U + i*n + j), *(B+i*lda+j));
		}
		for (size_t i=R;i<m; ++i) {
			for (size_t j=0; j<n; ++j)
				F.assign(*(U+i*n+j), zero);
		}
		/*  build L */
		for ( size_t i=0; i<m; ++i ){
			size_t j=0;
			for (; j< ((i<R)?i:R) ; ++j )
				F.assign( *(L+i*m+j), *(B+i*lda+j));
			for (; j<m; ++j )
				F.assign( *(L+i*m+j), zero);
		}

		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				m,0,(int)R, L, m, Q);

		for ( size_t i=0; i<m; ++i )
			F.assign(*(L+i*(m+1)), one);

		/*  reconstruct the diagonal */
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
		
		/*  Compute LQUP */
		FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
				m,0,(int) R, U, n, P);
		FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
				n,0,(int)R, U, n, Q);
		FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
			      m,n,m, F.one, L,m, U,n, F.zero, C,n);

		// write_field(F,cout<<"L := ",L,m,m,m,true);
		// write_field(F,cerr<<"U := ",U,m,n,n,true);
		// write_field(F,cerr<<"A := ",A,m,n,lda,true);
		// write_field(F,cerr<<"B := ",B,m,n,lda,true);
		// write_field(F,cerr<<"C := ",C,m,n,n,true);


	}
	else { /*  trans == FFLAS::FflasTrans */

		L = FFLAS::fflas_new(F,m,n);
		U = FFLAS::fflas_new(F,n,n);

		/*  build L */
		for (size_t i=0; i<R; ++i){
			for (size_t j=0; j<i; ++j)
				F.assign ( *(L + i + j*n), F.zero);
			for (size_t j=i+1; j<m; ++j)
				F.assign (*(L + i + j*n), *(B+ i+j*lda));
		}
		for (size_t i=R;i<n; ++i) {
			for (size_t j=0; j<m; ++j)
				F.assign(*(L+i+j*n), F.zero);
		}
		/*  build U */
		for ( size_t i=0; i<n; ++i ){
			size_t j=0;
			for (;  j< ((i<R)?i:R) ; ++j )
				F.assign( *(U + i+j*n), *(B+i+j*lda));
			for (; j<n; ++j )
				F.assign( *(U+i+j*n), F.zero);
		}
		//write_field(F,cerr<<"L = "<<endl,L,m,n,n);
		//write_field(F,cerr<<"U = "<<endl,U,n,n,n);

		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans,
				n,0,(int)R, U, n, Q);

		for (size_t i=0; i<n; ++i)
			F.assign (*(U+i*(n+1)),F.one);

		/*  reconstruct the diagonal */
		if (diag == FFLAS::FflasNonUnit) {
			for ( size_t i=0; i<R; ++i )
				F.assign (*(L+i*(n+1)), *(B+i*(lda+1)));
		}
		else{ // diag == FFLAS::FflasUnit
			for ( size_t i=0; i<R; ++i ){
				*(U+Q[i]*(n+1)) = *(B+Q[i]+i*lda);
				F.assign (*(L+i*(n+1)),F.one);
			}
		}
		// write_field(F,cerr<<"L = "<<endl,L,m,n,n);
		// write_field(F,cerr<<"U = "<<endl,U,n,n,n);

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

	FFLAS::fflas_delete(L);
	FFLAS::fflas_delete(U);
	FFLAS::fflas_delete(P);
	FFLAS::fflas_delete(Q);
	FFLAS::fflas_delete(B);
	FFLAS::fflas_delete(C);

	cout<<" LQUP MP: "<<(fail?"FAILED":"PASSED")<<" [rank="<<R<<"] ("<<chrono.usertime()<<")"<<endl;
	
	return fail;
}


//#define __ALL
template<class Field, FFLAS::FFLAS_DIAG diag, FFLAS::FFLAS_TRANSPOSE trans>
bool launch_test(const Field & F,
		 size_t r,
		 size_t m, size_t n)
{
	//typedef typename Field::Element Element ;
	typedef typename Field::Element_ptr Element_ptr ;
	bool fail = false ;

	{ /*  user given and lda bigger */
		size_t lda = n;//+10 ;
		Element_ptr A = FFLAS::fflas_new(F,m,lda);
		RandomMatrixWithRank(F,A,r,m,n,lda);
		fail |= test_lu<Field,diag,trans>(F,A,r,m,n,lda);
		if (fail) std::cout << "failed" << std::endl;
		FFLAS::fflas_delete(A);
	}
#ifdef __ALL
	{ /*  user given and lda bigger. Rank is max */ 
		size_t lda = n+10 ;
		size_t R = std::min(m,n);
		Element_ptr A = FFLAS::fflas_new(F,m,lda);
		RandomMatrixWithRank(F,A,R,m,n,lda);
		fail |= test_lu<Field,diag,trans>(F,A,R,m,n,lda);  
		if (fail) std::cout << "failed" << std::endl;
		FFLAS::fflas_delete(A);
	}
	{ /*  user given and lda bigger. Rank is min */
		size_t lda = n+10 ;
		size_t R = 0;
		Element_ptr A = FFLAS::fflas_new(F,m,lda);
		RandomMatrixWithRank(F,A,R,m,n,lda);
		fail |= test_lu<Field,diag,trans>(F,A,R,m,n,lda);
		if (fail) std::cout << "failed" << std::endl;
		FFLAS::fflas_delete(A);
	}

	{ /*  square  */
		size_t M = std::max(m,n);
		size_t N = M ;
		size_t R = M/2 ;
		size_t lda = N+10 ;
		Element_ptr A = FFLAS::fflas_new(F,M,lda);
		RandomMatrixWithRank(F,A,R,M,N,lda);
		fail |= test_lu<Field,diag,trans>(F,A,R,M,N,lda);
		if (fail) std::cout << "failed" << std::endl;
		FFLAS::fflas_delete(A);
	}
	{ /*  wide  */
		size_t M = std::max(m,n);
		size_t N = 2*M ;
		size_t R = 3*M/4 ;
		size_t lda = N+5 ;
		Element_ptr A = FFLAS::fflas_new(F,M,lda);
		RandomMatrixWithRank(F,A,R,M,N,lda);
		fail |= test_lu<Field,diag,trans>(F,A,R,M,N,lda);
		if (fail) std::cout << "failed" << std::endl;
		FFLAS::fflas_delete(A);
	}
	{ /*  narrow  */
		size_t M = std::max(m,n);
		size_t N = M/2 ;
		size_t R = 3*M/8 ;
		size_t lda = N+5 ;
		Element_ptr A = FFLAS::fflas_new(F,M,lda);
		RandomMatrixWithRank(F,A,R,M,N,lda);
		fail |= test_lu<Field,diag,trans>(F,A,R,M,N,lda);
		if (fail) std::cout << "failed" << std::endl;
		FFLAS::fflas_delete(A);
	}
#endif
	return fail;
}


int main(int argc, char** argv)
{
	cerr<<setprecision(20);
	int b = 100;
	size_t m = 50;
	size_t n = 50;
	size_t r = 20;
	int iter = 1 ;
	bool fail = false;
	double dr=0.8;
	static Argument as[] = {
		{ 'b', "-b B", "Set the bitsize of the field characteristic.",         TYPE_INT , &b },
		{ 'n', "-n N", "Set the number of cols in the matrix.", TYPE_INT , &n },
		{ 'm', "-m M", "Set the number of rows in the matrix.", TYPE_INT , &m },
		{ 'i', "-i I", "Set number of repetitions.",            TYPE_INT , &iter },
		{ 'r', "-r R", "Set defaut rank (%).",            TYPE_DOUBLE , &dr },
		// { 'f', "-f file", "Set input file", TYPE_STR, &file },
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(argc,argv,as);
	r=dr*std::min(m,n);
	
	typedef Modular<integer> Field;
	FFPACK::Integer p;
	FFPACK::Integer::random_exact_2exp(p, b);			
	nextprime(p,p);
	Field F(p);
	cout<<"p:="<<p<<";"<<endl;
	for (int i = 0 ; i < iter ; ++i) {
		 fail |= launch_test<Field,FFLAS::FflasNonUnit,FFLAS::FflasNoTrans>(F,r,m,n);
		 fail |= launch_test<Field,FFLAS::FflasUnit,FFLAS::FflasNoTrans>(F,r,m,n);
		 fail |= launch_test<Field,FFLAS::FflasUnit,FFLAS::FflasTrans>(F,r,m,n);
		 fail |= launch_test<Field,FFLAS::FflasNonUnit,FFLAS::FflasTrans>(F,r,m,n); 
	 }
	
	cout<<(fail?"FAILED":"PASSED")<<endl;

	return fail ;
 
}
