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

//--------------------------------------------------------------------------
//          Test for the computations of rank profiles
//--------------------------------------------------------------------------
#define  __FFLASFFPACK_SEQUENTIAL
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"

#include <iostream>
#include <iomanip>
#include <givaro/modular.h>

#include "test-utils.h"
#include "Matio.h"

using namespace FFPACK;


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

		size_t lda = n;
		typename Field::Element_ptr A=FFLAS::fflas_new (*F, m,lda);
		typename Field::Element_ptr B=FFLAS::fflas_new (*F, m,lda);
		RandomMatrixWithRankandRandomRPM(*F,A,lda,r,m,n);
		FFLAS::fassign (*F, m, n, A, lda, B, lda); 

		{
				// Testing if LUdivine and PLUQ return the same result
			size_t* RP1, * RP2;
			FFPACK::RowRankProfile (*F, m, n, A, lda, RP1, FFPACK::FfpackSlabRecursive);
			FFLAS::fassign (*F, m, n, B, lda, A, lda); 
			FFPACK::RowRankProfile (*F, m, n, A, lda, RP2, FFPACK::FfpackTileRecursive);
			for (size_t i=0; i<r; i++)
				ok &= (RP1[i] == RP2[i]);
			FFLAS::fflas_delete(RP1);
			FFLAS::fflas_delete(RP2);

			FFLAS::fassign (*F, m, n, B, lda, A, lda); 
			FFPACK::ColumnRankProfile (*F, m, n, A, lda, RP1, FFPACK::FfpackSlabRecursive);
			FFLAS::fassign (*F, m, n, B, lda, A, lda); 
			FFPACK::ColumnRankProfile (*F, m, n, A, lda, RP2, FFPACK::FfpackTileRecursive);
			for (size_t i=0; i<r; i++)
				ok &= (RP1[i] == RP2[i]);
			FFLAS::fflas_delete(RP1);
			FFLAS::fflas_delete(RP2);
		}
		{
			// Testing if 1 PLUQ computes the rank profiles of all leading submatrices 
			size_t* RP1, * RP2;
			size_t * P = FFLAS::fflas_new<size_t>(m);
			size_t * Q = FFLAS::fflas_new<size_t>(n);
			FFLAS::fassign (*F, m, n, B, lda, A, lda); 
			PLUQ(*F, FFLAS::FflasNonUnit, m, n, A, lda, P, Q);
			
			for (size_t i=0; i<1;i++){
				size_t mm = 1 + (rand() % m);
				size_t nn = 1 + (rand() % n);
				FFLAS::fassign (*F, m, n, B, lda, A, lda); 
				size_t rr = FFPACK::ColumnRankProfile (*F, mm, nn, A, lda, RP1, FFPACK::FfpackSlabRecursive);
				FFLAS::fassign (*F, m, n, B, lda, A, lda); 
				FFPACK::RowRankProfile (*F, mm, nn, A, lda, RP2, FFPACK::FfpackSlabRecursive);
				size_t* RRP = FFLAS::fflas_new<size_t>(r);
				size_t* CRP = FFLAS::fflas_new<size_t>(r);
				
				LeadingSubmatrixRankProfiles (m,n,r,mm,nn,P,Q,RRP,CRP);
				for (size_t ii=0; ii<rr; ii++)
					ok &= (RP1[ii] == CRP[ii]) && (RP2[ii] == RRP[ii]);
				
				FFLAS::fflas_delete(CRP);
				FFLAS::fflas_delete(RRP);
				FFLAS::fflas_delete(RP1);
				FFLAS::fflas_delete(RP2);
				
			}
			FFLAS::fflas_delete(P);
			FFLAS::fflas_delete(Q);
		}
		{
            // Testing PLUQ and LUDivine return a specified rank profile 
			size_t* RRP = FFLAS::fflas_new<size_t>(r);
			size_t* CRP = FFLAS::fflas_new<size_t>(r);
			size_t* RRPLUD, * RRPPLUQ, *CRPLUD, *CRPPLUQ;
			RandomRankProfile (m, r, RRP);
			RandomRankProfile (n, r, CRP);
			
			RandomMatrixWithRankandRPM(*F,A,lda,r,m,n, RRP, CRP);
			FFLAS::fassign (*F, m, n, A, lda, B, lda); 
			size_t cs = FFPACK::ColumnRankProfile (*F, m, n, A, lda, CRPLUD, FFPACK::FfpackSlabRecursive);
			FFLAS::fassign (*F, m, n, B, lda, A, lda); 
			size_t ct = FFPACK::ColumnRankProfile (*F, m, n, A, lda, CRPPLUQ, FFPACK::FfpackTileRecursive);
			FFLAS::fassign (*F, m, n, B, lda, A, lda); 
			size_t rs = FFPACK::RowRankProfile (*F, m, n, A, lda, RRPLUD, FFPACK::FfpackSlabRecursive);
			FFLAS::fassign (*F, m, n, B, lda, A, lda); 
			size_t rt = FFPACK::RowRankProfile (*F, m, n, A, lda, RRPPLUQ, FFPACK::FfpackTileRecursive);
			// write_perm (std::cout<<"RRP     = ", RRP, r);
			// write_perm (std::cout<<"CRP     = ", CRP, r);
			std::sort(CRP,CRP+r);
			std::sort(RRP,RRP+r);	
			ok &= (cs==ct)&(cs==rs)&(cs==rt)&(cs==r);
			for (size_t i=0; i<r; i++)
				ok &= (CRPLUD[i] == CRP[i]) && (CRPPLUQ[i] == CRP[i]) && 
					(RRPLUD[i] == RRP[i]) && (RRPPLUQ[i] == RRP[i]);
			FFLAS::fflas_delete(CRP);
			FFLAS::fflas_delete(RRP);
			FFLAS::fflas_delete(CRPLUD);
			FFLAS::fflas_delete(RRPLUD);
			FFLAS::fflas_delete(CRPPLUQ);
			FFLAS::fflas_delete(RRPPLUQ);
		}

		
		FFLAS::fflas_delete(A);
		FFLAS::fflas_delete(B);
		delete F;

		nbit--;
		if (!ok)
				//std::cout << "\033[1;31mFAILED\033[0m "<<std::endl;		
			std::cout << "FAILED "<<std::endl;
		else
				//std::cout << "\033[1;32mPASSED\033[0m "<<std::endl;
			std::cout << "PASSED "<<std::endl;
	}
	return ok;
	}

int main(int argc, char** argv){
	std::cerr<<std::setprecision(20);

	Givaro::Integer q = -1;
	size_t b = 0;
	size_t m = 150;
	size_t n = 280;
	size_t r = 85;
	size_t iters = 6 ;
	bool loop=false;
	static Argument as[] = {
		{ 'q', "-q Q", "Set the field cardinality.",         TYPE_INTEGER , &q },
		{ 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
		{ 'n', "-n N", "Set the number of cols in the matrix.", TYPE_INT , &n },
		{ 'm', "-m N", "Set the number of rows in the matrix.", TYPE_INT , &m },
		{ 'r', "-r r", "Set the rank of the matrix."          , TYPE_INT , &r },
		{ 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
		{ 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
		    // { 'f', "-f file", "Set input file", TYPE_STR, &file },
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
		ok&=run_with_field<Givaro::ModularBalanced<double> >   (q,b,m,n,r,iters);
		ok&=run_with_field<Givaro::Modular<int32_t> >   (q,b,m,n,r,iters);
		ok&=run_with_field<Givaro::ModularBalanced<int32_t> >   (q,b,m,n,r,iters);
		ok&=run_with_field<Givaro::Modular<int64_t> >   (q,b,m,n,r,iters);
		ok&=run_with_field<Givaro::ModularBalanced<int64_t> >   (q,b,m,n,r,iters);
		ok&=run_with_field<Givaro::Modular<Givaro::Integer> >(q,(b?b:128),m/4+1,n/4+1,r/4+1,iters); 
	} while (loop && ok);

	return !ok;
}
