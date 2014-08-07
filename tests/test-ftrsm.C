/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Ziad Sultan
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
//                        Test for pftrsm
//
//--------------------------------------------------------------------------
// Ziad Sultan
//-------------------------------------------------------------------------


#define __FFLASFFPACK_USE_OPENMP
//#define __FFLASFFPACK_USE_KAAPI

//#define __FFLASFFPACK_FORCE_SEQ

#define __FFLAS__TRSM_READONLY

#include <iomanip>
#include <iostream>
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "time.h"
#include "fflas-ffpack/utils/args-parser.h"


using namespace std;
using namespace FFPACK;

//typedef Modular<double> Field;
typedef ModularBalanced<double> Field;
//typedef ModularBalanced<float> Field;



Field::Element_ptr makemat(const Field& F, int m, int n){
	Field::Element_ptr res = FFLAS::fflas_new(F,m,n);
	Field::RandIter RF(F);
        for (long i = 0; i < m; ++i)
                for (long j = 0; j < n; ++j) {
                        RF.random(res[j+i*n]);
                }
        return res;
}

Field::Element_ptr maketriangmat(const Field& F, int n){
	Field::Element_ptr res = FFLAS::fflas_new(F, n,n);
	Field::RandIter RF(F);
	Field::NonZeroRandIter NRF(F,RF);
        for (long i = 0; i < n; ++i){
                for (long j = 0; j < n; ++j)
                        RF.random(res[j+i*n]);
		NRF.random(res[i*(n+1)]);
	}
        return res;
}



template<class Field>
bool check_TRSM(const Field                   & F,
		enum FFLAS::FFLAS_SIDE side,
		enum FFLAS::FFLAS_UPLO uplo,
		enum FFLAS::FFLAS_TRANSPOSE   & trans,
		enum FFLAS::FFLAS_DIAG  & diag,
		const size_t                    m,
		const size_t                    n,
		const size_t                    k,
		const typename Field::Element & alpha,
		typename Field::Element * A,
		typename Field::Element * B,
		const typename Field::Element * B2 // res
		)
{

	typename Field::Element invalpha;
	F.inv(invalpha, alpha);

	FFLAS::ftrmm (F, side, uplo, trans, diag, m, n, invalpha, A, k, B, n);
	bool wrong = false;

	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<n;++j)
			if ( !F.areEqual(*(B2+i*n+j), *(B+i*n+j))){
				cerr<<"B2 ["<<i<<", "<<j<<"] = "<<(*(B2+i*n+j))
				    <<" ; B ["<<i<<", "<<j<<"] = "<<(*(B+i*n+j))
				    <<endl;
				wrong = true;
			}

	if ( wrong ){
		cerr<<"FAIL"<<endl;
	}
	//	else   cerr<<"PASS"<<endl;

	return !wrong;
}


//BEGIN_PARALLEL_MAIN(int argc, char** argv)
int main(int argc, char** argv)
{
        srand((int)time(NULL));
        srand48(time(NULL));
	size_t m,n, nbit, iters;
	iters=5;
	unsigned long q;
	q=65521;
        m = 20+(size_t)random()% 100;
        n = 20+(size_t)random()% 100;
	bool p = false;
	int s = 1; int u = 0; int t = 1; int d = 0;


        static Argument as[] = {
                { 'q', "-q Q", "Set the field characteristic.",                TYPE_INT , &q },
                { 'm', "-m M", "Set the row dimension of unknown matrix.",      TYPE_INT , &m },
                { 'n', "-n N", "Set the column dimension of the unknown matrix.", TYPE_INT , &n },
                { 's', "-s Side", "Set the side of the matrix.",               TYPE_INT , &s },
                { 'u', "-u UpLow", "Set the triangular side of the matrix.",   TYPE_INT , &u },
                { 't', "-t Trans", "Set the transposition of the matrix.",     TYPE_INT , &t },
                { 'd', "-d Diag", "Set the Diag of the matrix.",                TYPE_INT , &d },
                { 'i', "-i i", "Set number of repetitions.",               TYPE_INT , &iters },
		{ 'p', "-par Y/N", "run the parallel ftrsm.", TYPE_BOOL , &p },
                END_OF_ARGUMENTS
        };

	FFLAS::parseArguments(argc,argv,as);


	FFLAS::FFLAS_SIDE side =  (s) ? FFLAS::FflasRight :  FFLAS::FflasLeft;
	FFLAS::FFLAS_UPLO uplo =  u ? FFLAS::FflasLower : FFLAS::FflasUpper;
	FFLAS::FFLAS_TRANSPOSE trans = t ? FFLAS::FflasTrans :  FFLAS::FflasNoTrans;
	FFLAS::FFLAS_DIAG diag = d ? FFLAS::FflasUnit : FFLAS::FflasNonUnit;

	srand48(BaseTimer::seed());

	Field F(q);

	nbit= iters; // number of times the product is performed

	size_t k = m;
	if (side == FFLAS::FflasRight)
		k = n;
	Field::Element_ptr A = maketriangmat (F,k);
	Field::Element_ptr B = makemat(F,m,n);
	Field::Element_ptr B2 = FFLAS::fflas_new (F, m,n);
	FFLAS::fcopy(F, m, n, B, n, B2, n);

	Field::Element alpha;
	F.init (alpha, 1.0);

	if (   ((side == FFLAS::FflasRight) &&(k != n))
	    || ((side == FFLAS::FflasLeft)&&(k != m))) {
		cerr<<"Error in the dimensions of the input matrices"<<endl;
		exit(-1);
	}


        const FFLAS::CuttingStrategy Strategy = FFLAS::BLOCK_THREADS;

	for(size_t i = 0;i<nbit;++i){

		if (p){
			PAR_REGION{
				FFLAS::TRSMHelper<FFLAS::StructureHelper::Iterative,
						  FFLAS::ParSeqHelper::Parallel> PH (FFLAS::ParSeqHelper::Parallel(MAX_THREADS,Strategy));
				FFLAS::pftrsm (F, side, uplo, trans, diag, m, n, alpha, A, k, B, n, PH);
		}
			BARRIER;
		} else {
			FFLAS::ftrsm (F, side, uplo, trans, diag, m, n, alpha, A, k, B, n);
		}
		if (i+1<nbit)
                        for (size_t j=0; j<m*n;++j)
                                F.assign(*(B+j),*(B2+j));

	}


	bool ok = true ;
	ok = ok && check_TRSM(F, side, uplo, trans, diag, m, n, k, alpha, A, B, B2);

	FFLAS::fflas_delete(A);
	FFLAS::fflas_delete(B);
	FFLAS::fflas_delete(B2);
	return ok ? 0 : -1;
}
//END_PARALLEL_MAIN()
