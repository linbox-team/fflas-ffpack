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
//                        Test for fgemm : 1 computation
//
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------

#ifndef NEWINO
#define NEWWINO
#endif
#define TIME 1

#include <iomanip>
#include <iostream>

//#include "fflas-ffpack/field/modular-positive.h"
//#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/field/modular-int32.h"
#include "utils/timer.h"
#include "Matio.h"
#include "fflas-ffpack/fflas/fflas.h"

#include "fflas-ffpack/utils/args-parser.h"
#include "test-utils.h"


using namespace std;
using namespace FFPACK;


//typedef Modular<double> Field;
//typedef Modular<float> Field;
// typedef ModularBalanced<double> Field;
//typedef ModularBalanced<float> Field;
//typedef Modular<int> Field;

// checks that D = alpha . C + beta . A ^ta * B ^tb
template<class Field>
bool check_MM(const Field                   & F,
	      const typename Field::Element * Cd, // c0
	      enum FFLAS::FFLAS_TRANSPOSE   & ta,
	      enum FFLAS::FFLAS_TRANSPOSE   & tb,
	      const size_t                    m,
	      const size_t                    n,
	      const size_t                    k,
	      const typename Field::Element & alpha,
	      const typename Field::Element * A,
	      const size_t                    lda,
	      const typename Field::Element * B,
	      const size_t                    ldb,
	      const typename Field::Element & beta,
	      const typename Field::Element * C, // res
	      const size_t                    ldc
	      )
{
	bool wrong = false;

	typedef typename Field::Element Element;
	Element aij, bij,  tmp;

	Element * D  = new Element[m*n];
	FFLAS::fcopy(F,m,n,D,n,Cd,n);

	for (size_t i = 0; i < m; ++i)
		for (size_t j = 0; j < n; ++j){
			F.mulin(*(D+i*n+j),beta);
			F.assign (tmp, F.zero);
			for ( size_t l = 0; l < k ; ++l ){
				if ( ta == FFLAS::FflasNoTrans )
					aij = *(A+i*lda+l);
				else
					aij = *(A+l*lda+i);
				if ( tb == FFLAS::FflasNoTrans )
					bij = *(B+l*ldb+j);
				else
					bij = *(B+j*ldb+l);
				F.axpyin (tmp, aij, bij);
			}
			F.axpyin (*(D+i*n+j), alpha, tmp);
			if ( !F.areEqual( *(D+i*n+j), *(C+i*ldc+j) ) ) {
				wrong = true;
			}
		}
	if ( wrong ){
		std::cerr<<"FAIL"<<std::endl;
		std::cerr << "m   :" << m   << ", n   : " <<  n  << ", k   : " << k << std::endl;
		std::cerr << "ldA :" << lda << ", ldB : " << ldb << ", ldC : " << ldc << std::endl;
		for (size_t i=0; i<m; ++i){
			for (size_t j =0; j<n; ++j)
				if (!F.areEqual( *(C+i*ldc+j), *(D+i*n+j) ) )
					std::cerr<<"Erreur C["<<i<<","<<j<<"]="
					<<(*(C+i*ldc+j))<<" C[d"<<i<<","<<j<<"]="
					<<(*(D+i*n+j))<<std::endl;
		}
	}

	delete[] D;

	return !wrong ;

}


template<class Field>
bool launch_MM(const Field & F,
	      const size_t   m,
	      const size_t   n,
	      const size_t   k,
	      const typename Field::Element alpha,
	      const typename Field::Element beta,
	      const size_t ldc,
	      const size_t lda,
	      enum FFLAS::FFLAS_TRANSPOSE    ta,
	      const size_t ldb,
	      enum FFLAS::FFLAS_TRANSPOSE    tb,
	      size_t iters,
	      int nbw )
{
	bool ok = true;

	typedef typename Field::Element Element;
	Element * A = new Element[m*lda];
	Element * B = new Element[k*ldb];
	Element * C = new Element[m*ldc];
	Element * D = new Element[m*n];
	for(size_t i = 0;i<iters;++i){
		RandomMatrix(F,A,m,k,lda);
		RandomMatrix(F,B,k,n,ldb);
		RandomMatrix(F,C,m,n,ldc);
		FFLAS::fcopy(F,m,n,D,n,C,ldc);
		FFLAS::fgemm (F, ta, tb,m,n,k,alpha, A,lda, B,ldb,
			      beta,C,ldc,nbw);
		ok &= check_MM(F, D, ta, tb,m,n,k,alpha, A,lda, B,ldb,
			      beta,C,ldc);

		if (!ok)
			break;

	}
	delete[] A;
	delete[] B;
	delete[] C;
	delete[] D;

	return ok ;
}


template<class Field>
bool launch_MM_dispatch(const Field &F,
			const size_t nn,
			const typename Field::Element alpha,
			const typename Field::Element beta,
			const size_t iters,
			const int nbw)
{
	bool ok = true;
	size_t m,n,k;
	size_t lda,ldb,ldc;
		//!@bug test for ldX equal
		//!@bug test for transpo
		//!@todo does nbw actually do nbw recursive calls and then call blas (check ?) ?
	size_t ld = 13 ;

	{
		FFLAS::FFLAS_TRANSPOSE ta = FFLAS::FflasNoTrans ;
		FFLAS::FFLAS_TRANSPOSE tb = FFLAS::FflasNoTrans ;
		if (random()%2) ta = FFLAS::FflasTrans ;
		if (random()%2) tb = FFLAS::FflasTrans ;

		m = 10+(size_t)random()%nn;
		n = 10+(size_t)random()%nn;
		k = 10+(size_t)random()%nn;
		lda = std::max(k,m)+(size_t)random()%ld;
		ldb = std::max(n,k)+(size_t)random()%ld;
		ldc = n+(size_t)random()%ld;
		ok &= launch_MM<Field>(F,m,n,k,
				       alpha,beta,
				       ldc,
				       lda, ta,
				       ldb, tb,
				       iters,nbw);
	}

	{
		FFLAS::FFLAS_TRANSPOSE ta = FFLAS::FflasNoTrans ;
		FFLAS::FFLAS_TRANSPOSE tb = FFLAS::FflasNoTrans ;
		// if (random()%2) ta = FFLAS::FflasTrans ;
		// if (random()%2) tb = FFLAS::FflasTrans ;

		m = 0;
		n = 10+(size_t)random()%nn;
		k = 10+(size_t)random()%nn;
		lda = std::max(k,m)+(size_t)random()%ld;
		ldb = std::max(n,k)+(size_t)random()%ld;
		ldc = n+(size_t)random()%ld;
		ok &= launch_MM<Field>(F,m,n,k,
				       alpha,beta,
				       ldc,
				       lda, ta,
				       ldb, tb,
				       iters,nbw);
	}

	{
	FFLAS::FFLAS_TRANSPOSE ta = FFLAS::FflasNoTrans ;
		FFLAS::FFLAS_TRANSPOSE tb = FFLAS::FflasNoTrans ;
		// if (random()%2) ta = FFLAS::FflasTrans ;
		// if (random()%2) tb = FFLAS::FflasTrans ;

		m = 10+(size_t)random()%nn;
		n = 0 ;
		k = 10+(size_t)random()%nn;
		lda = std::max(k,m)+(size_t)random()%ld;
		ldb = 1+std::max(n,k)+(size_t)random()%ld;
		ldc = 1+n+(size_t)random()%ld;
		ok &= launch_MM<Field>(F,m,n,k,
				       alpha,beta,
				       ldc,
				       lda, ta,
				       ldb, tb,
				       iters,nbw);
	}

	{
	FFLAS::FFLAS_TRANSPOSE ta = FFLAS::FflasNoTrans ;
		FFLAS::FFLAS_TRANSPOSE tb = FFLAS::FflasNoTrans ;
		// if (random()%2) ta = FFLAS::FflasTrans ;
		// if (random()%2) tb = FFLAS::FflasTrans ;

		m = 10+(size_t)random()%nn;
		n = 10+(size_t)random()%nn;
		k = 0;
		lda = 1+std::max(k,m)+(size_t)random()%ld;
		ldb = std::max(n,k)+(size_t)random()%ld;
		ldc = n+(size_t)random()%ld;
		ok &= launch_MM<Field>(F,m,n,k,
				       alpha,beta,
				       ldc,
				       lda, ta,
				       ldb, tb,
				       iters,nbw);
	}

	{
	FFLAS::FFLAS_TRANSPOSE ta = FFLAS::FflasNoTrans ;
		FFLAS::FFLAS_TRANSPOSE tb = FFLAS::FflasNoTrans ;
		// if (random()%2) ta = FFLAS::FflasTrans ;
		// if (random()%2) tb = FFLAS::FflasTrans ;

		m = 1;
		n = 10+(size_t)random()%nn;
		k = 10+(size_t)random()%nn;
		lda = std::max(k,m)+(size_t)random()%ld;
		ldb = std::max(n,k)+(size_t)random()%ld;
		ldc = n+(size_t)random()%ld;
		ok &= launch_MM<Field>(F,m,n,k,
				       alpha,beta,
				       ldc,
				       lda, ta,
				       ldb, tb,
				       iters,nbw);
	}


	{
		FFLAS::FFLAS_TRANSPOSE ta = FFLAS::FflasNoTrans ;
		FFLAS::FFLAS_TRANSPOSE tb = FFLAS::FflasNoTrans ;
		// if (random()%2) ta = FFLAS::FflasTrans ;
		// if (random()%2) tb = FFLAS::FflasTrans ;
		m = 10+(size_t)random()%nn;
		n = 1 ;
		k = 10+(size_t)random()%nn;
		lda = std::max(k,m)+(size_t)random()%ld;
		ldb = std::max(n,k)+(size_t)random()%ld;
		ldc = n+(size_t)random()%ld;
		// std::cout << m << ','
		// << n << ','
		// << k << ','
		// << lda << ','
		// << ldb << ','
		// << ldc << std::endl;
		ok &= launch_MM<Field>(F,m,n,k,
				       alpha,beta,
				       ldc,
				       lda, ta ,
				       ldb, tb ,
				       iters,nbw);
	}

	{
		FFLAS::FFLAS_TRANSPOSE ta = FFLAS::FflasNoTrans ;
		FFLAS::FFLAS_TRANSPOSE tb = FFLAS::FflasNoTrans ;
		// if (random()%2) ta = FFLAS::FflasTrans ;
		// if (random()%2) tb = FFLAS::FflasTrans ;

		m = 10+(size_t)random()%nn;
		n = 10+(size_t)random()%nn;
		k = 1;
		lda = std::max(k,m)+(size_t)random()%ld;
		ldb = std::max(n,k)+(size_t)random()%ld;
		ldc = n+(size_t)random()%ld;
		ok &= launch_MM<Field>(F,m,n,k,
				       alpha,beta,
				       ldc,
				       lda, ta,
				       ldb, tb,
				       iters,nbw);
	}


	return ok ;
}

int main(int argc, char** argv)
{


	static size_t iters =10 ;
	static unsigned long p = 65521 ;
	static size_t n = 100 ;
	static int nbw = 2 ;

	static Argument as[] = {
		{ 'p', "-p P", "Set the field characteristic.",         TYPE_INT , &p },
		{ 'n', "-n N", "Set the dimension of the matrix.",      TYPE_INT , &n },
		{ 'w', "-w N", "Set the number of winograd levels.",    TYPE_INT , &nbw },
		{ 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(argc,argv,as);

	bool ok = true ;

	{
		typedef Modular<double> Field ;
		typedef Field::Element Element ;
		typedef Field::RandIter Randiter ;
		typedef Field::Element  Element ;

		unsigned long q = std::min((unsigned long)Field::getMaxModulus(),p);
		Field F(q);

		Randiter R1(F);
		NonzeroRandIter<Field,Randiter> R(F,R1);

		ok &= launch_MM_dispatch<Field>(F,n,F.one ,F.zero,iters,nbw);
		ok &= launch_MM_dispatch<Field>(F,n,F.zero,F.zero,iters,nbw);
		ok &= launch_MM_dispatch<Field>(F,n,F.mOne,F.zero,iters,nbw);

		ok &= launch_MM_dispatch<Field>(F,n,F.one ,F.one,iters,nbw);
		ok &= launch_MM_dispatch<Field>(F,n,F.zero,F.one,iters,nbw);
		ok &= launch_MM_dispatch<Field>(F,n,F.mOne,F.one,iters,nbw);

		ok &= launch_MM_dispatch<Field>(F,n,F.one ,F.mOne,iters,nbw);
		ok &= launch_MM_dispatch<Field>(F,n,F.zero,F.mOne,iters,nbw);
		ok &= launch_MM_dispatch<Field>(F,n,F.mOne,F.mOne,iters,nbw);

		Element alpha,beta ;
		R.random(alpha);

		ok &= launch_MM_dispatch<Field>(F,n,F.one ,alpha,iters,nbw);
		ok &= launch_MM_dispatch<Field>(F,n,F.zero,alpha,iters,nbw);
		ok &= launch_MM_dispatch<Field>(F,n,F.mOne,alpha,iters,nbw);

		ok &= launch_MM_dispatch<Field>(F,n,alpha,F.one ,iters,nbw);
		ok &= launch_MM_dispatch<Field>(F,n,alpha,F.zero,iters,nbw);
		ok &= launch_MM_dispatch<Field>(F,n,alpha,F.mOne,iters,nbw);

		for (size_t j = 0 ; j < 9 ; ++j) {
			R.random(alpha);
			R.random(beta);
			ok &= launch_MM_dispatch<Field>(F,n,alpha,beta,iters,nbw);
		}
	}

	return !ok ;
}


