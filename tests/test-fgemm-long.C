/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/*
 * Copyright (C) the FFLAS-FFPACK group
 * Written by Clément Pernet
 *            Brice Boyer <bbboyer@ncsu.edu>
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

// #ifndef NEWINO
// #define NEWWINO
// #endif

// #define WINOTHRESHOLD 100

#define TIME 1

#include <iomanip>
#include <iostream>

//#include "fflas-ffpack/field/modular-positive.h"
//#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/field/modular-int32.h"
#include "fflas-ffpack/utils/timer.h"
#include "Matio.h"
#include "fflas-ffpack/fflas/fflas.h"

#include "fflas-ffpack/utils/args-parser.h"
#include "test-utils.h"
#include "givaro/givintprime.h"


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
	typedef typename Field::Element_ptr Element_ptr;
	typedef typename Field::ConstElement_ptr ConstElement_ptr;
	Element tmp;
	ConstElement_ptr ail,blj;
	Element_ptr D  = FFLAS::fflas_new (F,m,n);
	FFLAS::fcopy(F,m,n,Cd,n,D,n);

	for (size_t i = 0; i < m; ++i)
		for (size_t j = 0; j < n; ++j){
			F.mulin(*(D+i*n+j),beta);
			F.assign (tmp, F.zero);
			for ( size_t l = 0; l < k ; ++l ){
				if ( ta == FFLAS::FflasNoTrans )
					ail = A+i*lda+l;
				else
					ail = A+l*lda+i;
				if ( tb == FFLAS::FflasNoTrans )
					blj = B+l*ldb+j;
				else
					blj = B+j*ldb+l;
				F.axpyin (tmp, *ail, *blj);
			}
			F.axpyin (*(D+i*n+j), alpha, tmp);
			if ( !F.areEqual( *(D+i*n+j), *(C+i*ldc+j) ) ) {
				wrong = true;
			}
		}
	if ( wrong ){
		size_t ici = 20 ;
		std::cout<<"FAIL"<<std::endl;
		std::cout << "a   :" << alpha<<", b   : " << beta << std::endl;
		std::cout << "m   :" << m   << ", n   : " <<  n  << ", k   : " << k << std::endl;
		std::cout << "ldA :" << lda << ", ldB : " << ldb << ", ldC : " << ldc << std::endl;
		for (size_t i=0; i<m && ici; ++i){
			for (size_t j =0; j<n && ici; ++j)
				if (!F.areEqual( *(C+i*ldc+j), *(D+i*n+j) ) ) {
					std::cout<<"Error C["<<i<<","<<j<<"]="
					<<(*(C+i*ldc+j))<<" D["<<i<<","<<j<<"]="
					<<(*(D+i*n+j))<<std::endl;
					ici--;
				}
		}
		if (m<80 && n<80) {
			for (size_t i=0; i<m ; ++i){
				for (size_t j =0; j<n ; ++j) {
					if ( !F.areEqual( *(C+i*ldc+j), *(D+i*n+j) ) )
						std::cout << 'X' ;
					else
						std::cout << '.' ;
				}
				std::cout << std::endl;
			}
		}
	}
	// else std::cout<<"COOL"<<std::endl;

	FFLAS::fflas_delete (D);

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
	       int nbw,
	       bool par)
{
	bool ok = true;

	typedef typename Field::Element_ptr Element_ptr;
	Element_ptr A ;
	Element_ptr B ;
	Element_ptr C = FFLAS::fflas_new (F,m,ldc);
	FFLASFFPACK_check(ldc >= n);
	FFLAS::fzero(F,m,n,C,ldc);
	Element_ptr D = FFLAS::fflas_new (F, m, n);
	for(size_t i = 0;i<iters;++i){
		if (ta == FFLAS::FflasNoTrans) {
			FFLASFFPACK_check(lda >= k);
			A = FFLAS::fflas_new (F, m, lda);
			FFLAS::fzero(F,m,lda,A,lda);
			RandomMatrix(F,A,m,k,lda);
		}
		else {
			FFLASFFPACK_check(lda >= m);
			A = FFLAS::fflas_new (F, k, lda);
			FFLAS::fzero(F,k,lda,A,lda);
			RandomMatrix(F,A,k,m,lda);
		}
		if (tb == FFLAS::FflasNoTrans) {
			FFLASFFPACK_check(ldb >= n);
			B = FFLAS::fflas_new (F,k,ldb);
			FFLAS::fzero(F,k,ldb,B,ldb);
			RandomMatrix(F,B,k,n,ldb);
		}
		else {
			FFLASFFPACK_check(ldb >= k);
			B = FFLAS::fflas_new (F,n,ldb);
			FFLAS::fzero(F,n,ldb,B,ldb);
			RandomMatrix(F,B,n,k,ldb);
		}
		RandomMatrix(F,C,m,n,ldc);
		FFLAS::fcopy(F,m,n,C,ldc,D,n);
		if (par){
			FFLAS::MMHelper<Field,FFLAS::MMHelperAlgo::Winograd> WH(F,nbw,FFLAS::ParSeqHelper::Parallel());
			PAR_REGION{
				FFLAS::fgemm (F, ta, tb,m,n,k,alpha, A,lda, B,ldb, beta,C,ldc,WH);
			}
		}else{
			FFLAS::MMHelper<Field,FFLAS::MMHelperAlgo::Winograd> WH(F,nbw,FFLAS::ParSeqHelper::Sequential());
			FFLAS::fgemm (F, ta, tb,m,n,k,alpha, A,lda, B,ldb, beta,C,ldc,WH);
		}
		ok &= check_MM(F, D, ta, tb,m,n,k,alpha, A,lda, B,ldb, beta,C,ldc);

		FFLAS::fflas_delete(A);
		FFLAS::fflas_delete(B);

		if (!ok)
			break;


	}
	FFLAS::fflas_delete (C);
	FFLAS::fflas_delete (D);

	return ok ;
}


template<class Field>
bool launch_MM_dispatch(const Field &F,
			const size_t nn,
			const typename Field::Element alpha,
			const typename Field::Element beta,
			const size_t iters,
			const int nbw,
			const bool par)
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

		m = 1+(size_t)random()%nn;
		n = 1+(size_t)random()%nn;
		k = 1+(size_t)random()%nn;

		int logdim = (int)floor(log2(std::min(std::min(m,k),n)));
		int nw = std::min (logdim,nbw);

		lda = std::max(k,m)+(size_t)random()%ld;
		ldb = std::max(n,k)+(size_t)random()%ld;
		ldc = n+(size_t)random()%ld;
		std::cout <<"q = "<<F.characteristic()<<" nw = "<<nw<<" m,k,n = "<<m<<", "<<k<<", "<<n<<" C := "
			  <<alpha<<".A"<<((ta==FFLAS::FflasTrans)?"^T":"")
			  <<" * B"<<((tb==FFLAS::FflasTrans)?"^T":"");
		if (!F.isZero(beta))
			cout<<" + "<<beta<<" C";
		ok &= launch_MM<Field>(F,m,n,k,
				       alpha,beta,
				       ldc,
				       lda, ta,
				       ldb, tb,
				       iters,nw, par);
		std::cout<<(ok?" -> ok ":" -> KO")<<std::endl;
	}
	return ok ;
}
template <class Field>
bool run_with_field (int q, unsigned long b, size_t n, int nbw, size_t iters, bool par ){
	bool ok = true ;
	unsigned long p=q;
	int nbit=(int)iters;
	while (ok &&  nbit){
		typedef typename  Field::Element Element ;
		typedef typename Field::RandIter Randiter ;
		typedef typename Field::Element  Element ;
		if (nbw<0)
			nbw = (int) random() % 7;
		if (q==-1){
			b = 2 + (rand() % (int)floor(log2(Field::getMaxModulus())+1));
		}
		Givaro::IntPrimeDom IPD;
		Givaro::Integer tmp;
		if (b > 1){
			    // Choose characteristic as a random prime of b bits
			do{
				Givaro::Integer _p;
				Givaro::Integer::random_exact_2exp(_p,b);//max % (2<<30);
				IPD.prevprime( tmp, _p+1 );
				p =  tmp;
			}while( (p < 2) );
		}
		if (p > (unsigned long)Field::getMaxModulus()){
			IPD.prevprime( tmp, Field::getMaxModulus()+1 );
			p=tmp;
		}
		p = (int)std::max((unsigned long) Field::getMinModulus(),(unsigned long)p);
		Field F((int)p);

		Randiter R1(F);
		NonzeroRandIter<Field,Randiter> R(F,R1);

		    //size_t k = 0 ;
		    //std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(F,n,F.one ,F.zero,iters,nbw, par);
		    //std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(F,n,F.zero,F.zero,iters,nbw, par);
		    //std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(F,n,F.mOne,F.zero,iters,nbw, par);
		    //std::cout << k << "/24" << std::endl; ++k;

		ok &= launch_MM_dispatch<Field>(F,n,F.one ,F.one,iters,nbw, par);
		    //std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(F,n,F.zero,F.one,iters,nbw, par);
		    //std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(F,n,F.mOne,F.one,iters,nbw, par);
		    //std::cout << k << "/24" << std::endl; ++k;

		ok &= launch_MM_dispatch<Field>(F,n,F.one ,F.mOne,iters,nbw, par);
		    //std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(F,n,F.zero,F.mOne,iters,nbw, par);
		//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(F,n,F.mOne,F.mOne,iters,nbw, par);
		//std::cout << k << "/24" << std::endl; ++k;

		Element alpha,beta ;
		R.random(alpha);

		ok &= launch_MM_dispatch<Field>(F,n,F.one ,alpha,iters,nbw, par);
		//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(F,n,F.zero,alpha,iters,nbw, par);
		//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(F,n,F.mOne,alpha,iters,nbw, par);
		//std::cout << k << "/24" << std::endl; ++k;

		ok &= launch_MM_dispatch<Field>(F,n,alpha,F.one ,iters,nbw, par);
		//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(F,n,alpha,F.zero,iters,nbw, par);
		//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(F,n,alpha,F.mOne,iters,nbw, par);
		//std::cout << k << "/24" << std::endl; ++k;

		for (size_t j = 0 ; j < 3 ; ++j) {
			R.random(alpha);
			R.random(beta);
			ok &= launch_MM_dispatch<Field>(F,n,alpha,beta,iters,nbw, par);
			//std::cout << k << "/24" << std::endl; ++k;
		}
		    //std::cout<<std::endl;
		nbit--;
	}
	return ok;
}
int main(int argc, char** argv)
{
	std::cout<<setprecision(17);
	srand((int)time(NULL));
	srand48(time(NULL));

	static size_t iters = 3 ;
	static int p = -1 ;
	static unsigned long b = 0 ;
	static size_t n = 50 ;
	static int nbw = -1 ;
	static bool loop = false;
	static bool par = false;
	static Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INT , &p },
		{ 'b', "-b B", "Set the bitsize of the random characteristic.",         TYPE_INT , &b },
		{ 'n', "-n N", "Set the dimension of the matrix.",      TYPE_INT , &n },
		{ 'w', "-w N", "Set the number of winograd levels (-1 for random).",    TYPE_INT , &nbw },
		{ 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
		{ 'l', "-loop Y/N", "run the test in an infinte loop.", TYPE_BOOL , &loop },
		{ 'p', "-par Y/N", "run the parallel fgemm.", TYPE_BOOL , &par },
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(argc,argv,as);


	bool ok = true;
	do{
		std::cout<<"Modular<double>"<<std::endl;
		ok &= run_with_field<Modular<double> >(p,b,n,nbw,iters,par);
		std::cout<<"ModularBalanced<double>"<<std::endl;
		ok &= run_with_field<ModularBalanced<double> >(p,b,n,nbw,iters,par);
		std::cout<<"Modular<float>"<<std::endl;
		ok &= run_with_field<Modular<float> >(p,b,n,nbw,iters,par);
		std::cout<<"ModularBalanced<float>"<<std::endl;
		ok &= run_with_field<ModularBalanced<float> >(p,b,n,nbw,iters,par);
		std::cout<<"Modular<int32_t>"<<std::endl;
		ok &= run_with_field<Modular<int32_t> >(p,b,n,nbw,iters,par);
		std::cout<<"ModularBalanced<int32_t>"<<std::endl;
		ok &= run_with_field<ModularBalanced<int32_t> >(p,b,n,nbw,iters,par);
		    // std::cout<<"Modular<int64_t>"<<std::endl;
		// ok &= run_with_field<Modular<int64_t> >(p,b,n,nbw,iters);
		// std::cout<<"ModularBalanced<int64_t>"<<std::endl;
		// ok &= run_with_field<ModularBalanced<int64_t> >(p,b,n,nbw,iters);
	} while (loop && ok);

	return !ok ;
}


