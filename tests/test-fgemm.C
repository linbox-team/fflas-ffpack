/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/*
 * Copyright (C) the FFLAS-FFPACK group
 * Written by Clément Pernet
 *            Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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


// #ifndef NEWINO
// #define NEWWINO
// #endif

// #define WINOTHRESHOLD 100
// #define OLD_DYNAMIC_PEELING
//#define DEBUG 1

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iomanip>
#include <iostream>
#include <givaro/modular.h>
#include <givaro/udl.h>
#include <recint/rint.h>

#include <givaro/givintprime.h>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"

#include "fflas-ffpack/utils/args-parser.h"
#include "test-utils.h"
#include "fflas-ffpack/utils/Matio.h"



using namespace std;
using namespace FFPACK;
using Givaro::Modular;
using Givaro::ModularBalanced;


// checks that D = alpha . C + beta . A ^ta * B ^tb
template<class Field>
bool check_MM(const Field                   & F,
	      const typename Field::Element_ptr  Cd, // c0
	      enum FFLAS::FFLAS_TRANSPOSE   & ta,
	      enum FFLAS::FFLAS_TRANSPOSE   & tb,
	      const size_t                    m,
	      const size_t                    n,
	      const size_t                    k,
	      const typename Field::Element & alpha,
	      const typename Field::Element_ptr  A,
	      const size_t                    lda,
	      const typename Field::Element_ptr  B,
	      const size_t                    ldb,
	      const typename Field::Element & beta,
	      const typename Field::Element_ptr  C, // res
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
	FFLAS::fassign(F,m,n,Cd,n,D,n);

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
		std::cerr<<"FAIL"<<std::endl;
		std::cerr << "a   :" << alpha<<", b   : " << beta << std::endl;
		std::cerr << "m   :" << m   << ", n   : " <<  n  << ", k   : " << k << std::endl;
		std::cerr << "ldA :" << lda << ", ldB : " << ldb << ", ldC : " << ldc << std::endl;
		for (size_t i=0; i<m && ici; ++i){
			for (size_t j =0; j<n && ici; ++j)
				if (!F.areEqual( *(C+i*ldc+j), *(D+i*n+j) ) ) {
					std::cerr<<"Error C["<<i<<","<<j<<"]="
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
			   bool par, 
			   size_t b)
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
			RandomMatrix(F,A,m,k,lda,b);
		}
		else {
			FFLASFFPACK_check(lda >= m);
			A = FFLAS::fflas_new (F, k, lda);
			FFLAS::fzero(F,k,lda,A,lda);
			RandomMatrix(F,A,k,m,lda,b);
		}
		if (tb == FFLAS::FflasNoTrans) {
			FFLASFFPACK_check(ldb >= n);
			B = FFLAS::fflas_new (F,k,ldb);
			FFLAS::fzero(F,k,ldb,B,ldb);
			RandomMatrix(F,B,k,n,ldb,b);
		}
		else {
			FFLASFFPACK_check(ldb >= k);
			B = FFLAS::fflas_new (F,n,ldb);
			FFLAS::fzero(F,n,ldb,B,ldb);
			RandomMatrix(F,B,n,k,ldb,b);
		}
		RandomMatrix(F,C,m,n,ldc,b);
		FFLAS::fassign(F,m,n,C,ldc,D,n);
		if (par){
			FFLAS::MMHelper<Field,FFLAS::MMHelperAlgo::Auto, typename FFLAS::ModeTraits<Field>::value, FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::ThreeDAdaptive> > WH (F, nbw);
			PAR_BLOCK{
				FFLAS::fgemm (F, ta, tb,m,n,k,alpha, A,lda, B,ldb, beta,C,ldc,WH);
			}
		}else{
			FFLAS::MMHelper<Field,FFLAS::MMHelperAlgo::Auto> WH(F,nbw,FFLAS::ParSeqHelper::Sequential());
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
						const int mm,
						const int nn,
						const int kk,
						const typename Field::Element alpha,
						const typename Field::Element beta,
						const size_t iters,
						const int nbw,
						const bool par,
						const size_t b)
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
		if (! par) {
			if (random()%2) ta = FFLAS::FflasTrans ;
			if (random()%2) tb = FFLAS::FflasTrans ;
		}

		if (mm<0)
			m = 1+(size_t)random() % -mm;
		else m = mm;
		if (nn<0)
			n = 1+(size_t)random() % -nn;
		else n = nn;
		if (kk<0)
			k = 1+(size_t)random() % -kk;
		else k = kk;

		int logdim = (int)floor(log(std::min(std::min(m,k),n))/log(2.));
		int nw = std::min (logdim,nbw);

		lda = std::max(k,m)+(size_t)random()%ld;
		ldb = std::max(n,k)+(size_t)random()%ld;
		ldc = n+(size_t)random()%ld;
#ifdef DEBUG
		std::cerr <<"q = "<<F.characteristic()<<" nw = "<<nw<<" m,k,n = "<<m<<", "<<k<<", "<<n<<" C := "
			  <<alpha<<".A"<<((ta==FFLAS::FflasTrans)?"^T":"")
			  <<" * B"<<((tb==FFLAS::FflasTrans)?"^T":"");
		if (!F.isZero(beta))
			cerr<<" + "<<beta<<" C";
#endif
		ok &= launch_MM<Field>(F,m,n,k,
							   alpha,beta,
							   ldc,
							   lda, ta,
							   ldb, tb,
							   iters,nw, par, b);
#ifdef DEBUG
		std::cerr<<(ok?" -> ok ":" -> KO")<<std::endl;
#endif
	}
	return ok ;
}

template <class Field>
bool run_with_field (Givaro::Integer q, uint64_t b, int m, int n, int k, int nbw, size_t iters, bool par ){
	bool ok = true ;

	int nbit=(int)iters;

	while (ok &&  nbit){
		typedef typename Field::Element Element ;
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

		if (nbw<0)
			nbw = (int) random() % 7;
#ifdef DEBUG
		F->write(std::cerr) << std::endl;
#endif
		typedef typename Field::RandIter Randiter ;
		typedef typename Field::Element  Element ;
		Randiter R1(*F,b);
		Givaro::GeneralRingNonZeroRandIter<Field,Randiter> R(R1);

			//size_t k = 0 ;
			//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(*F,m,n,k,F->one,F->zero,iters,nbw, par, b);
			//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(*F,m,n,k,F->zero,F->zero,iters,nbw, par, b);
			//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(*F,m,n,k,F->mOne,F->zero,iters,nbw, par, b);
			//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(*F,m,n,k,F->one ,F->one,iters,nbw, par, b);
			//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(*F,m,n,k,F->zero,F->one,iters,nbw, par, b);
			//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(*F,m,n,k,F->mOne,F->one,iters,nbw, par, b);
			//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(*F,m,n,k,F->one ,F->mOne,iters,nbw, par, b);
			//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(*F,m,n,k,F->zero,F->mOne,iters,nbw, par, b);
		 	//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(*F,m,n,k,F->mOne,F->mOne,iters,nbw, par, b);
			//std::cout << k << "/24" << std::endl; ++k;

		 Element alpha,beta ;
		 R.random(alpha);
		ok &= launch_MM_dispatch<Field>(*F,m,n,k,F->one ,alpha,iters,nbw, par, b);
		//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(*F,m,n,k,F->zero,alpha,iters,nbw, par, b);
		//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(*F,m,n,k,F->mOne,alpha,iters,nbw, par, b);
		//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(*F,m,n,k,alpha,F->one ,iters,nbw, par, b);
		//std::cout << k << "/24" << std::endl; ++k;
		 ok &= launch_MM_dispatch<Field>(*F,m,n,k,alpha,F->zero,iters,nbw, par, b);
		//std::cout << k << "/24" << std::endl; ++k;
		ok &= launch_MM_dispatch<Field>(*F,m,n,k,alpha,F->mOne,iters,nbw, par, b);
		//std::cout << k << "/24" << std::endl; ++k;

		for (size_t j = 0 ; j < 3 ; ++j) {
			R.random(alpha);
			R.random(beta);
			ok &= launch_MM_dispatch<Field>(*F,m,n,k,alpha,beta,iters,nbw, par, b);
			//std::cout << k << "/24" << std::endl; ++k;
		}
			//std::cout<<std::endl;
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
	std::cout<<setprecision(17);
	std::cerr<<setprecision(17);
	srand((int)time(NULL));
	srand48(time(NULL));

	static size_t iters = 3 ;
	static Givaro::Integer q = -1 ;
	static uint64_t b = 0 ;
	static int m = -50 ;
	static int n = -50 ;
	static int k = -50 ;
	static int nbw = -1 ;
	static bool loop = false;
	static bool p = false;
	static Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
		{ 'b', "-b B", "Set the bitsize of the random characteristic.",         TYPE_INT , &b },
		{ 'm', "-m M", "Set the dimension of the matrix (negative values, mean, any random value between 0 and |n|).",      TYPE_INT , &m },
		{ 'n', "-n N", "Set the dimension of the matrix (negative values, mean, any random value between 0 and |n|).",      TYPE_INT , &n },
		{ 'k', "-k K", "Set the dimension of the matrix (negative values, mean, any random value between 0 and |k|).",      TYPE_INT , &k },
		{ 'w', "-w N", "Set the number of winograd levels (-1 for random).",    TYPE_INT , &nbw },
		{ 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
		{ 'l', "-l Y/N", "run the test in an infinte loop.", TYPE_BOOL , &loop },
		{ 'p', "-p Y/N", "run the parallel fgemm.", TYPE_BOOL , &p },
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(argc,argv,as);

	bool ok = true;
	do{
		ok &= run_with_field<Modular<double> >(q,b,m,n,k,nbw,iters,p);
		ok &= run_with_field<ModularBalanced<double> >(q,b,m,n,k,nbw,iters,p);
		ok &= run_with_field<Modular<float> >(q,b,m,n,k,nbw,iters,p);
		ok &= run_with_field<ModularBalanced<float> >(q,b,m,n,k,nbw,iters,p);
		ok &= run_with_field<Modular<int32_t> >(q,b,m,n,k,nbw,iters,p);
		ok &= run_with_field<ModularBalanced<int32_t> >(q,b,m,n,k,nbw,iters,p);
		ok &= run_with_field<Modular<RecInt::rint<7> > >(q,b?b:63_ui64,m,n,k,nbw,iters, p);
		ok &= run_with_field<Modular<RecInt::rint<8> > >(q,b?b:127_ui64,m,n,k,nbw,iters, p);
		ok &= run_with_field<Modular<int64_t> >(q,b,m,n,k,nbw,iters, p);
		ok &= run_with_field<ModularBalanced<int64_t> >(q,b,m,n,k,nbw,iters, p);
		ok &= run_with_field<Modular<Givaro::Integer> >(q,(b?b:512_ui64),m,n,k,nbw,iters,p);
		ok &= run_with_field<Givaro::ZRing<Givaro::Integer> >(0,(b?b:512_ui64),m,n,k,nbw,iters,p);

	} while (loop && ok);

	return !ok ;
}
