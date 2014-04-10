/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
//--------------------------------------------------------------------------
//                        Test for the  fgemm winograd
//
//--------------------------------------------------------------------------

/*
 * Copyright (c) FFLAS-FFPACK
 * Written by Clement Pernet <clement.pernet@imag.fr>
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
 */

#define NEWWINO

#include <iostream>
#include <iomanip>
using namespace std;
//#include "fflas-ffpack/fieldmodular-int.h"
//#include "fflas-ffpack/fieldmodular-positive.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "utils/timer.h"
#include "Matio.h"
#include "fflas-ffpack/fflas/fflas.h"

#ifndef __FFLASFFPACK_HAVE_GIVARO
#error you need givaro (and gmp) here
#endif


#include "givaro/givintprime.h"



using namespace FFPACK;

//typedef ModularBalanced<float> Field;
typedef ModularBalanced<double> Field;
//typedef Modular<double> Field;
//typedef Modular<float> Field;
//typedef Modular<int> Field;
//typedef GivaroZpz<Std32> Field;
//typedef GivaroGfq Field;

int main(int argc, char** argv){
	Timer tim;
	Givaro::IntPrimeDom IPD;
	Field::Element alpha, beta;
	long p;
	size_t M, K, N, Wino;
	bool keepon = true;
	Givaro::Integer _p,tmp;
	cerr<<setprecision(10);
	size_t TMAXM = 100, TMAXK = 100, TMAXN = 100;
	size_t PRIMESIZE = 23;
	size_t WINOMAX = 8;

	if (argc > 1 )
		TMAXM = atoi(argv[1]);
	if (argc > 2 )
		PRIMESIZE = atoi(argv[2]);
	if (argc > 3 )
		WINOMAX = atoi(argv[3]);
	if (argc > 4 )
		TMAXK = atoi(argv[4]);
    else
        TMAXK = TMAXM;
	if (argc > 5 )
		TMAXN = atoi(argv[5]);
    else
        TMAXN = TMAXM;


	enum FFLAS::FFLAS_TRANSPOSE ta, tb;
	size_t lda,ldb;
	Field::Element * A;
	Field::Element * B;
	Field::Element * C, *Cbis, *Cter;

	while (keepon){
		srandom(_p);
		do{
			//		max = Integer::random(2);
			_p = random();//max % (2<<30);
			IPD.prevprime( tmp, (_p% (1<<PRIMESIZE)) );
			p =  tmp;

		}while( (p <= 2) );

		Field F( p );
		Field::RandIter RValue( F );
		//NonzeroRandIter<Field> RnValue( F, RValue );


		do{
			M = (size_t)  random() % TMAXM;
			K = (size_t)  random() % TMAXK;
			N = (size_t)  random() % TMAXN;
			Wino = random() % WINOMAX;
		} while (!( (K>>Wino > 0) && (M>>Wino > 0) && (N>>Wino > 0) ));

		if (random()%2){
			ta = FFLAS::FflasTrans;
			lda = M;
		}
		else{
			ta = FFLAS::FflasNoTrans;
			lda = K;
		}
		if (random()%2){
			tb = FFLAS::FflasTrans;
			ldb = K;
		}
		else{
			tb = FFLAS::FflasNoTrans;
			ldb = N;
		}

		A = new Field::Element[M*K];
		B = new Field::Element[K*N];
		C = new Field::Element[M*N];
		Cbis = new Field::Element[M*N];
		Cter = new Field::Element[M*N];

		for( size_t i = 0; i < M*K; ++i )
			RValue.random( *(A+i) );
		for( size_t i = 0; i < K*N; ++i )
			RValue.random( *(B+i) );
		for( size_t i = 0; i < M*N; ++i )
			*(Cter+i) = *(Cbis+i)= RValue.random( *(C+i) );

		RValue.random( alpha );
		RValue.random( beta );

		cout <<"p = "<<(size_t)p<<" M = "<<M
		     <<" K = "<<K
             <<" N = "<<N
             <<" Winolevel = "<<Wino<<" "
		     <<alpha
		     <<((ta==FFLAS::FflasNoTrans)?".Ax":".A^Tx")
		     <<((tb==FFLAS::FflasNoTrans)?"B + ":"B^T + ")
		     <<beta<<".C"
		     <<"....";

		tim.clear();
		tim.start();
		FFLAS::fgemm (F, ta, tb, M, N, K, alpha, A, lda, B, ldb, beta, C, N, Wino);
		tim.stop();
// 		for (int j = 0; j < n; ++j ){
// 			FFLAS::fgemv( F, FFLAS::FflasNoTrans, m, k, alpha, A, k, B+j, n, beta, Cbis+j, n);
// 			for (int i=0; i<m; ++i)
// 				if ( !F.areEqual( *(Cbis+i*n+j), *(C+i*n+j) ) )
// 					keepon = false;
// 		}
		Field::Element aij, bij, temp;
		// Field::Element boa ;
		//F.div(boa, beta, alpha);
		for (size_t i = 0; i < M; ++i )
			for ( size_t j = 0; j < N; ++j ){
				//				F.mulin(*(Cbis+i*N+j),boa);
				F.mulin(*(Cbis+i*N+j),beta);
				for ( size_t l = 0; l < K ; ++l ){
					if ( ta == FFLAS::FflasNoTrans )
						aij = *(A+i*lda+l);
					else
						aij = *(A+l*lda+i);
					if ( tb == FFLAS::FflasNoTrans )
						bij = *(B+l*ldb+j);
					else
						bij = *(B+j*ldb+l);
					F.mul(temp,aij,bij);
					F.axpyin( *(Cbis+i*N+j), alpha, temp);
					//F.axpyin( *(Cbis+i*N+j), aij, bij );
				}
				//F.mulin( *(Cbis+i*N+j),alpha );
				if ( !F.areEqual( *(Cbis+i*N+j), *(C+i*N+j) ) ) {
					cerr<<"error for i,j="<<i<<" "<<j<<" "<<*(C+i*N+j)<<" "<<*(Cbis+i*N+j)<<endl;
					keepon = false;
				}
			}

		if (keepon){
			cout<<"Passed "
			    <<(2*double(M*N)/1000.0*double(K)/tim.usertime()/1000.0)<<"Mfops"<<endl;
			delete[] A;
			delete[] B;
			delete[] C;
			delete[] Cbis;
			delete[] Cter;
		}
		else{
			// cerr<<"C="<<endl;
// 			write_field( F, cerr, C, M, N, N );
// 			cerr<<"Cbis="<<endl;
// 			write_field( F, cerr, Cbis, M, N, N );
		}
	}
	cout<<endl;
	cerr<<"FAILED with p = "<<(size_t)p<<" M = "<<M<<" N = "<<N<<" K = "<<K
	    <<" Winolevel = "<<Wino
	    <<" alpha = "<<(size_t)alpha<<" beta = "<<(size_t)beta<<endl;
	cerr<<"A:"<<endl;
	if ( ta ==FFLAS::FflasNoTrans ){
		cerr<<M<<" "<<K<<" M"<<endl;
		for (size_t i=0; i<M; ++i)
			for (size_t j=0; j<K; ++j)
				cerr<<i+1<<" "<<j+1<<" "<<((size_t) *(A+i*lda+j) )<<endl;
	}
	else{
		cerr<<K<<" "<<M<<" M"<<endl;
		for (size_t i=0; i<K; ++i)
			for (size_t j=0; j<M; ++j)
				cerr<<i+1<<" "<<j+1<<" "<<((size_t) *(A+j*lda+i) )<<endl;

	}
	cerr<<"0 0 0"<<endl<<endl;
	cerr<<"B:"<<endl;
	if ( tb ==FFLAS::FflasNoTrans ){
		cerr<<K<<" "<<N<<" M"<<endl;
		for (size_t i=0; i<K; ++i)
			for (size_t j=0; j<N; ++j)
				cerr<<i+1<<" "<<j+1<<" "<<((size_t) *(B+i*ldb+j) )<<endl;
	}
	else{
		cerr<<N<<" "<<K<<" M"<<endl;
		for (size_t i=0; i<N; ++i)
			for (size_t j=0; j<K; ++j)
				cerr<<i+1<<" "<<j+1<<" "<<((size_t) *(B+i+j*ldb) )<<endl;
	}
	cerr<<"0 0 0"<<endl<<endl;
	cerr<<"C:"<<endl
	    <<M<<" "<<N<<" M"<<endl;
	for (size_t i=0; i<M; ++i)
		for (size_t j=0; j<N; ++j)
			cerr<<i+1<<" "<<j+1<<" "<<((size_t) *(Cter+i*N+j) )<<endl;
	cerr<<"0 0 0"<<endl;

	delete[] A;
	delete[] B;
	delete[] C;
	delete[] Cbis;
	delete[] Cter;
}














