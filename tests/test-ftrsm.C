/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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
#define  __FFLASFFPACK_SEQUENTIAL

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <givaro/modular-integer.h>

#include <iomanip>
#include <iostream>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "test-utils.h"
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>


using namespace std;
using namespace FFPACK;
using Givaro::Modular;
using Givaro::ModularBalanced;

template<typename T>
void write_matrix(Givaro::Integer p, size_t m, size_t n, T* C, size_t ldc){

	size_t www=(p.bitsize()*log(2.))/log(10.);
	for (size_t i=0;i<m;++i){
		cout<<"[ ";
		cout.width(www+1);
		cout<<std::right<<C[i*ldc];
		for (size_t j=1;j<n;++j){
			cout<<" ";
			cout.width(www);
			cout<<std::right<<C[i*ldc+j];
		}
		cout<<"]"<<endl;
	}
	cout<<endl;

}


template<typename Field>
bool check_ftrsm (const Field &F, size_t m, size_t n, const typename Field::Element &alpha, FFLAS::FFLAS_SIDE side, FFLAS::FFLAS_UPLO uplo, FFLAS::FFLAS_TRANSPOSE trans, FFLAS::FFLAS_DIAG diag){

	typedef typename Field::Element Element;
	Element * A, *B, *B2, *C, tmp;
	size_t k = (side==FFLAS::FflasLeft?m:n);
	size_t lda,ldb,ldc;
	lda=k+13;
	ldb=n+14;
	ldc=n+15;
	A  = FFLAS::fflas_new(F,k,lda);
	B  = FFLAS::fflas_new(F,m,ldb);
	B2 = FFLAS::fflas_new(F,m,ldb);
	C  = FFLAS::fflas_new(F,m,ldc);

	typename Field::RandIter Rand(F);
	typename Field::NonZeroRandIter NZRand(Rand);

	for (size_t i=0;i<k;++i){
		for (size_t j=0;j<i;++j)
			A[i*lda+j]= (uplo == FFLAS::FflasLower)? Rand.random(tmp) : F.zero;
		A[i*lda+i]= (diag == FFLAS::FflasNonUnit)? NZRand.random(tmp) : F.one;
		for (size_t j=i+1;j<k;++j)
			A[i*lda+j]= (uplo == FFLAS::FflasUpper)? Rand.random(tmp) : F.zero;
	}
	for (size_t i=0;i<m;++i){
		for(size_t j=0; j<n; ++j){
			B[i*ldb+j]= Rand.random(tmp);
			B2[i*ldb+j]=B[i*ldb+j];
		}
	}

	string ss=string((uplo == FFLAS::FflasLower)?"Lower_":"Upper_")+string((side == FFLAS::FflasLeft)?"Left_":"Right_")+string((trans == FFLAS::FflasTrans)?"Trans_":"NoTrans_")+string((diag == FFLAS::FflasUnit)?"Unit":"NonUnit");

	cout<<std::left<<"Checking FTRSM_";
	cout.fill('.');
	cout.width(35);
	cout<<ss;


	FFLAS::Timer t; t.clear();
	double time=0.0;
	t.clear();
	t.start();
	FFLAS::ftrsm (F, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
	t.stop();
	time+=t.usertime();

	Element invalpha;
	F.init(invalpha);
	F.inv(invalpha, alpha);

	//FFLAS::ftrmm (F, side, uplo, trans, diag, m, n, invalpha, A, k, B, n);

	if (side == FFLAS::FflasLeft)
		FFLAS::fgemm(F, trans, FFLAS::FflasNoTrans, m, n, m, invalpha, A, lda, B, ldb, F.zero, C, ldc);
	else
		FFLAS::fgemm(F, FFLAS::FflasNoTrans, trans, m, n, n, invalpha, B, ldb, A, lda, F.zero, C, ldc);


	bool wrong = false;
	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<n;++j)
			if ( !F.areEqual(*(B2+i*ldb+j), *(C+i*ldc+j))){
				wrong = true;
			}
	if ( wrong ){
		    //cout << "\033[1;31mFAILED\033[0m ("<<time<<")"<<endl;
		cout << "FAILED ("<<time<<")"<<endl;
		//cerr<<"FAILED ("<<time<<")"<<endl;

	} else
		    //cout << "\033[1;32mPASSED\033[0m ("<<time<<")"<<endl;
		cout << "PASSED ("<<time<<")"<<endl;
	    //cerr<<"PASSED ("<<time<<")"<<endl;

	F.mulin(invalpha,alpha);
	if (!F.isOne(invalpha)){
		cerr<<"invalpha is wrong !!!"<<endl;;
	}

	FFLAS::fflas_delete(A);
	FFLAS::fflas_delete(B);
	FFLAS::fflas_delete(B2);
	FFLAS::fflas_delete(C);
	return !wrong;
}
template <class Field>
bool run_with_field (Givaro::Integer q, size_t b, size_t m, size_t n, int s, size_t iters){
	bool ok = true ;
	int nbit=(int)iters;

	while (ok &&  nbit){
		//typedef typename Field::Element Element ;
		// choose Field
		Field* F= chooseField<Field>(q,b);
		if (F==nullptr)
			return true;

		typename Field::Element alpha;
		F->init (alpha, (typename Field::Element)s);
		cout<<"Checking with ";F->write(cout)<<endl;

		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit);
		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit);
		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit);
		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit);
		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit);
		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit);
		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit);
		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit);
		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit);
		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit);
		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit);
		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit);
		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit);
		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit);
		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit);
		ok = ok && check_ftrsm(*F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit);
		nbit--;
		delete F;
	}
	return ok;
}

int main(int argc, char** argv)
{
	cerr<<setprecision(10);
	static Givaro::Integer q=-1;
	static size_t b=0;
	static size_t m=128;
	static size_t n=128;
	static size_t s=1;
	static size_t iters=1;
	static bool loop=false;
	static Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
                { 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
                { 'm', "-m M", "Set the row dimension of unknown matrix.",      TYPE_INT , &m },
                { 'n', "-n N", "Set the column dimension of the unknown matrix.", TYPE_INT , &n },
                { 's', "-s S", "Set the scaling of trsm",                         TYPE_INT , &s },
		{ 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
		{ 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
                END_OF_ARGUMENTS
        };

	FFLAS::parseArguments(argc,argv,as);

	bool ok = true;
	do{
	 	ok &= run_with_field<Modular<double> >(q,b,m,n,s,iters);
		ok &= run_with_field<ModularBalanced<double> >(q,b,m,n,s,iters);
		ok &= run_with_field<Modular<float> >(q,b,m,n,s,iters);
		ok &= run_with_field<ModularBalanced<float> >(q,b,m,n,s,iters);
		ok &= run_with_field<Modular<int32_t> >(q,b,m,n,s,iters);
		ok &= run_with_field<ModularBalanced<int32_t> >(q,b,m,n,s,iters);
		ok &= run_with_field<Modular<int64_t> >(q,b,m,n,s,iters);
		ok &= run_with_field<ModularBalanced<int64_t> >(q,b,m,n,s,iters);
		ok &= run_with_field<Modular<Givaro::Integer> >(q,(b?b:512),m/4+1,n/4+1,s,iters); 
	} while (loop && ok);

	return !ok ;
}
