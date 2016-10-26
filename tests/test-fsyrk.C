/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2016 FFLAS-FFPACK
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

#define ENABLE_ALL_CHECKINGS 1

#include "fflas-ffpack/fflas-ffpack-config.h"

#include <iomanip>
#include <iostream>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "test-utils.h"
#include <givaro/modular.h>


using namespace std;
using namespace FFLAS;
using Givaro::Modular;
using Givaro::ModularBalanced;

template<typename Field, class RandIter>
bool check_fsyrk (const Field &F, size_t n, size_t k,
				  const typename Field::Element &alpha, const typename Field::Element &beta,
				  FFLAS::FFLAS_UPLO uplo, FFLAS::FFLAS_TRANSPOSE trans, RandIter& Rand){

	typedef typename Field::Element Element;
	Element * A, *C, *C2;
	size_t lda = k+13;
	size_t ldc = n+15;
	A  = FFLAS::fflas_new(F,n,lda);
	C  = FFLAS::fflas_new(F,n,ldc);
	C2  = FFLAS::fflas_new(F,n,ldc);

	FFPACK::RandomTriangularMatrix (F, n, n, uplo, FflasNonUnit, true, C, ldc, Rand);
	FFPACK::RandomMatrix (F, n, k, A, lda, Rand);
	FFLAS::fassign (F, n, n, C, ldc, C2, ldc);

	string ss=string((uplo == FFLAS::FflasLower)?"Lower_":"Upper_")+string((trans == FFLAS::FflasTrans)?"Trans_":"NoTrans_");

	cout<<std::left<<"Checking FSYRK_";
	cout.fill('.');
	cout.width(35);
	cout<<ss;

	write_field(F, std::cerr<<"A = "<<std::endl,A, n,k,lda);
	write_field(F, std::cerr<<"C = "<<std::endl,C, n,n,ldc);
	FFLAS::Timer t; t.clear();
	double time=0.0;
	t.clear(); t.start();
   
	fsyrk (F, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);

	t.stop();
	time+=t.usertime();

	write_field(F, std::cerr<<"Apres fsyrk C = "<<std::endl,C, n,n,ldc);

	fgemm (F, trans, (trans==FflasNoTrans)?FflasTrans:FflasNoTrans, n, n, k, alpha, A, lda, A, lda, beta, C2, ldc);

	write_field(F, std::cerr<<"Apres fgemm C2 = "<<std::endl,C2, n,n,ldc);

	bool ok = true;
	if (uplo == FflasUpper){
		for (size_t i=0; i<n; i++)
			for (size_t j=i; j<n; j++)
				ok &= F.areEqual(C2[i*ldc+j], C[i*ldc+j]);
	} else {
		for (size_t i=0; i<n; i++)
			for (size_t j=0; j<=i; j++)
				ok &= F.areEqual(C2[i*ldc+j], C[i*ldc+j]);
	}
	if (ok)
	    //cout << "\033[1;32mPASSED\033[0m ("<<time<<")"<<endl;
		cout << "PASSED ("<<time<<")"<<endl;
		//cerr<<"PASSED ("<<time<<")"<<endl;
	else
	    //cout << "\033[1;31mFAILED\033[0m ("<<time<<")"<<endl;
		cout << "FAILED ("<<time<<")"<<endl;
		//cerr<<"FAILED ("<<time<<")"<<endl;
	

	FFLAS::fflas_delete(A);
	FFLAS::fflas_delete(C2);
	FFLAS::fflas_delete(C);
	return ok;
}

template <class Field>
bool run_with_field (Givaro::Integer q, size_t b, size_t m, size_t n, size_t a, size_t c, size_t iters, uint64_t seed){
	bool ok = true ;
	int nbit=(int)iters;

	while (ok &&  nbit){
		//typedef typename Field::Element Element ;
		// choose Field
		Field* F= FFPACK::chooseField<Field>(q,b);
		typename Field::RandIter G(*F,0,seed);
		if (F==nullptr)
			return true;

		typename Field::Element alpha, beta;
		F->init (alpha, (typename Field::Element)a);
		F->init (beta, (typename Field::Element)c);
		cout<<"Checking with ";F->write(cout)<<endl;

		ok = ok && check_fsyrk(*F,m,n,alpha,beta,FflasUpper,FflasNoTrans,G);
		// ok = ok && check_fsyrk(*F,m,n,alpha,FflasLower,FflasNoTrans,G);
		// ok = ok && check_fsyrk(*F,m,n,alpha,FflasLower,FflasTrans,G);
		// ok = ok && check_fsyrk(*F,m,n,alpha,FflasUpper,FflasTrans,G);
		nbit--;
		delete F;
	}
	return ok;
}

int main(int argc, char** argv)
{
	cerr<<setprecision(10);
	Givaro::Integer q=-1;
	size_t b=0;
	size_t m=128;
	size_t n=128;
	size_t a=1;
	size_t c=1;
	size_t iters=1;
	bool loop=false;
	uint64_t seed = time(NULL);
	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
		{ 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
		{ 'm', "-m M", "Set the row dimension",      TYPE_INT , &m },
		{ 'n', "-n N", "Set the column dimension.", TYPE_INT , &n },
		{ 'a', "-a A", "Set the scaling alpha",                         TYPE_INT , &a },
		{ 'c', "-c C", "Set the scaling beta",                         TYPE_INT , &c },
		{ 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
		{ 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
		{ 's', "-s seed", "Set seed for the random generator", TYPE_INT, &seed },
                END_OF_ARGUMENTS
        };

	parseArguments(argc,argv,as);

	bool ok = true;
	do{
		ok &= run_with_field<Modular<double> >(q,b,m,n,a,c,iters,seed);
		// ok &= run_with_field<ModularBalanced<double> >(q,b,m,n,a,c,iters,seed);
		// ok &= run_with_field<Modular<float> >(q,b,m,n,a,c,iters,seed);
		// ok &= run_with_field<ModularBalanced<float> >(q,b,m,n,a,c,iters,seed);
		// ok &= run_with_field<Modular<int32_t> >(q,b,m,n,a,c,iters,seed);
		// ok &= run_with_field<ModularBalanced<int32_t> >(q,b,m,n,a,c,iters,seed);
		// ok &= run_with_field<Modular<int64_t> >(q,b,m,n,a,c,iters,seed);
		// ok &= run_with_field<ModularBalanced<int64_t> >(q,b,m,n,a,c,iters,seed);
//		ok &= run_with_field<Modular<Givaro::Integer> >(q,5,m/4+1,n/4+1,a,c,iters,seed);
//		ok &= run_with_field<Modular<Givaro::Integer> >(q,(b?b:512),m/4+1,n/4+1,a,c,iters,seed);
	} while (loop && ok);

	return !ok ;
}
