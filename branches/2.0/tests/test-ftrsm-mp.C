/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

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


#include <iomanip>
#include <iostream>
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/field/modular-integer.h"

using namespace std;


typedef FFPACK::Modular<FFPACK::integer> Field;
//typedef FFPACK::Modular<double> Field;


template<typename T>
void write_matrix(FFPACK::Integer p, size_t m, size_t n, T* C, size_t ldc){

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
void check_ftrsm (const Field &F, size_t m, size_t n, const typename Field::Element &alpha, FFLAS::FFLAS_SIDE side, FFLAS::FFLAS_UPLO uplo, FFLAS::FFLAS_TRANSPOSE trans, FFLAS::FFLAS_DIAG diag){

	typedef typename Field::Element Element;
	Element * A, *B, *B2, *C, tmp;
	size_t k = (side==FFLAS::FflasLeft?m:n);
	A  = new Element[k*k];
	B  = new Element[m*n];
	B2 = new Element[m*n];
	C  = new Element[m*n]; 
	
	typename Field::RandIter Rand(F);
	
	for (size_t i=0;i<k;++i){
		for (size_t j=0;j<i;++j) 
			A[i*k+j]= (uplo == FFLAS::FflasLower)? Rand.random(tmp) : F.zero;
		A[i*k+i]= (diag == FFLAS::FflasNonUnit)? Rand.random(tmp) : F.one;
		for (size_t j=i+1;j<k;++j) 
			A[i*k+j]= (uplo == FFLAS::FflasUpper)? Rand.random(tmp) : F.zero;
	}
	for (size_t i=0;i<m;++i){
		for(int j=0; j<n; ++j){
			B[i*n+j]= Rand.random(tmp);
			B2[i*n+j]=B[i*n+j];
		}  
	}
	
	string ss=string((uplo == FFLAS::FflasLower)?"Lower_":"Upper_")+string((side == FFLAS::FflasLeft)?"Left_":"Right_")+string((trans == FFLAS::FflasTrans)?"Trans_":"NoTrans_")+string((diag == FFLAS::FflasUnit)?"Unit":"NonUnit");

	cerr<<std::left<<"Checking FTRSM_";
	cerr.fill('.');
	cerr.width(35);
	cerr<<ss;

	::Timer t; t.clear();
	double time=0.0;	
	t.clear();
	t.start();
	FFLAS::ftrsm (F, side, uplo, trans, diag, m, n, alpha, A, k, B, n);
	t.stop();
	time+=t.usertime();
	
	
	Element invalpha;
	F.init(invalpha);
	F.inv(invalpha, alpha);	
		
	FFLAS::ftrmm (F, side, uplo, trans, diag, m, n, invalpha, A, k, B, n);
	
	// if (side == FFLAS::FflasLeft)
	// 	FFLAS::fgemm(F, trans, FFLAS::FflasNoTrans, m, n, m, invalpha, A, k, B, n, F.zero, C, n);
	// else
	// 	FFLAS::fgemm(F, FFLAS::FflasNoTrans, trans, m, n, n, invalpha, B, n, A, k, F.zero, C, n);


	F.mulin(invalpha,alpha);
	if (!F.isOne(invalpha)){
		cout<<"invalpha is wrong !!!";
	}
	bool wrong = false;

	for (int i=0;i<m;++i)
		for (int j=0;j<n;++j)
			if ( !F.areEqual(*(B2+i*n+j), *(B+i*n+j))){
				wrong = true;
			}
	
	if ( wrong ){
		cerr<<"FAILED ("<<time<<")"<<endl;
		
	} else
		cerr<<"PASSED ("<<time<<")"<<endl;
	
	delete[] A;
	delete[] B;
	delete[] B2;
	delete[] C;
	
}


int main(int argc, char** argv)
{
	cerr<<setprecision(10);
	
	if (argc != 5)	{ 
		cerr<<"Usage : test-ftrsm <pbits> <m> <n> <alpha>"
		    <<endl;
		exit(-1);
	}
	size_t b=atoi(argv[1]);
	size_t m=atoi(argv[2]);
	size_t n=atoi(argv[3]);

	FFPACK::Integer p;
	FFPACK::Integer::random_exact_2exp(p, b);			
	nextprime(p,p);
	Field F((Field::Element)p);
	cout<<F<<endl;
	Field::Element alpha;
	F.init (alpha, (Field::Element)FFPACK::Integer(argv[4])); 

	
	check_ftrsm(F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit);
	check_ftrsm(F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit);
	check_ftrsm(F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit);
	check_ftrsm(F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit);

	check_ftrsm(F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit);
	check_ftrsm(F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit);	
	check_ftrsm(F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit);
	check_ftrsm(F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit);

	check_ftrsm(F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit);
	check_ftrsm(F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit);
	check_ftrsm(F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit);
	check_ftrsm(F,m,n,alpha,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit);

	check_ftrsm(F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit);
	check_ftrsm(F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit);		
	check_ftrsm(F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit);
	check_ftrsm(F,m,n,alpha,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit);
	
	return 0;
}
