/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 *
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



#include <iostream>
#include <vector>
#include <string>
using namespace std; 

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"


#ifdef	BENCH_FLINT
#define __GMP_BITS_PER_MP_LIMB 64
extern "C" {
#include "flint/longlong.h"
#include "flint/long_extras.h" 
#include "flint/fmpz_mat.h"
#include "flint/fmpz.h"
#include "flint/flint.h"
} 
 

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


string check_res(size_t m, size_t n, FFPACK::Integer* A, size_t lda, const fmpz_mat_t & B){
	fmpz_t tmp; fmpz_init(tmp);
	bool allz=true;
	for (size_t i=0;i<m;i++) 
		for (size_t j=0;j<n;j++){
			fmpz_set_mpz(tmp, *(reinterpret_cast<const mpz_t*>(A+i*lda+j)));
			if (fmpz_cmp(fmpz_mat_entry(B,i,j),tmp) != 0)
				{
					return " (result is wrong)";	
				}
			if (A[i*lda+j]!=FFPACK::Integer(0))
				allz=false;
		}
	return (allz ? "(result: allzero) " :" (result is correct)");
}

#endif




int main(int argc, char** argv){
	if ((argc <2) || (argc > 6))
		std::cout<<"usage is : <entry bitsize> <m,k,n>  [seed]\n"; 

	size_t b = atoi(argv[1]);
	size_t m = atoi(argv[2]);
	size_t k=m,n=m;
	size_t seed= time(NULL);
	 
	if (argc==4)
		k=atoi(argv[3]);
	if (argc==5)
		n=atoi(argv[4]);
	if (argc==6)
		seed=atoi(argv[5]);
	 
	typedef FFPACK::Modular<FFPACK::integer> Field;
	
	FFPACK::Integer p;
	FFPACK::Integer::random_exact_2exp(p, b);			
	nextprime(p,p);
	
	Field F(p);
	size_t lda,ldb,ldc;
	lda=k;
	ldb=n;
	ldc=n;;

	typename Field::RandIter Rand(F,seed);
	FFPACK::Integer *A,*B,*C;
	C= new FFPACK::Integer[m*n];
	A= new FFPACK::Integer[m*k];
	B= new FFPACK::Integer[k*n];
	
	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<k;++j)
			Rand.random(A[i*k+j]);			
	for (size_t i=0;i<k;++i)
		for (size_t j=0;j<n;++j)
			Rand.random(B[i*n+j]);				
	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<n;++j)
			Rand.random(C[i*n+j]);			
	
	
	FFPACK::Integer alpha,beta;
	alpha=1;
	beta=0;
	

	Timer chrono;
#ifdef	BENCH_FLINT	
	// FLINT MUL //
	fmpz_t modp,tmp;
	fmpz_init(modp);
        fmpz_init(tmp);
        fmpz_set_mpz(modp, *(reinterpret_cast<const mpz_t*>(&p)));
	fmpz_mat_t AA,BB,CC,DD;
	fmpz_mat_init (AA, m, k); 
	fmpz_mat_init (BB, k, n); 
	fmpz_mat_init (CC, m, n); 
	fmpz_mat_init (DD, m, n);
	fmpz_t aalpha, bbeta;
	fmpz_set_mpz(aalpha,*(reinterpret_cast<const mpz_t*>(&alpha)));
	fmpz_set_mpz(bbeta,*(reinterpret_cast<const mpz_t*>(&beta)));

 	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<k;++j)
			fmpz_set_mpz(fmpz_mat_entry(AA,i,j),*(reinterpret_cast<const mpz_t*>(A+i*k+j)));
	for (size_t i=0;i<k;++i)
		for (size_t j=0;j<n;++j)
			fmpz_set_mpz(fmpz_mat_entry(BB,i,j),*(reinterpret_cast<const mpz_t*>(B+i*n+j)));
	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<n;++j)
			fmpz_set_mpz(fmpz_mat_entry(CC,i,j),*(reinterpret_cast<const mpz_t*>(C+i*n+j)));				
	chrono.clear();chrono.start();
	// DD= A.B
	fmpz_mat_mul(DD,AA,BB);	
	// CC = beta.C 
	fmpz_mat_scalar_mul_fmpz(CC,CC,bbeta);
	// CC = CC + DD.alpha
	fmpz_mat_scalar_addmul_fmpz(CC,DD,aalpha);
	// CC = CC mod p
	for (size_t i=0;i<m;++i)
                for (size_t j=0;j<n;++j)
                        fmpz_mod(fmpz_mat_entry(CC,i,j),fmpz_mat_entry(CC,i,j),modp);
	
	chrono.stop();
	fmpz_mat_clear(AA);
	fmpz_mat_clear(BB);
	cout<<" FLINT MUL: "<<chrono.usertime()<<endl;	 
#endif
	//END FLINT CODE //

	// RNS MUL_LA
	chrono.clear();chrono.start();	
	FFLAS::fgemm(F,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
	chrono.stop();
#ifdef BENCH_FLINT
	cout<<" RNS MUL LA: "<<chrono.usertime()<<check_res(m,n,C,n,CC)<<" "<<endl;	 
#else
	cout<<" RNS MUL LA: "<<chrono.usertime()<<endl;
#endif
	delete[] A;
	delete[] B;
	delete[] C;
	return 0;
}
 
