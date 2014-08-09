/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 *
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
/** @file fflas/fflas_ftrsm_mp.inl
 * @brief triangular system with matrix right hand side over multiprecision domain (either over Z or over Z/pZ)
 */
#ifndef __FFPACK_ftrsm_mp_INL
#define __FFPACK_ftrsm_mp_INL

#include <math.h>
#include "fflas-ffpack/field/rns-integer-mod.h"
#include "fflas-ffpack/field/modular-integer.h"
#include "fflas-ffpack/field/integer.h"
#include "fflas-ffpack/fflas/fflas.h"

#ifdef __FFLASFFPACK_HAVE_INTEGER

namespace FFLAS {
 
	
	void ftrsm (const FFPACK::Modular<FFPACK::integer> & F, 
		    const FFLAS_SIDE Side,
		    const FFLAS_UPLO Uplo,
		    const FFLAS_TRANSPOSE TransA,
		    const FFLAS_DIAG Diag, 
		    const size_t M, const size_t N,
		    const FFPACK::Integer alpha,
		    FFPACK::Integer * A, const size_t lda,
		    FFPACK::Integer * B, const size_t ldb){
		
#ifdef BENCH_PERF_TRSM
		double t_init=0, t_trsm=0, t_mod=0, t_rec=0;
		Timer chrono;
		chrono.start(); 
#endif
		FFPACK::Integer p;
		F.cardinality(p);
		size_t logp=p.bitsize();
		size_t K;
		if (Side == FFLAS::FflasLeft)		 	
			K=M;
		else
			K=N;			
		
		// compute bit size of feasible prime 
		size_t _k=max(K,logp/20), lk=0;
		while ( _k ) {_k>>=1; ++lk;}    
		size_t prime_bitsize= (53-lk)>>1;	
	
		// construct rns basis
		FFPACK::Integer maxC= (p-1)*(p-1)*(p-1)*K; 
		size_t n_pr =maxC.bitsize()/prime_bitsize;				
		maxC=(p-1)*(p-1)*K*(1<<prime_bitsize)*n_pr; 		

		FFPACK::rns_double RNS(maxC, prime_bitsize, true);
		FFPACK::RNSIntegerMod<FFPACK::rns_double> Zp(p, RNS);			
#ifdef BENCH_PERF_TRSM
		chrono.stop();
		t_init+=chrono.usertime();
		chrono.clear();chrono.start();
#endif	
		// compute A and B in RNS
		FFPACK::rns_double::Element_ptr Ap,Bp;		
		Ap._ptr = new double[K*K*RNS._size];
		Ap._stride = K*K;
		Bp._ptr = new double[M*N*RNS._size];
		Bp._stride = M*N;
		finit(Zp,K,K,(logp/16)+(logp%16?1:0),A,lda,Ap); 
		finit(Zp,M,N,(logp/16)+(logp%16?1:0),B,ldb,Bp);		
#ifdef BENCH_PERF_TRSM
		chrono.stop();
		t_mod+=chrono.usertime();
		chrono.clear();chrono.start();
#endif		
		// call ftrsm in rns		
		ftrsm(Zp, Side, Uplo, TransA, Diag, M, N, Zp.one, Ap, K, Bp, N); 
		
#ifdef BENCH_PERF_TRSM
		chrono.stop();
		t_trsm+=chrono.usertime();
		chrono.clear();chrono.start();
#endif			
		// reconstruct the result 
		fconvert(Zp,M,N,F.zero,B,ldb,Bp);	
		// reduce it modulo p 		
		finit(F,M,N,B,ldb);
		// scale it with alpha
		if (!F.isOne(alpha))
			fscalin(F,M,N,alpha,B,ldb);
		
#ifdef BENCH_PERF_TRSM
		chrono.stop();
		t_rec+=chrono.usertime();		
		cout<<"FTRSM RNS PERF:"<<endl;
		cout<<"  ***      init  : "<<t_init<<endl;
		cout<<"  ***  rns  mod  : "<<t_mod<<endl;
		cout<<"  ***  rns trsm  : "<<t_trsm<<endl;
		cout<<"  ***  rns  rec  : "<<t_rec<<endl;
#endif	

		delete[] Ap._ptr;
		delete[] Bp._ptr;
	} 


} // END OF NAMESPACE FFLAS

#endif
#endif
