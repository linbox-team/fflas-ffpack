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

#ifdef __FFLASFFPACK_HAVE_INTEGER

#include <math.h>

#include "fflas-ffpack/fflas/fflas_bounds.inl"
#include "fflas-ffpack/fflas/fflas_level3.inl"
#include "fflas-ffpack/field/rns-integer-mod.h"
#include "fflas-ffpack/field/rns-integer.h"
#include "fflas-ffpack/field/modular-integer.h"
#include "fflas-ffpack/field/integer.h"

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

#ifdef BENCH_PERF_TRSM_MP
		double t_init=0, t_trsm=0, t_mod=0, t_rec=0;
		FFLAS::Timer chrono;
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
#ifdef BENCH_PERF_TRSM_MP
		chrono.stop();
		t_init+=chrono.usertime();
		chrono.clear();chrono.start();
#endif
		// compute A and B in RNS
		FFPACK::rns_double::Element_ptr Ap,Bp;
		Ap = FFLAS::fflas_new(Zp,K,K);
		Bp = FFLAS::fflas_new(Zp,M,N);

		if (Side == FFLAS::FflasLeft){
			finit_rns(Zp,K,K,(logp/16)+(logp%16?1:0),A,lda,Ap);
			finit_rns(Zp,M,N,(logp/16)+(logp%16?1:0),B,ldb,Bp);
		}
		else {
			finit_trans_rns(Zp,K,K,(logp/16)+(logp%16?1:0),A,lda,Ap);
			finit_trans_rns(Zp,M,N,(logp/16)+(logp%16?1:0),B,ldb,Bp);					
		}		
#ifdef BENCH_PERF_TRSM_MP
		chrono.stop();
		t_mod+=chrono.usertime();
		chrono.clear();chrono.start();
#endif
		// call ftrsm in rns
		//ftrsm(Zp, Side, Uplo, TransA, Diag, M, N, Zp.one, Ap, K, Bp, N);
		if (Side == FFLAS::FflasLeft)
			ftrsm(Zp, Side, Uplo, TransA, Diag, M, N, Zp.one, Ap, K, Bp, N);
		else {
			if (Uplo == FFLAS::FflasUpper)
				ftrsm(Zp, FFLAS::FflasLeft, FFLAS::FflasLower, TransA, Diag, N, M, Zp.one, Ap, K, Bp, M);
			else
				ftrsm(Zp, FFLAS::FflasLeft, FFLAS::FflasUpper, TransA, Diag, N, M, Zp.one, Ap, K, Bp, M);
		}

#ifdef BENCH_PERF_TRSM_MP
		chrono.stop();
		t_trsm+=chrono.usertime();
		chrono.clear();chrono.start();
#endif
		// reconstruct the result
		if (Side == FFLAS::FflasLeft)
			fconvert_rns(Zp,M,N,F.zero,B,ldb,Bp);
		else{
			fconvert_trans_rns(Zp,M,N,F.zero,B,ldb,Bp);
		}	
	
		// reduce it modulo p
		finit(F,M,N,B,ldb);
		// scale it with alpha
		if (!F.isOne(alpha))
			fscalin(F,M,N,alpha,B,ldb);

#ifdef BENCH_PERF_TRSM_MP
		chrono.stop();
		t_rec+=chrono.usertime();
		cout<<"FTRSM RNS PERF:"<<endl;
		cout<<"  ***      init  : "<<t_init<<endl;
		cout<<"  ***  rns  mod  : "<<t_mod<<endl;
		cout<<"  ***  rns trsm  : "<<t_trsm<<" ( igemm="<<Zp.t_igemm<<" scal="<<Zp.t_scal<<" modp="<<Zp.t_modp<<endl;;
		cout<<"  ***  rns  rec  : "<<t_rec<<endl;
#endif

		FFLAS::fflas_delete(Ap);
		FFLAS::fflas_delete(Bp);
	}


	inline void cblas_imptrsm(const enum CBLAS_ORDER Order, 
				  const enum CBLAS_SIDE Side, 
				  const enum CBLAS_UPLO Uplo, 
				  const enum CBLAS_TRANSPOSE TransA,
				  const enum CBLAS_DIAG Diag, 
				  const int M, const int N, const FFPACK::rns_double_elt alpha, 
				  const FFPACK::rns_double_elt_ptr A, const int lda,
				  FFPACK::rns_double_elt_ptr B, const int ldb) {}

#ifndef DOXYGEN_SHOULD_SKIP_THIS	
	namespace Protected {

		template<>
		inline size_t TRSMBound (const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F)
		{
			return 1;
		}
	
		template <>
		inline size_t DotProdBoundClassic (const FFPACK::RNSIntegerMod<FFPACK::rns_double>& F,
						   const FFPACK::rns_double_elt& beta)
		{
			FFPACK::Integer p,b,M;
			F.cardinality(p);
			p--;
			F.convert(b,beta);
			M=F.rns()._M;
			size_t kmax= (M-b*p)/(p*p);
			return kmax;
		}
		
#ifndef __FTRSM_MP_FAST
#define __FFLAS_MULTIPRECISION

#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT



#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT
#endif // #ifdef __FTRSM_MP_FAST	
		
	} // end of namespace protected
#endif // #ifndef DOXYGEN_SHOULD_SKIP_THIS	
} // END OF NAMESPACE FFLAS
		
#endif
#endif
