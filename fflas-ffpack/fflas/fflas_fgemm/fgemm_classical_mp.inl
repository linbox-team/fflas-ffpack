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
/** @file fflas_fgemm/fgemm_classical_mp.inl
 * @brief matrix multiplication with multiprecision input (either over Z or over Z/pZ)
 */


#ifndef __FFPACK_fgemm_classical_INL
#define __FFPACK_fgemm_classical_INL

#include "fflas-ffpack/field/unparametric.h"
#include "fflas-ffpack/field/rns-integer.h"
#include "fflas-ffpack/field/modular-integer.h"
#include "fflas-ffpack/field/rns-integer-mod.h"

#ifdef __FFLASFFPACK_HAVE_INTEGER

namespace FFLAS {

	template<typename Field,
		 typename AlgoTrait,
		 typename ParSeqTrait>
	struct MMHelper<Field, AlgoTrait,FieldCategories::MultiPrecisionTag, ParSeqTrait> {
		FFPACK::Integer normA,normB;
		MMHelper() : normA(0), normB(0) {}
		MMHelper(FFPACK::Integer Amax, FFPACK::Integer Bmax) : normA(Amax), normB(Bmax) {}
		MMHelper(const Field& F, size_t m=0, size_t n=0, size_t k=0, ParSeqTrait PS=ParSeqTrait()) {F.characteristic(normA);F.characteristic(normB);}
	};


	template<typename RNS>
	struct FieldTraits<FFPACK::RNSInteger<RNS> > {typedef FieldCategories::MultiPrecisionTag value;};
	template<typename RNS>
	struct FieldTraits<FFPACK::RNSIntegerMod<RNS> > {typedef FieldCategories::MultiPrecisionTag value;};
	template<>
	struct FieldTraits<FFPACK::Modular<FFPACK::Integer> > {typedef FieldCategories::MultiPrecisionTag value;};
	template<>
	struct FieldTraits<FFPACK::UnparametricField<FFPACK::Integer> > {typedef FieldCategories::MultiPrecisionTag value;};


	/***********************************
	 *** MULTIPRECISION FGEMM OVER Z ***
	 ***********************************/

	// fgemm for RnsInteger with Winograd Helper
	template<typename RNS>
	inline  typename FFPACK::RNSInteger<RNS>::Element_ptr fgemm (const FFPACK::RNSInteger<RNS> &F,
								     const FFLAS_TRANSPOSE ta,
								     const FFLAS_TRANSPOSE tb,
								     const size_t m, const size_t n,const size_t k,
								     const typename FFPACK::RNSInteger<RNS>::Element alpha,
								     typename FFPACK::RNSInteger<RNS>::Element_ptr Ad, const size_t lda,
								     typename FFPACK::RNSInteger<RNS>::Element_ptr Bd, const size_t ldb,
								     const typename FFPACK::RNSInteger<RNS>::Element beta,
								     typename FFPACK::RNSInteger<RNS>::Element_ptr Cd, const size_t ldc,
								     MMHelper<FFPACK::RNSInteger<RNS>, MMHelperAlgo::Winograd, FieldCategories::MultiPrecisionTag, ParSeqHelper::Sequential> & H)
	{
		// compute each fgemm componentwise
		for(size_t i=0;i<F.size();i++){
			FFLAS::fgemm(F.rns()._field_rns[i],ta,tb,//FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
				     m, n, k, alpha._ptr[i], Ad._ptr+i*Ad._stride, lda, Bd._ptr+i*Bd._stride, ldb, beta._ptr[i], Cd._ptr+i*Cd._stride, ldc);
		}
		return Cd;
	}


	// fgemm for UnparametricField<Integer> with Winograd Helper
	inline FFPACK::Integer* fgemm (const FFPACK::UnparametricField<FFPACK::Integer>& F,
				       const FFLAS_TRANSPOSE ta,
				       const FFLAS_TRANSPOSE tb,
				       const size_t m, const size_t n,const size_t k,
				       const FFPACK::Integer alpha,
				       FFPACK::Integer* A, const size_t lda,
				       FFPACK::Integer* B, const size_t ldb,
				       FFPACK::Integer beta,
				       FFPACK::Integer* C, const size_t ldc,
				       MMHelper<FFPACK::UnparametricField<FFPACK::Integer>, MMHelperAlgo::Winograd, FieldCategories::MultiPrecisionTag,ParSeqHelper::Sequential> & H)
	{
		// compute bit size of feasible prime for FFLAS
		size_t _k=k,lk=0;
		while ( _k ) {_k>>=1; ++lk;}
		size_t prime_bitsize= (53-lk)>>1;

		// compute bound on the output
		FFPACK::Integer  mA,mB,mC;
		size_t logA,logB;
		mA=H.normA;
		mB=H.normB;
		if (H.normA==0){
			logA=A[0].bitsize();
			for (size_t i=1;i<m*k;i++)
				logA = max(logA,A[i].bitsize());
			H.normA=1; H.normA<<=(logA);
		}
		else {
			logA=H.normA.bitsize();
		}
		if (H.normB==0){
			logB=B[0].bitsize();
			for (size_t i=1;i<k*n;i++)
				logB=max(logB,B[i].bitsize());
			H.normB=1; H.normB<<=(logB);
		}
		else {
			logB=H.normA.bitsize();
		}
		mC = 2*k*H.normA*H.normB*abs(alpha); // need to use 2x bound to reach both positive and negative

		// construct an RNS structure and its associated Domain
		FFPACK::rns_double RNS(mC, prime_bitsize);
		typedef FFPACK::RNSInteger<FFPACK::rns_double> RnsDomain;
		RnsDomain Zrns(RNS);

		size_t Arowd,Acold,Browd,Bcold;
		if (ta == FFLAS::FflasNoTrans){ Arowd = m; Acold = k; }
		else {Arowd = k; Acold = m;}
		if (tb == FFLAS::FflasNoTrans){ Browd = k; Bcold = n; }
		else {Browd = n; Bcold = k;}
		    // allocate data for RNS representation
		typename RnsDomain::Element_ptr Ap,Bp,Cp;
		Ap = FFLAS::fflas_new(Zrns,Arowd,Acold);
		Bp = FFLAS::fflas_new(Zrns,Browd,Bcold);
		Cp = FFLAS::fflas_new(Zrns,m,n);

		//Ap._ptr = new double[m*k*RNS._size];
		//Ap._stride = m*k;
		// Bp._ptr = new double[k*n*RNS._size];
		// Bp._stride = k*n;
		// Cp._ptr = new double[m*n*RNS._size];
		// Cp._stride = m*n;

		// convert the input matrices to RNS representation
		finit(Zrns,Arowd,Acold,(logA/16)+((logA%16)?1:0),A,lda,Ap);
		finit(Zrns,Browd,Bcold,(logB/16)+((logB%16)?1:0),B,ldb,Bp);

		// perform the fgemm in RNS
		MMHelper<RnsDomain, MMHelperAlgo::Winograd> H2;// H2(Zrns,0,H.parseq);

		// compute alpha and beta in RNS
		typename RnsDomain::Element alphap, betap;
		Zrns.init(alphap, alpha);
		Zrns.init(betap, F.zero);

		// call  fgemm
		fgemm(Zrns,ta,tb,m,n,k,alphap,Ap,Acold,Bp,Bcold,betap,Cp,n,H2);

		// convert the RNS output to integer representation (C=beta.C+ RNS^(-1)(Cp) )
		fconvert(Zrns,m,n,beta,C,ldc,Cp);

		FFLAS::fflas_delete(Ap);
		FFLAS::fflas_delete(Bp);
		FFLAS::fflas_delete(Cp);

		return C;
	}


	/************************************
	 *** MULTIPRECISION FGEMM OVER Fp ***
	 ************************************/

	// fgemm for RNSIntegerMod  with Winograd Helper
	template<typename RNS>
	inline typename RNS::Element_ptr fgemm (const FFPACK::RNSIntegerMod<RNS> &F,
						 const FFLAS_TRANSPOSE ta,
						 const FFLAS_TRANSPOSE tb,
						 const size_t m, const size_t n,const size_t k,
						 const typename RNS::Element alpha,
						 typename RNS::Element_ptr Ad, const size_t lda,
						 typename RNS::Element_ptr Bd, const size_t ldb,
						 const typename RNS::Element beta,
						 typename RNS::Element_ptr Cd, const size_t ldc,
						 MMHelper<FFPACK::RNSIntegerMod<RNS>, MMHelperAlgo::Winograd, FieldCategories::MultiPrecisionTag, ParSeqHelper::Sequential> & H)
	{
		// compute the product over Z
		typedef FFPACK::RNSInteger<RNS> RnsDomain;
		MMHelper<RnsDomain, MMHelperAlgo::Winograd, FieldCategories::MultiPrecisionTag,ParSeqHelper::Sequential> H2;
		RnsDomain Zrns(F.rns());
		fgemm(Zrns,ta,tb,m,n,k,alpha,Ad,lda,Bd,ldb,beta,Cd,ldc,H2);

		// reduce the product mod p (note that entries are larger than p, due to RNS modulo reduction)
		finit(F,m,n,Cd,ldc);
		return Cd;
	}


	// fgemm for IntegerDomain with Winograd Helper
	inline FFPACK::integer* fgemm (const FFPACK::Modular<FFPACK::Integer>& F,
				       const FFLAS_TRANSPOSE ta,
				       const FFLAS_TRANSPOSE tb,
				       const size_t m, const size_t n,const size_t k,
				       const FFPACK::Integer alpha,
				       FFPACK::Integer *A, const size_t lda,
				       FFPACK::Integer *B, const size_t ldb,
				       const FFPACK::Integer beta,
				       FFPACK::Integer* C, const size_t ldc,
				       MMHelper<FFPACK::Modular<FFPACK::integer>, MMHelperAlgo::Winograd, FieldCategories::MultiPrecisionTag,ParSeqHelper::Sequential> & H)
	{
		// compute the product over Z
		typedef FFPACK::UnparametricField<FFPACK::Integer> IntegerDomain;
		FFPACK::Integer p;
		F.cardinality(p);
		IntegerDomain Z;
		MMHelper<IntegerDomain, MMHelperAlgo::Winograd, FieldCategories::MultiPrecisionTag,ParSeqHelper::Sequential> H2(p,p);//H2(Z,0,H.parseq);
		fgemm(Z,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H2);

		// reduce the product mod p
		finit(F,m,n,C,ldc);

		return C;
	}



}// END of namespace FFLAS

#endif
#endif
