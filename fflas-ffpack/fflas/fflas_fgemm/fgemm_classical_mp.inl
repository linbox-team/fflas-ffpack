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

#include <givaro/modular-integer.h>
#include <givaro/zring.h>

#include "fflas-ffpack/field/rns-double.h"
#include "fflas-ffpack/field/rns-integer.h"
#include "fflas-ffpack/field/rns-integer-mod.h"
#include "fflas-ffpack/field/field-traits.h"
#include "fflas-ffpack/fflas/fflas_helpers.inl" 

namespace FFLAS {
 
	template<typename Field,
		 typename AlgoTrait,
		 typename ParSeqTrait>
	struct MMHelper<Field, AlgoTrait,ModeCategories::ConvertTo<ElementCategories::RNSElementTag>, ParSeqTrait> {
		Givaro::Integer normA,normB;
		int recLevel;
		ParSeqTrait parseq;
		MMHelper() : normA(0), normB(0), recLevel(-1) {}
		template <class F2, class A2, class M2, class PS2>
		MMHelper(MMHelper<F2, A2, M2, PS2> H2) : 
				normA(H2.normA), normB(H2.normB), recLevel(H2.recLevel), parseq(H2.parseq) {}
		MMHelper(Givaro::Integer Amax, Givaro::Integer Bmax) : normA(Amax), normB(Bmax), recLevel(-1) {}
		MMHelper(const Field& F, size_t m, size_t n, size_t k, ParSeqTrait PS=ParSeqTrait())
			: recLevel(-1), parseq(PS)
		{F.characteristic(normA);F.characteristic(normB);}
		MMHelper(const Field& F, int wino, ParSeqTrait PS=ParSeqTrait()) : recLevel(wino), parseq(PS)
		{F.characteristic(normA);F.characteristic(normB);}
		void setNorm(Givaro::Integer p){normA=normB=p;}
	};
	template<typename E,
		 typename AlgoTrait,
		 typename ParSeqTrait>
	struct MMHelper<FFPACK::RNSInteger<E>, AlgoTrait,ModeCategories::DefaultTag, ParSeqTrait> {
		Givaro::Integer normA,normB;
		int recLevel;
		ParSeqTrait parseq;
		MMHelper() : normA(0), normB(0), recLevel(-1) {}
		MMHelper(Givaro::Integer Amax, Givaro::Integer Bmax) : normA(Amax), normB(Bmax), recLevel(-1) {}
		MMHelper(const FFPACK::RNSInteger<E>& F, size_t m, size_t n, size_t k, ParSeqTrait PS=ParSeqTrait())
			: recLevel(-1), parseq(PS)
		{F.characteristic(normA);F.characteristic(normB);}
		MMHelper(const FFPACK::RNSInteger<E>& F, int wino, ParSeqTrait PS=ParSeqTrait()) : recLevel(wino), parseq(PS)
		{F.characteristic(normA);F.characteristic(normB);}
		void setNorm(Givaro::Integer p){normA=normB=p;}
	};
	template<typename E,
		 typename AlgoTrait,
		 typename ParSeqTrait>
	struct MMHelper<FFPACK::RNSIntegerMod<E>, AlgoTrait,ModeCategories::DefaultTag, ParSeqTrait> {
		Givaro::Integer normA,normB;
		int recLevel;
		ParSeqTrait parseq;
		MMHelper() : normA(0), normB(0), recLevel(-1) {}
		MMHelper(Givaro::Integer Amax, Givaro::Integer Bmax) : normA(Amax), normB(Bmax), recLevel(-1) {}
		MMHelper(const FFPACK::RNSIntegerMod<E>& F, size_t m, size_t n, size_t k, ParSeqTrait PS=ParSeqTrait())
			: recLevel(-1), parseq(PS)
		{F.characteristic(normA);F.characteristic(normB);}
		MMHelper(const FFPACK::RNSIntegerMod<E>& F, int wino, ParSeqTrait PS=ParSeqTrait()) : recLevel(wino), parseq(PS)
		{F.characteristic(normA);F.characteristic(normB);}
		void setNorm(Givaro::Integer p){normA=normB=p;}
	};

	/***********************************
	 *** MULTIPRECISION FGEMM OVER Z ***
	 ***********************************/

	// fgemm for RnsInteger with Winograd Helper
	template<typename RNS>
	inline  typename FFPACK::RNSInteger<RNS>::Element_ptr 
	fgemm (const FFPACK::RNSInteger<RNS> &F,
	       const FFLAS_TRANSPOSE ta,
	       const FFLAS_TRANSPOSE tb,
	       const size_t m, const size_t n,const size_t k,
	       const typename FFPACK::RNSInteger<RNS>::Element alpha,
	       typename FFPACK::RNSInteger<RNS>::ConstElement_ptr Ad, const size_t lda,
	       typename FFPACK::RNSInteger<RNS>::ConstElement_ptr Bd, const size_t ldb,
	       const typename FFPACK::RNSInteger<RNS>::Element beta,
	       typename FFPACK::RNSInteger<RNS>::Element_ptr Cd, const size_t ldc,
	       MMHelper<FFPACK::RNSInteger<RNS>, MMHelperAlgo::Winograd> & H)
	{		

		// compute each fgemm componentwise
		for(size_t i=0;i<F.size();i++){
			MMHelper<typename RNS::ModField,MMHelperAlgo::Winograd> H2(F.rns()._field_rns[i], H.recLevel, H.parseq);
			FFLAS::fgemm(F.rns()._field_rns[i],ta,tb,//FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
				     m, n, k, alpha._ptr[i*alpha._stride],
				     Ad._ptr+i*Ad._stride, lda, 
				     Bd._ptr+i*Bd._stride, ldb, 
				     beta._ptr[i*beta._stride], 
				     Cd._ptr+i*Cd._stride, ldc);			
		}
		
		return Cd;
	} 


	// fgemm for UnparametricField<Integer> with Winograd Helper (bb: file is classical ??)
	inline Givaro::Integer* 
	fgemm (const Givaro::ZRing<Givaro::Integer>& F,
	       const FFLAS_TRANSPOSE ta,
	       const FFLAS_TRANSPOSE tb,
	       const size_t m, const size_t n,const size_t k,
	       const Givaro::Integer alpha,
	       const Givaro::Integer* A, const size_t lda,
	       const Givaro::Integer* B, const size_t ldb,
	       Givaro::Integer beta,
	       Givaro::Integer* C, const size_t ldc,
	       MMHelper<Givaro::ZRing<Givaro::Integer>, MMHelperAlgo::Winograd, ModeCategories::ConvertTo<ElementCategories::RNSElementTag> >  & H)
	{
		if (alpha == 0){
			fscalin(F,m,n,beta,C,ldc);
			return C;
		}
		if (k==0) return C;
		// compute bit size of feasible prime for FFLAS
		size_t _k=k,lk=0;
		while ( _k ) {_k>>=1; ++lk;}
		size_t prime_bitsize= (53-lk)>>1;

		// compute bound on the output
		Givaro::Integer  mA,mB,mC;
		size_t logA,logB;
		mA=H.normA;
		mB=H.normB;
		if (H.normA==0){
			logA=A[0].bitsize();
			for (size_t i=1;i<m*k;i++)
				logA = std::max(logA,A[i].bitsize());
			H.normA=1; H.normA<<=uint64_t(logA);
		}
		else {
			logA=H.normA.bitsize();
		}
		if (H.normB==0){
			logB=B[0].bitsize();
			for (size_t i=1;i<k*n;i++)
				logB=std::max(logB,B[i].bitsize());
			H.normB=1; H.normB<<=uint64_t(logB);
		}
		else {
			logB=H.normA.bitsize();
		}
		mC = 2*uint64_t(k)*H.normA*H.normB*abs(alpha); // need to use 2x bound to reach both positive and negative

		// construct an RNS structure and its associated Domain
		FFPACK::rns_double RNS(mC, prime_bitsize);
		typedef FFPACK::RNSInteger<FFPACK::rns_double> RnsDomain;
		RnsDomain Zrns(RNS);
		
		size_t Acold,Arowd,Bcold,Browd;
		if (ta == FFLAS::FflasNoTrans){Arowd=m; Acold = k; }
		else { Arowd=k; Acold = m;}
		if (tb == FFLAS::FflasNoTrans){Browd=k; Bcold = n; }
		else { Browd=n; Bcold = k;}
		
		// allocate data for RNS representation
		typename RnsDomain::Element_ptr Ap,Bp,Cp;
		Ap = FFLAS::fflas_new(Zrns,Arowd,Acold);
		Bp = FFLAS::fflas_new(Zrns,Browd,Bcold);
		Cp = FFLAS::fflas_new(Zrns,m,n);

		// convert the input matrices to RNS representation
		finit_rns(Zrns,Arowd,Acold,(logA/16)+((logA%16)?1:0),A,lda,Ap);
		finit_rns(Zrns,Browd,Bcold,(logB/16)+((logB%16)?1:0),B,ldb,Bp);

		// perform the fgemm in RNS
		MMHelper<RnsDomain, MMHelperAlgo::Winograd, ModeCategories::DefaultTag>  H2(Zrns,H.recLevel,H.parseq);

		// compute alpha and beta in RNS
		typename RnsDomain::Element alphap, betap;
		Zrns.init(alphap, alpha);
		Zrns.init(betap, F.zero);

		// call  fgemm
		fgemm(Zrns,ta,tb,m,n,k,alphap,Ap,Acold,Bp,Bcold,betap,Cp,n,H2);
		
		// convert the RNS output to integer representation (C=beta.C+ RNS^(-1)(Cp) )
		fconvert_rns(Zrns,m,n,beta,C,ldc,Cp);

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
						typename RNS::ConstElement_ptr Ad, const size_t lda,
						typename RNS::ConstElement_ptr Bd, const size_t ldb,
						const typename RNS::Element beta,
						typename RNS::Element_ptr Cd, const size_t ldc,
						MMHelper<FFPACK::RNSIntegerMod<RNS>, MMHelperAlgo::Winograd> & H)
	{
		// compute the product over Z
		typedef FFPACK::RNSInteger<RNS> RnsDomain;
		RnsDomain Zrns(F.rns());
		MMHelper<RnsDomain, MMHelperAlgo::Winograd> H2(Zrns, H.recLevel,H.parseq);
#ifdef BENCH_PERF_FGEMM_MP
		FFLAS::Timer chrono;chrono.start();
#endif
		fgemm(Zrns,ta,tb,m,n,k,alpha,Ad,lda,Bd,ldb,beta,Cd,ldc,H2);
		// reduce the product mod p (note that entries are larger than p, due to RNS modulo reduction)
		freduce (F, m, n, Cd, ldc);
#ifdef BENCH_PERF_FGEMM_MP
		chrono.stop();
		F.t_igemm+=chrono.usertime();
#endif

		return Cd;
	}


	// fgemm for IntegerDomain with Winograd Helper
	inline Givaro::Integer* fgemm (const Givaro::Modular<Givaro::Integer>& F,
				       const FFLAS_TRANSPOSE ta,
				       const FFLAS_TRANSPOSE tb,
				       const size_t m, const size_t n,const size_t k,
				       const Givaro::Integer alpha,
				       const Givaro::Integer *A, const size_t lda,
				       const Givaro::Integer *B, const size_t ldb,
				       const Givaro::Integer beta,
				       Givaro::Integer* C, const size_t ldc,
				       MMHelper<Givaro::Modular<Givaro::Integer>, MMHelperAlgo::Winograd, ModeCategories::ConvertTo<ElementCategories::RNSElementTag> > & H)
	//MMHelper<Givaro::Modular<Givaro::Integer>, MMHelperAlgo::Winograd, FieldCategories::MultiPrecisionTag,ParSeqHelper::Sequential> & H)

	{
		// compute the product over Z
		typedef Givaro::ZRing<Givaro::Integer> IntegerDomain;
		Givaro::Integer p;
		F.cardinality(p);
		IntegerDomain Z;
		MMHelper<IntegerDomain,MMHelperAlgo::Winograd, ModeCategories::ConvertTo<ElementCategories::RNSElementTag> > H2(Z,H.recLevel,H.parseq);
		H2.setNorm(p);
		
		fgemm(Z,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H2);
		
		// reduce the product mod p
		freduce (F, m, n, C, ldc);

		return C;
	}

	// fgemm for IntegerDomain with Winograd Helper
	template<class Cut, class Param>
	inline Givaro::Integer* fgemm (const Givaro::Modular<Givaro::Integer>& F,
				       const FFLAS_TRANSPOSE ta,
				       const FFLAS_TRANSPOSE tb,
				       const size_t m, const size_t n,const size_t k,
				       const Givaro::Integer alpha,
				       const Givaro::Integer *A, const size_t lda,
				       const Givaro::Integer *B, const size_t ldb,
				       const Givaro::Integer beta,
				       Givaro::Integer* C, const size_t ldc,
				       MMHelper<Givaro::Modular<Givaro::Integer>, MMHelperAlgo::Winograd, ModeCategories::ConvertTo<ElementCategories::RNSElementTag>,ParSeqHelper::Parallel<Cut, Param> > & H)

	{
		// compute the product over Z
		typedef Givaro::ZRing<Givaro::Integer> IntegerDomain;
		Givaro::Integer p;
		F.cardinality(p);
		IntegerDomain Z;
		MMHelper<IntegerDomain,MMHelperAlgo::Winograd,  ModeCategories::ConvertTo<ElementCategories::RNSElementTag>, ParSeqHelper::Sequential> H2(Z,H.recLevel);
		H2.setNorm(p);
		
		fgemm(Z,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H2);
		
		// reduce the product mod p
		freduce (F, m, n, C, ldc);

		return C;
	}

	// PARALLEL VERSION (NOT PARALLEL YET)
	template<class Cut, class Param>
	inline Givaro::Integer* fgemm (const Givaro::ZRing<Givaro::Integer>& F,
				       const FFLAS_TRANSPOSE ta,
				       const FFLAS_TRANSPOSE tb,
				       const size_t m, const size_t n,const size_t k,
				       const Givaro::Integer alpha,
				       const Givaro::Integer* A, const size_t lda,
				       const Givaro::Integer* B, const size_t ldb,
				       Givaro::Integer beta,
				       Givaro::Integer* C, const size_t ldc,
				       MMHelper<Givaro::ZRing<Givaro::Integer>,MMHelperAlgo::Winograd,FieldCategories::UnparametricTag,ParSeqHelper::Parallel<Cut,Param> > & H){
		MMHelper<Givaro::ZRing<Givaro::Integer>,MMHelperAlgo::Winograd> H2(F,H.recLevel);
		return fgemm(F,ta,tb,m,n,k,alpha,A,lda,B,lda,beta,C,ldc,H2);
	}

}// END of namespace FFLAS

#endif

