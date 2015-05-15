/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_fgemm.inl
 * Copyright (C) 2005 Clement Pernet
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by Clement Pernet  < Clement.Pernet@imag.fr >
 *            BB <bbboyer@ncsu.edu>
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

#ifndef __FFLASFFPACK_fgemm_INL
#define __FFLASFFPACK_fgemm_INL

#include <givaro/modular.h>
#include <givaro/modular-balanced.h>

#include "fflas-ffpack/utils/debug.h"

namespace FFLAS { namespace Protected{

	template <typename FloatElement, class Field, class FieldMode>
	inline typename Field::Element_ptr
	fgemm_convert (const Field& F,
		       const FFLAS_TRANSPOSE ta,
		       const FFLAS_TRANSPOSE tb,
		       const size_t m, const size_t n, const size_t k,
		       const typename Field::Element alpha,
		       typename Field::ConstElement_ptr A,const size_t lda,
		       typename Field::ConstElement_ptr B,const size_t ldb,
		       const typename Field::Element beta,
		       typename Field::Element_ptr C, const size_t ldc,
		       MMHelper<Field, MMHelperAlgo::Winograd, FieldMode> & H)
	{
		// CP: lda, ldb, ldc can be zero (if m,n or k is 0) and since  this may have not 
		// been checked by the caller at this point.
		// FFLASFFPACK_check(lda);
		// FFLASFFPACK_check(ldb);
		// FFLASFFPACK_check(ldc);

		Givaro::ModularBalanced<FloatElement> G((FloatElement) F.characteristic());
		FloatElement tmp,alphaf, betaf;
		    // This conversion is quite tricky, but convert and init are required
		    // in sequence e.g. for when F is a ModularBalanced field and alpha == -1
		F.convert (tmp, beta);
		G.init(betaf, tmp);
		F.convert (tmp, alpha);
		G.init(alphaf, tmp);

		FloatElement* Af = FFLAS::fflas_new(G, m, k);
		FloatElement* Bf = FFLAS::fflas_new(G, k, n);
		FloatElement* Cf = FFLAS::fflas_new(G, m, n);

		size_t ma, ka, kb, nb; //mb, na
		if (ta == FflasTrans) { ma = k; ka = m; }
		else { ma = m; ka = k; }
		if (tb == FflasTrans) { kb = n; nb = k; }
		else {  kb = k; nb = n; }
		size_t ldaf = ka, ldbf = nb, ldcf= n;

		fconvert(F, ma, ka, Af, ka, A, lda);
		freduce(G, ma, ka, Af, ka);
		fconvert(F, kb, nb, Bf, nb, B, ldb);
		freduce(G, kb, nb, Bf, nb);

		if (!F.isZero(beta)){
			fconvert(F, m, n, Cf, n, C, ldc);
			freduce (G, m, n, Cf, n);
		}
		MMHelper<Givaro::ModularBalanced<FloatElement>, 
			 MMHelperAlgo::Winograd> 
			HG(G,H.recLevel, ParSeqHelper::Sequential());
		fgemm (G, ta, tb, m, n, k, alphaf, Af, ldaf, Bf, ldbf, betaf, Cf, ldcf, HG);

		finit (F, m, n, Cf, n, C, ldc);

		fflas_delete (Af);
		fflas_delete (Bf);
		fflas_delete (Cf);
		return C;
	}
	}//Protected
}//FFLAS

namespace FFLAS{ namespace Protected{
	template <class Field, class Element, class AlgoT>
	inline bool NeedPreAddReduction (Element& Outmin, Element& Outmax,
					 Element& Op1min, Element& Op1max,
					 Element& Op2min, Element& Op2max,
					 MMHelper<Field, AlgoT, ModeCategories::LazyTag >& WH)
	{
		Outmin = Op1min + Op2min;
		Outmax = Op1max + Op2max;
		if (std::max(-Outmin, Outmax) > WH.MaxStorableValue){
			// Reducing both Op1 and Op2
			Op1min = Op2min = WH.FieldMin;
			Op1max = Op2max = WH.FieldMax;
			Outmin = 2*WH.FieldMin;
			Outmax = 2*WH.FieldMax;
			return true;
		} else return false;
	}

	template <class Field, class Element, class AlgoT, class ModeT>
	inline bool NeedPreAddReduction (Element& Outmin, Element& Outmax,
					 Element& Op1min, Element& Op1max,
					 Element& Op2min, Element& Op2max,
					 MMHelper<Field, AlgoT, ModeT >& WH)
	{
		Outmin = WH.FieldMin;
		Outmax = WH.FieldMax;
		return false;
	}

	template <class Field, class Element, class AlgoT>
	inline bool NeedPreSubReduction (Element& Outmin, Element& Outmax,
					 Element& Op1min, Element& Op1max,
					 Element& Op2min, Element& Op2max,
					 MMHelper<Field, AlgoT, ModeCategories::LazyTag >& WH)
	{
		Outmin = Op1min - Op2max;
		Outmax = Op1max - Op2min;
		if (std::max(-Outmin, Outmax) > WH.MaxStorableValue){
			// Reducing both Op1 and Op2
			Op1min = Op2min = WH.FieldMin;
			Op1max = Op2max = WH.FieldMax;
			Outmin = WH.FieldMin-WH.FieldMax;
			Outmax = -Outmin;
			return true;
		} else return false;
	}

	template <class Field, class Element, class AlgoT, class ModeT>
	inline bool NeedPreSubReduction (Element& Outmin, Element& Outmax,
					 Element& Op1min, Element& Op1max,
					 Element& Op2min, Element& Op2max,
					 MMHelper<Field, AlgoT, ModeT >& WH)
	{
		    // Necessary? -> CP: Yes, for generic Mode of op
		Outmin = WH.FieldMin;
		Outmax = WH.FieldMax;
		return false;
	}


	template<class Field, class Element, class AlgoT>
	inline bool NeedDoublePreAddReduction (Element& Out1min, Element& Out1max,
					       Element& Out2min, Element& Out2max,
					       Element& Op1min, Element& Op1max,
					       Element& Op2min, Element& Op2max, Element beta,
					       MMHelper<Field, AlgoT, ModeCategories::LazyTag >& WH)
	{
		// Testing if P5 need to be reduced
		Out2min = Op1min + std::min(beta*Op2min,beta*Op2max);
		Out2max = Op1max + std::max(beta*Op2min,beta*Op2max);
		Out1min = Op1min + std::min(beta*Op2min,beta*Op2max);
		Out1max = Op1max + std::max(beta*Op2min,beta*Op2max);
		if (std::max(-Out2min, std::max(Out2max, std::max(-Out1min, Out1max))) > WH.MaxStorableValue){
			Out2min = Op1min + std::min(beta*Op2min,beta*Op2max);
			Out2max = Op1max + std::max(beta*Op2min,beta*Op2max);
			Out1min = Op1min + std::min(beta*Op2min,beta*Op2max);
			Out1max = Op1max + std::max(beta*Op2min,beta*Op2max);
			return true;
		} else return false;
	}

	template<class Field, class Element, class AlgoT, class ModeT>
	inline bool NeedDoublePreAddReduction (Element& Out1min, Element& Out1max,
					       Element& Out2min, Element& Out2max,
					       Element& Op1min, Element& Op1max,
					       Element& Op2min, Element& Op2max, Element beta,
					       MMHelper<Field, AlgoT, ModeT>& WH)
	{
		Out2min = WH.FieldMin;
		Out2max = WH.FieldMax;
		Out1min = WH.FieldMin;
		Out1max = WH.FieldMax;
		return false;
	}

	template <class Field, class AlgoT>
	inline void ScalAndReduce (const Field& F, const size_t N,
				   const typename Field::Element alpha,
				   typename Field::Element_ptr X, const size_t incX,
				   const MMHelper<Field, AlgoT, ModeCategories::LazyTag >& H)
	{
		if (!F.isOne(alpha) && !F.isMOne(alpha)){
			typename MMHelper<Field, AlgoT, ModeCategories::LazyTag >::DFElt al; 
			F.convert(al, alpha);
			if (al < 0) al = -al;
			if (al*std::max(-H.Outmin, H.Outmax) > H.MaxStorableValue){
				freduce (F, N, X, incX);
				fscalin (F, N, alpha, X, incX);
			} else {
				fscalin (H.delayedField, N, alpha, X, incX);
				freduce (F, N, X, incX);
			}
		} else
			freduce (F, N, X, incX);
	}

	template <class Field, class AlgoT>
	inline void ScalAndReduce (const Field& F, const size_t M, const size_t N,
				   const typename Field::Element alpha,
				   typename Field::Element_ptr A, const size_t lda,
				   const MMHelper<Field, AlgoT, ModeCategories::LazyTag >& H)
	{
		if (!F.isOne(alpha) && !F.isMOne(alpha)){
			typename MMHelper<Field, AlgoT, ModeCategories::LazyTag >::DFElt al; 
			F.convert(al, alpha);
			if (al<0) al = -al;
			if (al*std::max(-H.Outmin, H.Outmax) > H.MaxStorableValue){
				freduce (F, M, N, A, lda);
				fscalin (F, M, N, alpha, A, lda);
			} else {
				fscalin (H.delayedField, M, N, alpha, (typename MMHelper<Field, AlgoT, ModeCategories::LazyTag >::DFElt*)A, lda);
				freduce (F, M, N, A, lda);
			}
		} else
			freduce (F, M, N, A, lda);
	}

} // Protected
} // FFLAS

namespace FFLAS {

	template<class Field>
	inline  typename Field::Element_ptr
	fgemm (const Field& F,
	       const FFLAS_TRANSPOSE ta,
	       const FFLAS_TRANSPOSE tb,
	       const size_t m, const size_t n, const size_t k,
	       const typename Field::Element alpha,
	       typename Field::ConstElement_ptr A, const size_t lda,
	       typename Field::ConstElement_ptr B, const size_t ldb,
	       const typename Field::Element beta,
	       typename Field::Element_ptr C, const size_t ldc,
	       MMHelper<Field, MMHelperAlgo::Winograd, ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> > & H)
	{
		if (F.cardinality() < DOUBLE_TO_FLOAT_CROSSOVER)
			return Protected::fgemm_convert<float,Field>(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H);
		else
			return Protected::fgemm_convert<double,Field>(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H);
	}
}// FFLAS


// fgemm
namespace FFLAS {

	template<typename Field>
	inline typename Field::Element_ptr
	fgemm( const Field& F,
	       const FFLAS_TRANSPOSE ta,
	       const FFLAS_TRANSPOSE tb,
	       const size_t m,
	       const size_t n,
	       const size_t k,
	       const typename Field::Element alpha,
	       typename Field::ConstElement_ptr A, const size_t lda,
	       typename Field::ConstElement_ptr B, const size_t ldb,
	       const typename Field::Element beta,
	       typename Field::Element_ptr C, const size_t ldc,
	       const ParSeqHelper::Sequential seq)
	{
		    // The entry point to fgemm.
		    // Place where the algorithm is chosen. Winograd's alg. is now the default.
		MMHelper<Field, MMHelperAlgo::Winograd, typename FFLAS::ModeTraits<Field>::value, ParSeqHelper::Sequential > HW (F, m, k, n, seq);
		return 	fgemm (F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, HW);
	}

	template<typename Field>
	inline typename Field::Element_ptr
	fgemm( const Field& F,
	       const FFLAS_TRANSPOSE ta,
	       const FFLAS_TRANSPOSE tb,
	       const size_t m,
	       const size_t n,
	       const size_t k,
	       const typename Field::Element alpha,
	       typename Field::ConstElement_ptr A, const size_t lda,
	       typename Field::ConstElement_ptr B, const size_t ldb,
	       const typename Field::Element beta,
	       typename Field::Element_ptr C, const size_t ldc,
	       const ParSeqHelper::Parallel par)
	{

		MMHelper<Field, MMHelperAlgo::Winograd, typename FFLAS::ModeTraits<Field>::value, ParSeqHelper::Parallel > HW (F, m, k, n, par);
		return 	fgemm (F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, HW);
	}

	template<typename Field>
	inline typename Field::Element_ptr
	fgemm( const Field& F,
	       const FFLAS_TRANSPOSE ta,
	       const FFLAS_TRANSPOSE tb,
	       const size_t m,
	       const size_t n,
	       const size_t k,
	       const typename Field::Element alpha,
	       typename Field::ConstElement_ptr A, const size_t lda,
	       typename Field::ConstElement_ptr B, const size_t ldb,
	       const typename Field::Element beta,
	       typename Field::Element_ptr C, const size_t ldc)
	{
		if (!m || !n) {return C;}
		
		if (!k || F.isZero (alpha)){
			fscalin(F, m, n, beta, C, ldc);
		 	return C;
		}
		return fgemm(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,FFLAS::ParSeqHelper::Sequential());
	}


	template<class Field>
	inline typename Field::Element_ptr
	fgemm( const Field& F,
	       const FFLAS_TRANSPOSE ta,
	       const FFLAS_TRANSPOSE tb,
	       const size_t m,
	       const size_t n,
	       const size_t k,
	       const typename Field::Element alpha,
	       typename Field::ConstElement_ptr A, const size_t lda,
	       typename Field::ConstElement_ptr B, const size_t ldb,
	       const typename Field::Element beta,
	       typename Field::Element_ptr C, const size_t ldc,
	       MMHelper<Field, MMHelperAlgo::Winograd, ModeCategories::DelayedTag> & H)
	{		
		if (!m || !n) {return C;}

		if (!k || F.isZero (alpha)){
			fscalin(F, m, n, beta, C, ldc);
			return C;
		}
#ifndef NDEBUG
		/*  check if alpha is invertible.
		 *  XXX do it in F.isInvertible(Element&) ?
		 *  XXX do it in return status of F.inv(Element&,Element&)
		 */
		typename Field::Element e ;
		F.init(e,beta);
		F.divin(e,alpha);
		F.mulin(e,alpha);
		FFLASFFPACK_check(F.areEqual(e,beta));
#endif

#if 0
		// detect fgemv
		if (n == 1 and ...) {}
		// detect fger
		if (k==1 and ...) {}
#endif
		if (Protected::AreEqual<Field, Givaro::Modular<double> >::value ||
		    Protected::AreEqual<Field, Givaro::ModularBalanced<double> >::value){
			    //Givaro::Modular<double> need to switch to float if p too small
			if (F.characteristic() < DOUBLE_TO_FLOAT_CROSSOVER)
				return Protected::fgemm_convert<float,Field>(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H);
		}
		typename Field::Element alpha_,beta_;
		if ( !F.isOne(alpha) && !F.isMOne(alpha)){
			F.assign (alpha_, F.one);
			F.div (beta_, beta, alpha);
		} else {
			F.assign (alpha_,alpha);
			F.assign (beta_,beta);
		}

		MMHelper<Field, MMHelperAlgo::Winograd, ModeCategories::LazyTag>  HD(H);

		fgemm (F, ta, tb, m, n, k, alpha_, A, lda, B, ldb, beta_, C, ldc, HD);

		Protected::ScalAndReduce (F, m, n, alpha, C, ldc, HD);

		H.initOut();

		return C;
	}
} // FFLAS

// #include "fflas_fgemm/matmul_algos.inl"
#include "fflas_fgemm/fgemm_classical.inl"
#include "fflas_fgemm/fgemm_winograd.inl"
// #include "fflas_fgemm/gemm_bini.inl"

// fsquare
namespace FFLAS {
	template < class Field >
	inline typename Field::Element_ptr
	fsquare (const Field& F,
		 const FFLAS_TRANSPOSE ta,
		 const size_t n, const typename Field::Element alpha,
		  typename Field::ConstElement_ptr A, const size_t lda,
		 const typename Field::Element beta,
		 typename Field::Element_ptr C, const size_t ldc)
	{

		double alphad, betad;
		F.convert (alphad, alpha);
		if (F.isMOne (beta))
			betad = -1.0;
		else
			F.convert (betad, beta);

		//! @bug why double ?
		// Double  matrices initialisation
		Givaro::DoubleDomain::Element_ptr Ad = fflas_new (Givaro::DoubleDomain(),n,n);
		Givaro::DoubleDomain::Element_ptr Cd = fflas_new (Givaro::DoubleDomain(),n,n);
		// Conversion finite Field = >  double
		fconvert (F, n, n, Ad, n, A, lda);
		if (!F.isZero(beta)) fconvert(F, n, n, Cd, n, C, ldc);

		// Call to the blas Multiplication
		FFLASFFPACK_check(n);
		cblas_dgemm (CblasRowMajor, (CBLAS_TRANSPOSE)ta,
			     (CBLAS_TRANSPOSE)ta, (int)n, (int)n, (int)n,
			     (Givaro::DoubleDomain::Element) alphad, Ad, (int)n, Ad, (int)n,
			     (Givaro::DoubleDomain::Element) betad, Cd, (int)n);
		// Conversion double = >  Finite Field
		fflas_delete (Ad);
		finit (F,n,n, Cd, n, C, ldc);
		fflas_delete (Cd);
		return C;
	}

	namespace Protected {

		// F is Modular(Balanced)<float/double>
		template < class Field >
		inline typename Field::Element_ptr
		fsquareCommon (const Field& F,
			       const FFLAS_TRANSPOSE ta,
			       const size_t n, const typename Field::Element alpha,
			       typename Field::ConstElement_ptr A, const size_t lda,
			       const typename Field::Element beta,
			       typename Field::Element_ptr C, const size_t ldc)
		{
			if (C==A) {
				typename Field::Element_ptr Ad = fflas_new (F, n, n);
				fassign(F,n,n,A,lda,Ad,n);
				fgemm (F, ta, ta, n, n, n, alpha, Ad, n, Ad, n, beta, C, ldc);
				fflas_delete (Ad);
			}
			else
				fgemm (F, ta, ta, n, n, n, alpha, A, lda, A, lda, beta, C, ldc);
			return C;

		}

	} // Protected

	template <>
	inline double* fsquare (const  Givaro::ModularBalanced<double> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const double alpha,
				const double* A, const size_t lda,
				const double beta,
				double* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

	template <>
	inline float * fsquare (const  Givaro::ModularBalanced<float> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const float alpha,
				const float* A, const size_t lda,
				const float beta,
				float* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

	template <>
	inline double* fsquare (const  Givaro::Modular<double> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const double alpha,
				const double* A, const size_t lda,
				const double beta,
				double* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

	template <>
	inline float * fsquare (const  Givaro::Modular<float> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const float alpha,
				const float* A, const size_t lda,
				const float beta,
				float* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

} // FFLAS

#endif // __FFLASFFPACK_fgemm_INL
