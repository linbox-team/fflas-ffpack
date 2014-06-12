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

//#include "fflas_fgemm/fflas_bounds_winograd.inl"

namespace FFLAS { namespace Protected{
		
	template <typename FloatElement, class Field>
	inline typename Field::Element*
	fgemm_convert (const Field& F,
		       const FFLAS_TRANSPOSE ta,
		       const FFLAS_TRANSPOSE tb,
		       const size_t m, const size_t n, const size_t k,
		       const typename Field::Element alpha,
		       const typename Field::Element* A,const size_t lda,
		       const typename Field::Element* B,const size_t ldb,
		       const typename Field::Element beta,
		       typename Field::Element * C, const size_t ldc, const int w)
	{
		FFLASFFPACK_check(lda);
		FFLASFFPACK_check(ldb);
		FFLASFFPACK_check(ldc);

		FFPACK::Modular<FloatElement> G((FloatElement) F.characteristic());
		FloatElement tmp,alphaf, betaf;
		    // This conversion is quite tricky, but convert and init are required
		    // in sequence e.g. for when F is a ModularBalanced field and alpha == -1
		F.convert (tmp, beta);
		G.init(betaf, tmp);
		F.convert (tmp, alpha);
		G.init(alphaf, tmp);
		
		FloatElement * Af = new FloatElement[m*k];
		FloatElement * Bf = new FloatElement[k*n];
		FloatElement * Cf = new FloatElement[m*n];

		size_t ma, ka, kb, nb; //mb, na
		if (ta == FflasTrans) { ma = k; ka = m; }
		else { ma = m; ka = k; }
		if (tb == FflasTrans) { kb = n; nb = k; }
		else {  kb = k; nb = n; }
		size_t ldaf = ka, ldbf = nb, ldcf= n;

		fconvert(F, ma, ka, Af, ka, A, lda);
		fconvert(F, kb, nb, Bf, nb, B, ldb);
		
		if (!F.isZero(beta))
			fconvert(F, m, n, Cf, n, C, ldc);

		fgemm (G, ta, tb, m, n, k, alphaf, Af, ldaf, Bf, ldbf, betaf, Cf, ldcf, w);

		finit(F, m, n, Cf, n, C, ldc);

		delete[] Af;
		delete[] Bf;
		delete[] Cf;
		return C;
	}
	}//Protected
}//FFLAS

namespace FFLAS{ namespace Protected{
	template <class AlgoT, class Field>
	inline bool NeedPreAddReduction (double& Outmin, double& Outmax, 
					 double& Op1min, double& Op1max, 
					 double& Op2min, double& Op2max, 
					 MMHelper<AlgoT, FieldCategories::ModularFloatingPointTag, Field >& WH)
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

	template <class AlgoT, class FieldT, class Field>
	inline bool NeedPreAddReduction (double& Outmin, double& Outmax, 
				     double& Op1min, double& Op1max, 
				     double& Op2min, double& Op2max, 
				     MMHelper<AlgoT, FieldT, Field >& WH)
	{
		Outmin = WH.FieldMin;
		Outmax = WH.FieldMax;
		return false;
	}

	template <class AlgoT, class Field>
	inline bool NeedPreSubReduction (double& Outmin, double& Outmax, 
					 double& Op1min, double& Op1max, 
					 double& Op2min, double& Op2max, 
					 MMHelper<AlgoT, FieldCategories::ModularFloatingPointTag, Field >& WH)
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

	template <class AlgoT, class FieldT, class Field>
	inline bool NeedPreSubReduction (double& Outmin, double& Outmax, 
					 double& Op1min, double& Op1max, 
					 double& Op2min, double& Op2max, 
					 MMHelper<AlgoT, FieldT, Field >& WH)
	{
		    // Necessary?
		Outmin = WH.FieldMin;
		Outmax = WH.FieldMax;
		return false;
	}

	template <class AlgoT, class FieldT, class Field>
	inline void ScalAndInit (const Field& F, const size_t N, 
				 const typename Field::Element alpha,
				 typename Field::Element * X, const size_t incX,
				 const MMHelper<AlgoT, FieldT, Field >& H)
	{
		finit(F, N, X, incX);
	}
	
	template <class AlgoT, class Field>
	inline void ScalAndInit (const Field& F, const size_t N,
				 const typename Field::Element alpha,
				 typename Field::Element * X, const size_t incX,
				 const MMHelper<AlgoT, FieldCategories::ModularFloatingPointTag, Field >& H)
	{
		if (!F.isOne(alpha) && !F.isMOne(alpha)){
			if (abs(alpha)*std::max(-H.Outmin, H.Outmax) > H.MaxStorableValue){
				finit (F, N, X, incX);
				fscalin (F, N, alpha, X, incX);
			} else {
				fscalin (H.delayedField, N, alpha, X, incX);
				finit(F, N, X, incX);
			}
		} else
			finit(F, N, X, incX);
	}

	template <class AlgoT, class FieldT, class Field>
	inline void ScalAndInit (const Field& F, const size_t M, const size_t N, 
				 const typename Field::Element alpha,
				 typename Field::Element * A, const size_t lda,
				 const MMHelper<AlgoT, FieldT, Field >& H)
	{
		finit(F, M, N, A, lda);
	}
	
	template <class AlgoT, class Field>
	inline void ScalAndInit (const Field& F, const size_t M, const size_t N,
				 const typename Field::Element alpha,
				 typename Field::Element * A, const size_t lda,
				 const MMHelper<AlgoT, FieldCategories::ModularFloatingPointTag, Field >& H)
	{
		if (!F.isOne(alpha) && !F.isMOne(alpha)){
			if (abs(alpha)*std::max(-H.Outmin, H.Outmax) > H.MaxStorableValue){
				finit (F, M, N, A, lda);
				fscalin (F, M, N, alpha, A, lda);
			} else {
				fscalin (H.delayedField, M, N, alpha, A, lda);
				finit(F, M, N, A, lda);
			}
		} else
			finit(F, M, N, A, lda);
	}

} // Protected
} // FFLAS

namespace FFLAS {
	
	template<class Field> 
	inline  typename Field::Element*
	fgemm (const Field& F, 
	       const FFLAS_TRANSPOSE ta, 
	       const FFLAS_TRANSPOSE tb, 
	       const size_t m, const size_t n, const size_t k, 
	       const typename Field::Element alpha, 
	       typename Field::Element * A, const size_t lda, 
	       typename Field::Element * B, const size_t ldb, 
	       const typename Field::Element beta, 
	       typename Field::Element * C, const size_t ldc, 
	       MMHelper<MMHelperCategories::Winograd, FieldCategories::FloatingPointConvertibleTag, Field> & H)
	{ 
		if (F.characteristic() < DOUBLE_TO_FLOAT_CROSSOVER)
			return Protected::fgemm_convert<float,Field>(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H.recLevel);
		else
			return Protected::fgemm_convert<double,Field>(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H.recLevel);
	}
}// FFLAS

// #include "fflas_fgemm/matmul_algos.inl"
#include "fflas_fgemm/fgemm_classical.inl"
#include "fflas_fgemm/fgemm_winograd.inl"
// #include "fflas_fgemm/gemm_bini.inl"

// fgemm
namespace FFLAS {

	template<class Field>
	inline typename Field::Element*
	fgemm( const Field& F,
	       const FFLAS_TRANSPOSE ta,
	       const FFLAS_TRANSPOSE tb,
	       const size_t m,
	       const size_t n,
	       const size_t k,
	       const typename Field::Element alpha,
	       const typename Field::Element* A, const size_t lda,
	       const typename Field::Element* B, const size_t ldb,
	       const typename Field::Element beta,
	       typename Field::Element* C, const size_t ldc, const int w=-1)
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
		typename Field::Element alpha_,beta_;
		F.assign (alpha_,alpha);
		F.assign (beta_,beta);
		if (Protected::AreEqual<Field, FFPACK::Modular<double> >::value ||
		    Protected::AreEqual<Field, FFPACK::ModularBalanced<double> >::value){
			    // Modular<double> need to switch to float if p too small	
			if (F.characteristic() < DOUBLE_TO_FLOAT_CROSSOVER)
				return Protected::fgemm_convert<float,Field>(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,w);
		}
		if (Protected::AreEqual<typename FieldTraits<Field>::value,FieldCategories::ModularFloatingPointTag>::value){
			if ( !F.isOne(alpha) && !F.isMOne(alpha)){
				F.assign (alpha_, F.one);
				F.div (beta_, beta, alpha);
			}
		}

		typedef MMHelper<MMHelperCategories::Winograd, typename FieldTraits<Field>::value, Field > MMH_t;
		MMH_t H;
		if (w==-1) // Auto set Winograd rec. level
			H = MMH_t (F, m, k, n);
		else
			H = MMH_t (F, w);

		fgemm (F, ta, tb, m, n, k, alpha_, 
		       const_cast<typename Field::Element*>(A), lda, 
		       const_cast<typename Field::Element*>(B), ldb, 
		       beta_, C, ldc, H);
		
		Protected::ScalAndInit (F, m, n, alpha, C, ldc, H);
		// if (Protected::AreEqual<typename FieldTraits<Field>::value,FieldCategories::ModularFloatingPointTag>::value){
		// 	if ( !F.isOne(alpha) && !F.isMOne(alpha)){
		// 		if (abs(alpha)*std::max(-H.Outmin, H.Outmax)>H.MaxStorableValue){
		// 			finit(F,m,n,C,ldc);
		// 			fscalin(F, m,n,alpha,C,ldc);
		// 		} else {
		// 			fscalin(H.delayedField, m,n,alpha,C,ldc);
		// 			finit(F,m,n,C,ldc);
		// 		}
		// 		return C;
		// 	}
		// }
		// finit(F,m,n,C,ldc);
		return C;
	}

} // FFLAS


// fsquare
namespace FFLAS {
	template < class Field >
	inline typename Field::Element*
	fsquare (const Field& F,
		 const FFLAS_TRANSPOSE ta,
		 const size_t n, const typename Field::Element alpha,
		  typename Field::Element* A, const size_t lda,
		 const typename Field::Element beta,
		 typename Field::Element* C, const size_t ldc)
	{

		double alphad, betad;
		F.convert (alphad, alpha);
		if (F.isMOne (beta))
			betad = -1.0;
		else
			F.convert (betad, beta);

		//! @bug why double ?
		// Double  matrices initialisation
		DoubleDomain::Element * Ad = new DoubleDomain::Element[n*n];
		DoubleDomain::Element * Cd = new DoubleDomain::Element[n*n];
		// Conversion finite Field = >  double
		fconvert (F, n, n, Ad, n, A, lda);
		if (!F.isZero(beta)) fconvert(F, n, n, Cd, n, C, ldc);

		// Call to the blas Multiplication
		FFLASFFPACK_check(n);
		cblas_dgemm (CblasRowMajor, (CBLAS_TRANSPOSE)ta,
			     (CBLAS_TRANSPOSE)ta, (int)n, (int)n, (int)n,
			     (DoubleDomain::Element) alphad, Ad, (int)n, Ad, (int)n,
			     (DoubleDomain::Element) betad, Cd, (int)n);
		// Conversion double = >  Finite Field
		delete[] Ad;
		finit (F,n,n, Cd, n, C, ldc);
		delete[] Cd;
		return C;
	}

	namespace Protected {

		// F is Modular(Balanced)<float/double>
		template < class Field >
		inline typename Field::Element*
		fsquareCommon (const Field& F,
			       const FFLAS_TRANSPOSE ta,
			       const size_t n, const typename Field::Element alpha,
			       typename Field::Element* A, const size_t lda,
			       const typename Field::Element beta,
			       typename Field::Element* C, const size_t ldc)
		{
			typedef typename Field::Element Element ; // double or float
			if (C==A) {
				Element * Ad = new Element[n*n];
				fcopy(F,n,n,A,lda,Ad,n);
				fgemm (F, ta, ta, n, n, n, alpha, Ad, n, Ad, n, beta, C, ldc);
				delete[] Ad;
			}
			else
				fgemm (F, ta, ta, n, n, n, alpha, A, lda, A, lda, beta, C, ldc);
			// Conversion double/float = >  Finite Field
			finit(F,n,n,C,ldc);
			return C;

		}

	} // Protected

	template <>
	inline double* fsquare (const  FFPACK:: ModularBalanced<double> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const double alpha,
				 double* A, const size_t lda,
				const double beta,
				double* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

	template <>
	inline float * fsquare (const  FFPACK:: ModularBalanced<float> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const float alpha,
				 float* A, const size_t lda,
				const float beta,
				float* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

	template <>
	inline double* fsquare (const  FFPACK:: Modular<double> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const double alpha,
				 double* A, const size_t lda,
				const double beta,
				double* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

	template <>
	inline float * fsquare (const  FFPACK:: Modular<float> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const float alpha,
				 float* A, const size_t lda,
				const float beta,
				float* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

} // FFLAS

#endif // __FFLASFFPACK_fgemm_INL
