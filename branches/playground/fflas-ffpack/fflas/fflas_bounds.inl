/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* fflas/fflas_bounds.inl
 * Copyright (C) 2008 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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

#ifndef __FFLASFFPACK_fflas_bounds_INL
#define __FFLASFFPACK_fflas_bounds_INL

#define FFLAS_INT_TYPE long unsigned int

#include "fflas-ffpack/fflas-ffpack-config.h"

namespace FFLAS { namespace Protected {

		template<class T>
		unsigned long Mantissa () {return DBL_MANT_DIG;}
		template<>
		unsigned long Mantissa<float> () {return FLT_MANT_DIG;}
		// unsigned long Mantissa (const FFLAS_BASE base)
		// {return (base == FflasDouble) ? DBL_MANT_DIG : FLT_MANT_DIG;}


} // Protected

} // FFLAS

namespace FFLAS { namespace Protected {

	template <class Field>
	inline double computeFactorClassic (const Field& F)
	{
		FFLAS_INT_TYPE p=0;
		F.characteristic(p);
		return (double) (p-1);
	}

	/*************************************************************************************
	 * Specializations for ModularPositive and ModularBalanced over double and float
	 *************************************************************************************/
	template <>
	inline double computeFactorClassic (const FFPACK:: ModularBalanced<double>& F)
	{
		FFLAS_INT_TYPE p;
		F.characteristic(p);
		return double((p-1) >> 1);
	}

	//BB: ajout, pourquoi pas ?
	template <>
	inline double computeFactorClassic (const FFPACK:: ModularBalanced<float>& F)
	{
		FFLAS_INT_TYPE p;
		F.characteristic(p);
		return double((p-1) >> 1);
	}

	template <class Field>
	inline size_t DotProdBoundClassic (const Field& F,
					   const typename Field::Element& beta,
					   const FFLAS_BASE base)
	{

		FFLAS_INT_TYPE p=0;
		F.characteristic(p);

		unsigned long mantissa = Protected::Mantissa<typename Field::Element>();

		//(base == FflasDouble) ? DBL_MANT_DIG : FLT_MANT_DIG;

		if (p == 0)
			return 1;

		double kmax;
		{

			double c = computeFactorClassic(F);

			double cplt=0;
			if (!F.isZero (beta)){
				if (F.isOne (beta) || F.areEqual (beta, F.mOne)) cplt = c;
				else{
					double be;
					F.convert(be, beta);
					cplt = fabs(be)*c;
				}
			}
			kmax = floor ( (double (double(1ULL << mantissa) - cplt)) / (c*c));
			if (kmax  <= 1) return 1;
		}
		    //kmax--; // we computed a strict upper bound
		return  (size_t) std::min ((unsigned long long)kmax, 1ULL << 31);
	}


} // FFLAS
} // Protected

namespace FFLAS { namespace Protected {


	/**
	 * TRSMBound
	 *
	 * \brief  computes the maximal size for delaying the modular reduction
	 *         in a triangular system resolution
	 *
	 * This is the default version over an arbitrary field.
	 * It is currently never used (the recursive algorithm is run until n=1 in this case)
	 *
	 * \param F Finite Field/Ring of the computation
	 *
	 */
	template <class Field>
	inline size_t TRSMBound (const Field& F)
	{
		return 1;
	}

	// /**
	//  * Specialization for positive modular representation over double
	//  * Computes nmax s.t. (p-1)/2*(p^{nmax-1} + (p-2)^{nmax-1}) < 2^53
	//  * See [Dumas Giorgi Pernet 06, arXiv:cs/0601133]
	//  */
	// template<>
	// inline size_t TRSMBound (const FFPACK:: Modular<double>& F)
	// {

	// 	FFLAS_INT_TYPE pi;
	// 	F.characteristic(pi);
	// 	unsigned long p = pi;
	// 	unsigned long long p1(1UL), p2(1UL);
	// 	size_t nmax = 0;
	// 	unsigned long long max = ( (1ULL << (DBL_MANT_DIG + 1) ) / ((unsigned long long)(p - 1)));
	// 	while ( (p1 + p2) < max ){
	// 		p1*=p;
	// 		p2*=p-2;
	// 		nmax++;
	// 	}
	// 	return nmax;
	// }


	/**
	 * Specialization for positive modular representation over float.
	 * Computes nmax s.t. (p-1)/2*(p^{nmax-1} + (p-2)^{nmax-1}) < 2^24
	 * @pbi
	 * See [Dumas Giorgi Pernet 06, arXiv:cs/0601133]
	 */
	template<class Element>
	inline size_t TRSMBound (const FFPACK:: Modular<Element>& F)
	{

		FFLAS_INT_TYPE pi;
		F.characteristic(pi);
		double p = pi;
		double p1 = 1.0, p2 = 1.0;
		double pm1 = (p - 1) / 2;
		size_t nmax = 0;
		unsigned long long max = 1ULL << Protected::Mantissa<Element>();
		while ( (p1 + p2)*pm1 < max ){
			p1*=p;
			p2*=p-2;
			nmax++;
		}
		return nmax;
	}

	/**
	 * Specialization for balanced modular representation over double.
	 * Computes nmax s.t. (p-1)/2*(((p+1)/2)^{nmax-1}) < 2^53
	 * @bib
	 * - Dumas Giorgi Pernet 06, arXiv:cs/0601133
	 */
	template<class Element>
	inline size_t TRSMBound (const FFPACK:: ModularBalanced<Element>& F)
	{

		FFLAS_INT_TYPE pi;
		F.characteristic (pi);
		double pp1 = (pi + 1) / 2;
		double pm1 = (pi - 1) / 2;
		double p1 = 1.0;
		size_t nmax = 0;
		double max = 1ULL << Protected::Mantissa<Element>();
		while (pm1*p1 < max){
			p1 *= pp1;
			nmax++;
		}
		return 1;//nmax;
	}

	// /**
	//  * Specialization for balanced modular representation over float
	//  * Computes nmax s.t. (p-1)/2*(((p+1)/2)^{nmax-1}) < 2^24
	//  * See [Dumas Giorgi Pernet 06, arXiv:cs/0601133]
	//  */
	// template<>
	// inline size_t TRSMBound (const FFPACK:: ModularBalanced<float>& F)
	// {

	// 	FFLAS_INT_TYPE pi;
	// 	F.characteristic (pi);
	// 	unsigned long p = (pi + 1) / 2;
	// 	unsigned long long p1(1UL);
	// 	size_t nmax = 0;
	// 	unsigned long long max = (1ULL << (FLT_MANT_DIG + 1)) ;
	// 	while ((pi-1)*p1 < max){
	// 		p1 *= p;
	// 		nmax++;
	// 	}
	// 	return nmax;

	// }
} // Protected
} // FFLAS

#endif // __FFLASFFPACK_fflas_bounds_INL
