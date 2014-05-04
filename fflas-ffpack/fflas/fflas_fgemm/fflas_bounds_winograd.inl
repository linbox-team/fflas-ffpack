/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* fflas/fflas_bounds.inl
 * Copyright (C) 2014 FFLAS-FFPACK group
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

#ifndef __FFLASFFPACK_fflas_bounds_winograd_INL
#define __FFLASFFPACK_fflas_bounds_winograd_INL

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas_bounds_classic.inl"

// WinogradSteps
namespace FFLAS {

	/** \brief Computes the number of recursive levels to perform.
	 *
	 * \param m the common dimension in the product AxB
	 */
	template<class Field>
	inline int WinogradSteps (const Field & F, const size_t & m)
	{
		int w = 0;
		size_t mt = m;
		while ( mt >= WINOTHRESHOLD ) {
			++w;
			mt >>= 1;
		}
		return w;
	}

	template<>
	inline int WinogradSteps (const FFPACK:: Modular<double> & F, const size_t & m)
	{
		int w = 0;
		size_t mt = m;
		while ( mt >= __FFLASFFPACK_WINOTHRESHOLD ) {
			++w;
			mt >>= 1;
		}
		return w;
	}

	template<>
	inline int WinogradSteps (const FFPACK:: ModularBalanced<double> & F, const size_t & m)
	{
		int w = 0;
		size_t mt = m;
		while ( mt >= __FFLASFFPACK_WINOTHRESHOLD_BAL ) {
			++w;
			mt >>= 1;
		}
		return w;
	}

	template<>
	inline int WinogradSteps (const FFPACK:: Modular<float> & F, const size_t & m)
	{
		int w = 0;
		size_t mt = m;
		while ( mt >= __FFLASFFPACK_WINOTHRESHOLD_FLT ) {
			++w;
			mt >>= 1;
		}
		return w;
	}

	template<>
	inline int WinogradSteps (const FFPACK:: ModularBalanced<float> & F, const size_t & m)
	{
		int w = 0;
		size_t mt = m;
		while ( mt >= __FFLASFFPACK_WINOTHRESHOLD_BAL_FLT ) {
			++w;
			mt >>= 1;
		}
		return w;
	}
} // FFLAS

// BaseCompute
namespace FFLAS { namespace Protected {
	/**  BaseCompute determines the type of floating point representation to
	 * convert to, for BLAS computations.
	 * \param F Finite Field/Ring of the computation
	 * \param w Number of recursive levels in Winograd's algorithm
	 *
	 */
	//! @bug these thresholds come from nowhere.
	template <class Field>
	inline FFLAS_BASE BaseCompute (const Field& F, const size_t w)
	{

		FFLAS_INT_TYPE pi=0;
		F.characteristic(pi);
		FFLAS_BASE base;
		switch (w) {
		case 0: base = (pi < FLOAT_DOUBLE_THRESHOLD_0)? FflasFloat : FflasDouble;
			break;
		case 1:  base = (pi < FLOAT_DOUBLE_THRESHOLD_1)? FflasFloat : FflasDouble;
			 break;
		case 2:  base = (pi < FLOAT_DOUBLE_THRESHOLD_2)? FflasFloat : FflasDouble;
			 break;
		default: base = FflasDouble;
			 break;
		}
		return base;
	}
#if 1
	template <>
	inline FFLAS_BASE BaseCompute (const FFPACK:: Modular<double>& ,
				       const size_t )
	{
		return FflasDouble;
	}

	template <>
	inline FFLAS_BASE BaseCompute (const FFPACK:: Modular<float>& ,
				       const size_t )
	{
		return FflasFloat;
	}

	template <>
	inline FFLAS_BASE BaseCompute (const FFPACK:: ModularBalanced<double>& ,
				       const size_t )
	{
		return FflasDouble;
	}

	template <>
	inline FFLAS_BASE BaseCompute (const FFPACK:: ModularBalanced<float>& ,
				       const size_t )
	{
		return FflasFloat;
	}
#endif

} // Protected
} // FFLAS


// computeFactorWinograd
namespace FFLAS { namespace Protected {
	/** @internal
	 * @brief Internal function for the bound computation
	 * Generic implementation for positive representations
	 */
	template <class Field>
	inline double computeFactorWinograd (const Field& F, const size_t w)
	{
		FFLAS_INT_TYPE p=0;
		F.characteristic(p);
		size_t ex=1;
		for (size_t i=0; i < w; ++i) 	ex *= 3;
		return double(p - 1) * double(1 + ex) / double(2);
	}

	/*************************************************************************************
	 * Specializations for ModularPositive and ModularBalanced over double and float
	 *************************************************************************************/

	template <>
	inline double computeFactorWinograd (const FFPACK:: ModularBalanced<double>& F, const size_t w)
	{
		FFLAS_INT_TYPE p;
		F.characteristic(p);
		size_t ex=1;
		for (size_t i=0; i < w; ++i) 	ex *= 3;
		return  double((p - 1) * ex / 2);
	}

	template <>
	inline double computeFactorWinograd (const FFPACK:: ModularBalanced<float>& F, const size_t w)
	{
		FFLAS_INT_TYPE p;
		F.characteristic(p);
		size_t ex=1;
		for (size_t i=0; i < w; ++i) 	ex *= 3;
		return  double((p - 1) * ex / 2);
	}


} // Protected
} // FFLAS

// DotProdBoundWinograd
namespace FFLAS { namespace Protected {
	/** DotProdBound computes the maximal size for delaying the modular reduction
	 * in a dotproduct.
	 *
	 * This is the default version assuming a conversion to a positive modular representation
	 *
	 * \param F Finite Field/Ring of the computation
	 * \param winogradRecLevel Number of recusrive Strassen-Winograd levels (if any, \p 0 otherwise)
	 * \param beta Computing <code>AB + beta C</code>
	 * \param base Type of floating point representation for delayed modular computations
	 *
	 */
	template <class Field>
	inline size_t DotProdBoundWinograd (const Field& F,
					    const size_t w,
					    const typename Field::Element& beta,
					    const FFLAS_BASE base)
	{

		double kmax;
		if (w > 0) {
			FFLAS_INT_TYPE p=0;
			F.characteristic(p);

			if (p == 0)
				return 1;

			unsigned long mantissa = Mantissa (F, base);
			//(base == FflasDouble) ? DBL_MANT_DIG : FLT_MANT_DIG;

			double c = computeFactorWinograd (F,w);

			double d = (double (1ULL << mantissa) /(c*c) + 1);
			if (d < 2)
				return 1;
			kmax = floor (d * double(1ULL << w));
			// if (kmax  <= 1) return 1;
		} else {
			// why call Winograd ?
			kmax = (double) DotProdBoundClassic(F,beta,base);
		}

		//kmax--; // we computed a strict upper bound
		return  (size_t) std::min ((unsigned long long)kmax, 1ULL << 31);
	}

} // Protected
} // FFLAS

// MatMulParameters (Winograd)
namespace FFLAS { namespace Protected {
	/**
	 * Computes the threshold parameters for the cascade
	 *        Matmul algorithm.
	 *
	 *
	 * \param F Finite Field/Ring of the computation.
	 * \param m row dim of A
	 * \param n col dim of B
	 * \param k Common dimension of A and B, in the product A x B
	 * \param beta Computing \f$AB + \beta C\f$
	 * \param delayedDim Returns the size of blocks that can be multiplied
	 *                   over Z with no overflow
	 * \param base Returns the type of BLAS representation to use
	 * \param winogradRecLevel Returns the number of recursion levels of
	 *                     Strassen-Winograd's algorithm to perform
	 * \param winogradLevelProvided tells whether the user forced the number of
	 *                          recursive level of Winograd's algorithm
	 *
	 * @bib
	 * - Dumas, Giorgi, Pernet, arXiv cs/0601133  <a href=http://arxiv.org/abs/cs.SC/0601133>here</a>
	 */
	template <class Field>
	inline void MatMulParametersWinograd (const Field& F,
				      const size_t &m,
				      const size_t &n,
				      const size_t &k,
				      const typename Field::Element& beta,
				      size_t& delayedDim,
				      FFLAS_BASE& base,
				      int& winogradRecLevel,
				      bool winogradLevelProvided=false)
	{

		// Strategy : determine Winograd's recursion first, then choose appropriate
		// floating point representation, and finally the blocking dimension.
		// Can be improved for some cases.

		if (!winogradLevelProvided)
			winogradRecLevel = WinogradSteps (F, min3(m,k,n));
		base = BaseCompute (F, winogradRecLevel);
		// std::cout << typeid(typename Field::Element).name() << "->" << ((base == FflasFloat)?"f":"d") << std::endl;
		delayedDim = DotProdBoundWinograd (F, winogradRecLevel, beta, base);

		size_t oldk = k;
		int winogradDel = winogradRecLevel;
		// Computes the delayedDim, only depending on the recursive levels
		// that must be performed over Z
		if (!winogradLevelProvided){ // Automatic mode
			while (winogradDel > 0 && delayedDim < oldk) {
				// If the topmost rec call is done with modular reductions
				// Then do one less rec step
				winogradDel--;
				// recompute the bound
				delayedDim = DotProdBoundWinograd (F, winogradDel, beta, base);
				if (winogradDel == 0) {
					typename Field::Element mbeta;
					F.neg (mbeta,beta);
					// beta only comes into play when winodel==0
					// Then -beta = p-beta is also used in some calls
					delayedDim = std::min(delayedDim, DotProdBoundClassic (F, mbeta, base));

				}
				oldk >>= 1;
			}
			winogradRecLevel = winogradDel;
		}
		delayedDim = std::min (k, delayedDim);

		//std::cerr<<"p = "<<F.characteristic()<<" k = "<<k<<" winogradRecLevel = "<<winogradRecLevel<<" -> delayedDim = "<<delayedDim<<std::endl;
	}

	template <>
	inline void MatMulParametersWinograd (const DoubleDomain & F,
				      const size_t &m,
				      const size_t &n,
				      const size_t &k,
				      const DoubleDomain::Element& beta,
				      size_t& delayedDim,
				      FFLAS_BASE& base,
				      int& winogradRecLevel,
				      bool winogradLevelProvided)
	{
		if (!winogradLevelProvided)
			winogradRecLevel = WinogradSteps (F, min3(m,k,n)) ;

		delayedDim = k+1;
		base = FflasDouble;
	}

	template <>
	inline void MatMulParametersWinograd (const FloatDomain & F,
				      const size_t &m,
				      const size_t &n,
				      const size_t &k,
				      const FloatDomain::Element& beta,
				      size_t& delayedDim,
				      FFLAS_BASE& base,
				      int& winogradRecLevel,
				      bool winogradLevelProvided)
	{
		if (!winogradLevelProvided)
			winogradRecLevel = WinogradSteps (F, min3(m,k,n)) ;

		delayedDim = k+1;
		base = FflasFloat;
	}


} // Protected
} // FFLAS

#endif // __FFLASFFPACK_fflas_bounds_winograd_INL
