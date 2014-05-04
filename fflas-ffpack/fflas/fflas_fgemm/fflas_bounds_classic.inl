/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
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

#ifndef __FFLASFFPACK_fflas_bounds_classic_INL
#define __FFLASFFPACK_fflas_bounds_classic_INL


//  computeFactorClassic
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


} // FFLAS
} // Protected

// DotProdBoundClassic
namespace FFLAS { namespace Protected {

	template <class Field>
		inline size_t DotProdBoundClassic (const Field& F,
				// const size_t w,
				const typename Field::Element& beta,
				const FFLAS_BASE base)
		{

			FFLAS_INT_TYPE p=0;
			F.characteristic(p);

			unsigned long mantissa = Mantissa (F, base);

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

// MatMulParameters (Classic)
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
		inline void MatMulParametersClassic (const Field& F,
				const size_t &m,
				const size_t &n,
				const size_t &k,
				const typename Field::Element& beta,
				size_t& delayedDim,
				FFLAS_BASE& base)
		{

			// Strategy : determine Winograd's recursion first, then choose appropriate
			// floating point representation, and finally the blocking dimension.
			// Can be improved for some cases.

			base = BaseCompute (F, 0);
			// std::cout << typeid(typename Field::Element).name() << "->" << ((base == FflasFloat)?"f":"d") << std::endl;
			delayedDim = DotProdBoundClassic (F, beta, base);

			size_t oldk = k;
			// Computes the delayedDim, only depending on the recursive levels
			// that must be performed over Z
			delayedDim = std::min (k, delayedDim);

			//std::cerr<<"p = "<<F.characteristic()<<" k = "<<k<<" winogradRecLevel = "<<winogradRecLevel<<" -> delayedDim = "<<delayedDim<<std::endl;
		}

	template <>
		inline void MatMulParametersClassic (const DoubleDomain & F,
				const size_t &m,
				const size_t &n,
				const size_t &k,
				const DoubleDomain::Element& beta,
				size_t& delayedDim,
				FFLAS_BASE& base)
		{

			delayedDim = k+1;
			base = FflasDouble;
		}

	template <>
		inline void MatMulParametersClassic (const FloatDomain & F,
				const size_t &m,
				const size_t &n,
				const size_t &k,
				const FloatDomain::Element& beta,
				size_t& delayedDim,
				FFLAS_BASE& base)
		{

			delayedDim = k+1;
			base = FflasFloat;
		}


} // Protected
} // FFLAS



#endif // __FFLASFFPACK_fflas_bounds_classic_INL


