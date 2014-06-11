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


