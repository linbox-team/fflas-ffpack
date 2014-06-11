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

namespace FFLAS {

	template <class T>
	const T& min3(const T & m, const T & n , const T & k)	{return std::min(m,std::min(n,k));}

	template <class T>
	const T& max3(const T & m, const T & n , const T & k)	{return std::max(m,std::min(n,k));}

	template <class T>
	const T& min4(const T & m, const T & n , const T & k, const T & l)
	{return std::min(std::min(m,n),std::min(k,l));}

	template <class T>
	const T& max4(const T & m, const T & n , const T & k, const T & l)
	{return std::max(std::max(m,n),std::max(k,l));}


} // FFLAS


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

	/**
	 * Specialization for positive modular representation over double
	 * Computes nmax s.t. (p-1)/2*(p^{nmax-1} + (p-2)^{nmax-1}) < 2^53
	 * See [Dumas Giorgi Pernet 06, arXiv:cs/0601133]
	 */
	template<>
	inline size_t TRSMBound (const FFPACK:: Modular<double>& F)
	{

		FFLAS_INT_TYPE pi;
		F.characteristic(pi);
		unsigned long p = pi;
		unsigned long long p1(1UL), p2(1UL);
		size_t nmax = 0;
		unsigned long long max = ( (1ULL << (DBL_MANT_DIG + 1) ) / ((unsigned long long)(p - 1)));
		while ( (p1 + p2) < max ){
			p1*=p;
			p2*=p-2;
			nmax++;
		}
		return nmax;
	}


	/**
	 * Specialization for positive modular representation over float.
	 * Computes nmax s.t. (p-1)/2*(p^{nmax-1} + (p-2)^{nmax-1}) < 2^24
	 * @pbi
	 * See [Dumas Giorgi Pernet 06, arXiv:cs/0601133]
	 */
	template<>
	inline size_t TRSMBound (const FFPACK:: Modular<float>& F)
	{

		FFLAS_INT_TYPE pi;
		F.characteristic(pi);
		unsigned long p = pi;
		unsigned long long p1(1UL), p2(1UL);
		size_t nmax = 0;
		unsigned long long max = ( (1ULL << (FLT_MANT_DIG + 1) ) / ((unsigned long long)(p - 1)));
		while ( (p1 + p2) < max ){
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
	template<>
	inline size_t TRSMBound (const FFPACK:: ModularBalanced<double>& F)
	{

		FFLAS_INT_TYPE pi;
		F.characteristic (pi);
		unsigned long p = (pi + 1) / 2;
		unsigned long long p1(1UL);
		size_t nmax = 0;
		unsigned long long max = ((1ULL << (DBL_MANT_DIG + 1)) / ((unsigned long long)(p - 1)));
		while (p1 < max){
			p1 *= p;
			nmax++;
		}
		return nmax;
	}

	/**
	 * Specialization for balanced modular representation over float
	 * Computes nmax s.t. (p-1)/2*(((p+1)/2)^{nmax-1}) < 2^24
	 * See [Dumas Giorgi Pernet 06, arXiv:cs/0601133]
	 */
	template<>
	inline size_t TRSMBound (const FFPACK:: ModularBalanced<float>& F)
	{

		FFLAS_INT_TYPE pi;
		F.characteristic (pi);
		unsigned long p = (pi + 1) / 2;
		unsigned long long p1(1UL);
		size_t nmax = 0;
		unsigned long long max = ((1ULL << (FLT_MANT_DIG + 1)) / ((unsigned long long) (p - 1)));
		while (p1 < max){
			p1 *= p;
			nmax++;
		}
		return nmax;

	}
} // Protected
} // FFLAS

#endif // __FFLASFFPACK_fflas_bounds_INL
