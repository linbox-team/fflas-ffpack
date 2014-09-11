/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* fflas/fflas_finit.inl
 * Copyright (C) 2014 Pascal Giorgi
 *
 * Written by Pascal Giorgi <Pascal.Giorgi@lirmm.fr>
 * BB<bboyer@ncsu.edu>
 * Bastien Vialla <Bastien.Vialla@lirmm.fr>
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

#ifndef __FFLASFFPACK_fflas_simd_functions_H
#define __FFLASFFPACK_fflas_simd_functions_H

#ifdef __FFLASFFPACK_HAVE_CXX11

#include <type_traits>

#if defined(__FFLASFFPACK_USE_SSE) or defined(__FFLASFFPACK_USE_AVX) or defined(__FFLASFFPACK_USE_AVX2)
#define __FFLASFFPACK_USE_SIMD // see configure...
#endif

#ifdef __FFLASFFPACK_USE_SIMD
#include "fflas-ffpack/utils/simd.h"
#endif

namespace FFLAS { namespace vectorised {

	template<class T>
	T monfmod(T A, T B)
	{
		return fmod(A,B);
	}

	template<>
	float monfmod(float A, float B)
	{
		return fmodf(A,B);
	}

	template<class T>
	T monrint(T A)
	{
		return rint(A);
	}

	template<>
	float monrint(float A)
	{
		return rintf(A);
	}

#ifdef __FFLASFFPACK_USE_SIMD

	template<class SimdT>
	inline typename std::enable_if<is_simd<SimdT>::value, void>::type
	VEC_MOD(SimdT & C, SimdT & Q, SimdT & T, const SimdT & P, const SimdT & NEGP, const SimdT & INVP, const SimdT & MIN,  const SimdT & MAX)
	{
		using simd = Simd<typename simdToType<SimdT>::type>;
		Q = simd::mul(C, INVP);
		Q = simd::floor(Q);
		C = simd::nmadd(C,Q,P);
		Q = simd::greater(C, MAX);
		T = simd::lesser(C, MIN);
		Q = simd::vand(Q, NEGP);
		T = simd::vand(T, P);
		Q = simd::vor(Q, T);
		C = simd::add(C, Q);
	}

	template<bool positive, bool round, class Element, class T1, class T2>
	inline typename std::enable_if<std::is_floating_point<Element>::value, void>::type
	modp(Element * T, const Element * U, size_t n, Element p, Element invp, T1 min_, T2 max_)
	{
		Element min = (Element)min_, max = (Element)max_;
		using simd = Simd<Element>;
		using vect_t = typename simd::vect_t;

		size_t i = 0;
		if (n < simd::vect_size)
		{
			for (; i < n ; i++)
			{
				if (round)
				{
					T[i] = monrint(U[i]);
					T[i] = monfmod(T[i],p);
				}
				else
				{
					T[i]=monfmod(U[i],p);
				}
				if (!positive)
				{
					T[i]-=(T[i]>max)?p:0;
				}
				T[i]+=(T[i]<min)?p:0;
			}
			return;
		}

		long st = long(T) % simd::alignment;

		// the array T is not 32 byte aligned (process few elements s.t. (T+i) is 32 bytes aligned)

		if (st)
		{
			for (size_t j = static_cast<size_t>(st) ; j < simd::alignment ; j += sizeof(Element), i++)
			{
				if (round)
				{
					T[i] = monrint(U[i]);
					T[i] = monfmod(T[i],p);
				}
				else
				{
					T[i] = monfmod(U[i],p);
				}
				if (!positive)
				{
					T[i] -= (T[i] > max) ? p : 0;
				}
				T[i] += (T[i] < min) ? p : 0;
			}
		}

		FFLASFFPACK_check((long(T+i) % simd::alignment == 0));

		vect_t C, Q, TMP;
		vect_t P = simd::set1(p);
		vect_t NEGP = simd::set1(-p);
		vect_t INVP = simd::set1(invp);
		vect_t MIN = simd::set1(min);
		vect_t MAX = simd::set1(max);

		if((long(U+i) % simd::alignment == 0))
		{
			// perform the loop using 256 bits SIMD
			for (; i<= n - simd::vect_size ; i += simd::vect_size)
			{
				C = simd::load(U + i);

				if (round)
				{
					C = simd::round(C);
				}

				VEC_MOD(C,Q,TMP, P, NEGP,INVP,MIN,MAX);
				simd::store(T+i, C);
			}
		}

		// perform the last elt from T without SIMD
		for (;i<n;i++)
		{
			if (round)
			{
				T[i] = monrint(U[i]);
				T[i] = monfmod(T[i],p);
			}
			else
			{
				T[i] = monfmod(U[i],p);
			}
			if (!positive)
			{
				T[i] -= (T[i] > max) ? p : 0;
			}
			T[i] += (T[i] < min) ? p : 0;
		}
	}

	template<class SimdT>
	inline typename std::enable_if<is_simd<SimdT>::value, void>::type
	VEC_ADD(SimdT & C, SimdT & A, SimdT & B, SimdT & Q, SimdT & T, SimdT & P, SimdT & NEGP, SimdT & MIN, SimdT & MAX)
	{
		using simd = Simd<typename simdToType<SimdT>::type>;
		C = simd::add(A, B);
		Q = simd::greater(C, MAX);
		T = simd::lesser(C, MIN);
		Q = simd::vand(Q, NEGP);
		T = simd::vand(T, P);
		Q = simd::vor(Q, T);
		C = simd::add(C, Q);
	}

	template<bool positive, class Element, class T1, class T2>
	inline typename std::enable_if<std::is_floating_point<Element>::value, void>::type
	addp(Element * T, const Element * TA, const Element * TB,  size_t n,  Element p,  T1 min_,  T2 max_)
	{
		Element min= (Element)min_, max= (Element)max_;
		using simd = Simd<Element>;
		using vect_t = typename simd::vect_t;

		size_t i = 0;

		if (n < simd::vect_size)
		{
			for (; i < n ; i++)
			{
				T[i] = TA[i] + TB[i];
				T[i] -= (T[i] > max) ? p : 0;
				if (!positive)
				{
					T[i] += (T[i] < min) ? p : 0;
				}
			}
			return;

		}
		vect_t A,B,C,Q,P,NEGP,TMP,MIN,MAX;
		P   = simd::set1(p);
		NEGP= simd::set1(-p);
		MIN = simd::set1(min);
		MAX = simd::set1(max);
		long st = long(T)%simd::alignment;
		if (st)
		{ // the array T is not 32 byte aligned (process few elements s.t. (T+i) is 32 bytes aligned)
			for (size_t j=static_cast<size_t>(st) ; j < simd::alignment ; j += sizeof(Element), i++)
			{
				T[i] = TA[i] + TB[i];
				T[i] -= (T[i] > max) ? p : 0;
				if (!positive)
					T[i] += (T[i] < min) ? p : 0;
			}
		}
		FFLASFFPACK_check((long(T+i) % simd::alignment == 0));
		if ( (long(TA+i)%simd::alignment==0) && (long(TB+i)%simd::alignment==0))
		{
			// perform the loop using 256 bits SIMD
			for (; i <= n - simd::vect_size ; i += simd::vect_size)
			{
				C = simd::load(T+i);
				A = simd::load(TA+i);
				B = simd::load(TB+i);
				VEC_ADD(C, A, B, Q, TMP, P, NEGP, MIN, MAX);
				simd::store(T+i, C);
			}
		}
		// perform the last elt from T without SIMD
		for (; i < n ; i++)
		{
			T[i] = TA[i] + TB[i];
			T[i] -= (T[i] > max) ? p : 0;
			if (!positive)
				T[i] += (T[i] < min) ? p : 0;
		}
	}

	template<class SimdT>
	inline typename std::enable_if<is_simd<SimdT>::value, void>::type
	VEC_SUB(SimdT & C, SimdT & A, SimdT & B, SimdT & Q, SimdT & T, SimdT & P, SimdT & NEGP, SimdT & MIN, SimdT & MAX)
	{
		using simd = Simd<typename simdToType<SimdT>::type>;
		C = simd::sub(A, B);
		Q = simd::greater(C, MAX);
		T = simd::lesser(C, MIN);
		Q = simd::vand(Q, NEGP);
		T = simd::vand(T, P);
		Q = simd::vor(Q, T);
		C = simd::add(C, Q);
	}

	template<bool positive, class Element, class T1, class T2>
	inline typename std::enable_if<std::is_floating_point<Element>::value, void>::type
	subp(Element * T, const Element * TA, const Element * TB, const size_t n, const Element p, const T1 min_, const T2 max_)
	{
		Element min = (Element)min_, max = (Element)max_;
		using simd = Simd<Element>;
		using vect_t = typename simd::vect_t;

		size_t i = 0;

		if (n < simd::vect_size)
		{
			for (; i < n ; i++)
			{
				T[i] = TA[i] - TB[i];
				if (!positive)
					T[i] -= (T[i] > max) ? p : 0;
				T[i] += (T[i] < min) ? p : 0;
			}
			return;

		}
		vect_t A,B,C,Q,P,NEGP,TMP,MIN,MAX;
		P   = simd::set1(p);
		NEGP= simd::set1(-p);
		MIN = simd::set1(min);
		MAX = simd::set1(max);
		long st = long(T) % simd::alignment;
		if (st)
		{ // the array T is not 32 byte aligned (process few elements s.t. (T+i) is 32 bytes aligned)
			for (size_t j = static_cast<size_t>(st) ; j < simd::alignment ; j += sizeof(Element), i++)
			{
				T[i] = TA[i] - TB[i];
				if (!positive)
					T[i] -= (T[i] > max) ? p : 0;
				T[i] += (T[i] < min) ? p : 0;
			}
		}
		FFLASFFPACK_check((long(T+i) % simd::alignment == 0));
		if ( (long(TA+i) % simd::alignment == 0) && (long(TB+i) % simd::alignment == 0))
		{
			// perform the loop using 256 bits SIMD
			for (; i <= n - simd::vect_size ; i += simd::vect_size)
			{
				C = simd::load(T+i);
				A = simd::load(TA+i);
				B = simd::load(TB+i);
				VEC_SUB(C, A, B, Q, TMP, P, NEGP, MIN, MAX);
				simd::store(T+i, C);
			}
		}

		// perform the last elt from T without SIMD
		for (; i < n ; i++)
		{
			T[i] = TA[i] - TB[i];
			if (!positive)
				T[i] -= (T[i] > max) ? p : 0;
			T[i] += (T[i] < min) ? p : 0;
		}
	}

	template<class SimdT>
	inline typename std::enable_if<is_simd<SimdT>::value, void>::type
	VEC_SCAL(SimdT & C, SimdT & ALPHA, SimdT & Q, SimdT & T, SimdT & P, SimdT & NEGP, SimdT & INVP, SimdT & MIN, SimdT & MAX)
	{
		using simd = Simd<typename simdToType<SimdT>::type>;
		Q = simd::mul(C,INVP);
		C = simd::mul(C,ALPHA);
		Q = simd::floor(Q);
		C = simd::nmadd(C,Q,P);
		Q = simd::greater(C,MAX);
		T = simd::lesser(C,MIN);
		Q = simd::vand(Q,NEGP);
		T = simd::vand(T,P);
		Q = simd::vor(Q,T);
		C = simd::add(C,Q);
	}

	template<class Element, class T1, class T2>
	inline typename std::enable_if<std::is_floating_point<Element>::value, void>::type
	scalp(Element *T, const Element alpha, const Element * U, size_t n, Element p, Element invp, T1 min_, T2 max_)
	{
		Element min = (Element)min_, max=(Element)max_;
		using simd = Simd<Element>;
		using vect_t = typename simd::vect_t;

		size_t i = 0;

		if (n < simd::vect_size)
		{
			for (; i < n ; i++)
			{
				T[i]=monfmod(alpha*U[i], p);
				T[i] -= (T[i] > max) ? p : 0;
				T[i] += (T[i] < min) ? p : 0;
			}
			return;

		}
		vect_t C,Q,P,NEGP,INVP,TMP,MIN,MAX,ALPHA;
		ALPHA = simd::set1(alpha);
		P   = simd::set1(p);
		NEGP = simd::set1(-p);
		INVP = simd::set1(invp);
		MIN = simd::set1(min);
		MAX = simd::set1(max);
		long st = long(T) % simd::alignment;
		if (st)
		{ // the array T is not 32 byte aligned (process few elements s.t. (T+i) is 32 bytes aligned)
			for (size_t j = static_cast<size_t>(st) ; j < simd::alignment ; j+=sizeof(Element), i++)
			{
				T[i] = monfmod(alpha*U[i], p);
				T[i] -= (T[i] > max) ? p : 0;
				T[i] += (T[i] < min) ? p : 0;
			}
		}
		FFLASFFPACK_check((long(T+i) % simd::alignment == 0));
		if ((long(U+i)%simd::alignment==0))
		{
			// perform the loop using 256 bits SIMD
			for (;i <= n - simd::vect_size ; i += simd::vect_size)
			{
				C = simd::load(U+i);
				VEC_SCAL(C, ALPHA, Q, TMP, P, NEGP, INVP, MIN, MAX);
				simd::store(T+i,C);
			}
		}
		// perform the last elt from T without SIMD
		for (; i < n ; i++)
		{
			T[i] = monfmod(alpha*U[i],p);
			T[i] -= (T[i] > max) ? p : 0;
			T[i] += (T[i] < min) ? p : 0;
		}
	}
#else // not AVX

	template<bool positive, bool round, class Element, class T1, class T2>
	void
	modp(Element * T, const Element * U, size_t n, Element p, Element invp, T1 min_, T2 max_)
	{
		Element min = (Element)min_, max = (Element)max_;

		size_t i = 0;
		for (; i < n ; i++)
		{
			if (round)
			{
				T[i] = monrint(U[i]);
				T[i] = monfmod(T[i],p);
			}
			else
			{
				T[i]=monfmod(U[i],p);
			}
			if (!positive)
			{
				T[i]-=(T[i]>max)?p:0;
			}
			T[i]+=(T[i]<min)?p:0;
		}
	}

	template<bool positive, class Element, class T1, class T2>
	void
	addp(Element * T, const Element * TA, const Element * TB,  size_t n,  Element p,  T1 min_,  T2 max_)
	{
		Element min= (Element)min_, max= (Element)max_;
		size_t i = 0;
		for (; i < n ; i++)
		{
			T[i] = TA[i] + TB[i];
			T[i] -= (T[i] > max) ? p : 0;
			if (!positive)
			{
				T[i] += (T[i] < min) ? p : 0;
			}
		}
		return;

	}

	template<bool positive, class Element, class T1, class T2>
	void
	subp(Element * T, const Element * TA, const Element * TB, const size_t n, const Element p, const T1 min_, const T2 max_)
	{
		Element min = (Element)min_, max = (Element)max_;
		size_t i = 0;
		for (; i < n ; i++)
		{
			T[i] = TA[i] - TB[i];
			if (!positive)
				T[i] -= (T[i] > max) ? p : 0;
			T[i] += (T[i] < min) ? p : 0;
		}
		return;

	}

	template<class Element, class T1, class T2>
	void
	scalp(Element *T, const Element alpha, const Element * U, size_t n, Element p, Element invp, T1 min_, T2 max_)
	{
		Element min = (Element)min_, max=(Element)max_;

		size_t i = 0;

		{
			for (; i < n ; i++)
			{
				T[i]=monfmod(alpha*U[i], p);
				T[i] -= (T[i] > max) ? p : 0;
				T[i] += (T[i] < min) ? p : 0;
			}
			return;

		}

	}

#endif // AVX
} /* vectorised */
} /* FFLAS */

#else /* C++11 */
#error "You need a c++11 compiler."
#endif /* c++11 */

#endif /* __FFLASFFPACK_fflas_simd_functions_H */
