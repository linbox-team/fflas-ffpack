/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* field/modular-randiter.h
 * Copyright (C) 2008 Clement Pernet
 *
 * Written by Clement Pernet <clement.pernet@gmail.com>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __FFLAFLAS_modular_randiter_H
#define __FFLAFLAS_modular_randiter_H

#include <sys/time.h>
#include <stdlib.h>
#include <limits>

namespace FFPACK {
	template< class Element >
	class Modular;

	template <class Element>
	class ModularRandIter
	{
	public:
		ModularRandIter (const Modular<Element> &F) :
			_F(F)
		{
			struct timeval tp;
			gettimeofday(&tp, 0) ;
			long _seed = (long)(tp.tv_usec);
			srand48(_seed);
		}

		ModularRandIter (const ModularRandIter<Element> &R) :
			_F (R._F)
		{}

		/*! @bug not so random... (at all)  */
		Element &random (Element &a) const
		{
			return _F.init(a,(double)lrand48()-std::numeric_limits<long>::max());
		}

	private:
		Modular<Element> _F;

	};

	template <class T>
	class ModularBalanced;

	template <class Element>
	class ModularBalancedRandIter
	{
	public:
		ModularBalancedRandIter (const ModularBalanced<Element> &F):_F(F){}
		ModularBalancedRandIter (const ModularBalancedRandIter<Element> &R)
		: _F (R._F)
		{}
		Element &random (Element &a) const
		{
			return _F.init(a,(double)rand());
		}
	private:
		ModularBalanced<Element> _F;

	};

} // FFPACK

#endif // __FFLAFLAS_modular_randiter_H
