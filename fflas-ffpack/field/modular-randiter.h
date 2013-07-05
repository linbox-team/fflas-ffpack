/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* field/modular-randiter.h
 * Copyright (C) 2008 Clement Pernet
 *
 * Written by Clement Pernet <clement.pernet@gmail.com>
 *
 * ------------------------------------
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

#ifndef __FFLASFFPACK_modular_randiter_H
#define __FFLASFFPACK_modular_randiter_H

#include <sys/time.h>
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE
#endif
#include <stdlib.h>
#include <limits>

namespace FFPACK {
	template< class Element >
	class Modular;

	template <class Element>
	class ModularRandIter {
	public:
                ModularRandIter (const Modular<Element> &F, size_t seed=0) :
                        _F(F)
                {
                        if (seed==0) {
                                struct timeval tp;
                                gettimeofday(&tp, 0) ;
                                seed = (size_t) tp.tv_usec;
                        }
                        srand48((long)seed);
                }

		ModularRandIter (const ModularRandIter<Element> &R) :
			_F (R._F)
		{}

		/*! @bug not so random... (at all)  */
		Element &random (Element &a) const
		{
			return _F.init(a,(double)(lrand48()-std::numeric_limits<long>::max()));
		}

	private:
		Modular<Element> _F;

	};

	template <class T>
	class ModularBalanced;

	template <class Element>
	class ModularBalancedRandIter {
	public:
		ModularBalancedRandIter (const ModularBalanced<Element> &F):_F(F){}
		ModularBalancedRandIter (const ModularBalancedRandIter<Element> &R) :
			_F (R._F)
		{}
		Element &random (Element &a) const
		{
			return _F.init(a,(double)rand());
		}
	private:
		ModularBalanced<Element> _F;

	};

} // FFPACK

#endif // __FFLASFFPACK_modular_randiter_H
