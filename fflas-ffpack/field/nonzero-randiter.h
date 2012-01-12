/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas-ffpack/nonzero-randiter.h
 * Copyright (C) 2001-2002 Bradford Hovinen
 *               2008 Clement Pernet
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *            Clement Pernet <clement.pernet@gmail.com>
 *
 * taken for LinBox
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

#ifndef __NONZERO_RANDITER_H
#define __NONZERO_RANDITER_H

#include <sys/time.h>
#include <stdlib.h>

#include <string>

namespace FFPACK {
/** Random iterator for nonzero random numbers
 *
 * Wraps around an existing random iterator and ensures that the output
 * is entirely nonzero numbers.
 **/
template <class Field, class RandIter = typename Field::RandIter>
class NonzeroRandIter
{
public:

	typedef typename Field::Element Element;

	NonzeroRandIter (const Field &F, const RandIter &r)
		: _F (F), _r (r)
	{}

	NonzeroRandIter (const NonzeroRandIter& R)
		: _F (R._F), _r (R._r) {}

	~NonzeroRandIter()
	{}

	NonzeroRandIter& operator=(const NonzeroRandIter& R)
	{
		if (this != &R) { // guard against self-assignment
			_F = R._F;
			_r = R._r;
		}

		return *this;
	}

	Element &random (Element &a)  const
	{
		do _r.random (a); while (_F.isZero (a));
		return a;
	}

private:

	Field    _F;
	RandIter _r;

}; // class NonzeroRandIter

} // FFPACK

#endif // __NONZERO_RANDITER_H
