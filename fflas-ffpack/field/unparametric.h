/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* field/unparametric.h
 * Copyright (C) 1999-2005 William J Turner,
 *               2001 Bradford Hovinen
 *               2005 Clement Pernet
 *
 * Written by W. J. Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 * Modified By C. Pernet and inserted into Fflas_Ffpack
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __FIELD_UNPARAMETRIC_H
#define __FIELD_UNPARAMETRIC_H

#include <string>
#include <algorithm>

template< class K>
class UnparametricField;

template<class K>
class UnparametricField
{
public:

	typedef K Element;


	UnparametricField(){}
	~UnparametricField () {}
    	UnparametricField &operator=(const UnparametricField &F) { return *this; }

	Element &init (Element &x, const unsigned long &y=0) const
	{ return x = y; }

	Element &init (Element &x, const double &y=0) const
	{ return x = y; }
	Element &init (Element &x, const float &y=0) const
	{ return x = y; }


	unsigned long convert (unsigned long &x, const Element &y) const
	{
		return x = (unsigned long) y;
	}


	K &convert (K &x, const Element &y) const { return x = y; }

	Element &assign (Element &x, const Element &y) const { return x = y; }

	unsigned long &cardinality (unsigned long &c) const { return c = 0; }

        unsigned long &characteristic (unsigned long &c) const { return c = 0; }

	//  x == y
	bool areEqual (const Element &x, const Element &y) const { return x == y; }

	//  x == 0
	bool isZero (const Element &x) const { return x == Element (0); }

	///  x == 1
	bool isOne (const Element &x) const { return x == Element (1); }

	// x := y + z
	Element &add (Element &x, const Element &y, const Element &z) const
	{ return x = y + z; }
    	// x := y - z
	Element &sub (Element &x, const Element &y, const Element &z) const
	{ return x = y - z; }
    	// x := y*z
	Element &mul (Element &x, const Element &y, const Element &z) const
	{ return x = y * z; }
	/// x := y/z
	Element &div (Element &x, const Element &y, const Element &z) const
	{ return x = y / z; }
    	/// x := -y
	Element &neg (Element &x, const Element &y) const { return x = - y; }
    	/// x := 1/y
	Element &inv (Element &x, const Element &y) const
	{ return x = Element (1) / y; }
    	/// z := a*x + y
	Element &axpy (Element &z, const Element &a, const Element &x,
		       const Element &y) const
	{ return z = a * x + y; }
 	// x := x + y
	Element &addin (Element &x, const Element &y) const { return x += y; }
 	// x := x - y
	Element &subin (Element &x, const Element &y) const { return x -= y; }
	/// x := x*y
	Element &mulin (Element &x, const Element &y) const { return x *= y; }
	/// x := x/y
	Element &divin (Element &x, const Element &y) const { return x /= y; }
	/// x := -x
	Element &negin (Element &x) const { return x = - x; }
	/// x := 1/x
	Element &invin (Element &x) const { return x = Element (1) / x; }
	/// y := a*x + y
	Element &axpyin (Element &y, const Element &a, const Element &x) const
	{ return y += a * x; }


	std::ostream &write (std::ostream &os) const { return os << "unparamterized field"; }
    	std::istream &read (std::istream &is) const { return is; }

	std::ostream &write (std::ostream &os, const Element &x) const
	{ return os << x; }
	std::istream &read (std::istream &is, Element &x) const { return is >> x; }

};
#endif // __FIELD_UNPARAMETRIC_H_
