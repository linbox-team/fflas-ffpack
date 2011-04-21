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

#ifndef __FFLAFLAS_field_unparametric_H
#define __FFLAFLAS_field_unparametric_H

#include <iostream> // std::cout
#include <string>
#include <algorithm>


namespace FFPACK
{

	template <typename Target, typename Source>
	Target& Caster (Target& t, const Source& s)
	{
		return t = static_cast<Target>(s);
	}

	/** \brief Unparameterized field adapter.
	 * \ingroup field
	 * \defgroup UnparametricField UnparametricField
	 *
	 * A field having an interface similar to that of floats is adapted to LinBox.
	 *
	 *  Used to generate efficient field classes for unparameterized fields (or hidden parameter fields).
	 *
	 *  Some fields are implemented by definition of the C++ arithmetic operators, such as z = x*y,
	 *  for z, y, z instances of a type K.   The LinBox field
	 *  Unparametric<K> is the adaptation to LinBox.
	 *
	 *  For a typical unparametric field, some of the methods must be defined in a specialization.
	 */
	template <class K>
	class UnparametricField {
	protected:
		long int _p ; long int _card ;
	public:

		/** The field's element type.
		 * Type K must provide a default constructor,
		 * a copy constructor, a destructor, and an assignment operator.
		 */

		typedef K Element;

		/** @name Field Object Basics.
		*/
		//@{

		/** Builds this field to have characteristic q and cardinality q<sup>e</sup>.
		 *  This constructor must be defined in a specialization.
		 */
		UnparametricField(long int q = 0, size_t e = 1) :
			_p(q), _card(q == 0 ? -1 : pow(q, e) )
			{}  // assuming q is a prime or zero.

		/// construct this field as copy of F.
		UnparametricField (const UnparametricField &F) :
			_p(F._p), _card(F._card)
		{}

		// UnparametricField(){}
		///
		~UnparametricField () {}

		/* Assignment operator.
		 * Assigns UnparametricField object F to field.
		 * @param  F UnparametricField object.
		 */
		// I believe this should be virtual -bds
		///
		const UnparametricField &operator=(const UnparametricField &F) const
		{
		       	return *this;
		}
		//@} Field Object Basics.

		/** @name Data Object Management.
		 * first argument is set and the value is also returned.
		 */
		//@{

		Element& init (Element& x) const
		{
			return x;

		}

		/// x := y.  Caution: it is via cast to long.  Good candidate for specialization.
		template <typename Src>
		Element& init (Element& x, const Src& s) const
		{
			return Caster (x, s);

		}


		/// x :=  y.  Caution: it is via cast to long.  Good candidate for specialization. --dpritcha

		template <typename T>
		T& convert (T &x, const Element &y) const
		{
			return Caster (x,y);
		}

		///
		Element &assign (Element &x, const Element &y) const
		{
		       	return x = y;
		}


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

} // FFPACK
#endif // __FIELD_UNPARAMETRIC_H_
