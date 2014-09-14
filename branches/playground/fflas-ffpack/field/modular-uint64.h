/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2014 LinBox
 * Written by BB <brice.boyer@lip6.fr>
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
 *
 */

/*! @file field/modular-uint64.h
 * @ingroup field
 * @brief  representation of <code>Z/mZ</code> over \c uint64_t .
 */

#ifndef __FFLASFFPACK_modular_uint64_H
#define __FFLASFFPACK_modular_uint64_H

#include <math.h>
#include <sys/time.h>
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"
#include "fflas-ffpack/utils/debug.h"


// Namespace in which all LinBox code resides
namespace FFPACK
{
	template< class Element >
	class Modular;


	/** \brief Specialization of Modular to uint64_t element type with efficient dot product.
	 *
	 * Efficient element operations for dot product, mul, axpy, by using floating point
	 * inverse of modulus (borrowed from NTL) and some use of non-normalized intermediate values.
	 *
	 * For some uses this is the most efficient field for primes in the range from half word
	 * to 2^30.
	 *
	 * Requires: Modulus < 2^30.
	 * Intended use: 2^15 < prime modulus < 2^30.
	 * \ingroup field

	 * @todo what about this _two64  not so usefull here ?? (but in linbox)
	 */
	template <>
	class Modular<uint64_t> {
	public:

		typedef uint64_t Element;
		typedef Element* Element_ptr ;
		typedef const Element* ConstElement_ptr;

		const Element zero,one,mOne ;
		int64_t _modulus;

		Modular () :
			zero(0),one(1),mOne(0)
		{}

		Modular (uint64_t value)  :
			zero(0),one(1),mOne(value-1)
			,_modulus(value)
		{
			init_two_64 ();
		}

		const Modular &operator=(const Modular &F)
		{
			_modulus = F._modulus;
			_two_64 = F._two_64;
			F.assign(const_cast<Element&>(one),F.one);
			F.assign(const_cast<Element&>(zero),F.zero);
			F.assign(const_cast<Element&>(mOne),F.mOne);

			return *this;
		}

		Element &init (Element &x, const long int &y ) const
		{
			x = abs (y) % (long int) (_modulus);
			if (y < 0) x = _modulus - x;
			return x;
		}

		Element &init (Element &x, const int &y ) const
		{
			x = abs (y) % (long int) (_modulus);
			if (y < 0) x = _modulus - x;
			return x;
		}

		Element &init (Element &x, const long unsigned int &y ) const
		{
			x = Element(y %  (_modulus));
			return x;
		}

		Element &init (Element &x, const unsigned int &y ) const
		{
			x = Element(y %  (_modulus));
			return x;
		}

		Element &init (Element &x, const double &y) const
		{
			double z = fmod(y, (double)_modulus);
			if (z < 0) z += (double) _modulus;
			return x = (Element) (z);
		}

		template< class XXX>
		Element& init(Element & x, const XXX & y) const
		{
			return init(x,double(y));
		}

		Element &init (Element &x) const
		{
			return x = zero ;
		}

		Element &convert (Element &x, const Element &y) const
		{
			return x = y;
		}

		double &convert (double &x, const Element &y) const
		{
			return x = (double)y;
		}

		float &convert (float &x, const Element &y) const
		{
			return x = (float)y;
		}


		Element &add (Element &x, const Element &y, const Element &z) const
		{
			x = y + z;
			if ((uint64_t) x >= (uint64_t) _modulus) x -= _modulus;
			return x;
		}

		Element &sub (Element &x, const Element &y, const Element &z) const
		{
			x = y - z;
			if ((int64_t) x < 0) x += _modulus;
			return x;
		}

		Element &mul (Element &x, const Element &y, const Element &z) const
		{
			return x = Element( ((uint64_t) y * (uint64_t) z) % (uint64_t) _modulus);
		}

		Element &div (Element &x, const Element &y, const Element &z) const
		{
			Element temp;
			inv (temp, z);
			return mul (x, y, temp);
		}

		Element &neg (Element &x, const Element &y) const
		{
			if (y == 0)
				return x = y;
			else
				return x = _modulus - y;
		}

		Element &inv (Element &x, const Element &y) const
		{
			// The extended Euclidean algoritm
			int64_t x_int, y_int, q, tx, ty, temp;
			x_int = _modulus;
			y_int = y;
			tx = 0;
			ty = 1;

			while (y_int != 0) {
				// always: gcd (modulus,residue) = gcd (x_int,y_int)
				//         sx*modulus + tx*residue = x_int
				//         sy*modulus + ty*residue = y_int
				q = x_int / y_int; // integer quotient
				temp = y_int;  y_int  = x_int  - q * y_int;
				x_int  = temp;
				temp = ty; ty = tx - q * ty;
				tx = temp;
			}

			if (tx < 0) tx += _modulus;

			// now x_int = gcd (modulus,residue)
			return x = Element(tx);
		}

		Element &axpy (Element &r,
			       const Element &a,
			       const Element &x,
			       const Element &y) const
		{
			r = Element( ((uint64_t) a * (uint64_t) x + (uint64_t) y) % (uint64_t) _modulus );
			if ((int64_t) r < 0) r += _modulus;
			return r;
		}

		Element &addin (Element &x, const Element &y) const
		{
			x += y;
			if ((uint64_t) x >= (uint64_t) _modulus) x -= _modulus;
			return x;
		}

		Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			if ((int64_t) x < 0) x += _modulus;
			return x;
		}

		Element &mulin (Element &x, const Element &y) const
		{
			x = Element( ((uint64_t) x * (uint64_t) y) % (uint64_t) _modulus );
			return x;
		}

		Element &divin (Element &x, const Element &y) const
		{
			Element temp;
			inv (temp, y);
			return mulin (x, temp);
		}

		Element &negin (Element &x) const
		{
			if (x == 0)
				return x;
			else
				return x = _modulus - x;
		}

		Element &invin (Element &x) const
		{
			return inv (x, x);
		}

		Element &axpyin (Element &r, const Element &a, const Element &x) const
		{
			r = Element( ((uint64_t) r + (uint64_t) a * (uint64_t) x) % (uint64_t) _modulus );
			if ((int64_t) r < 0) r += _modulus;
			return r;
		}


		static Element getMaxModulus()
		{
			return 1073741824;// 2^30 (ou plus ?)
		}

		Element &assign (Element &x, const Element &y) const
		{
			return x = y;
		}

		unsigned long &cardinality (unsigned long &c) const
		{
			return c = _modulus;
		}
		unsigned long cardinality () const
		{
			return  _modulus;
		}

		unsigned long &characteristic (unsigned long &c) const
		{
			return c = _modulus;
		}

		unsigned long characteristic () const
		{
			return  _modulus;
		}


		bool isZero (const Element &x) const
		{
			return x == 0;
		}

		bool isOne (const Element &x) const
		{
			return x == 1;
		}

		inline bool isMOne (const Element &x) const
		{
			return x == mOne ;
		}

		bool areEqual (const Element &x, const Element &y) const
		{
			return x == y;
		}

	std::ostream &write (std::ostream &os) const
		{
			return os << "uint64_t mod " << _modulus;
		}

		std::istream &read (std::istream &is)
		{
			is >> _modulus;
			// modulusinv = 1 /((double) _modulus );
#ifdef DEBUG
			// if(modulus <= 1) throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			// if(modulus > getMaxModulus()) throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
#endif
			// _two64 = (int64_t) ((uint64_t) (-1) % (uint64_t) modulus);
			// _two64 += 1;
			// if (_two64 >= modulus) _two64 -= modulus;

			return is;
		}

		std::ostream &write (std::ostream &os, const Element &x) const
		{
			return os << x;
		}

		std::istream &read (std::istream &is, Element &x) const
		{
			int64_t tmp;
			is >> tmp;
			init(x,tmp);
			return is;
		}

		Element minElement() const
		{
			return zero ;
		}

		Element maxElement() const
		{
			return mOne ;
		}

	private:
		void init_two_64 ()
		{
			uint64_t two_64 = 2;

			for (int i = 0; i < 6; ++i)
				two_64 = (two_64 * two_64) % _modulus;

			_two_64 = (Element) two_64;
		}
	public: /*  ?? */

		Element _two_64;


	};


}

// const int64_t FFPACK::Modular<int64_t>::one  =  1UL;
// const int64_t FFPACK::Modular<int64_t>::zero =  0UL;




#include "field-general.h"

#endif //__LINBOX_modular_uint64_H

