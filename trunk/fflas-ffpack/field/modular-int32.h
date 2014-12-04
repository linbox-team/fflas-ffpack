/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2010,2014 LinBox
 * Adapted by B Boyer <brice.boyer@imag.fr>
 * (from other modular-balanced* files)
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

/*! @file field/modular-int32.h
 * @ingroup field
 * @brief  representation of <code>Z/mZ</code> over \c int32_t .
 */

#ifndef __FFLASFFPACK_modular_int32_H
#define __FFLASFFPACK_modular_int32_H

#include <math.h>
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"
#include "fflas-ffpack/utils/debug.h"

#ifndef LINBOX_MAX_INT
#define LINBOX_MAX_INT INT32_MAX
#endif

#define NORMALISE_LO(x) \
{ \
			if (x < 0.) x += modulus; \
}

#define NORMALISE_HI(x) \
{ \
			if (x >= modulus) x -= modulus; \
}


namespace FFPACK {

	template< class Element >
	class Modular;

	/** \brief Specialization of Modular to int32_t element type with efficient dot product.
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
	class Modular<int32_t> {

	protected:

		int32_t modulus;
		double modulusinv;
		uint32_t lmodulus;
		int32_t _two64;

	public:

		typedef int32_t Element;
		typedef int32_t* Element_ptr;
		typedef const int32_t* ConstElement_ptr;

		const Element one   ;
		const Element zero  ;
		const Element mOne ;

	public:

		static const bool balanced = false ;

		typedef ModularRandIter<Element> RandIter;
		typedef NonzeroRandIter<Modular<Element>, RandIter> NonZeroRandIter;

		//default modular field,taking 65521 as default modulus
		Modular () :
			modulus(65521),lmodulus((uint32_t) modulus)
			,one(1),zero(0),mOne(modulus -1)
		{
			modulusinv=1/(double)65521;

			_two64 = (Element) ((uint64_t) (-1) % (uint64_t) 65521);
			_two64 += 1;
			if (_two64 >= 65521) _two64 -= 65521;
		}

		Modular (Element value, int32_t exp = 1) :
			modulus(value),lmodulus((uint32_t)value)
			,one(1),zero(0),mOne(modulus -1)
		{
			modulusinv = 1 / ((double) value);
#ifdef DEBUG
			if(exp != 1) throw Failure(__func__,__FILE__,__LINE__,"exponent must be 1");
			if(value<=1) throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			if(value>getMaxModulus())  {
				std::cerr << value << '>' << getMaxModulus() << std::endl;
				throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
			}
#endif
			_two64 = (Element) ((uint64_t) (-1) % (uint64_t) value);
			_two64 += 1;
			if (_two64 >= value) _two64 -= value;
		}

		Modular (const Modular<Element>& mf) :
			modulus(mf.modulus),modulusinv(mf.modulusinv)
			,lmodulus(mf.lmodulus),_two64(mf._two64)
			,one(mf.one),zero(mf.zero),mOne(mf.mOne)
		{
		}

		Modular <Element>& assign(const Modular<Element> &F)
		{
			modulus = F.modulus;
			modulusinv = F.modulusinv;
			lmodulus   = F.lmodulus;
			_two64     = F._two64;
			//inv_modulus = F.inv_modulus;
			F.assign(const_cast<Element&>(one),F.one);
			F.assign(const_cast<Element&>(zero),F.zero);
			F.assign(const_cast<Element&>(mOne),F.mOne);
			return *this;
		}

#if 1
		const Modular &operator=(const Modular<Element> &F)
		{
			modulus = F.modulus;
			modulusinv = F.modulusinv;
			lmodulus   = F.lmodulus;
			_two64     = F._two64;
			//inv_modulus = F.inv_modulus;
			F.assign(const_cast<Element&>(one),F.one);
			F.assign(const_cast<Element&>(zero),F.zero);
			F.assign(const_cast<Element&>(mOne),F.mOne);
			return *this;
		}
#endif

		inline uint32_t &cardinality ( uint32_t &c) const
		{
			return c = lmodulus;
		}

		inline uint32_t cardinality () const
		{
			return lmodulus;
		}

		inline uint32_t &characteristic (uint32_t &c) const
		{
			return c = lmodulus;
		}

		inline uint32_t characteristic () const
		{
			return lmodulus;
		}

			inline Element &convert (Element &x, const Element &y) const
		{
			return x = y;
		}

		inline double & convert(double &x, const Element &y) const
		{
			return x = (double) y;
		}

		inline float & convert(float &x, const Element &y) const
		{
			return x = (float) y;
		}

		inline std::ostream &write (std::ostream &os) const
		{
			return os << "int32_t mod " << modulus;
		}

		inline std::istream &read (std::istream &is)
		{
			is >> modulus;
			modulusinv = 1 /((double) modulus );
#ifdef DEBUG
			if(modulus <= 1) throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			if(modulus > getMaxModulus()) throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
#endif
			_two64 = (Element) ((uint64_t) (-1) % (uint64_t) modulus);
			_two64 += 1;
			if (_two64 >= modulus) _two64 -= modulus;

			return is;
		}

		inline std::ostream &write (std::ostream &os, const Element &x) const
		{
			return os << x;
		}

		inline std::istream &read (std::istream &is, Element &x) const
		{
			Element tmp;
			is >> tmp;
			init(x,tmp);
			return is;
		}

		inline Element &init (Element & x, const double &y) const
		{
			x = (Element) fmod(y,(double)modulus);
			NORMALISE_LO(x);
			return x;
		}

		inline Element &init (Element & x, const float &y) const
		{
			return init(x , (double) y);
		}

		template<class Element1>
		inline Element &init (Element & x, const Element1 &y) const
		{
			return reduce (x, Element(y));
		}

		inline Element &reduce (Element & x, const Element &y) const
		{
			x = (y % modulus);
			NORMALISE_LO(x);
			return x;
		}

		inline Element &reduce (Element & x) const
		{
			x %= modulus;
			NORMALISE_LO(x);
			return x;
		}

		inline Element& init(Element&x) const
		{
			return x = 0;
		}

		inline Element& init(Element& x, int64_t y ) const
		{
			x = Element(y % (int64_t)modulus);
			NORMALISE_LO(x);
			return x;
		}

		inline Element& init(Element& x, uint32_t y ) const
		{
			x = (Element)(y % lmodulus);
			NORMALISE_LO(x);
			return x;
		}

		inline Element& init (Element& x, uint64_t y) const
		{
			x = (Element)(y % lmodulus);
			NORMALISE_LO(x);
			return x;
		}

		inline Element& assign(Element& x, const Element& y) const
		{
			return x = y;
		}

		inline bool areEqual (const Element &x, const Element &y) const
		{
			return x == y;
		}

		inline  bool isZero (const Element &x) const
		{
			return x == 0;
		}

		inline bool isOne (const Element &x) const
		{
			return x == 1;
		}

		inline bool isMOne (const Element &x) const
		{
			return x == mOne ;
		}

		inline Element &add (Element &x, const Element &y, const Element &z) const
		{
			x = y + z;
			NORMALISE_HI(x);
			return x;
		}

		inline Element &sub (Element &x, const Element &y, const Element &z) const
		{
			x = y - z;
			NORMALISE_LO(x);
			return x;
		}

		inline Element &mul (Element &x, const Element &y, const Element &z) const
		{
			Element q;

			q  = (Element) ((((double) y) * ((double) z)) * modulusinv);  // q could be off by (+/-) 1
			x = (Element) (y*z - q*modulus);

			NORMALISE_LO(x);
			NORMALISE_HI(x);
			return x;
		}

		inline Element &div (Element &x, const Element &y, const Element &z) const
		{
			Element temp;
			return mul (x, y, inv(temp,z));
		}

		inline Element &neg (Element &x, const Element &y) const
		{
			if(y == 0) return x=0;
			else return x = modulus-y;
		}

		inline Element &inv (Element &x, const Element &y) const
		{
			FFLASFFPACK_check(!isZero(y));
			Element d, t;
			XGCD(d, x, t, y, modulus);
			if (d != 1)
			{
#ifdef DEBUG
				throw Failure(__func__,__FILE__,__LINE__,"InvMod: Input is not invertible ");
#endif
			}
			NORMALISE_LO(x);
			return x;

		}

		inline Element &axpy (Element &r,
			       const Element &a,
			       const Element &x,
			       const Element &y) const
		{
			Element q;

			q  = (Element) (((((double) a) * ((double) x)) + (double)y) * modulusinv);  // q could be off by (+/-) 1
			r = (Element) (a * x + y - q*modulus);

			NORMALISE_LO(r);
			NORMALISE_HI(r);

			return r;

		}

		inline Element &axmy (Element &r,
			       const Element &a,
			       const Element &x,
			       const Element &y) const
		{
			Element q;

			q  = (Element) (((((double) a) * ((double) x)) - (double)y) * modulusinv);  // q could be off by (+/-) 1
			r = (Element) (a * x - y - q*modulus);

			NORMALISE_LO(r);
			NORMALISE_HI(r);

			return r;

		}

		inline Element &maxpy (Element &r,
			       const Element &a,
			       const Element &x,
			       const Element &y) const
		{
			return negin(axmy(r,a,x,y));
		}

		inline Element &addin (Element &x, const Element &y) const
		{
			x += y;
			NORMALISE_HI(x);
			return x;
		}

		inline Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			NORMALISE_LO(x);
			return x;
		}

		inline Element &mulin (Element &x, const Element &y) const
		{
			return mul(x,x,y);
		}

		inline Element &divin (Element &x, const Element &y) const
		{
			return div(x,x,y);
		}

		inline Element &negin (Element &x) const
		{
			if (x == 0) return x;
			else return x = modulus - x;
		}

		inline Element &invin (Element &x) const
		{
			return inv (x, x);
		}

		inline Element &axpyin (Element &r, const Element &a, const Element &x) const
		{
			Element q;

			q  = (Element) (((((double) a) * ((double) x)) + (double) r) * modulusinv);  // q could be off by (+/-) 1
			r = (Element) (a * x + r - q*modulus);

			NORMALISE_LO(r);
			NORMALISE_HI(r);

			return r;
		}

		inline Element &maxpyin (Element &r, const Element &a, const Element &x) const
		{
			Element q;
            maxpy(q,a,x,r);
            return assign(r,q);
		}

		static inline Element getMaxModulus()
		{
			return 46341 ;
		}

		static  Element getMinModulus()	{return 2.0;}

		Element minElement() const
		{
			return zero ;
		}

		Element maxElement() const
		{
			return mOne ;
		}

	private:

		static void XGCD(Element& d, Element& s, Element& t, Element a, Element b)
		{
			Element  u, v, u0, v0, u1, v1, u2, v2, q, r;
			Element aneg = 0, bneg = 0;

			if (a < 0) {
#ifdef DEBUG
				if (a < -LINBOX_MAX_INT)
					throw Failure(__func__,__FILE__,__LINE__,"XGCD: integer overflow");
#endif
				a = -a;
				aneg = 1;
			}

			if (b < 0) {
#ifdef DEBUG
				if (b < -LINBOX_MAX_INT)
					throw Failure(__func__,__FILE__,__LINE__,"XGCD: integer overflow");
#endif
				b = -b;
				bneg = 1;
			}

			u1 = 1; v1 = 0;
			u2 = 0; v2 = 1;
			u = a; v = b;

			while (v != 0) {
				q = u / v;
				r = u % v;
				u = v;
				v = r;
				u0 = u2;
				v0 = v2;
				u2 =  u1 - q*u2;
				v2 = v1- q*v2;
				u1 = u0;
				v1 = v0;
			}

			if (aneg)
				u1 = -u1;

			if (bneg)
				v1 = -v1;

			d = u;
			s = u1;
			t = v1;
		}

	};

}

#undef LINBOX_MAX_INT
#undef NORMALISE_LO
#undef NORMALISE_HI

#include "field-general.h"

#endif //__LINBOX_modular_int32_H

