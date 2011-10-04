/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2010 LinBox
 *
 *
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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

/*! @file field/modular-int32_t.h
 * @ingroup field
 * @brief  representation of <code>Z/mZ</code> over \c int32_t .
 */
#ifndef __FFLASFFPACK_modular_int32_H
#define __FFLASFFPACK_modular_int32_H

#include <math.h>
#include <sys/time.h>
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"
#include "fflas-ffpack/utils/debug.h"

#ifndef LINBOX_MAX_INT
#define LINBOX_MAX_INT INT32_MAX
#endif

// Namespace in which all LinBox code resides
namespace FFPACK
{

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
		unsigned long lmodulus;
		int32_t _two64;


	public :
		typedef int32_t Element;
		static const Element one  = 1 ;
		static const Element zero = 0 ;
		Element mone ; // can't be const because of operator=

	public:


		static const bool balanced = false ;
		typedef ModularRandIter<Element> RandIter;
		typedef NonzeroRandIter<Modular<Element>, ModularRandIter<Element> > NonZeroRandIter;

		//default modular field,taking 65521 as default modulus
		Modular () :
			modulus(65521),lmodulus(modulus)
			,mone(modulus -1)
		{
			modulusinv=1/(double)65521;

			_two64 = (int32_t) ((uint64_t) (-1) % (uint64_t) 65521);
			_two64 += 1;
			if (_two64 >= 65521) _two64 -= 65521;
		}

		Modular (int32_t value, int32_t exp = 1) :
			modulus(value),lmodulus(value)
			,mone(modulus -1)
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
			_two64 = (int32_t) ((uint64_t) (-1) % (uint64_t) value);
			_two64 += 1;
			if (_two64 >= value) _two64 -= value;
		}

		Modular (unsigned long int value) :
			modulus((Element) value),lmodulus(value)
			,mone(modulus -1)
		{
			modulusinv = 1 / ((double) value);
#ifdef DEBUG
			if(value<=1) throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			if (value>INT32_MAX)  // stupidly big ?
				throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
			if((Element)value>getMaxModulus()) // we can cast now
				throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
#endif
			_two64 = (int32_t) ((uint64_t) (-1) % (uint64_t) value);
			_two64 += 1;
			if ((unsigned long)_two64 >= value)
				_two64 = _two64 - (int32_t) value;
		}

		Modular (long int value) :
			modulus((Element) value),lmodulus((long int)value)
			,mone(modulus -1)
		{
			modulusinv = 1 / ((double) value);
#ifdef DEBUG
			if(value<=1) throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			if (value>INT32_MAX)  // stupidly big ?
				throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
			if((Element)value>getMaxModulus()) // we can cast now
				throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
#endif
			_two64 = (int32_t) ((uint64_t) (-1) % (uint64_t) value);
			_two64 += 1;
			if ((long int)_two64 >= value)
				_two64 = _two64 - (int32_t) value;
		}


		Modular(const Modular<int32_t>& mf) :
			modulus(mf.modulus),modulusinv(mf.modulusinv)
			,lmodulus(mf.lmodulus),_two64(mf._two64)
			,mone(modulus -1)
		{}

		const Modular &operator=(const Modular<int32_t> &F)
		{
			modulus    = F.modulus;
			modulusinv = F.modulusinv;
			lmodulus   = F.lmodulus ;
			_two64     = F._two64;
			mone       = F.mone ;
			return *this;
		}


		unsigned long &cardinality (unsigned long &c) const
		{
			return c = lmodulus;
		}

		unsigned long &characteristic (unsigned long &c) const
		{
			return c = lmodulus;
		}

		unsigned long characteristic () const
		{
			return lmodulus;
		}

		unsigned long cardinality () const
		{
			return lmodulus;
		}


		int32_t &convert (int32_t &x, const Element &y) const
		{
			return x = y;
		}

		double &convert (double &x, const Element &y) const
		{
			return x = (double) y;
		}

		float &convert (float &x, const Element &y) const
		{
			return x = (float) y;
		}

		std::ostream &write (std::ostream &os) const
		{
			return os << "int32_t mod " << modulus;
		}

		std::istream &read (std::istream &is)
		{
			is >> modulus;
			modulusinv = 1 /((double) modulus );
#ifdef DEBUG
			if(modulus <= 1) throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			if(modulus > getMaxModulus()) throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
#endif
			_two64 = (int32_t) ((uint64_t) (-1) % (uint64_t) modulus);
			_two64 += 1;
			if (_two64 >= modulus) _two64 -= modulus;

			return is;
		}

		std::ostream &write (std::ostream &os, const Element &x) const
		{
			return os << x;
		}

		std::istream &read (std::istream &is, Element &x) const
		{
			long int tmp;
			is >> tmp;
			init(x,tmp);
			return is;
		}

		Element &init (Element & x, const double &y) const
		{
			double z = fmod(y, (double)modulus);
			if (z < 0)
				z += (double)modulus;
			//z += 0.5; // C Pernet Sounds nasty and not necessary
			return x = static_cast<Element>(z); //rounds towards 0
		}

		Element &init (Element & x, const float &y) const
		{
			return init(x , (double) y);
		}

		template<class Element1>
		Element &init (Element & x, const Element1 &y) const
		{
			x = Element(y) % modulus;
			if (x < 0) x += modulus;
			return x;
		}

		Element& init(Element& x, int y =0) const
		{
			x = y % modulus;
			if ( x < 0 ) x += modulus;
			return x;
		}

		Element& init(Element& x, long y) const
		{
			x = Element(y % (long)modulus);
			if ( x < 0 ) x += modulus;
			return x;
		}

		Element& assign(Element& x, const Element& y) const
		{
			return x = y;
		}


		bool areEqual (const Element &x, const Element &y) const
		{
			return x == y;
		}

		bool isZero (const Element &x) const
		{
			return x == 0;
		}

		bool isOne (const Element &x) const
		{
			return x == 1;
		}

		Element &add (Element &x, const Element &y, const Element &z) const
		{
			x = y + z;
			if ( x >= modulus ) x -= modulus;
			return x;
		}

		Element &sub (Element &x, const Element &y, const Element &z) const
		{
			x = y - z;
			if (x < 0)
				x += (Element) modulus;
			return x;
		}

		Element &mul (Element &x, const Element &y, const Element &z) const
		{
			int32_t q;

			q  = (int32_t) ((((double) y)*((double) z)) * modulusinv);  // q could be off by (+/-) 1
			x = (int32_t) (y*z - q*modulus);


			if (x >= modulus)
				x -= (Element) modulus;
			else if (x < 0)
				x += (Element) modulus;

			return x;
		}

		Element &div (Element &x, const Element &y, const Element &z) const
		{
			FFLASFFPACK_check(!isZero(z));
			Element temp;
			inv (temp, z);
			return mul (x, y, temp);
		}

		Element &neg (Element &x, const Element &y) const
		{
			if(y == 0) return x=0;
			else return x = modulus-y;
		}

		Element &inv (Element &x, const Element &y) const
		{
			FFLASFFPACK_check(!isZero(y));
			int32_t d, t;
			XGCD(d, x, t, y, modulus);
			if (d != 1)
			{
#ifdef DEBUG
				throw Failure(__func__,__FILE__,__LINE__,"InvMod: Input is not invertible ");
#endif
			}
			if (x < 0)
				x += modulus;
			return x;

		}

		Element &axpy (Element &r,
			       const Element &a,
			       const Element &x,
			       const Element &y) const
		{
			int32_t q;

			q  = (int32_t) (((((double) a) * ((double) x)) + (double)y) * modulusinv);  // q could be off by (+/-) 1
			r = (int32_t) (a * x + y - q*modulus);


			if (r >= modulus)
				r -= modulus;
			else if (r < 0)
				r += modulus;

			return r;

		}

		Element &addin (Element &x, const Element &y) const
		{
			x += y;
			if (  x >= modulus ) x -= modulus;
			return x;
		}

		Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			if (x < 0) x += modulus;
			return x;
		}

		Element &mulin (Element &x, const Element &y) const
		{
			return mul(x,x,y);
		}

		Element &divin (Element &x, const Element &y) const
		{
			return div(x,x,y);
		}

		Element &negin (Element &x) const
		{
			if (x == 0) return x;
			else return x = modulus - x;
		}

		Element &invin (Element &x) const
		{
			FFLASFFPACK_check(!isZero(x));
			return inv (x, x);
		}

		Element &axpyin (Element &r, const Element &a, const Element &x) const
		{
			int32_t q;

			q  = (int32_t) (((((double) a) * ((double) x)) + (double) r) * modulusinv);  // q could be off by (+/-) 1
			r = (int32_t) (a * x + r - q*modulus);


			if (r >= modulus)
				r -= modulus;
			else if (r < 0)
				r += modulus;

			return r;
		}

		unsigned long AccBound(const Element&r) const
		{
			// Element one, zero ;
			// init(one,1UL) ;
			// init(zero,0UL);
			double max_double = (double) (INT32_MAX) - modulus ;
			double p = modulus-1 ;
			if (areEqual(zero,r))
				return (unsigned long) (max_double/p) ;
			else if (areEqual(one,r))
			{
				if (modulus>= getMaxModulus())
					return 0 ;
				else
					return (unsigned long) max_double/(modulus*modulus) ;
			} else
				throw "Bad input, expecting 0 or 1";
			return 0;
		}


		static  int32_t getMaxModulus()
		{
			// return INT32_MAX ; // 2^31-1
			return 1073741824;// 2^30
			// return 46341 ;
		}

	private:

		static void XGCD(int32_t& d, int32_t& s, int32_t& t, int32_t a, int32_t b)
		{
			int32_t  u, v, u0, v0, u1, v1, u2, v2, q, r;

			int32_t aneg = 0, bneg = 0;

			if (a < 0)
			{
#ifdef DEBUG
				if (a < -LINBOX_MAX_INT) throw Failure(__func__,__FILE__,__LINE__,"XGCD: integer overflow");
#endif
				a = -a;
				aneg = 1;
			}

			if (b < 0)
			{
#ifdef DEBUG
				if (b < -LINBOX_MAX_INT) throw Failure(__func__,__FILE__,__LINE__,"XGCD: integer overflow");
#endif
				b = -b;
				bneg = 1;
			}

			u1 = 1; v1 = 0;
			u2 = 0; v2 = 1;
			u = a; v = b;

			while (v != 0)
			{
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

#include "field-general.h"

#endif //__LINBOX_modular_int32_H

