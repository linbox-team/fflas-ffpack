/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2009 LinBox
 * Written by C Pernet
 * updated to compilable condition by <brice.boyer@imag.fr>
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


/*! @file field/modular-balanced-int32.h
 * @ingroup field
 * @brief Balanced representation of <code>Z/mZ</code> over \c int32_t .
 */

#ifndef __FFLASFFPACK_modular_balanced_int32_H
#define __FFLASFFPACK_modular_balanced_int32_H

#include <math.h>
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"
#include "fflas-ffpack/utils/debug.h"

#ifndef LINBOX_MAX_INT
#define LINBOX_MAX_INT INT32_MAX
#endif


namespace FFPACK
{


	/// \ingroup field
	template <>
	class ModularBalanced<int32_t> {
	protected:
		int32_t modulus;
		int32_t half_mod;
		int32_t mhalf_mod;
		double modulusinv;

	public:

		typedef int32_t Element;
		typedef ModularBalancedRandIter<int32_t> RandIter;
		typedef NonzeroRandIter<ModularBalanced<int32_t>,RandIter> NonZeroRandIter;

		const Element one  ;
		const Element zero ;
		const Element mOne ;

		static const bool balanced = true;

		//default modular field,taking 65521 as default modulus
		ModularBalanced () :
			modulus(65521)
			,one(1),zero(0),mOne(-1)
		{
			modulusinv = 1/(double)65521;
			half_mod = (65521 >> 1);
			mhalf_mod = half_mod-65520;
		}

		ModularBalanced (int32_t value, int exp = 1)  :
			modulus(value)
			,one(1),zero(0),mOne(-1)
		{
			half_mod = (modulus >> 1);
			mhalf_mod = half_mod-modulus+1;
			modulusinv = 1 / ((double) value);
#ifdef DEBUG
			if(exp != 1) throw Failure(__func__,__FILE__,__LINE__,"exponent must be 1");
			if(value <= 1) throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			// std::cout << value << '<' << getMaxModulus() << std::endl;
			if(value > getMaxModulus() ) throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
			if( ! (value % 2) ) throw Failure(__func__,__FILE__,__LINE__,"modulus must be odd");
#endif

		}

		ModularBalanced (const ModularBalanced<int32_t>& mf) :
			modulus(mf.modulus),
			half_mod(mf.half_mod),
			mhalf_mod(mf.mhalf_mod),
			modulusinv(mf.modulusinv)
			,one(mf.one),zero(mf.zero),mOne(mf.mOne)
		{ }

		ModularBalanced<Element> & assign(const ModularBalanced<Element> &F)
		{
			modulus = F.modulus;
			half_mod  = F.half_mod;
			mhalf_mod = F.mhalf_mod;
			// lmodulus   = F.lmodulus;
			modulusinv = F.modulusinv;
			F.assign(const_cast<Element&>(one),F.one);
			F.assign(const_cast<Element&>(zero),F.zero);
			F.assign(const_cast<Element&>(mOne),F.mOne);
			return *this;
		}

#if 1
		const ModularBalanced &operator=(const ModularBalanced<int32_t> &F)
		{
			modulus = F.modulus;
			half_mod  = F.half_mod;
			mhalf_mod = F.mhalf_mod;
			// lmodulus   = F.lmodulus;
			modulusinv = F.modulusinv;
			F.assign(const_cast<Element&>(one),F.one);
			F.assign(const_cast<Element&>(zero),F.zero);
			F.assign(const_cast<Element&>(mOne),F.mOne);
			return *this;
		}
#endif

		size_t characteristic () const
		{
			return (size_t) modulus;
		}

		unsigned long &characteristic (unsigned long & c) const
		{
			return c=(unsigned long)modulus;
		}

		size_t cardinality () const
		{
			return (size_t) modulus;
		}


		double & convert(double &x, const Element &y) const
		{
			return x = (double) y;
		}

		float & convert(float &x, const Element &y) const
		{
			return x = (float) y;
		}



		std::ostream &write (std::ostream &os) const
		{
			return os << "balanced int32_t mod " << modulus;
		}

		std::istream &read (std::istream &is)
		{
			is >> modulus;
			half_mod = modulus/2;
			mhalf_mod = half_mod-modulus+1;
			modulusinv = 1 /((double) modulus );
#ifdef DEBUG
			if(modulus <= 1) throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");

			if(modulus > getMaxModulus() ) throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
			if( ! (modulus % 2) ) throw Failure(__func__,__FILE__,__LINE__,"modulus must be oddd");
#endif

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


		Element &init (Element &x, const double &y) const
		{
			x = (Element) fmod(y,(double)modulus);
			if (x < mhalf_mod) x += modulus;
			else if (x > half_mod) x -= modulus;
			return x;
		}


		Element &init (Element &x, const size_t &y) const
		{
			x = Element((Element)y % (Element)modulus);
			if (x < mhalf_mod)
				x += modulus;
			else if (x > half_mod)
				x -= modulus;
			return x;
		}


		inline Element& init(Element& x, int y =0) const
		{
			x = Element(y % modulus);

			if ( x < mhalf_mod )
				x += modulus;
			else if (x > half_mod )
				x -= modulus;

			return x;
		}

		inline Element& init(Element& x, long y) const
		{
			x = Element(y % modulus);
			if ( x < mhalf_mod )
				x += modulus;
			else if ( x > half_mod )
				x -= modulus;

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

		inline Element &add (Element &x, const Element &y, const Element &z) const
		{
			x = y + z;
			if ( x > half_mod ) x -= modulus;
			else if ( x < mhalf_mod ) x += modulus;

			return x;
		}

		inline Element &sub (Element &x, const Element &y, const Element &z) const
		{
			x = y - z;
			if (x > half_mod) x -= modulus;
			else if (x < mhalf_mod) x += modulus;
			return x;
		}

		inline Element &mul (Element &x, const Element &y, const Element &z) const
		{
			int32_t q;

			q  = (int32_t) ((((double) y) * ((double) z)) * modulusinv);  // q could be off by (+/-) 1
			x = (int32_t) (y*z - q*modulus);

			if (x > half_mod)
				x -= modulus;
			else if (x < mhalf_mod)
				x += modulus;

			return x;
		}

		inline Element &div (Element &x, const Element &y, const Element &z) const
		{
			Element temp;
			inv (temp, z);
			return mul (x, y, temp);
		}

		inline Element &neg (Element &x, const Element &y) const
		{
			return x = -y;
		}

		inline Element &inv (Element &x, const Element &y) const
		{
			int32_t d, t;
			XGCD(d, x, t, y, modulus);
#ifdef DEBUG
			if (d != 1)
				throw Failure(__func__,__FILE__,__LINE__,"InvMod: inverse undefined");
#endif
			if (x > half_mod)
				x -= modulus;
			else if (x < mhalf_mod)
				x += modulus;

			return x;

		}

		inline Element &axpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			int32_t q;

			q  = (int32_t) (((((double) a) * ((double) x)) + (double)y) * modulusinv);  // q could be off by (+/-) 1
			r = (int32_t) (a * x + y - q*modulus);


			if (r > half_mod)
				r -= modulus;
			else if (r < mhalf_mod)
				r += modulus;

			return r;

		}

		inline Element &addin (Element &x, const Element &y) const
		{
			x += y;
			if ( x > half_mod ) x -= modulus;
			else if (x < mhalf_mod) x += modulus;

			return x;
		}

		inline Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			if (x > half_mod)
				x -= modulus;
			else if (x < mhalf_mod)
				x += modulus;

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
			return x = -x;
		}

		inline Element &invin (Element &x) const
		{
			return inv (x, x);
		}

		inline Element &axpyin (Element &r, const Element &a, const Element &x) const
		{
			int32_t q;

			q  = (int32_t) (((((double) a)*((double) x)) + (double)r) * modulusinv);  // q could be off by (+/-) 1
			r = (int32_t) (a * x + r - q*modulus);


			if (r > half_mod)
				r -= modulus;
			else if (r < mhalf_mod)
				r += modulus;

			return r;
		}

		static inline int32_t getMaxModulus()
		{
			// return 1073741824; // 2^30
			return 2147483647; // 2^31 -1
			// return 1 << 15; // 2^15
			// FFLASFFPACK_check(46340LL*46341LL<2147483647LL);
			// FFLASFFPACK_check(46342LL*46341LL>2147483647LL);
			// return  92681 ;
		}

	private:

		inline static void XGCD(int32_t& d, int32_t& s, int32_t& t, int32_t a, int32_t b) {
			int32_t  u, v, u0, v0, u1, v1, u2, v2, q, r;

			int32_t aneg = 0, bneg = 0;

			if (a < 0) {
#ifdef DEBUG
				if (a < -LINBOX_MAX_INT) throw Failure(__func__,__FILE__,__LINE__,"XGCD: integer overflow");
#endif
				a = -a;
				aneg = 1;
			}

			if (b < 0) {
#ifdef DEBUG
				if (b < -LINBOX_MAX_INT) throw Failure(__func__,__FILE__,__LINE__,"XGCD: integer overflow");
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


#include "field-general.h"
#endif // __FFLASFFPACK_modular_balanced_int32_H

