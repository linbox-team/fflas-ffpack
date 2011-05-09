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
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


/*! @file field/modular-balanced-int32.h
 * @ingroup field
 * @brief Balanced representation of <code>Z/mZ</code> over \c int32_t .
 */

#ifndef __FFLAFLAS_modular_balanced_int32_H
#define __FFLAFLAS_modular_balanced_int32_H

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
		int32_t halfmodulus;
		int32_t nhalfmodulus;
		double modulusinv;

	public:

		typedef int32_t Element;
		typedef ModularBalancedRandIter<int32_t> RandIter;
		typedef NonzeroRandIter<ModularBalanced<int32_t>,RandIter> NonZeroRandIter;

		const bool balanced;

		//default modular field,taking 65521 as default modulus
		ModularBalanced () :
			modulus(65521),balanced(true)
		{
			modulusinv = 1/(double)65521;
			halfmodulus = (65521 >> 1);
			nhalfmodulus = halfmodulus-65520;
		}

		ModularBalanced (int32_t value, int exp = 1)  :
			modulus(value),balanced(true)
		{
			halfmodulus = (modulus >> 1);
			nhalfmodulus = halfmodulus-modulus+1;
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
			halfmodulus(mf.halfmodulus),
			nhalfmodulus(mf.nhalfmodulus),
			modulusinv(mf.modulusinv),
			balanced(true)
		{ }

		const ModularBalanced &operator=(const ModularBalanced<int32_t> &F)
		{
			modulus      = F.modulus;
			halfmodulus  = F.halfmodulus;
			nhalfmodulus = F.nhalfmodulus;
			modulusinv   = F.modulusinv;

			return *this;
		}

		size_t characteristic () const
		{
			return modulus;
		}

		size_t cardinality () const
		{
			return modulus;
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
			halfmodulus = modulus/2;
			nhalfmodulus = halfmodulus-modulus+1;
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
			double tmp;
			is >> tmp;
			init(x,tmp);
			return is;
		}


		Element &init (Element &x, const double &y) const
		{
			x = (Element) fmod(y,(double)modulus);
			if (x < nhalfmodulus) x += modulus;
			else if (x > halfmodulus) x -= modulus;
			return x;
		}


		Element &init (Element &x, const size_t &y) const
		{
			x = y % (modulus);
			if (x < nhalfmodulus) x += modulus;
			else if (x > halfmodulus) x -= modulus;
			return x;
		}


		inline Element& init(Element& x, int y =0) const
		{
			x = y % modulus;

			if ( x < nhalfmodulus ) x += modulus;
			else if (x > halfmodulus ) x -= modulus;

			return x;
		}

		inline Element& init(Element& x, long y) const
		{
			x = y % modulus;
			if ( x < nhalfmodulus ) x += modulus;
			else if ( x > halfmodulus ) x -= modulus;

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
			if ( x > halfmodulus ) x -= modulus;
			else if ( x < nhalfmodulus ) x += modulus;

			return x;
		}

		inline Element &sub (Element &x, const Element &y, const Element &z) const
		{
			x = y - z;
			if (x > halfmodulus) x -= modulus;
			else if (x < nhalfmodulus) x += modulus;
			return x;
		}

		inline Element &mul (Element &x, const Element &y, const Element &z) const
		{
			int32_t q;

			q  = (int32_t) ((((double) y) * ((double) z)) * modulusinv);  // q could be off by (+/-) 1
			x = (int32_t) (y*z - q*modulus);

			if (x > halfmodulus)
				x -= modulus;
			else if (x < nhalfmodulus)
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
			if (x > halfmodulus)
				x -= modulus;
			else if (x < nhalfmodulus)
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


			if (r > halfmodulus)
				r -= modulus;
			else if (r < nhalfmodulus)
				r += modulus;

			return r;

		}

		inline Element &addin (Element &x, const Element &y) const
		{
			x += y;
			if ( x > halfmodulus ) x -= modulus;
			else if (x < nhalfmodulus) x += modulus;

			return x;
		}

		inline Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			if (x > halfmodulus)
				x -= modulus;
			else if (x < nhalfmodulus)
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


			if (r > halfmodulus)
				r -= modulus;
			else if (r < nhalfmodulus)
				r += modulus;

			return r;
		}

		static inline int32_t getMaxModulus()
		{
			// return 1073741824; // 2^30
			return 1 << 15; // 2^15
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
#endif // __FFLAFLAS_modular_balanced_int32_H

