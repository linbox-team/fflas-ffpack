/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2010 LinBox
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


/*! @file field/modular-balanced-int64.h
 * @ingroup field
 * @brief Balanced representation of <code>Z/mZ</code> over \c int64_t .
 * @warning NOT DEFINED for EVEN modulus
 */

#ifndef __FFLASFFPACK_modular_balanced_int64_H
#define __FFLASFFPACK_modular_balanced_int64_H

#include <math.h>
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"
#include "fflas-ffpack/utils/debug.h"


#ifndef LINBOX_MAX_INT64
#ifdef __x86_64__
#define LINBOX_MAX_INT64 INT64_MAX
#else
#define LINBOX_MAX_INT64 INT64_MAX
#endif
#endif

// todo INT64_MAX
namespace FFPACK
{


	/// \ingroup field
	template <>
	class ModularBalanced<int64_t>  {
	protected:
		int64_t modulus;
		int64_t half_mod;
		int64_t mhalf_mod;
		double modulusinv;

	public:


		typedef int64_t Element;
		typedef int64_t* Element_ptr;
		typedef const int64_t* ConstElement_ptr;
		typedef ModularBalancedRandIter<int64_t> RandIter;


		const Element one  ;
		const Element zero ;
		const Element mOne ;


		static const bool balanced = true ;

		//default modular field,taking 65521 as default modulus
		ModularBalanced () :
			modulus(65521)
			,one(1),zero(0),mOne(-1)
		{
			modulusinv = 1/(double)65521;
			half_mod = (65521 >> 1);
			mhalf_mod = half_mod-65520;
			FFLASFFPACK_check(isOdd(modulus));
		}

		ModularBalanced (int64_t value, int exp = 1)  :
			modulus(value)
			,one(1),zero(0),mOne(-1)
		{
			half_mod = (modulus >> 1);
			mhalf_mod = half_mod-modulus+1;
			modulusinv = 1 / ((double) value);
			FFLASFFPACK_check(isOdd(modulus));
#ifdef DEBUG
			if(exp != 1) throw Failure(__func__,__FILE__,__LINE__,"exponent must be 1");
			if(value <= 1) throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			if(value > getMaxModulus() ) throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
			if( ! (value % 2) ) throw Failure(__func__,__FILE__,__LINE__,"modulus must be odd");
#endif

		}

		ModularBalanced (const ModularBalanced<int64_t>& mf) :
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
		const ModularBalanced &operator=(const ModularBalanced<int64_t> &F)
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

		uint64_t characteristic () const
		{
		       	return (uint64_t)modulus;
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
			return os << "balanced int64_t mod " << modulus;
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
			if( ! (modulus % 2) ) throw Failure(__func__,__FILE__,__LINE__,"modulus must be odd");
#endif

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

		Element &init (Element & x, const double &y) const
		{
			x = (Element) fmod(y,(double)modulus);
			if (x < mhalf_mod) return x += modulus;
			else if (x > half_mod) return x -= modulus;
			else return x;
		}

		inline Element &init (Element & x, const float &y) const
		{
			return init(x , (double) y);
		}

		template<class Element1>
		Element &init (Element & x, const Element1 &y) const
		{
			x = y % modulus;

			if ( x < mhalf_mod ) return x += modulus;
			else if (x > half_mod ) return x -= modulus;
                        else
				return x;
		}

		Element &init (Element &x, const size_t &y) const
		{
			x = (Element)y % Element(modulus);
			if (x < mhalf_mod) return x += modulus;
			else if (x > half_mod) return x -= modulus;
			else
				return x;
		}


		inline Element& init(Element& x, int y =0) const
		{
			x = y % modulus;

			if ( x < mhalf_mod ) return x += modulus;
			else if (x > half_mod ) return x -= modulus;
                        else
				return x;
		}

		inline Element& init(Element& x, long y) const
		{
			x = y % modulus;
			if ( x < mhalf_mod ) return x += modulus;
			else if ( x > half_mod ) return x -= modulus;
                        else
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
			if ( x > half_mod ) return x -= modulus;
			else if ( x < mhalf_mod ) return x += modulus;
                        else return x;
		}

		inline Element &sub (Element &x, const Element &y, const Element &z) const
		{
			x = y - z;
			if (x > half_mod) return x -= modulus;
			else if (x < mhalf_mod) return x += modulus;
			else return x;
		}

		inline Element &mul (Element &x, const Element &y, const Element &z) const
		{
			int64_t q;

			q  = (int64_t) ((((double) y) * ((double) z)) * modulusinv);  // q could be off by (+/-) 1
			x = (int64_t) (y*z - q*modulus);

			if (x > half_mod) return x -= modulus;
			else if (x < mhalf_mod) return x += modulus;
                        else return x;
		}

		inline Element &div (Element &x, const Element &y, const Element &z) const
		{
			Element temp;
			return mul (x, y, inv (temp, z) ) ;
		}

		inline Element &neg (Element &x, const Element &y) const
		{
			return x = -y;
		}

		inline Element &inv (Element &x, const Element &y) const
		{
			int64_t d;
			XINV(d, x, y, modulus);
#ifdef DEBUG
			if (d != 1)
				throw Failure(__func__,__FILE__,__LINE__,"InvMod: inverse undefined");
#endif
			if (x > half_mod) return x -= modulus;
			else if (x < mhalf_mod) return x += modulus;
                        else return x;

		}

		inline Element &axpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			int64_t q;

			q  = (int64_t) (((((double) a) * ((double) x)) + (double)y) * modulusinv);  // q could be off by (+/-) 1
			r = (int64_t) (a * x + y - q*modulus);


			if (r > half_mod) return r -= modulus;
			else if (r < mhalf_mod) return r += modulus;
                        else return r;

		}

		inline Element &axmy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			int64_t q;

			q  = (int64_t) (((((double) a) * ((double) x)) - (double)y) * modulusinv);  // q could be off by (+/-) 1
			r = (int64_t) (a * x - y - q*modulus);


			if (r > half_mod) return r -= modulus;
			else if (r < mhalf_mod) return r += modulus;
                        else return r;

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
			if ( x > half_mod ) return x -= modulus;
			else if (x < mhalf_mod) return x += modulus;
                        else return x;
		}

		inline Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			if (x > half_mod) return x -= modulus;
			else if (x < mhalf_mod) return x += modulus;
                        else return x;
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
			int64_t q;

			q  = (int64_t) (((((double) a)*((double) x)) + (double)r) * modulusinv);  // q could be off by (+/-) 1
			r = (int64_t) (a * x + r - q*modulus);


			if (r > half_mod) return r -= modulus;
			else if (r < mhalf_mod) return r += modulus;
                        else return r;
		}

		inline Element &maxpyin (Element &r, const Element &a, const Element &x) const
		{
			Element q;
                        maxpy(q,a,x,r);
                        return assign(r,q);
		}

		static inline int64_t getMaxModulus()
		{
// #if 1
// #ifdef __x86_64__
// 			return 4611686018427387904L; // 2^62
// 			// return 8589934591L;
// #else
// 			return 4611686018427387904LL; // 2^62
// 			// return 8589934591LL;
// #endif
// #endif
			// return 1 << 31;

            // (p-1)(p+1)/4 < 2^{63}
#ifdef __x86_64__
			return 6074001000L;
#else
			return 6074001000LL;
#endif
		}
		static  Element getMinModulus()	{return 3.0;}

	private:

		inline static int64_t& XINV(int64_t& d, int64_t& s, int64_t a, int64_t b)
		{
			int64_t  v, u2;
			int64_t aneg = 0;

			if (a < 0) {
#ifdef DEBUG
                            if (a < -LINBOX_MAX_INT64) throw Failure(__func__,__FILE__,__LINE__,"XINV: integer overflow");
#endif
                            v = -a;
                            aneg = 1;
			} else {
                            v = a;
                        }

			s = 0;
			u2 = 1;
			d = b;

			while (v != 0) {
				int64_t  q = d / v;
				int64_t  r = d % v;
				d = v;
				v = r;
				r = u2;
				u2 = s - q*u2;
				s = r;
			}

			if (aneg) return s = -s;
                        else return s;
		}

	};

}

#undef LINBOX_MAX_INT64

#include "field-general.h"

#endif //__FFLASFFPACK_modular_balanced_int64_H

