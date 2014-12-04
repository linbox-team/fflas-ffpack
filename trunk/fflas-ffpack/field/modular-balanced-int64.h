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

#define NORMALISE(x) \
{ \
	if (x < mhalf_mod) return x += modulus; \
	else if (x > half_mod) return x -= modulus; \
}

#define NORMALISE_HI(x) \
{ \
			if (x > half_mod) x -= modulus; \
}

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

		const Element one  ;
		const Element zero ;
		const Element mOne ;

	public:

		static const bool balanced = true ;

		typedef ModularBalancedRandIter<Element> RandIter;
		typedef NonzeroRandIter<ModularBalanced<Element>, RandIter> NonZeroRandIter;

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

		ModularBalanced (Element value, int32_t exp = 1) :
			modulus(value)
			,one(1),zero(0),mOne(-1)
		{
			half_mod = (modulus >> 1);
			mhalf_mod = half_mod-modulus+1;
			modulusinv = 1 / ((double) value);
			// FFLASFFPACK_check(isOdd(modulus));
#ifdef DEBUG
			if(exp != 1) throw Failure(__func__,__FILE__,__LINE__,"exponent must be 1");
			if(value <= 1) throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			if(value > getMaxModulus() ) throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
			// if( ! (value % 2) ) throw Failure(__func__,__FILE__,__LINE__,"modulus must be odd");
#endif

		}

		ModularBalanced (const ModularBalanced<Element>& mf) :
			modulus(mf.modulus),
			half_mod(mf.half_mod),
			mhalf_mod(mf.mhalf_mod),
			modulusinv(mf.modulusinv)
			,one(mf.one),zero(mf.zero),mOne(mf.mOne)
		{
		}

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
		const ModularBalanced &operator=(const ModularBalanced<Element> &F)
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

		inline uint64_t &cardinality ( uint64_t &c) const
		{
			return c = (uint64_t)modulus;
		}

		inline uint64_t cardinality () const
		{
			return (uint64_t) modulus;
		}

		inline uint64_t &characteristic (uint64_t &c) const
		{
			return c = (uint64_t) modulus;
		}

		inline uint64_t characteristic () const
		{
		       	return (uint64_t)modulus;
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
			return os << "balanced int64_t mod " << modulus;
		}

		inline std::istream &read (std::istream &is)
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
			NORMALISE(x);
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
			NORMALISE(x);
			return x;
		}

		inline Element& reduce (Element & x, const Element &y) const
		{
			x = (y % modulus);
			NORMALISE(x);
			return x;
		}

		inline Element& reduce (Element & x) const
		{
			x %= modulus;
			NORMALISE(x);
			return x;
		}

		inline Element& init(Element&x) const
		{
			return x = 0;
		}

		inline Element& init(Element& x, Element y) const
		{
			x = (y % modulus);
			NORMALISE(x);
			return x;
		}

		inline Element& init(Element& x, int32_t y ) const
		{
			x = (Element)y % modulus;
			NORMALISE(x);
			return x;
		}

		inline Element& init(Element& x, uint32_t y ) const
		{
			x = y % modulus;
			NORMALISE_HI(x);
			return x;
		}

		inline Element& init (Element& x, uint64_t y) const
		{
			x = (Element)y % (uint64_t)modulus;
			NORMALISE_HI(x);
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
			NORMALISE(x);
                        return x;
		}

		inline Element &sub (Element &x, const Element &y, const Element &z) const
		{
			x = y - z;
			NORMALISE(x);
			return x;
		}

		inline Element &mul (Element &x, const Element &y, const Element &z) const
		{
			Element q;

			q  = (Element) ((((double) y) * ((double) z)) * modulusinv);  // q could be off by (+/-) 1
			x = (Element) (y*z - q*modulus);

			NORMALISE(x);
                        return x;
		}

		inline Element &div (Element &x, const Element &y, const Element &z) const
		{
			Element temp;
			return mul (x, y, inv(temp,z));
		}

		inline Element &neg (Element &x, const Element &y) const
		{
			return x = -y;
		}

		inline Element &inv (Element &x, const Element &y) const
		{
			Element d;
			XINV(d, x, y, modulus);
#ifdef DEBUG
			if (d != 1)
				throw Failure(__func__,__FILE__,__LINE__,"InvMod: inverse undefined");
#endif
			NORMALISE(x);
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


			NORMALISE(r);

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


			NORMALISE(r);
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
			NORMALISE(x);
                        return x;
		}

		inline Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			NORMALISE(x);
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
			Element q;

			q  = (Element) (((((double) a) * ((double) x)) + (double) r) * modulusinv);  // q could be off by (+/-) 1
			r = (Element) (a * x + r - q*modulus);


			NORMALISE(r);

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
#ifdef __x86_64__
			return 6074001000L;
#else
			return 6074001000LL;
#endif
		}

		static  Element getMinModulus()	{return 3.0;}

		Element minElement() const
		{
			return mhalf_mod ;
		}

		Element maxElement() const
		{
			return half_mod ;
		}

	private:
            inline static Element& XINV(Element& q, Element& u, Element a, Element b) {
                int64_t u3;
                int64_t v1,v3;
                u = 1; u3 = a;
                v1 = 0; v3 = b;
                while (v3 != 0)
                {
                    int64_t t1, t3;
                    q = u3 / v3;
                    t1 = u - q * v1; t3 = u3 - q * v3;
                    u = v1; u3 = v3; v1 = t1; v3 = t3;
                }
                return (u3<0?u=-u:u);
            }

	};

}

#undef LINBOX_MAX_INT64
#undef NORMALISE
#undef NORMALISE_HI

#include "field-general.h"

#endif //__FFLASFFPACK_modular_balanced_int64_H

