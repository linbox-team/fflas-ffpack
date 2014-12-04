/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* field/modular-balanced-float.h
 * Copyright (C) 2003 Pascal Giorgi
 *               2005,2008 Clement Pernet
 * Written by Clement Pernet <clement.pernet@gmail.com>
 *            Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 * Modified   Brice Boyer <bboyer@imag.fr>
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
 *
 */

/*! @file field/modular-balanced-float.h
 * @ingroup field
 * @brief Balanced representation of <code>Z/mZ</code> over \c float .
 * @warning NOT DEFINED for EVEN modulus
 */

#ifndef __FFLASFFPACK_modular_balanced_float_H
#define __FFLASFFPACK_modular_balanced_float_H

#include <math.h>
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"
#include "fflas-ffpack/utils/debug.h"
#include <float.h>

#define NORMALISE(x) \
{ \
	if (x < mhalf_mod) return x += modulus; \
	else if (x > half_mod) return x -= modulus; \
}

#define NORMALISE_HI(x) \
{ \
			if (x > half_mod) x -= modulus; \
}


namespace FFPACK {

	//! spec for float
	template<>
	class ModularBalanced<float> {
	protected:
		float modulus;
		float half_mod;
		float mhalf_mod;
		unsigned long   lmodulus;

	public:

		typedef float  Element;
		typedef float* Element_ptr;
		typedef const float* ConstElement_ptr;

		const Element one  ;
		const Element zero ;
		const Element mOne ;

	public:

		static const bool balanced = true;

		typedef unsigned long FieldInt;
		typedef ModularBalancedRandIter<Element> RandIter;
		typedef NonzeroRandIter<ModularBalanced<Element>, RandIter> NonZeroRandIter;


		ModularBalanced (int32_t p, int exp = 1) :
			modulus((Element)p),
			half_mod(Element((p-1)/2)),
			mhalf_mod( (Element) half_mod-modulus+1),
			lmodulus ((unsigned int)p)
			,one(1),zero(0),mOne(-1.f)
		{
			// FFLASFFPACK_check(isOdd(modulus));
#ifdef DEBUG
			if(modulus <= 1)
				throw Failure(__func__,__FILE__,
					      __LINE__,
					      "modulus must be > 1");
			if( exp != 1 ) throw Failure(__func__,__FILE__,
						     __LINE__,
						     "exponent must be 1");
			if ((Element) modulus > (Element) getMaxModulus())
				throw Failure (__func__,__FILE__,
					       __LINE__,
					       "modulus is too big");
#endif
		}

		ModularBalanced (Element p) :
			modulus (p),
			half_mod( Element(floor((p-1)/2))),
			mhalf_mod( half_mod-p+1),
			lmodulus ((unsigned long)p)
			,one(1),zero(0),mOne(-1.f)
		{
			FFLASFFPACK_check(isOdd(modulus));
#ifdef DEBUG
			if (modulus <= 1)
				throw Failure(__func__,__FILE__,
					      __LINE__,
					      "modulus must be > 1");
			if ((Element) modulus > (Element) getMaxModulus())
				throw Failure (__func__,__FILE__,
					       __LINE__,
					       "modulus is too big");
#endif
		}

		ModularBalanced (double p) :
			modulus (Element(p)),
			half_mod( Element(floor((p-1)/2))),
			mhalf_mod( half_mod-Element(p)+1),
			lmodulus ((unsigned long)p)
			,one(1),zero(0),mOne(-1.f)
		{
			FFLASFFPACK_check(isOdd(modulus));
#ifdef DEBUG
			if (modulus <= 1)
				throw Failure(__func__,__FILE__,
					      __LINE__,
					      "modulus must be > 1");
			if ((Element) modulus > (Element) getMaxModulus())
				throw Failure (__func__,__FILE__,
					       __LINE__,
					       "modulus is too big");
#endif
		}


		ModularBalanced (FieldInt p) :
			modulus((Element)p),
			half_mod( Element((p-1)/2)),
			mhalf_mod( (Element) half_mod-modulus+1),
			lmodulus(p)
			,one(1),zero(0),mOne(-1.f)
		{
			FFLASFFPACK_check(isOdd(modulus));
#ifdef DEBUG
			if ((Element) modulus <= 1)
				throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			if ((Element) modulus > (Element) getMaxModulus())
				throw Failure (__func__,__FILE__,
					       __LINE__,
					       "modulus is too big");
#endif
		}

		ModularBalanced<float>(const ModularBalanced<float>& mf) :
			modulus(mf.modulus),
			half_mod(mf.half_mod),
			mhalf_mod(mf.mhalf_mod),
			lmodulus (mf.lmodulus)
			,one(mf.one),zero(mf.zero),mOne(mf.mOne)
		{}

		ModularBalanced<Element> & assign(const ModularBalanced<Element> &F)
		{
			modulus = F.modulus;
			half_mod  = F.half_mod;
			mhalf_mod = F.mhalf_mod;
			lmodulus   = F.lmodulus;
			F.assign(const_cast<Element&>(one),F.one);
			F.assign(const_cast<Element&>(zero),F.zero);
			F.assign(const_cast<Element&>(mOne),F.mOne);
			return *this;
		}

#if 1
		const ModularBalanced<Element> &operator=(const ModularBalanced<Element> &F)
		{
			modulus = F.modulus;
			half_mod  = F.half_mod;
			mhalf_mod = F.mhalf_mod;
			lmodulus   = F.lmodulus;
			F.assign(const_cast<Element&>(one),F.one);
			F.assign(const_cast<Element&>(zero),F.zero);
			F.assign(const_cast<Element&>(mOne),F.mOne);
			return *this;
		}
#endif

		inline FieldInt &cardinality (FieldInt &c) const
		{
			return c = (FieldInt) modulus;
		}

		inline FieldInt cardinality () const
		{
			return  (FieldInt) modulus;
		}

		inline FieldInt &characteristic (FieldInt &c) const
		{
			return c = (FieldInt) modulus;
		}

		inline FieldInt characteristic() const
		{
			return (FieldInt) modulus;
		}

		inline FieldInt &convert (FieldInt &x, const Element &y) const
		{
			return x = (FieldInt)y;
		}

		inline Element &convert (Element &x, const Element& y) const
		{
			return x = y;
		}

		inline double &convert (double &x, const Element& y) const
		{
			return x = y;
		}

		inline std::ostream &write (std::ostream &os) const
		{
			return os << "balanced float mod " << int(modulus);
		}

		inline std::istream &read (std::istream &is)
		{
			is >> modulus;
#ifdef DEBUG
			if(modulus <= 1)
				throw Failure (__func__,
					       __LINE__,
					       "modulus must be > 1");
			if(modulus > getMaxModulus())
				throw Failure (__func__,
					       __LINE__,
					       "modulus is too big");
#endif

			return is;
		}

		inline std::ostream &write (std::ostream &os, const Element &x) const
		{
			return os  << x ;
		}

		inline std::istream &read (std::istream &is, Element &x) const
		{
			Element tmp;
			is >> tmp;
			init(x,tmp);
			return is;
		}

		inline Element &init (Element &x, const unsigned long &y) const
		{
			x = Element(y % lmodulus);
			NORMALISE_HI(x);
			return x;
		}

		inline Element &init (Element &x, const long &y) const
		{
			// no pb : float is small !
			x  = Element(y % (long) lmodulus);
			NORMALISE(x);
			return x;
		}

		inline Element &init (Element &x, const int &y) const
		{
			// no pb : float is small !
			x  = Element(y % (long) lmodulus);
			NORMALISE(x);
			return x;
		}

		inline Element& init(Element& x, const double y ) const
		{
			x = (Element) fmod (y, double(modulus));
			NORMALISE(x);
			return x ;
		}

		inline Element& init (Element& x, const Element y) const
		{
			return reduce (x,y);
		}

		inline Element& reduce (Element& x, const Element y) const
		{
			x = fmodf (y, modulus);
			NORMALISE(x);
			return x ;
		}

		inline Element& reduce (Element& x) const
		{
			x = fmodf (x, modulus);
			NORMALISE(x);
			return x ;
		}

		inline Element& init(Element& x) const
		{
			return x=0.f ;
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
			return x == 0.f;
		}

		inline bool isOne (const Element &x) const
		{
			return x == 1.f;
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
			x = y * z;
			return reduce (x);
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
			// The extended Euclidean algoritm
			int x_int, y_int, tx, ty;
			x_int = int (modulus);
			y_int = (y < 0.) ? int(y + modulus) : int(y);
			tx = 0;
			ty = 1;

			while (y_int != 0) {
				int q, temp;
				// always: gcd (modulus,residue) = gcd (x_int,y_int)
				//         sx*modulus + tx*residue = x_int
				//         sy*modulus + ty*residue = y_int
				q = x_int / y_int; // integer quotient
				temp = y_int; y_int = x_int - q * y_int;
				x_int = temp;
				temp = ty; ty = tx - q * ty;
				tx = temp;
			}
			x = (Element)tx;
			NORMALISE(x);
			return x;

		}

		inline Element &axpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			r = a * x + y;
			return reduce (r);
		}

		inline Element &axmy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			r = a * x - y;
			return reduce (r);
		}

		inline Element &maxpy (Element &r,
				       const Element &a,
				       const Element &x,
				       const Element &y) const
		{
			r = y - a * x;
			return reduce (r);
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
			r += a * x;
			return reduce (r);
		}

		inline Element &maxpyin (Element &r, const Element &a, const Element &x) const
		{
			r -= a * x;
			return reduce (r);
		}

		static inline Element getMaxModulus()
		{
			return 5792.0f;  // floor( 2^12.5 ) (A.B+betaC must fit in the float mantissa)
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

	};

} // FFPACK

#undef NORMALISE
#undef NORMALISE_HI

#include "field-general.h"
#endif // __FFLASFFPACK_modular_balanced_double_H

