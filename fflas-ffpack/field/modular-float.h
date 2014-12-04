/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas-ffpack/modular-positive.h
 * Copyright (C) 2003 Pascal Giorgi
 *               2008 Clement Pernet
 * Written by Clement Pernet <clement.pernet@gmail.com>
 *            Pascal Giorgi  <pascal.giorgi@ens-lyon.fr>
 *
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

/*! @file field/modular-float.h
 * @ingroup field
 * @brief  representation of <code>Z/mZ</code> over \c float.
 */

#ifndef __FFLASFFPACK_modular_float_H
#define __FFLASFFPACK_modular_float_H

#include <math.h>
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"
#include "fflas-ffpack/utils/debug.h"
#include <float.h>

#define NORMALISE_LO(x) \
{ \
			if (x < 0.) x += modulus; \
}

#define NORMALISE_HI(x) \
{ \
			if (x >= modulus) x -= modulus; \
}


namespace FFPACK {

	template <>
	class Modular<float> {

	protected:

		float modulus;
		unsigned long   lmodulus;
		//Element inv_modulus;

	public:

		typedef float Element;
		typedef float* Element_ptr;
		typedef const float* ConstElement_ptr;

		const Element one  ;
		const Element zero ;
		const Element mOne ;

	public:

		static const bool balanced = false ;

		typedef unsigned long FieldInt;
		typedef ModularRandIter<Element> RandIter;
		typedef NonzeroRandIter<Modular<Element>, RandIter> NonZeroRandIter;



		Modular () :
			modulus(0),lmodulus(0)
			,one(0),zero(0),mOne(0)
		{}

		Modular (int32_t p, int exp = 1)  :
			modulus((Element)p), lmodulus((unsigned long)p)//, inv_modulus(1./(Element)
			,one(1.f),zero(0.f),mOne(p==2 ? 1.f : modulus -1.f)
		{
#ifdef DEBUG
			if(modulus <= 1.f)
				throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			if( exp != 1 ) throw Failure(__func__,__FILE__,__LINE__,"exponent must be 1");
			if(modulus > getMaxModulus() )
				throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
#endif

		}
		Modular (Element p) :
			modulus(p),  lmodulus((unsigned long)p)
			,one(1.f),zero(0.f),mOne(p==2.f ? 1.f : modulus -1.f)
		{
#ifdef DEBUG
			if( modulus <= 1.f )
				throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			if( modulus > getMaxModulus())
				throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
#endif
		}

		Modular (unsigned long int p) :
			modulus((Element)p), lmodulus(p)
			,one(1.f),zero(0.f),mOne(p==2 ? 1.f : modulus -1.f)
		{
#ifdef DEBUG
			if( (Element) modulus <= 1.f )
				throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			if( (Element) modulus > getMaxModulus())
				throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
#endif
		}

		Modular(const Modular<Element>& mf) :
			modulus(mf.modulus), lmodulus(mf.lmodulus)//inv_modulus(mf.inv_modulus)
			,one(mf.one),zero(mf.zero),mOne(mf.mOne)
		{}

		Modular <Element>& assign(const Modular<Element> &F)
		{
			modulus = F.modulus;
			lmodulus= F.lmodulus;
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
			lmodulus= F.lmodulus;
			//inv_modulus = F.inv_modulus;
			F.assign(const_cast<Element&>(one),F.one);
			F.assign(const_cast<Element&>(zero),F.zero);
			F.assign(const_cast<Element&>(mOne),F.mOne);
			return *this;
		}
#endif

		inline FieldInt &cardinality (FieldInt &c) const
		{
			return c = lmodulus;
		}

		inline FieldInt  cardinality () const
		{
			return lmodulus ;
		}

		inline FieldInt &characteristic (FieldInt &c) const
		{
			return c = lmodulus ;
		}

		inline FieldInt characteristic() const
		{
			return lmodulus;
		}

		inline FieldInt &convert (FieldInt &x, const Element &y) const
		{
			return x = (FieldInt)y;
		}

		inline Element  &convert (Element &x, const Element &y) const
		{
			return x = y;
		}

		inline double &convert (double &x, const Element& y) const
		{
			return x = (double)y;
		}

		inline std::ostream &write (std::ostream &os) const
		{
			return os << "float mod " << (int)modulus;
		}

		inline std::istream &read (std::istream &is)
		{
			is >> modulus;
#ifdef DEBUG
			if(modulus <= 1.f)
				throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			if(modulus > 94906265)
				throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
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

		//!@warning possibly buggy. use % instead ?
		inline Element& init(Element& x, unsigned long int y) const
		{

			x = (Element)(y % lmodulus);

			NORMALISE_LO(x);
			return x;
		}

		inline Element& init(Element& x, long int y) const
		{
			x = (Element)(y % (long)lmodulus);

			NORMALISE_LO(x);
			return x;
		}

		inline Element& init(Element& x, Element y ) const
		{
			return reduce (x,y);
		}

		inline Element& reduce (Element& x, Element y ) const
		{
			x = fmodf (y, modulus);
			NORMALISE_LO(x);
			return x;
		}

		inline Element& reduce (Element& x) const
		{
			x = fmodf (x, modulus);
			NORMALISE_LO(x);
			return x;
		}

		inline Element& init(Element& x, double y ) const
		{
			x = (Element)fmod (y, (double)modulus);
			NORMALISE_LO(x);
			return x;
		}

		inline Element& init(Element& x, int y ) const
		{
			x = (Element) (y%(int) lmodulus);
			NORMALISE_LO(x);
			return x;
		}

		inline Element& init(Element& x, unsigned int y ) const
		{
			x = (Element)(y % (unsigned int)lmodulus);
			NORMALISE_LO(x);
			return x;
		}

		inline Element &init(Element& x) const
		{
			return x = 0.f;
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
			Element tmp= y*z;
			x= fmodf(tmp, modulus);
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
			// The extended Euclidean algoritm
			int x_int, y_int, tx, ty;
			x_int = int (modulus);
			y_int = int (y);
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

			x = (Element) tx ;
			NORMALISE_LO(x);

			// now x_int = gcd (modulus,residue)
			return x ;


		}

		inline Element &axpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			Element tmp = a * x + y;
			return r= fmodf(tmp, modulus);
			//return r= tmp- floor(tmp*inv_modulus)*modulus;

		}

		inline Element &axmy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			Element tmp = a * x - y;
			r= fmodf(tmp, modulus);
			NORMALISE_LO(r);
            return r;

		}

		inline Element &maxpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			Element tmp = y - a * x;
			r=fmodf(tmp, modulus);
			NORMALISE_LO(r);
			return r;

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
			if (x == 0.) return x;
			else return x = modulus - x;
		}

		inline Element &invin (Element &x) const
		{
			return inv (x, x);
		}

		inline Element &axpyin (Element &r, const Element &a, const Element &x) const
		{
			r += a * x;
			return r = fmodf(r, modulus);

			//return r= tmp- floor(tmp*inv_modulus)*modulus;
		}

		inline Element &maxpyin (Element &r, const Element &a, const Element &x) const
		{
			r -= a * x;
			r= fmodf(r, modulus);
			NORMALISE_LO(r);
            return r;
		}

		static inline Element getMaxModulus()
		{
			return 2896.0f;  // floor( 2^11.5 ) (A.B+betaC must fit in the float mantissa)
		}

		static  Element getMinModulus()	{return 2.0;}

		Element minElement() const
		{
			return zero ;
		}

		Element maxElement() const
		{
			return mOne;
		}

	};

} // FFPACK




#undef NORMALISE_LO
#undef NORMALISE_HI

#include "field-general.h"

#endif // __FFLASFFPACK_modular_float_H
