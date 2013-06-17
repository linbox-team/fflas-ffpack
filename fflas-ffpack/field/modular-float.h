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
 *.
 */

#ifndef __FFLASFFPACK_modular_float_H
#define __FFLASFFPACK_modular_float_H

#include <math.h>
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"
#include "fflas-ffpack/utils/debug.h"
#include <float.h>

namespace FFPACK {

	template <>
	class Modular<float> {

	public :
		typedef float Element;

	protected:

		Element         modulus;
		unsigned long   lmodulus;

		//Element inv_modulus;

	public:
		typedef unsigned long FieldInt;
		typedef ModularRandIter<Element> RandIter;
		typedef NonzeroRandIter<Modular<Element>, ModularRandIter<Element> > NonZeroRandIter;

		const Element one  ;
		const Element zero ;
		const Element mOne ;


		static const bool balanced = false;

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


		unsigned long &cardinality (unsigned long &c) const
		{
			return c = lmodulus;
		}

		unsigned long  cardinality() const
		{
			return lmodulus ;
		}

		unsigned long int & characteristic (long unsigned int& c) const
		{
			return c = lmodulus ;
		}

		unsigned long characteristic () const
		{
			return lmodulus;
		}


		unsigned long &convert (unsigned long &x, const Element &y) const
		{
			return x = (unsigned long)(y);
		}

		double  &convert (double &x, const Element &y) const
		{
			return x = y;
		}

		Element &convert (Element &x, const Element& y) const
		{
			return x=y;
		}

		std::ostream &write (std::ostream &os) const
		{
			return os << "float mod " << (int)modulus;
		}

		std::istream &read (std::istream &is)
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

		std::ostream &write (std::ostream &os, const Element &x) const
		{
			return os << x;
		}

		std::istream &read (std::istream &is, Element &x) const
		{
			float tmp = 0.f;
			is >> tmp;
			init(x,tmp);
			return is;
		}

		//!@warning possibly buggy. use % instead ?
		inline Element& init(Element& x, unsigned long int y) const
		{

			x = (Element)(y % lmodulus);

			if (x < 0) x += modulus;
			return x;
		}

		inline Element& init(Element& x, long int y) const
		{

			x = (Element)(y % (long)lmodulus);

			if (x < 0) x += modulus;
			return x;
		}

		inline Element& init(Element& x, Element y ) const
		{

			x = fmodf (y, modulus);
			if (x < 0) x += modulus;
			return x;
		}

		inline Element& init(Element& x, double y ) const
		{

			x = (Element)fmod (y, (double)modulus);
			if (x < 0) x += modulus;
			return x;
		}

		inline Element& init(Element& x, int y ) const
		{

			x = (Element) (y%(int) lmodulus);
			if (x < 0) x += modulus;
			return x;
		}

		inline Element& init(Element& x, unsigned int y ) const
		{

			x = (Element)(y % (unsigned int)lmodulus);
			if (x < 0) x += modulus;
			return x;
		}

		inline Element &init(Element& x) const
		{
			return x = 0. ;
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

		inline Element &add (Element &x, const Element &y, const Element &z) const
		{
			x = y + z;
			if ( x >= modulus ) x -= modulus;
			return x;
		}

		inline Element &sub (Element &x, const Element &y, const Element &z) const
		{
			x = y - z;
			if (x < 0) x += modulus;
			return x;
		}

		inline Element &mul (Element &x, const Element &y, const Element &z) const
		{
			Element tmp= y*z;
			x= fmodf(tmp, modulus);
			//x= tmp - floor(tmp*inv_modulus)*modulus;

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
			if(y == 0) return x = 0;
			else return x = modulus - y;
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

			if (tx < 0) tx += (int)modulus;

			// now x_int = gcd (modulus,residue)
			return x = (Element)tx;


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
 			if (r < 0.) r += modulus;
            return r;

		}

		inline Element &maxpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			Element tmp = y - a * x;
			r=fmodf(tmp, modulus);
			if (r < 0.) r += modulus;
            return r;

		}

		inline Element &addin (Element &x, const Element &y) const
		{
			x += y;
			if (  x >= modulus ) x -= modulus;
			return x;
		}

		inline Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			if (x < 0.) x += modulus;
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

		inline Element &maxpyin (Element &r,
				      const Element &a,
				      const Element &x) const
		{
			r -= a * x;
			r= fmodf(r, modulus);
			if (r < 0) r += modulus;
            return r;
		}

		static inline Element getMaxModulus()
		{
			return 4096.0f;  // floor( 2^12 )
			// return  1 << (FLT_MANT_DIG >> 1);  // 2^(DBL_MANT_DIG/2)
		}

	};

} // FFPACK


// const float FFPACK::Modular<float>::one  =  1UL;
// const float FFPACK::Modular<float>::zero =  0UL;



#include "field-general.h"

#endif // __FFLASFFPACK_modular_float_H
