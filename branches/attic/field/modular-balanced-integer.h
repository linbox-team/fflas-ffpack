/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas-ffpack/modular-positive.h
 * Copyright (C) 2014 The FFLAS-FFPACK group
 *
 * Written by BB <brice.boyer@lip6.fr>
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

/*! @file field/modular-balanced-integer.h
 * @ingroup field
 * @brief  balanced representation of <code>Z/mZ</code> over \c multiprecision integer.
 */


#ifndef __FFLASFFPACK_modular_balanced_integer_H
#define __FFLASFFPACK_modular_balanced_integer_H

#include <math.h>
#include "fflas-ffpack/utils/debug.h"
#include "fflas-ffpack/field/field-general.h"
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/integer.h"

// activate only if FFLAS-FFPACK haves multiprecision integer
#ifdef __FFLASFFPACK_HAVE_INTEGER

namespace FFPACK {

	template <>
	class ModularBalanced<Integer> {
	public:
		typedef Integer        	      	Element;
		typedef Element*        	Element_ptr;
		typedef const Element* 		ConstElement_ptr;

	protected:
		Element         modulus;


	public:
		const Element one  ;
		const Element zero ;
		const Element mOne ;

	protected:
		Element half_mod ;
		Element mhalf_mod ;

	public:

		typedef ModularBalancedRandIter<Element> RandIter;

		ModularBalanced () : modulus(0),one(0),zero(0),mOne(0)
				     ,half_mod(0),mhalf_mod(0)
		{}


		ModularBalanced (Element p) :
			modulus(p)
			,one(1),zero(0),mOne(-1)
				     ,half_mod((p-1)/2),mhalf_mod(half_mod-modulus+1)
		{
#ifdef DEBUG
			if( modulus <= 1 )
				throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
#endif
		}

		ModularBalanced (unsigned long int p) :
			modulus((Element)p),one(1),zero(0),mOne(-1)
				     ,half_mod((p-1)/2),mhalf_mod(half_mod-modulus+1)
		{
#ifdef DEBUG
			if( (Element) modulus <= 1 )
				throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
#endif
		}

		ModularBalanced(const ModularBalanced<Element>& mf) :
			modulus(mf.modulus),one(mf.one),zero(mf.zero),mOne(mf.mOne)
				     ,half_mod(mf.half_mod),mhalf_mod(mf.mhalf_mod)
		{}

		ModularBalanced<Element> & assign(const ModularBalanced<Element> &F)
		{
			modulus = F.modulus;
			F.assign(const_cast<Element&>(one),F.one);
			F.assign(const_cast<Element&>(zero),F.zero);
			F.assign(const_cast<Element&>(mOne),F.mOne);
			F.assign(const_cast<Element&>(half_mod),F.half_mod);
			F.assign(const_cast<Element&>(mhalf_mod),F.mhalf_mod);
			return *this;
		}

		const ModularBalanced &operator=(const ModularBalanced<Element> &F)
		{
			modulus = F.modulus;
			F.assign(const_cast<Element&>(one),F.one);
			F.assign(const_cast<Element&>(zero),F.zero);
			F.assign(const_cast<Element&>(mOne),F.mOne);
			F.assign(const_cast<Element&>(mhalf_mod),F.half_mod);
			F.assign(const_cast<Element&>(mhalf_mod),F.mhalf_mod);
			return *this;
		}

		Integer characteristic() const
		{
			return modulus ;
		}

		Integer &characteristic(Integer &c) const
		{
			return c = modulus ;
		}

		Integer &cardinality (Integer &c) const
		{
			return c = modulus ;
		}

		Integer cardinality () const
		{
			return modulus;
		}

		std::ostream &write (std::ostream &os) const
		{
			return os << "Integer balanced mod " << modulus;
		}

		std::istream &read (std::istream &is)
		{
			is >> modulus;
#ifdef DEBUG
			if(modulus <= 1)
				throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
#endif

			return is;
		}

		std::ostream &write (std::ostream &os, const Element &x) const
		{
			return os << x;
		}

		std::istream &read (std::istream &is, Element &x) const
		{
			Integer tmp;;
			is >> tmp;
			init(x,tmp);
			return is;
		}

		template<typename T>
		Element &init (Element &x, const T &y) const
		{
			x = y ;
			x%= modulus;
			if (x > half_mod) x -= modulus;

			return x;
		}

		Element& init(Element& x) const
		{
			return x=zero;
		}

		Element & convert(Element&x, const Element& y) const
		{
			return x = y ;
		}

		int64_t & convert(int64_t &x, const Element& y) const
		{
			return x = (int64_t) y ;
		}

		uint64_t & convert(uint64_t &x, const Element& y) const
		{
			return x = (uint64_t) y ;
		}


		double & convert(double &x, const Element& y) const
		{
			return x = (double) y ;
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
			return x == zero;
		}

		 bool isOne (const Element &x) const
		{
			return x == one;
		}

		 inline bool isMOne (const Element &x) const
		{
			return x == mOne;
		}

		 Element &add (Element &x, const Element &y, const Element &z) const
		{
			x = y + z;
			if ( x > half_mod ) x -= modulus;
			else if ( x < mhalf_mod ) return x += modulus;
			return x;
		}

		 Element &sub (Element &x, const Element &y, const Element &z) const
		{
			x = y - z;
			if ( x > half_mod ) x -= modulus;
			else if ( x < mhalf_mod ) return x += modulus;
			return x;
		}

		 Element &mul (Element &x, const Element &y, const Element &z) const
		{
			x = y*z;
			return init(x,x);
		}

		 Element &div (Element &x, const Element &y, const Element &z) const
		{
			Element temp;
			inv (temp, z);
			return mul (x, y, temp);
		}

		 Element &neg (Element &x, const Element &y) const
		{
			if(y == zero) return x = zero;
			else return x = - y; // hopefully
		}

		 Element &inv (Element &x, const Element &y) const
		 {
			 // The extended Euclidean algoritm
			 Integer x_int, y_int, tx, ty,q,temp;
			 y_int = y;
			 x_int = modulus;
			 tx = 0;
			 ty = 1;
			 while (y_int != zero) {
				 // always: gcd (modulus,residue) = gcd (x_int,y_int)
				 //         sx*modulus + tx*residue = x_int
				 //         sy*modulus + ty*residue = y_int
				 q = x_int / y_int; // Integer quotient
				 temp = y_int;  y_int  = x_int  - q*y_int;  x_int  = temp;
				 temp = ty; ty = tx - q*ty; tx = temp;
			 }
			 // now x_int = gcd (modulus,residue)
			 x = tx;

			 if (x > half_mod ) return x -= modulus;
			 else if ( x < mhalf_mod ) return x += modulus;

			 return x;
		}

		 Element &axpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			r = a * x + y;
			return init(r,r);
		}

		 Element &axmy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			r = a * x - y;
			return init(r,r);
		}

		 Element &maxpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			r = y - a * x;
			return init(r,r);
		}

		 Element &addin (Element &x, const Element &y) const
		{
			x += y;
			if (x > half_mod ) return x -= modulus;
			else if ( x < mhalf_mod ) return x += modulus;

			return x;
		}

		 Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			if (x > half_mod ) return x -= modulus;
			else if ( x < mhalf_mod ) return x += modulus;

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
			if (x == zero) return x;
			else return x =  - x;
		}

		 Element &invin (Element &x) const
		{
			return inv (x, x);
		}

		 Element &axpyin (Element &r, const Element &a, const Element &x) const
		{
			r = r + a * x;
			return init(r,r);

		}

		 Element &maxpyin (Element &r, const Element &a, const Element &x) const
		{
			r = r - a * x;
			init(r,r);
			// if (x > half_mod ) return x -= modulus;
			// else if ( x < mhalf_mod ) return x += modulus;

			return r;
		}

		Element minElement() const
		{
			return mhalf_mod;
		}

		Element maxElement() const
		{
			return half_mod;
		}

	};

	template <>
        class ModularBalancedRandIter<Integer> {
        public:
		typedef Integer Element;
                ModularBalancedRandIter (const ModularBalanced<Element> &F, size_t seed=0) :
                        _F(F)
                {
                        if (seed==0) {
                                struct timeval tp;
                                gettimeofday(&tp, 0) ;
                                seed = (size_t) tp.tv_usec;
                        }
			Integer::seeding(seed);
                }

                ModularBalancedRandIter (const ModularBalancedRandIter<Element> &R) :
                        _F (R._F)
                {}

                /*! @bug not so random... (at all)  */
                Element &random (Element &a) const
                {
                        return _F.init(a,Integer::random(a,_F.cardinality()));
                }

        private:
                ModularBalanced<Integer> _F;

        }; // end of class ModularBalanced<Integer>

} // FFPACK

#endif // __FFLASFFPACK_HAVE_INTEGER

#endif // __FFLASFFPACK_modular_balanced_integer_H
