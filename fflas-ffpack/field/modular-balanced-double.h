/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* field/modular-balanced-double.h
 * Copyright (C) 2003 Pascal Giorgi
 *               2005 Clement Pernet
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 * and Clement Pernet <Clement.Pernet@imag.fr>
 * and Brice Boyer <bboyer@imag.fr>
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

/*! @file field/modular-balanced-double.h
 * @ingroup field
 * @brief Balanced representation of <code>Z/mZ</code> over \c double .
 * @warning NOT DEFINED for EVEN modulus
 */

#ifndef __FFLASFFPACK_modular_balanced_double_H
#define __FFLASFFPACK_modular_balanced_double_H

#include <math.h>
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"
#include "fflas-ffpack/utils/debug.h"
#include <float.h>

namespace FFPACK
{

	//! it is forbiden to have char 2
	template<>
	class ModularBalanced <double>{

	protected:
		double modulus;
		double half_mod;
		double mhalf_mod;
		unsigned long lmodulus;

	public:
		typedef double Element;
		typedef double* Element_ptr;
		typedef const double* ConstElement_ptr;

		typedef unsigned long FieldInt;
		typedef ModularBalancedRandIter<double> RandIter;
		typedef NonzeroRandIter<ModularBalanced<double>, RandIter >  NonZeroRandIter;

		const Element one  ;
		const Element zero ;
		const Element mOne ; // except in char 2

		static const bool balanced = true;

		ModularBalanced (int32_t p, int exp = 1) :
			modulus((double)p),
			half_mod (double((p-1)/2)),
			mhalf_mod(half_mod-modulus+1),
			lmodulus ((unsigned int)p)
			,one(1),zero(0),mOne(-1)
		{
			FFLASFFPACK_check(isOdd(modulus));
#ifdef DEBUG
			if(modulus <= 1)
				throw FFPACK::Failure(__func__,__FILE__,
					      __LINE__,
					      "modulus must be > 1");
			if( exp != 1 ) throw Failure(__func__,__FILE__,
						     __LINE__,
						     "exponent must be 1");
			if (modulus > getMaxModulus())
				throw Failure (__func__,__FILE__,
					       __LINE__,
					       "modulus is too big");
#endif
		}

		ModularBalanced (double p) :
			modulus (p),
			half_mod (double((int)(p-1)/2)),
			mhalf_mod(half_mod-modulus+1),
			lmodulus ((unsigned long)p)
			,one(1),zero(0),mOne(-1)
		{
			FFLASFFPACK_check(isOdd(modulus));
#ifdef DEBUG
			if (modulus <= 1)
				throw Failure(__func__,__FILE__,
					      __LINE__,
					      "modulus must be > 1");
			if (modulus > getMaxModulus())
				throw Failure (__func__,__FILE__,
					       __LINE__,
					       "modulus is too big");
#endif
		}

		ModularBalanced (unsigned long p) :
			modulus ((double)p),
			half_mod (double((unsigned long)(p-1)/2)),
			mhalf_mod(half_mod-modulus+1),
			lmodulus (p)
			,one(1),zero(0),mOne(-1)
		{
			FFLASFFPACK_check(isOdd(modulus));
#ifdef DEBUG
			if (modulus <= 1)
				throw Failure(__func__,__FILE__,
					      __LINE__,
					      "modulus must be > 1");
			if (modulus > getMaxModulus())
				throw Failure (__func__,__FILE__,
					       __LINE__,
					       "modulus is too big");
#endif
		}

		ModularBalanced (long p) :
			modulus ((double)p),
			half_mod (double((unsigned long)(p-1)/2)),
			mhalf_mod(half_mod-modulus+1),
			lmodulus ((unsigned int)p)
			,one(1),zero(0),mOne(-1)
		{
			FFLASFFPACK_check(isOdd(modulus));
#ifdef DEBUG
			if (modulus <= 1)
				throw Failure(__func__,__FILE__,
					      __LINE__,
					      "modulus must be > 1");
			if (modulus > getMaxModulus())
				throw Failure (__func__,__FILE__,
					       __LINE__,
					       "modulus is too big");
#endif
		}


		ModularBalanced<double>(const ModularBalanced<double>& mf) :
			modulus(mf.modulus), half_mod(mf.half_mod)
			,mhalf_mod(mf.mhalf_mod), lmodulus(mf.lmodulus)
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
		const ModularBalanced<double> &operator=(const ModularBalanced<double> &F)
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


		FieldInt &cardinality (FieldInt &c) const
		{
			return c = (FieldInt) modulus;
		}

		FieldInt cardinality () const
		{
			return (FieldInt) modulus;
		}

		unsigned long characteristic() const
		{
			return lmodulus ;
		}

		FieldInt &characteristic (FieldInt &c) const
		{
			return c = (FieldInt) modulus;
		}

		unsigned long &convert (unsigned long &x, const Element &y) const
		{
			return x = (unsigned long)y;
		}

		double &convert (double &x, const Element& y) const
		{
			return x=y;
		}

		float &convert (float &x, const Element& y) const
		{
			return x=float(y);
		}

		std::ostream &write (std::ostream &os) const
		{
			return os << "balanced double mod " << (long int)modulus;
		}

		std::istream &read (std::istream &is)
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

		std::ostream &write (std::ostream &os, const Element &x) const
		{
			return os << (long) x;
		}

		std::istream &read (std::istream &is, Element &x) const
		{
			double tmp = 0.;
			is >> tmp;
			init(x,tmp);
			return is;
		}

		Element &init (Element &x, const unsigned long &y) const
		{
			x  = Element(y % lmodulus);
			if (x > half_mod) return x -=  modulus;
// 			else if (x<mhalf_mod) return x += modulus;
			else return x;
		}

		Element& init(Element& x, const double y) const
		{

			x = fmod (y, modulus);
			// x =(Element) ((long int)y%(long int)(modulus));
			if (x < mhalf_mod) return x +=  modulus;
			else if (x > half_mod) return x -=  modulus;
			else return x;
		}

		Element& init(Element& x) const
		{
			return x = 0.;
		}

		template<class T>
		Element& init(Element& x, const T y)  const
		{
			return init(x,double(y));
		}


		Element& assign(Element& x, const Element& y) const
		{
			return x = y;
		}

		/*! Tests equality.
		 * @param x element
		 * @param y element
		 * @warning \c x and \c y are supposed to be reduced.
		 */
		inline bool areEqual (const Element &x, const Element &y) const
		{
			return x == y;
		}

		inline  bool isZero (const Element &x) const
		{
			return x == 0.;
		}

		inline bool isOne (const Element &x) const
		{
			return x == 1.;
		}


		inline bool isMOne (const Element &x) const
		{
			return x == mOne;
		}

		inline Element &add (Element &x,
				     const Element &y,
				     const Element &z) const
		{
			x = y + z;
			if ( x > half_mod ) return x -= modulus;
			else if ( x < mhalf_mod ) return x += modulus;
			else return x;
		}

		inline Element &sub (Element &x,
				     const Element &y,
				     const Element &z) const
		{
			x = y - z;
			if (x > half_mod) return x -= modulus;
			else if (x < mhalf_mod) return x += modulus;
			else return x;
		}

		inline Element &mul (Element &x,
				     const Element &y, const Element &z) const
		{
			x = y * z;
			return init (x,x);
		}

		inline Element &div (Element &x,
				     const Element &y, const Element &z) const
		{
			Element temp;
			return mul (x, y, inv (temp, z) );
		}

		inline Element &neg (Element &x,
				     const Element &y) const
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
				int q,temp;
				// always: gcd (modulus,residue) = gcd (x_int,y_int)
				//         sx*modulus + tx*residue = x_int
				//         sy*modulus + ty*residue = y_int
				q = x_int / y_int; // integer quotient
				temp = y_int; y_int = x_int - q * y_int;
				x_int = temp;
				temp = ty; ty = tx - q * ty;
				tx = temp;
			}
                        x = (double)tx;
			if (x > half_mod ) return x -= modulus;
			else if ( x < mhalf_mod ) return x += modulus;
			else return x;
		}

		inline Element &axpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			r = a * x + y;
			return init (r, r);
		}

		inline Element &axmy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			r = a * x - y;
			return init (r, r);
		}

		inline Element &maxpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			r = y - a * x;
			return init (r, r);
		}

		inline Element &addin (Element &x, const Element &y) const
		{
			x += y;
			if ( x > half_mod ) return x -= modulus;
			else if ( x < mhalf_mod ) return x += modulus;
			else return x;
		}

		inline Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			if ( x > half_mod ) return x -= modulus;
			else if ( x < mhalf_mod ) return x += modulus;
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
			return x = - x;
		}

		inline Element &invin (Element &x) const
		{
			return inv (x, x);
		}

		inline Element &axpyin (Element &r, const Element &a, const Element &x) const
		{
			r += a * x;
			return init (r, r);
		}

		inline Element &maxpyin (Element &r, const Element &a, const Element &x) const
		{
			r -= a * x;
			return init (r, r);
		}

		static inline double getMaxModulus()
		{
                        // return 67108864.0;  // 2^26
			// return  1 << (DBL_MANT_DIG >> 1);  // 2^(DBL_MANT_DIG/2)
			// FFLASFFPACK_check(94906266LL*94906267LL>9007199254740991LL);
			// FFLASFFPACK_check(94906265LL*94906266LL<9007199254740991LL);
                        // (p-1)*(p+1) < 2^{53+2}
			return 189812531.0 ;
		}

	};


} // FFPACK

// const double FFPACK::ModularBalanced<double>::one  =  1;
// const double FFPACK::ModularBalanced<double>::mOne =  -1;
// const double FFPACK::ModularBalanced<double>::zero =  0;

#include "field-general.h"

#endif // __FFLASFFPACK_modular_balanced_double_H

