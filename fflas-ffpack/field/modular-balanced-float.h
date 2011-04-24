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
 * See COPYING for license information.
 */

/*! @file field/modular-balanced-float.h
 * @ingroup field
 * @brief Balanced representation of <code>Z/mZ</code> over \c float .
 */

#ifndef __FFLAFLAS_modular_balanced_float_H
#define __FFLAFLAS_modular_balanced_float_H

#include <math.h>
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"

namespace FFPACK {

	template<>
	class ModularBalanced <float>{

	public:
		typedef float Element;

	protected:
		Element modulus;
		Element half_mod;
		Element mhalf_mod;
		unsigned long   lmodulus;


	public:
		typedef unsigned long FieldInt;
		typedef ModularBalancedRandIter<float> RandIter;
		typedef NonzeroRandIter<ModularBalanced<float>, RandIter> NonZeroRandIter;

		const bool balanced;

		ModularBalanced (int32_t p, int exp = 1) :
			modulus((Element)p),
			half_mod( Element((p-1)/2)),
			mhalf_mod( half_mod-p+1),
			lmodulus (p),
			balanced(true)
		{
#ifdef DEBUG
			if(modulus <= 1)
				throw Failure(__func__,__FILE__,
					      __LINE__,
					      "modulus must be > 1");
			if( exp != 1 ) throw Failure(__func__,__FILE__,
						     __LINE__,
						     "exponent must be 1");
			integer max;
			if (modulus > (Element) FieldTraits<ModularBalanced<Element> >::maxModulus(max))
				throw Failure (__func__,__FILE__,
					       __LINE__,
					       "modulus is too big");
#endif
		}

		ModularBalanced (Element p) :
			modulus (p),
			half_mod( Element((p-1)/2)),
			mhalf_mod( half_mod-p+1),
			lmodulus ((unsigned long)p),
			balanced(true)

		{
#ifdef DEBUG
			if (modulus <= 1)
				throw Failure(__func__,__FILE__,
					      __LINE__,
					      "modulus must be > 1");
			integer max;
			if (modulus > (Element) FieldTraits<ModularBalanced<Element> >::maxModulus(max))
				throw Failure (__func__,__FILE__,
					       __LINE__,
					       "modulus is too big");
#endif
		}

		ModularBalanced (FieldInt p) :
			modulus((Element)p),
			half_mod( Element((p-1)/2)),
			mhalf_mod( half_mod-p+1),
			lmodulus(p),
			balanced(true)

		{
#ifdef DEBUG
			if ((Element) modulus <= 1)
				throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			integer max;
			if ((Element) modulus > (Element) FieldTraits<ModularBalanced<Element> >::maxModulus(max))
				throw Failure (__func__,__FILE__,
					       __LINE__,
					       "modulus is too big");
#endif
		}

		ModularBalanced<float>(const ModularBalanced<float>& mf) :
			modulus(mf.modulus),
			half_mod(mf.half_mod),
			mhalf_mod(mf.mhalf_mod),
			lmodulus (mf.lmodulus),
			balanced(true)
		{}

		const ModularBalanced<float> &operator=(const ModularBalanced<float> &F) {
			modulus   = F.modulus;
			half_mod  = F.half_mod;
			mhalf_mod = F.mhalf_mod;
			lmodulus  = F.lmodulus;
			return *this;
		}

		FieldInt &cardinality (FieldInt &c) const
		{
			return c = (FieldInt) modulus;
		}

		FieldInt cardinality () const
		{
			return  (FieldInt) modulus;
		}


		long unsigned int &characteristic (long unsigned int &c) const
		{
			return c = (FieldInt) modulus;
		}

		FieldInt characteristic () const
		{
			return (FieldInt) modulus;
		}

		unsigned long &convert (unsigned long &x, const Element &y) const
		{
			if ( y < 0. ) return x= (unsigned long) (y + modulus) ;
			else return x = (unsigned long)y;
		}

		float &convert (float &x, const Element& y) const
		{
			return x=y;
		}
		double &convert (double &x, const Element& y) const
		{
			return x=y;
		}

		std::ostream &write (std::ostream &os) const
		{
			return os << "balanced Element mod " << int(modulus);
		}

		std::istream &read (std::istream &is) {
			is >> modulus;
#ifdef DEBUG
			if(modulus <= 1)
				throw PreconditionFailed (__func__,
							  __LINE__,
							  "modulus must be > 1");
			if(modulus > getMaxModulus())
				throw PreconditionFailed (__func__,
							  __LINE__,
							  "modulus is too big");
#endif
			return is;
		}

		std::ostream &write (std::ostream &os, const Element &x) const
		{
			return os << int(x);
		}

		std::istream &read (std::istream &is, Element &x) const
		{
			unsigned long tmp;
			is >> tmp;
			init(x,tmp);
			return is;
		}

		Element &init (Element &x, const unsigned long &y) const  {
			Element tmp  = Element(y % lmodulus);
			if (tmp > half_mod) return x =  tmp-modulus;
			else if (tmp<mhalf_mod) return x = tmp+modulus;
			else return x=tmp;
		}

		inline Element& init(Element& x, const double y =0) const
		{
			x = (Element) fmod (y, double(modulus));
			if ( x > half_mod ) return x -= modulus;
			else if ( x < mhalf_mod ) return x +=  modulus;
			else return x ;
		}

		inline Element& init(Element& x, const Element y =0) const
		{

			x = fmodf (y, modulus);

			if ( x > half_mod ) return x -= modulus;
			else if ( x < mhalf_mod ) return x +=  modulus;
			else return x ;
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
			return x == 0.;
		}

		inline bool isOne (const Element &x) const
		{
			return x == 1.;
		}

		inline Element &add (Element &x,
				     const Element &y, const Element &z) const
		{
			x = y + z;
			if ( x > half_mod ) return x -= modulus;
			if ( x < mhalf_mod ) return x += modulus;
			else return x;
		}

		inline Element &sub (Element &x,
				     const Element &y, const Element &z) const
		{
			x = y - z;
			if (x > half_mod) return x -= modulus;
			if (x < mhalf_mod) return x += modulus;
			else return x;
		}

		inline Element &mul (Element &x,
				     const Element &y, const Element &z) const
		{
			x = y * z;
			return init (x,x);
		}

		inline Element &div (Element &x, const Element &y, const Element &z) const
		{
			inv (x, z);
			return mulin (x, y);
		}

		inline Element &neg (Element &x, const Element &y) const
		{
			return x = -y;
		}

		inline Element &inv (Element &x, const Element &y) const
		{
			// The extended Euclidean algoritm
			int x_int, y_int, q, tx, ty, temp;
			x_int = int (modulus);
			y_int = (y < 0.) ? int(y + modulus) : int(y);
			tx = 0;
			ty = 1;

			while (y_int != 0) {
				// always: gcd (modulus,residue) = gcd (x_int,y_int)
				//         sx*modulus + tx*residue = x_int
				//         sy*modulus + ty*residue = y_int
				q = x_int / y_int; // integer quotient
				temp = y_int; y_int = x_int - q * y_int;
				x_int = temp;
				temp = ty; ty = tx - q * ty;
				tx = temp;
			}

			if (tx > half_mod ) return x = tx - modulus;
			else if ( tx < mhalf_mod ) return x = tx + modulus;
			return x = (Element) tx;
		}

		inline Element &axpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			r = a * x + y;
			return init (r, r);
		}

		inline Element &addin (Element &x, const Element &y) const
		{
			x += y;
			if ( x > half_mod ) return x -= modulus;
			else if ( x < mhalf_mod ) return x += modulus;
			return x;
		}

		inline Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			if ( x > half_mod ) return x -= modulus;
			else if ( x < mhalf_mod ) return x += modulus;
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
			return init (r, r);
		}

		static inline Element getMaxModulus()
		{
			return 2048.;  // 2^11
		}

	};

} // FFPACK

#endif // __FFLAFLAS_modular_balanced_double_H

