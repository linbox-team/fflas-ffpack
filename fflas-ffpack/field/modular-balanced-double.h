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
 * See COPYING for license information.
 */

/*! @file field/modular-balanced-double.h
 * @ingroup field
 * @brief Balanced representation of <code>Z/mZ</code> over \c double .
 */

#ifndef __FFLAFLAS_modular_balanced_double_H
#define __FFLAFLAS_modular_balanced_double_H

#include <math.h>
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"
#include "fflas-ffpack/utils/debug.h"

namespace FFPACK
{

	template<>
	class ModularBalanced <double>{

	protected:
		double modulus;
		double half_mod;
		double mhalf_mod;
		unsigned long lmodulus;

	public:
		typedef double Element;
		typedef unsigned long FieldInt;
		typedef ModularBalancedRandIter<double> RandIter;
		typedef NonzeroRandIter<ModularBalanced<double>, RandIter >  NonZeroRandIter;


		const bool balanced;

		ModularBalanced (int32_t p, int exp = 1) :
			modulus((double)p),
			half_mod (double((p-1)/2)),
			mhalf_mod(half_mod-modulus+1),
			lmodulus (p),
			balanced(true)
		{
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
			lmodulus ((unsigned long)p),
			balanced(true)
		{
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
			lmodulus (p),
			balanced(true)
		{
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
			,mhalf_mod(mf.mhalf_mod), lmodulus(mf.lmodulus), balanced(true)
		{}

		const ModularBalanced<double> &operator=(const ModularBalanced<double> &F)
		{
			modulus   = F.modulus;
			half_mod  = F.half_mod;
			mhalf_mod = F.mhalf_mod;
			lmodulus  = F.lmodulus;
			// balanced  = F.balanced;
			return *this;
		}


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
			if ( y < 0. ) return x= (unsigned long) (y + modulus) ;
			else return x = (unsigned long)y;
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
			unsigned long tmp;
			is >> tmp;
			init(x,tmp);
			return is;
		}

		Element &init (Element &x, const unsigned long &y) const
		{
			x  = Element(y % lmodulus);
			if (x > half_mod) return x -=  modulus;
			else if (x<mhalf_mod) return x += modulus;
			else return x;
		}

		Element& init(Element& x, const double y =0) const
		{

			x = fmod (y, modulus);
			if (x > half_mod) return x -=  modulus;
			if (x < mhalf_mod) return x +=  modulus;
			return x;
		}

		template<class T>
		Element& init(Element& x, const T y =0) const
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

		inline Element &add (Element &x,
				     const Element &y,
				     const Element &z) const
		{
			x = y + z;
			if ( x > half_mod ) return x -= modulus;
			if ( x < mhalf_mod ) return x += modulus;
			return x;
		}

		inline Element &sub (Element &x,
				     const Element &y,
				     const Element &z) const
		{
			x = y - z;
			if (x > half_mod) return x -= modulus;
			if (x < mhalf_mod) return x += modulus;
			return x;
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
			inv (temp, z);
			return mul (x, y, temp);
		}

		inline Element &neg (Element &x,
				     const Element &y) const
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
			else return x = (double) tx;
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
			if ( x < mhalf_mod ) return x += modulus;
			return x;
		}

		inline Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			if ( x > half_mod ) return x -= modulus;
			if ( x < mhalf_mod ) return x += modulus;
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

		static inline double getMaxModulus()
		{
			return 67108864.0;  // 2^26
		}

	};


} // FFPACK

#include "field-general.h"

#endif // __FFLAFLAS_modular_balanced_double_H

