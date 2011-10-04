/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2010 LinBox
 * Adapted by B Boyer <brice.boyer@imag.fr>
 * (from other modular-balanced* files)
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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

/*! @file field/modular-int64.h
 * @ingroup field
 * @brief  representation of <code>Z/mZ</code> over \c int64_t .
 */
#ifndef __FFLASFFPACK_modular_int32_H
#define __FFLASFFPACK_modular_int32_H

#include <math.h>
#include <sys/time.h>
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"


#ifndef LINBOX_MAX_INT64
#ifdef __x86_64__
#define LINBOX_MAX_INT64 INT64_MAX
#else
#define LINBOX_MAX_INT64 INT64_MAX
#endif
#endif

// Namespace in which all LinBox code resides
namespace FFPACK
{

	template< class Element >
	class Modular;

	/** \brief Specialization of Modular to int64_t element type with efficient dot product.
	 *
	 * Efficient element operations for dot product, mul, axpy, by using floating point
	 * inverse of modulus (borrowed from NTL) and some use of non-normalized intermediate values.
	 *
	 * For some uses this is the most efficient field for primes in the range from half word
	 * to 2^62.
	 *
	 * Requires: Modulus < 2^62.
	 * Intended use: 2^30 < prime modulus < 2^62.
	 \ingroup field
	 */
	template <>
	class Modular<int64_t> {

	protected:

		int64_t modulus;
		double modulusinv;
		unsigned long lmodulus ;
		int64_t _two64 ;

	public:

		typedef int64_t Element;
		static const Element one  = 1 ;
		static const Element zero = 0 ;
		const Element mone ;


		typedef ModularRandIter<int64_t> RandIter;

		static const bool balanced = false ;

		//default modular field,taking 65521 as default modulus
		Modular () :
			modulus(65521),lmodulus(modulus)
			,mone(modulus -1)
		{
			modulusinv=1/(double)65521;
			_two64 = (int64) ((uint64) (-1) % (uint64) 65521);
			_two64 += 1;
			if (_two64 >= 65521) _two64 -= 65521;

		}

		Modular (int64_t value, int64_t exp = 1) :
			modulus(value),lmodulus(modulus)
			,mone(modulus -1)
		{
			modulusinv = 1 / ((double) value);
#ifdef DEBUG
			if(exp != 1) throw Failure(__func__,__FILE__,__LINE__,"exponent must be 1");
			if(value<=1) throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			int64_t max;
			FieldTraits<Modular<int64_t > >::maxModulus((uint64_t&)max) ;
			if( value > max ) throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
#endif
			 _two64 = (int64) ((uint64) (-1) % (uint64) value);
			 _two64 += 1;
			 if (_two64 >= value) _two64 -= value;

		}

		Modular(const Modular<int64_t>& mf) :
			modulus(mf.modulus),modulusinv(mf.modulusinv),lmodulus(modulus),_two64(mf._two64)
			,mone(modulus -1)
		{}

		const Modular &operator=(const Modular<int64_t> &F)
		{
			modulus    = F.modulus;
			modulusinv = F.modulusinv;
			lmodulus   = F.lmodulus;
			_two64     = F._two64;
			mone       = F.mone ;
			return *this;
		}

		inline unsigned long &cardinality ( unsigned long &c) const
		{
			return c = lmodulus;
		}

		inline unsigned long &characteristic (unsigned long &c) const
		{
			return c = lmodulus;
		}

		inline unsigned long characteristic () const
		{
			return lmodulus;
		}

		inline unsigned long cardinality () const
		{
			return lmodulus;
		}


		inline int64_t &convert (int64_t &x, const Element &y) const
		{
			return x = y;
		}

		inline double &convert (double &x, const Element &y) const
		{
			return x = (double) y;
		}

		inline float &convert (float &x, const Element &y) const
		{
			return x = (float) y;
		}

		inline std::ostream &write (std::ostream &os) const
		{
			return os << "int64_t mod " << modulus;
		}

		inline std::istream &read (std::istream &is)
		{
			is >> modulus;
			modulusinv = 1 /((double) modulus );
#ifdef DEBUG
			if(modulus <= 1) throw Failure(__func__,__FILE__,__LINE__,"modulus must be > 1");
			int64_t max;
			FieldTraits< Modular<int64_t> >::maxModulus((uint64_t&)max) ;
			if(modulus > max ) throw Failure(__func__,__FILE__,__LINE__,"modulus is too big");
#endif

			return is;
		}

		inline std::ostream &write (std::ostream &os, const Element &x) const
		{
			return os << x;
		}

		inline std::istream &read (std::istream &is, Element &x) const
		{
			unsigned long tmp;
			is >> tmp;
			init(x,tmp);
			return is;
		}

		inline Element &init (Element & x, const double &y) const
		{
			double z = fmod(y, (double)modulus);
			if (z < 0) z += (double)modulus;
			//z += 0.5; // C Pernet Sounds nasty and not necessary
			return x = static_cast<long>(z); //rounds towards 0
		}

		inline Element &init (Element & x, const float &y) const
		{
			return init(x , (double) y);
		}

		template<class Element1>
		inline Element &init (Element & x, const Element1 &y) const
		{
			x = y % modulus;
			if (x < 0) x += modulus;
			return x;
		}


		inline Element& init(Element& x, int y =0) const
		{
			x = y % modulus;
			if ( x < 0 ) x += modulus;
			return x;
		}

		inline Element& init(Element& x, long y) const
		{
			x = y % modulus;
			if ( x < 0 ) x += modulus;
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
			int64_t q;

			q  = (int64_t) ((((double) y)*((double) z)) * modulusinv);  // q could be off by (+/-) 1
			x = (int64_t) (y*z - q*modulus);


			if (x >= modulus)
				x -= modulus;
			else if (x < 0)
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
			if(y == 0) return x=0;
			else return x = modulus-y;
		}

		inline Element &inv (Element &x, const Element &y) const
		{
			int64_t d, t;
			XGCD(d, x, t, y, modulus);
			if (d != 1)
			{
#ifdef DEBUG
				throw Failure(__func__,__FILE__,__LINE__,"InvMod: Input is not invertible ");
#endif
			}
			if (x < 0)
				x += modulus;
			return x;

		}

		inline Element &axpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			int64_t q;

			q  = (int64_t) (((((double) a) * ((double) x)) + (double)y) * modulusinv);  // q could be off by (+/-) 1
			r = (int64_t) (a * x + y - q*modulus);


			if (r >= modulus)
				r -= modulus;
			else if (r < 0)
				r += modulus;

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
			if (x < 0) x += modulus;
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
			if (x == 0) return x;
			else return x = modulus - x;
		}

		inline Element &invin (Element &x) const
		{
			return inv (x, x);
		}

		inline Element &axpyin (Element &r, const Element &a, const Element &x) const
		{
			int64_t q;

			q  = (int64_t) (((((double) a) * ((double) x)) + (double) r) * modulusinv);  // q could be off by (+/-) 1
			r = (int64_t) (a * x + r - q*modulus);


			if (r >= modulus)
				r -= modulus;
			else if (r < 0)
				r += modulus;

			return r;
		}

		static inline int64_t getMaxModulus()
		{
#if 1
#ifdef __x86_64__
			return 4611686018427387904L;  // 2^62 in long long
#else
			return 4611686018427387904LL;  // 2^62 in long
#endif
#endif
			// return 1 << 31 ;
			// return 4294967296 ;
		}

	private:

		static void XGCD(int64_t& d, int64_t& s, int64_t& t, int64_t a, int64_t b)
		{
			int64_t  u, v, u0, v0, u1, v1, u2, v2, q, r;

			int64_t aneg = 0, bneg = 0;

			if (a < 0)
			{
#ifdef DEBUG
				if (a < -LINBOX_MAX_INT)
					throw Failure(__func__,__FILE__,__LINE__,"XGCD: integer overflow");
#endif
				a = -a;
				aneg = 1;
			}

			if (b < 0)
			{
#ifdef DEBUG
				if (b < -LINBOX_MAX_INT) throw
					Failure(__func__,__FILE__,__LINE__,"XGCD: integer overflow");
#endif
				b = -b;
				bneg = 1;
			}

			u1 = 1; v1 = 0;
			u2 = 0; v2 = 1;
			u = a; v = b;

			while (v != 0)
			{
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

#undef LINBOX_MAX_INT64

#include "field-general.h"

#endif //__LINBOX_modular_int64_H

