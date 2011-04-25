/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas-ffpack/field/modular-int.h
 * Copyright (C) 2003 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __FFLAFLAS_modular_int_H
#define __FFLAFLAS_modular_int_H


#include <math.h>
#include <sys/time.h>
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"

namespace FFPACK {
template< class Element >
class Modular;

template <>
class Modular<int> {

	typedef unsigned long long uint64;
	typedef unsigned long FieldInt;

protected:

	int modulus;
	double modulusinv;
	int _two64;

	public:

	typedef int Element;
	typedef ModularRandIter<int> RandIter;
	typedef NonzeroRandIter<Modular<int>,RandIter> NonZeroRandIter;

	const bool balanced;


	Modular (int value)  : modulus(value), balanced(false) {
		modulusinv = 1 / ((double) value);
		_two64 = (int) ((uint64) (-1) % (uint64) value);
		_two64 += 1;
		if (_two64 >= value) _two64 -= value;
	}

	Modular(const Modular<int>& mf) : modulus(mf.modulus),modulusinv(mf.modulusinv),_two64(mf._two64), balanced(false){}

	const Modular &operator=(const Modular<int> &F) {
		modulus = F.modulus;
		modulusinv = F.modulusinv;
		_two64 = F._two64;
		return *this;
	}


	FieldInt &cardinality (FieldInt &c) const{
		return c = modulus;
	}

	FieldInt &characteristic (FieldInt &c) const {
		return c = modulus;
	}

	FieldInt &convert (FieldInt &x, const Element &y) const {
		return x = y;
	}
	Element &convert (Element &x, const Element &y) const {
		return x = y;
	}

	double & convert (double &x, const Element &y) const {
		return x = (double) y;
	}

	float & convert (float &x, const Element &y) const {
		return x = (float) y;
	}

	std::ostream &write (std::ostream &os) const {
		return os << "int mod " << modulus;
	}

	std::istream &read (std::istream &is) {
		is >> modulus;
		modulusinv = 1 /((double) modulus );
		_two64 = (int) ((uint64) (-1) % (uint64) modulus);
		_two64 += 1;
		if (_two64 >= modulus) _two64 -= modulus;

		return is;
	}

	std::ostream &write (std::ostream &os, const Element &x) const {
		return os << x;
	}

	std::istream &read (std::istream &is, Element &x) const {
		FieldInt tmp;
		is >> tmp;
		init(x,tmp);
		return is;
	}


	template<class Element1>
	Element &init (Element & x, const Element1 &y) const {
		x = y % modulus;
		if (x < 0) x += modulus;
		return x;
	}

	Element &init (Element &x, const double &y) const  {
		double z = fmod(y, (double)modulus);
		if (z < 0) z += (double)modulus;
		z += 0.5;
		return x = static_cast<long>(z); //rounds towards 0
	}
	Element &init (Element &x, const float  &y) const  {
		float z = fmod(y, (float)modulus);
		if (z < 0) z += (float)modulus;
		z += 0.5;
		return x = static_cast<long>(z); //rounds towards 0
	}

	Element &init (Element &x, const FieldInt &y) const  {
		x = y % modulus;
		if (x < 0) x += modulus;
		return x;
	}

	inline Element& init(Element& x, int y =0) const {
		x = y % modulus;
		if ( x < 0 ) x += modulus;
		return x;
	}

	inline Element& init(Element& x, long y) const {
		x = y % modulus;
		if ( x < 0 ) x += modulus;
		return x;
	}

	inline Element& assign(Element& x, const Element& y) const {
		return x = y;
	}


	inline bool areEqual (const Element &x, const Element &y) const {
		return x == y;
	}

	inline  bool isZero (const Element &x) const {
		return x == 0;
	}

	inline bool isOne (const Element &x) const {
		return x == 1;
	}

	inline Element &add (Element &x, const Element &y, const Element &z) const {
		x = y + z;
		if ( x >= modulus ) x -= modulus;
		return x;
	}

	inline Element &sub (Element &x, const Element &y, const Element &z) const {
		x = y - z;
		if (x < 0) x += modulus;
		return x;
	}

	inline Element &mul (Element &x, const Element &y, const Element &z) const {
		int q;

		q  = (int) ((((double) y)*((double) z)) * modulusinv);  // q could be off by (+/-) 1
		x = (int) (y*z - q*modulus);


		if (x >= modulus)
			x -= modulus;
		else if (x < 0)
			x += modulus;

		return x;
	}

	inline Element &div (Element &x, const Element &y, const Element &z) const {
		Element temp;
		inv (temp, z);
		return mul (x, y, temp);
	}

	inline Element &neg (Element &x, const Element &y) const {
		if(y == 0) return x=0;
		else return x = modulus-y;
	}

	inline Element &inv (Element &x, const Element &y) const {
		int d, t;
		XGCD(d, x, t, y, modulus);
		if (x < 0)
			x += modulus;
		return x;

	}

	inline Element &axpy (Element &r,
			      const Element &a,
			      const Element &x,
			      const Element &y) const {
		int q;

		q  = (int) (((((double) a) * ((double) x)) + (double)y) * modulusinv);  // q could be off by (+/-) 1
		r = (int) (a * x + y - q*modulus);


		if (r >= modulus)
			r -= modulus;
		else if (r < 0)
			r += modulus;

		return r;

	}

	inline Element &addin (Element &x, const Element &y) const {
		x += y;
		if (  x >= modulus ) x -= modulus;
		return x;
	}

	inline Element &subin (Element &x, const Element &y) const {
		x -= y;
		if (x < 0) x += modulus;
		return x;
	}

	inline Element &mulin (Element &x, const Element &y) const {
		return mul(x,x,y);
	}

	inline Element &divin (Element &x, const Element &y) const {
		return div(x,x,y);
	}

	inline Element &negin (Element &x) const {
		if (x == 0) return x;
		else return x = modulus - x;
	}

	inline Element &invin (Element &x) const {
		return inv (x, x);
	}

	inline Element &axpyin (Element &r, const Element &a, const Element &x) const {
		int q;

		q  = (int) (((((double) a) * ((double) x)) + (double) r) * modulusinv);  // q could be off by (+/-) 1
		r = (int) (a * x + r - q*modulus);


		if (r >= modulus)
			r -= modulus;
		else if (r < 0)
			r += modulus;

		return r;
	}

private:

	static void XGCD(int& d, int& s, int& t, int a, int b) {
		int  u, v, u0, v0, u1, v1, u2, v2, q, r;

		int aneg = 0, bneg = 0;

		if (a < 0) {
			a = -a;
			aneg = 1;
		}

		if (b < 0) {
			b = -b;
			bneg = 1;
		}

		u1 = 1; v1 = 0;
		u2 = 0; v2 = 1;
		u = a; v = b;

		while (v != 0) {
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
} // FFPACK


#include "field-general.h"

#endif
