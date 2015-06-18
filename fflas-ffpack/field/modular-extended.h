/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) FFLAS group
 *
 * Written by Bastien Vialla <bastien.vialla@lirmm.fr>
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

#ifndef __FFLASFFPACK_MODULAR_EXTENDED_H
#define __FFLASFFPACK_MODULAR_EXTENDED_H

#include "givaro/givranditer.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-general.h"

// namespace Givaro{
//  template<class T>
//  class ModularExtended;// : public RingInterface<double>{};
// } // Givaro

namespace Givaro{
/*
 *
 * Modular double/float allowing big moduli
 * !!: RandIter does not works, use your own random
 *
 */
template<class _Element>
class ModularExtended// : public RingInterface<double>
{
public:

	typedef double Element;
	typedef Element* Element_ptr ;
	typedef const Element ConstElement;
	typedef const Element* ConstElement_ptr;
	// ----- Exported Types and constantes
	typedef ModularExtended<Element> Self_t;
	typedef uint64_t Residu_t;
	enum { size_rep = sizeof(Residu_t) };

private:
	// Verkampt Split
	inline void split(const Element x, Element &x_h, Element &x_l) const {
    	Element c;
    	if(std::is_same<Element, double>::value){
    		c = (Element)((1 << 27)+1);	
    	}else if(std::is_same<Element, float>::value){
    		c = (Element)((1 << 13)+1);	
    	}
    	 
    	x_h = (c*x)+(x-(c*x));
    	x_l = x - x_h;
	}	

	// Dekker mult, a * b = s + t
	inline void mult(const Element a, const Element b, Element &s, Element &t) const{
    	s = a*b;
//#ifdef __FMA__
    	t = std::fma(-a, b, s);
//#else
    	Element ah, al, bh, bl;
    	split(a, ah, al);
    	split(b, bh, bl);
    	t = ((((-s+ah*bh)+(ah*bl))+(al*bh))+(al*bl));
//#endif
	}

public:
	// ----- Constantes
	const Element zero = 0.0;
	const Element one = 1.0;
	const Element mOne = -1.0;

	// ----- Constructors
	ModularExtended() = default;

	template<class XXX> ModularExtended(const XXX& p)
	: zero(0.0), one(1.0), mOne((Element)p - 1.0), _p((Element)p), _invp(1/_p), _negp(-_p), _lp((Residu_t)p)
	{
	    assert(_p >= getMinModulus());
	    assert(_p <= maxCardinality());
	}

	//ModularExtended(const Self_t& F) = default;
	//ModularExtended(Self_t&& F) = default;
	// : zero(F.zero), one(F.one), mOne(F.mOne), _p(F._p), _lp(F._lp) {}

	// ----- Accessors
	inline Element minElement() const  { return zero; }
	inline Element maxElement() const  { return mOne; }

	// ----- Access to the modulus
	inline Residu_t residu() const { return _lp; }
	inline Residu_t size() const { return _lp; }
	inline Residu_t characteristic() const { return _lp; }
	template<class T> inline T& characteristic(T& p) const { return p = _lp; }
	inline Residu_t cardinality() const { return _lp; }
	template<class T> inline T& cardinality(T& p) const { return p = _lp; }
	static inline Residu_t maxCardinality() { 
		if(std::is_same<Element, double>::value)
			return 4503599627370496;
		else if(std::is_same<Element, float>::value)
			return 8388608;
	}
	static inline Residu_t getMinModulus() { return 2; }

	// ----- Checkers
	inline bool isZero(const Element& a) const  { return a == zero; }
	inline bool isOne (const Element& a) const  { return a == one; }
	inline bool isMOne(const Element& a) const  { return a == mOne; }
	inline bool areEqual(const Element& a, const Element& b) const  { return a == b; }
	inline size_t length(const Element a) const { return size_rep; }
	
	// ----- Ring-wise operators
	inline bool operator==(const Self_t& F) const { return _p == F._p; }
	inline bool operator!=(const Self_t& F) const { return _p != F._p; }
	inline Self_t& operator=(const Self_t& F)
	{
		F.assign(const_cast<Element&>(one),  F.one);
		F.assign(const_cast<Element&>(zero), F.zero);
		F.assign(const_cast<Element&>(mOne), F.mOne);
		_p = F._p;
		_negp = F._negp;
		_invp = F._invp;
		_lp= F._lp;
		return *this;
	}

	// ----- Initialisation
	Element &init (Element &x) const{
		return x = zero;
	}

	template<class XXX> Element& init(Element & x, const XXX & y) const{
		x=Element(y);
		return reduce(x);
	}

	Element &assign (Element &x, const Element &y) const{
		return x = y;
	}

	// ----- Convert and reduce
	Integer& convert  (Integer &x, const Element &y) const{
		return x = (Integer)y;
	}
	Residu_t& convert (Residu_t &x, const Element &y) const{
		return x = (Residu_t)y;
	}
	Element& convert   (Element &x, const Element &y) const{
		return x = y;
	}
	float& convert    (float &x, const Element &y) const{
		return x = (float)y;
	}

	Element& reduce (Element& x, const Element& y) const{
		Element q = floor(y*_invp);
		Element pqh, pql;
		mult(_p, q, pqh, pql);
		x = (x-pqh)-pql;
		if(x >= _p)
			x -= _p;
		else if(x < 0)
			x += _p;
		return x;	
	}
	Element& reduce (Element& x) const{
		Element q = floor(x*_invp);
		Element pqh, pql;
		mult(_p, q, pqh, pql);
		x = (x-pqh)-pql;
		if(x >= _p)
			x -= _p;
		else if(x < zero)
			x += _p;
		return x;	
	}

	// ----- Classic arithmetic
	Element& mul(Element& r, const Element& a, const Element& b) const {
		Element abh, abl, pqh, pql;
		mult(a, b, abh, abl);
		Element q = floor(abh*_invp);
		mult(_p, q, pqh, pql);		
		r = (abh-pqh)+(abl-pql);
		if(r > _p)
			r-= _p;
		else if(r < 0)
			r += _p;
		return r;
	}

	
	Element& div(Element& r, const Element& a, const Element& b) const{
		return mulin(inv(r, a), b);
	}
	Element& add(Element& r, const Element& a, const Element& b) const {
		r = a + b;
		if(r >= _p)
			r += _negp;
		return r;
	}
	Element& sub(Element& r, const Element& a, const Element& b) const {
		r = a - b;
		if(r < 0)
			r += _p;
		return r;
	}
	Element& neg(Element& r, const Element& a) const {
		r = -a;
		if(r < 0)
			r += _p;
		return r;
	}
	Element& inv(Element& x, const Element& y) const{
		int64_t x_int, y_int, tx, ty;
		x_int = int64_t(_lp);
		y_int = int64_t(y);
		tx = 0;
		ty = 1;

		while (y_int != 0) {
			// always: gcd (modulus,residue) = gcd (x_int,y_int)
			//         sx*modulus + tx*residue = x_int
			//         sy*modulus + ty*residue = y_int
			int64_t q = x_int / y_int; // integer quotient
			int64_t temp = y_int;  y_int  = x_int  - q * y_int;
			x_int  = temp;
			temp = ty; ty = tx - q * ty;
			tx = temp;
		}

		if (tx < 0) tx += int64_t(_p);

		// now x_int = gcd (modulus,residue)
		return x = Element(tx);
	}

	Element& mulin(Element& r, const Element& a) const {
		return mul(r, r, a);
	}
	Element& divin(Element& r, const Element& y) const{
		Element iy;
		return mulin(r, inv(iy, y));
	}
	Element& addin(Element& r, const Element& a) const {
		return add(r, r, a);
	}
	Element& subin(Element& r, const Element& a) const {
		return sub(r, r, a);
	}
	Element& negin(Element& r) const {
		return neg(r, r);
	}
	Element& invin(Element& r) const {
	  return inv(r, r);
	}
	
	// -- axpy:   r <- a * x + y
	// -- axpyin: r <- a * x + r
	Element& axpy  (Element& r, const Element& a, const Element& x, const Element& y) const {
		Element tmp;
		mul(tmp, a, x);
		return add(r, tmp, y);
	}
	Element& axpyin(Element& r, const Element& a, const Element& x) const {
		Element tmp(r);
		return axpy(r, a, x, tmp);
	}

	// -- axmy:   r <- a * x - y
	// -- axmyin: r <- a * x - r
	Element& axmy  (Element& r, const Element& a, const Element& x, const Element& y) const {
		Element tmp;
		mul(tmp, a, x);
		return sub(r, tmp, y);
	}
	Element& axmyin(Element& r, const Element& a, const Element& x) const {
		return axmy(r, a, x, r);
	}

	// -- maxpy:   r <- y - a * x
	// -- maxpyin: r <- r - a * x
	Element& maxpy  (Element& r, const Element& a, const Element& x, const Element& y) const {
		Element tmp;
		mul(tmp, a, x);
		return sub(r, y, tmp);
	}
	Element& maxpyin(Element& r, const Element& a, const Element& x) const {
		return maxpy(r, a, x, r);
	}

	// ----- Random generators
	// typedef ModularRandIter<Self_t> RandIter;
	// typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
 //    template< class Random > Element& random(const Random& g, Element& r) const { return init(r, g()); }
 //    template< class Random > Element& nonzerorandom(const Random& g, Element& a) const
 //    	{ while (isZero(init(a, g())));
 //    	  return a; }
		
protected:
	double _p = 0;
	double _invp = 0;
	double _negp = 0;
	Residu_t _lp = 0;

};

}// Givaro

#endif //__FFLASFFPACK_MODULAR_EXTENDED_H
