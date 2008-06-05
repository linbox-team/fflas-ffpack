/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* field/modular-balanced.h
 * Copyright (C) 2003 Pascal Giorgi
 *               2005 Clement Pernet
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 * and Clement Pernet <Clement.Pernet@imag.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __MODULAR_BALANCED_H
#define __MODULAR_BALANCED_H

#include <math.h>
#include "fflas-ffpack/modular-randiter.h"
#include "fflas-ffpack/nonzero-randiter.h"

template <class Element>
class ModularBalanced;

template<>
class ModularBalanced <double>{

protected:
	double modulus;
	double half_mod;
	unsigned long lmodulus;
       
public:	       
	typedef double Element;
	typedef unsigned long FieldInt;
	typedef ModularBalancedRandIter<double> RandIter;
	typedef NonzeroRandIter<ModularBalanced<double>, ModularBalancedRandIter<double> >  NonZeroRandIter;

	const bool balanced;

	ModularBalanced<double> (int p)	: modulus((double)p),
		half_mod( (p-1.)/2),
		balanced(true) {}

	ModularBalanced<double>(const ModularBalanced<double>& mf)
	: modulus(mf.modulus), half_mod(mf.half_mod), balanced(true){}
	
	const ModularBalanced<double> &operator=(const ModularBalanced<double> &F) {
			modulus = F.modulus;
			half_mod = F.half_mod;
			return *this;
		}

	
	
	
	
	FieldInt &cardinality (FieldInt &c) const{ 
		return c = (FieldInt) modulus;
	}

	FieldInt &characteristic (FieldInt &c) const {
		return c = (FieldInt) modulus; 
	}

	unsigned long &convert (unsigned long &x, const Element &y) const { 
		if ( y < 0. ) return x= (unsigned long) (y + modulus) ;
		else return x = (unsigned long)y;
	}
	
	double &convert (double &x, const Element& y) const {
		return x=y;
	}
		
	std::ostream &write (std::ostream &os) const {
		return os << "double mod " << (int)modulus;
	}
	
	std::istream &read (std::istream &is) {
		is >> modulus; 
		return is;
	}
		
	std::ostream &write (std::ostream &os, const Element &x) const {
		return os << x;
	}

	std::istream &read (std::istream &is, Element &x) const {
		unsigned long tmp;
		is >> tmp;
		init(x,tmp); 
		return is;
	}
		
	Element &init (Element &x, const unsigned long &y) const  {
		Element tmp  = y % (unsigned long)  (modulus);
		if (tmp > half_mod) return x =  tmp-modulus;
		else if (tmp<-half_mod) return x = tmp+modulus;
		else return x=tmp;
	}

	inline Element& init(Element& x, const double y =0) const {		  
		double tmp;
		
		//tmp = floor (y + 0.5);
		//tmp = fmod (tmp, modulus);
		tmp = fmod (y, modulus);
		
		if ( tmp > half_mod ) return x = tmp - modulus;
		else if ( tmp <-half_mod ) return x = tmp + modulus;
		else return x = tmp;
	}

		
		
	inline Element& assign(Element& x, const Element& y) const {
		return x = y;
	}
									
		
	inline bool areEqual (const Element &x, const Element &y) const {
		return x == y;
	}

	inline  bool isZero (const Element &x) const {
		return x == 0.; 
	}
		
	inline bool isOne (const Element &x) const {
		return x == 1.; 
	}

	inline Element &add (Element &x, const Element &y, const Element &z) const {
		x = y + z;
		if ( x > half_mod ) return x -= modulus;
		if ( x < -half_mod ) return x += modulus; 
		return x;
	}
 
	inline Element &sub (Element &x, const Element &y, const Element &z) const {
		x = y - z;
		if (x > half_mod) return x -= modulus;
		if (x < -half_mod) return x += modulus;
		return x;
	}
		
	inline Element &mul (Element &x, const Element &y, const Element &z) const {		
		double tmp= y*z;
		return init (x, tmp);
	}
 
	inline Element &div (Element &x, const Element &y, const Element &z) const {
		Element temp;
		inv (temp, z);
		return mul (x, y, temp);
	}
 
	inline Element &neg (Element &x, const Element &y) const {
		//			if (y == 0) return x=0;
		return x = -y;
	}
 
	inline Element &inv (Element &x, const Element &y) const {
		// The extended Euclidean algoritm 
		int x_int, y_int, q, tx, ty, temp;
		x_int = (int) modulus;
		y_int = (int) ((y<0.)?y+modulus:y);
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
		else if ( tx < -half_mod ) return x = tx + modulus;
		else return x = (double) tx;
	}

	inline Element &axpy (Element &r, 
			      const Element &a, 
			      const Element &x, 
			      const Element &y) const {
		double tmp= a*x+y;
		return init( r, tmp);
	}

	inline Element &addin (Element &x, const Element &y) const {
		x += y;
		if ( x > half_mod ) return x -= modulus;
		if ( x < -half_mod ) return x += modulus; 
		return x;
	}
 
	inline Element &subin (Element &x, const Element &y) const {
		x -= y;
		if ( x > half_mod ) return x -= modulus;
		if ( x < -half_mod ) return x += modulus; 
		return x;
	}
 
	inline Element &mulin (Element &x, const Element &y) const {
		return mul(x,x,y);
	}
 
	inline Element &divin (Element &x, const Element &y) const {
		return div(x,x,y);
	}
 
	inline Element &negin (Element &x) const {
		//if (x == 0.) return x; 
		return x = - x; 
	}
		
	inline Element &invin (Element &x) const {
		return inv (x, x);
	}
		
	inline Element &axpyin (Element &r, const Element &a, const Element &x) const {
		r += a * x;
		return init( r, r);
	}

	static inline double getMaxModulus()
	{ return 67108864.0; } // 2^26
		
};
template<>
class ModularBalanced <float>{

protected:
	float modulus;
	float half_mod;

       
public:	       
	typedef float Element;
	typedef unsigned long FieldInt;
	typedef ModularBalancedRandIter<float> RandIter;
	typedef NonzeroRandIter<ModularBalanced<float>, RandIter> NonZeroRandIter;

	const bool balanced;

	ModularBalanced<float> (int p)  : modulus((float)p), half_mod( (p-1.)/2), balanced(true) {}

	ModularBalanced<float>(const ModularBalanced<float>& mf) : modulus(mf.modulus), half_mod(mf.half_mod), balanced(true){}
	
	const ModularBalanced<float> &operator=(const ModularBalanced<float> &F) {
			modulus = F.modulus;
			half_mod = F.half_mod;
			return *this;
		}

	FieldInt &cardinality (FieldInt &c) const{ 
		return c = (FieldInt) modulus;
	}

	FieldInt &characteristic (FieldInt &c) const {
		return c = (FieldInt) modulus; 
	}

	unsigned long &convert (unsigned long &x, const Element &y) const { 
		if ( y < 0. ) return x= (unsigned long) (y + modulus) ;
		else return x = (unsigned long)y;
	}
	
	float &convert (float &x, const Element& y) const {
		return x=y;
	}
		
	std::ostream &write (std::ostream &os) const {
		return os << "int mod " << (int)modulus;
	}
	
	std::istream &read (std::istream &is) {
		is >> modulus; 
		return is;
	}
		
	std::ostream &write (std::ostream &os, const Element &x) const {
		return os << x;
	}

	std::istream &read (std::istream &is, Element &x) const {
		unsigned long tmp;
		is >> tmp;
		init(x,tmp); 
		return is;
	}
		
	Element &init (Element &x, const unsigned long &y) const  {
		Element tmp  = y % (unsigned long)  (modulus);
		if (tmp > half_mod) return x =  tmp-modulus;
		else if (tmp<-half_mod) return x = tmp+modulus;
		else return x=tmp;
	}

	inline Element& init(Element& x, float y =0) const {		  
		float tmp=y;
		
		//tmp = floor (y + 0.5);
		tmp = fmod (tmp, modulus);
		//tmp = fmod (y, modulus);
		
		if ( tmp > half_mod ) return x = tmp - modulus;
		else if ( tmp <-half_mod ) return x = tmp + modulus;
		else return x = tmp;
	}

	
	inline Element& init(Element& x, double y =0) const {
		double tmp=y;
		
		//tmp = floor (y + 0.5);
		//tmp = fmod (tmp, modulus);
		tmp = fmod (y, double(modulus));
		
		if ( tmp > half_mod ) return x = tmp - modulus;
		else if ( tmp <-half_mod ) return x = tmp + modulus;
		else return x = tmp;
	}
		
	inline Element& assign(Element& x, const Element& y) const {
		return x = y;
	}
									
		
	inline bool areEqual (const Element &x, const Element &y) const {
		return x == y;
	}

	inline  bool isZero (const Element &x) const {
		return x == 0.; 
	}
		
	inline bool isOne (const Element &x) const {
		return x == 1.; 
	}

	inline Element &add (Element &x, const Element &y, const Element &z) const {
		x = y + z;
		if ( x > half_mod ) return x -= modulus;
		if ( x < -half_mod ) return x += modulus; 
		else return x;
	}
 
	inline Element &sub (Element &x, const Element &y, const Element &z) const {
		x = y - z;
		if (x > half_mod) return x -= modulus;
		else if (x < -half_mod) return x += modulus;
		else return x;
	}
		
	inline Element &mul (Element &x, const Element &y, const Element &z) const {		
		x = y*z;
		return init (x, x);
	}
 
	inline Element &div (Element &x, const Element &y, const Element &z) const {
		inv (x, z);
		return mulin (x, y);
	}
 
	inline Element &neg (Element &x, const Element &y) const {
		//			if (y == 0) return x=0;
		return x = -y;
	}
 
	inline Element &inv (Element &x, const Element &y) const {
		// The extended Euclidean algoritm 
		int x_int, y_int, q, tx, ty, temp;
		x_int = (int) modulus;
		y_int = (int) ((y<0.)?y+modulus:y);
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
		else if ( tx < -half_mod ) return x = tx + modulus;
		else return x = (float) tx;
	}

	inline Element &axpy (Element &r, 
			      const Element &a, 
			      const Element &x, 
			      const Element &y) const {
		r = a * x + y;
		return init( r, r);

	}

	inline Element &addin (Element &x, const Element &y) const {
		x += y;
		if ( x > half_mod ) return x -= modulus;
		if ( x < -half_mod ) return x += modulus; 
		else return x;
	}
 
	inline Element &subin (Element &x, const Element &y) const {
		x -= y;
		if ( x > half_mod ) return x -= modulus;
		if ( x < -half_mod ) return x += modulus; 
		else return x;
	}
 
	inline Element &mulin (Element &x, const Element &y) const {
		return mul(x,x,y);
	}
 
	inline Element &divin (Element &x, const Element &y) const {
		return div(x,x,y);
	}
 
	inline Element &negin (Element &x) const {
		//if (x == 0.) return x; 
		return x = - x; 
	}
		
	inline Element &invin (Element &x) const {
		return inv (x, x);
	}
		
	inline Element &axpyin (Element &r, const Element &a, const Element &x) const {
		r += a * x;
		return init( r, r);
	}
	static inline float getMaxModulus()
	{ return 2048; } // 2^11
		
};

#endif

