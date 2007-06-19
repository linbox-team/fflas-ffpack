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

#include <math.h>


#ifndef __MODULAR_DOUBLE_BALANCED_H
#define __MODULAR_DOUBLE_BALANCED_H

template <class T>
class Modular;

template <class Element>
class ModularRandIter
{
public:
	ModularRandIter (const Modular<Element> &F):_F(F){}
	ModularRandIter (const ModularRandIter<Element> &R) 
		: _F (R._F) {}
	Element &random (Element &a) const
	{ return _F.init(a,(double)rand()); }
private:
	Modular<Element> _F;
	
};

template<>
class Modular <double>{

protected:
	double modulus;
	double inv_modulus;
	double half_mod;

       
public:	       
	typedef double Element;
	typedef unsigned long FieldInt;
	typedef ModularRandIter<double> RandIter;

	const bool balanced;

	Modular<double> (int p)  : modulus((double)p), inv_modulus(1./(double)p), half_mod( (p-1.)/2), balanced(true) {}

	Modular<double>(const Modular<double>& mf) : modulus(mf.modulus), inv_modulus(mf.inv_modulus), half_mod(mf.half_mod), balanced(true){}
	
	const Modular<double> &operator=(const Modular<double> &F) {
			modulus = F.modulus;
			inv_modulus = F.inv_modulus;
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

	inline Element& init(Element& x, double y =0) const {		  
		double tmp=y;
		
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
		else return x;
	}
 
	inline Element &sub (Element &x, const Element &y, const Element &z) const {
		x = y - z;
		if (x > half_mod) return x -= modulus;
		else if (x < -half_mod) return x += modulus;
		else return x;
	}
		
	inline Element &mul (Element &x, const Element &y, const Element &z) const {		
		double tmp= y*z;
		//x= tmp - floor(tmp*inv_modulus)*modulus;
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
		//return r = tmp- floor(tmp*inv_modulus)*modulus; 
		return init( r, tmp);

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
		double tmp=r+a*x;
		//return r= tmp- floor(tmp*inv_modulus)*modulus; 
		return init( r, tmp );
	}
		
};
template<>
class Modular <float>{

protected:
	float modulus;
	float inv_modulus;
	float half_mod;

       
public:	       
	typedef float Element;
	typedef unsigned long FieldInt;
	typedef ModularRandIter<float> RandIter;

	const bool balanced;

	Modular<float> (int p)  : modulus((float)p), inv_modulus(1./(float)p), half_mod( (p-1.)/2), balanced(true) {}

	Modular<float>(const Modular<float>& mf) : modulus(mf.modulus), inv_modulus(mf.inv_modulus), half_mod(mf.half_mod), balanced(true){}
	
	const Modular<float> &operator=(const Modular<float> &F) {
			modulus = F.modulus;
			inv_modulus = F.inv_modulus;
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
		float tmp= y*z;
		//x= tmp - floor(tmp*inv_modulus)*modulus;
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
		else return x = (float) tx;
	}

	inline Element &axpy (Element &r, 
			      const Element &a, 
			      const Element &x, 
			      const Element &y) const {
		float tmp= a*x+y;
		//return r = tmp- floor(tmp*inv_modulus)*modulus; 
		return init( r, tmp);

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
		float tmp=r+a*x;
		//return r= tmp- floor(tmp*inv_modulus)*modulus; 
		return init( r, tmp );
	}
		
};

#endif

