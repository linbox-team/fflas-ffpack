/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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
/*! @file field/rns-integer.h
 * @ingroup field
 * @brief  representation of <code>Z</code> using RNS representation (note: fixed precision)
 */

#ifndef __FFPACK_unparametric_rns_double_H 
#define __FFPACK_unparametric_rns_double_H

#include "fflas-ffpack/field/integer.h"
#include "fflas-ffpack/field/rns-double.h"

// activate only if FFLAS-FFPACK haves multiprecision integer
#ifdef __FFLASFFPACK_HAVE_INTEGER


namespace FFPACK {
	
	
	template<typename RNS>
	class RNSInteger {
	protected:
		const RNS *_rns; // the rns structure
		typedef typename RNS::BasisElement BasisElement;
		
	public:
		typedef typename RNS::Element                   Element;
		typedef typename RNS::Element_ptr           Element_ptr;
		typedef typename RNS::ConstElement_ptr ConstElement_ptr;

		Element                one, mOne,zero;
		
		RNSInteger(const RNS& rns) : _rns(&rns) 
		{
			init(one,1);
			init(zero,0);
			init(mOne,-1);			
		}
		
		
		const RNS& rns() const {return *_rns;}

		size_t size() const {return _rns->_size;}
		
		bool isOne(const Element& x) const {
			bool isone=true;
			for (size_t i=0;i<_rns->_size;i++)
				isone&= (one._ptr[i]== x._ptr[i]);
			return isone;
		}

		bool isMOne(const Element& x) const {
			bool ismone=true;
			for (size_t i=0;i<_rns->_size;i++)
				ismone&= (mOne._ptr[i]== x._ptr[i]);
			return ismone;
		}

		bool isZero(const Element& x) const {
			bool iszero=true;
			for (size_t i=0;i<_rns->_size;i++)
				iszero&= (zero._ptr[i]== x._ptr[i]);
			return iszero;
		}

		integer characteristic(integer &p) const { return p=0;}
		
		integer cardinality(integer &p) const { return p=-1;}

		Element& init(Element& x) const{
			if (x._ptr == NULL){
				x._ptr = new BasisElement[_rns->_size];
				x._stride=1;
				x._alloc=true;
			}
			return x;
		}
		Element& init(Element& x, const FFPACK::Integer& y) const{
			init(x);
			size_t k =(y.bitsize())/16+((y.bitsize())%16?1:0);
			_rns->init(1,1,x._ptr,x._stride, &y,1,k);
			return x;
		}
		FFPACK::Integer convert(FFPACK::Integer& x, const Element& y)const {
			_rns->convert(1,1,integer(0),&x,1,y._ptr,y._stride);
			return x;
		}
		
		Element& assign(Element& x, const Element& y) const {
			for(size_t i=0;i<_rns->_size;i++)
				x._ptr[i*x._stride] = y._ptr[i*y._stride];
			return x;
		}
		ostream& write(ostream& os, const Element& y) const {
			os<<"[ "<< (long) (y._ptr)[0];
			for(size_t i=1;i<_rns->_size;i++)
				os<<" , "<< (long) (y._ptr)[i*y._stride];
			return os<<" ]";
		}


		ostream& write(ostream& os) const {
			os<<"M:=[ "<< (long) _rns->_basis[0];
			for(size_t i=1;i<_rns->_size;i++)
				os<<" , "<< (long) _rns->_basis[i];
			return os<<" ]"<<endl;
		}



	}; // end of class Unparametric<rns_double>
	

}  // end of namespace FFPACK

namespace FFLAS {
 
	// specialization for the fflas alloc function
	template<>
	FFPACK::rns_double_elt_ptr fflas_new(FFPACK::RNSInteger<FFPACK::rns_double> &F, const size_t m, const size_t n){
		double *ptr=new double[m*n*F.size()];
		return FFPACK::rns_double_elt_ptr(ptr,m*n);
	}
		
	// function to convert from integer to RNS (note: this is not the finit function from FFLAS, extra k)
	template<typename RNS>
	void finit(const FFPACK::RNSInteger<RNS> &F, const size_t m, const size_t n, size_t k,
		   const FFPACK::integer *B, const size_t ldb, typename FFPACK::RNSInteger<RNS>::Element_ptr A)
	{
		F.rns().init(m,n,A._ptr,A._stride, B,ldb,k);
	}
	// function to convert from RNS to integer (note: this is not the fconvert function from FFLAS, extra alpha)
	template<typename RNS>
	void fconvert(const FFPACK::RNSInteger<RNS> &F, const size_t m, const size_t n,
		      FFPACK::integer alpha, FFPACK::integer *B, const size_t ldb, typename FFPACK::RNSInteger<RNS>::ConstElement_ptr A)
	{		
		F.rns().convert(m,n,alpha,B,ldb,A._ptr,A._stride);
	}


}; // end of namespace FFLAS
 
#endif
#endif 
 
