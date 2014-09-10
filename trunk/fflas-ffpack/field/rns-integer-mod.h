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
/*! @file field/rns-integer-mod.h
 * @ingroup field
 * @brief  representation of <code>Z/pZ</code> using RNS representation (note: fixed precision)
 */


#ifndef __FFPACK_rns_integer_mod_H
#define __FFPACK_rns_integer_mod_H

#include <vector>
using namespace std;

#include "fflas-ffpack/field/integer.h"
#include "fflas-ffpack/field/modular-integer.h"
#include "fflas-ffpack/field/rns-double.h"
#include "fflas-ffpack/fflas/fflas.h"
// activate only if FFLAS-FFPACK haves multiprecision integer
#ifdef __FFLASFFPACK_HAVE_INTEGER


namespace FFPACK {

	template<typename RNS>
	class RNSIntegerMod {
	public:
		typedef typename RNS::Element                   Element;
		typedef typename RNS::Element_ptr           Element_ptr;
		typedef typename RNS::ConstElement_ptr ConstElement_ptr;

	protected:
		typedef typename RNS::BasisElement BasisElement;
		typedef FFPACK::Modular<BasisElement> ModField;
		integer                              _p;
		vector<BasisElement>       _Mi_modp_rns;
		vector<BasisElement>       _iM_modp_rns;
		const RNS                         *_rns;
		FFPACK::Modular<FFPACK::Integer>     _F;
	public:
		Element                one, mOne,zero;

		RNSIntegerMod(const integer& p, const RNS& myrns) : _p(p),
								  _Mi_modp_rns(myrns._size*myrns._size),
								  _iM_modp_rns(myrns._size*myrns._size),
								  _rns(&myrns),
								  _F(p){
			init(one,1);
			init(zero,0);
			init(mOne,-1);
			integer iM=0;
			size_t mysize=myrns._size;
			integer sum=0;
			for (size_t i=0;i<mysize;i++){
				integer Mi = myrns._Mi[i] % _p;
				for (size_t j=0;j<mysize;j++){
					_iM_modp_rns[i+j*mysize]=  iM % myrns._basis[j];
					_Mi_modp_rns[i+j*mysize]=  Mi % myrns._basis[j];
				}
				iM+=myrns._M;iM%=_p;
				sum+=myrns._basis[i];
			}
		}

		const rns_double& rns() const {return *_rns;}

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
			return isZero;
		}

		integer characteristic(integer &p) const { return p=_p;}

		integer cardinality(integer &p) const { return p=_p;}

		Element& init(Element& x) const{
			if (x._ptr == NULL){
				x._ptr = FFLAS::fflas_new<BasisElement>(_rns->_size);
				x._stride=1;
				x._alloc=true;
			}
			return x;
		}
		Element& init(Element& x, const FFPACK::Integer& y) const{
			init(x);
			size_t k =(_p.bitsize())/16+((_p.bitsize())%16?1:0);
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

		Element& inv(Element& x, const Element& y) const {
			FFPACK::Integer tmp;
			convert(tmp,y);
			_F.invin(tmp);
			init(x,tmp);
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


		void reduce_modp(size_t n, BasisElement* A, size_t rda) const{

			size_t _size= _rns->_size;
			BasisElement *Gamma, *alpha;
			Gamma = FFLAS::fflas_new<BasisElement>(n*_size);
			alpha = FFLAS::fflas_new<BasisElement>(n);

			// compute Gamma
			for(size_t i=0;i<_size;i++)
				FFLAS::fscal(_rns->_field_rns[i], n, _rns->_MMi[i], A+i*rda, 1, Gamma+i*n,1);

			UnparametricField<BasisElement> D;

			// compute A = _Mi_modp_rns.Gamma (note must be reduced mod m_i, but this is postpone to the end)
			FFLAS::fgemm(D,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans, _size, n, _size, D.one, _Mi_modp_rns.data(), _size, Gamma, n, D.zero, A, rda);

			// compute alpha = _invbase.Gamma
			FFLAS::fgemv(D,FFLAS::FflasTrans, _size, n, D.one, Gamma, n, _rns->_invbasis.data(), 1 , D.zero, alpha, 1);

			// compute ((z-(alpha.M mod p)) mod m_i (perform the subtraction over Z and reduce at the end)
			for(size_t i=0;i<_size;i++){
				for(size_t j=0;j<n;j++){
					//long aa=floor(alpha[j]+0.5);
					long aa= (long)rint(alpha[j]);
					A[j+i*rda]-=_iM_modp_rns[aa+i*_size];
				}
			}

			// reduce each row of A modulo m_i
			for (size_t i=0;i<_size;i++)
				FFLAS::finit(_rns->_field_rns[i], n, A+i*rda, 1);

			delete[] Gamma;
			delete[] alpha;

	 	}

		void reduce_modp(size_t m, size_t n, BasisElement* A, size_t lda, size_t rda) const{
			//cout<<"REDUCE MOD WITH LDA!=N"<<endl;
			size_t _size= _rns->_size;
			size_t mn=m*n;
			BasisElement *Gamma, *alpha, *z;
			Gamma = FFLAS::fflas_new<BasisElement>(mn*_size);
			alpha = FFLAS::fflas_new<BasisElement>(mn);
			z     = FFLAS::fflas_new<BasisElement>(mn*_size);

			// compute Gamma
			for(size_t i=0;i<_size;i++)
				FFLAS::fscal(_rns->_field_rns[i], m, n, _rns->_MMi[i], A+i*rda, lda, Gamma+i*mn,n);

			// compute Gamma = _Mi_modp_rns.Gamma (note must be reduced mod m_i, but this is postpone to the end)
			UnparametricField<BasisElement> D;
			FFLAS::fgemm(D,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,_size, mn, _size, D.one, _Mi_modp_rns.data(), _size, Gamma, mn, D.zero, z, mn);

			// compute alpha = _invbase.Gamma
			FFLAS::fgemv(D, FFLAS::FflasTrans, _size, mn, D.one, Gamma, mn, _rns->_invbasis.data(), 1 , D.zero, alpha, 1);

			// compute A=((Gamma--(alpha.M mod p)) mod m_i (perform the subtraction over Z and reduce at the end)
 			for(size_t k=0;k<_size;k++){
				for(size_t i=0;i<m;i++)
					for(size_t j=0;j<n;j++){
						long aa=(long)floor(alpha[j+i*n]+0.5);
						A[j+i*lda+k*rda]= z[j+i*n+k*mn]-_iM_modp_rns[aa+k*_size];
					}
			}

			// reduce each row of A modulo m_i
			for (size_t i=0;i<_size;i++)
				FFLAS::finit(_rns->_field_rns[i], m, n, A+i*rda, lda);

			delete[] Gamma;
			delete[] alpha;
			delete[] z;
		}

	}; // end of class RNSIntegerMod

} // end of namespace FFPACK


namespace FFLAS {

	// specialization for the fflas alloc function
	template<>
	inline FFPACK::rns_double_elt_ptr
	fflas_new(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F,  const size_t m,  const size_t n, const Alignment align){
		return FFPACK::rns_double_elt_ptr(FFLAS::fflas_new<double>(m*n*F.size(),align),m*n);
	}

	// function to convert from integer to RNS (note: this is not the finit function from FFLAS, extra k)
	template<typename RNS>
	void finit(const FFPACK::RNSIntegerMod<RNS> &F, const size_t m, const size_t n, size_t k,
		   const FFPACK::integer *B, const size_t ldb, typename RNS::Element_ptr A)
	{
		F.rns().init(m,n,A._ptr,A._stride, B,ldb,k);
	}
	// function to convert from RNS to integer (note: this is not the fconvert function from FFLAS, extra alpha)
	template<typename RNS>
	void fconvert(const FFPACK::RNSIntegerMod<RNS> &F, const size_t m, const size_t n,
		      FFPACK::integer alpha, FFPACK::integer *B, const size_t ldb, typename RNS::ConstElement_ptr A)
	{
		F.rns().convert(m,n,alpha,B,ldb,A._ptr,A._stride);
	}

	// specialization of the level1 finit function for the field RNSInteger<rns_double>
	template<>
	void finit(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F, const size_t n, FFPACK::rns_double::Element_ptr A, size_t inc)
	{
		if (inc==1)
			F.reduce_modp(n,A._ptr,A._stride);
		else
			F.reduce_modp(n,1,A._ptr,inc,A._stride);
       		//throw FFPACK::Failure(__func__,__FILE__,__LINE__,"finit RNSIntegerMod  -> (inc!=1) NOT SUPPORTED");
	}
	// specialization of the level2 finit function for the field RNSInteger<rns_double>
	template<>
	void finit(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F, const size_t m, const size_t n, FFPACK::rns_double::Element_ptr A, size_t lda)
	{
		if (lda == n)
			F.reduce_modp(m*n,A._ptr,A._stride);
		else
			F.reduce_modp(m,n,A._ptr,lda,A._stride); // seems to be buggy
	}
	// specialization of the level1 fscalin function for the field RNSInteger<rns_double>
	template<>
	void fscalin(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F,  const size_t n,
		     const FFPACK::rns_double::Element alpha,
		     FFPACK::rns_double::Element_ptr A, const size_t inc) {

		for (size_t i=0;i<F.size();i++)
			fscalin(F.rns()._field_rns[i], n, alpha._ptr[i], A._ptr+i*A._stride,inc);
		finit(F,n,A,inc);
	}
	// specialization of the level2 fscalin function for the field RNSInteger<rns_double>
	template<>
	void fscalin(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F,  const size_t m, const size_t n,
		     const FFPACK::rns_double::Element alpha,
		     FFPACK::rns_double::Element_ptr A, const size_t lda) {
		for (size_t i=0;i<F.size();i++)
			fscalin(F.rns()._field_rns[i], m, n, alpha._ptr[i], A._ptr+i*A._stride,lda);
		finit(F,m,n,A,lda);
	}

} // end of namespace FFLAS

#endif
#endif

