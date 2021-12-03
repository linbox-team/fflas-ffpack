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

#ifndef __FFPACK_unparametric_rns_integer_H
#define __FFPACK_unparametric_rns_integer_H

#include <givaro/givinteger.h>

#include "fflas-ffpack/field/rns-double.h"

namespace FFPACK {


    template<typename RNS>
    class RNSInteger {
    protected:
        const RNS *_rns; // the rns structure
        typedef typename RNS::BasisElement BasisElement;
        typedef Givaro::Integer integer;

    public:
        typedef typename RNS::Element                   Element;
        typedef typename RNS::Element_ptr           Element_ptr;
        typedef typename RNS::ConstElement_ptr ConstElement_ptr;

        class RandIter : public rnsRandIter<RNS> {
        public:
            RandIter(const RNSInteger<RNS> &F, uint64_t seed=0) : rnsRandIter<RNS>(*F._rns,seed) {}
        };


        Element                one, mOne,zero;

        RNSInteger(const RNS& myrns) : _rns(&myrns)
        {
            init(one,1);
            init(zero,0);
            init(mOne,-1);
        }
        template<typename T>
        RNSInteger(const T &F) : _rns(&(F.rns()))
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
                x._ptr = FFLAS::fflas_new<BasisElement>(_rns->_size);
                x._stride=1;
                x._alloc=true;
            }
            return x;
        }
        Element& init(Element& x, const Givaro::Integer& y) const{
            init(x);
            size_t k =(y.bitsize())/16+((y.bitsize())%16?1:0);
            _rns->init(1,1,x._ptr,x._stride, &y,1,k);
            return x;
        }
        Element& reduce (Element& x, const Element& y) const {return assign (x,y);}

        Element& reduce (Element& x) const {return x;}

        Givaro::Integer convert(Givaro::Integer& x, const Element& y)const {
            _rns->convert(1,1,integer(0),&x,1,y._ptr,y._stride);
            return x;
        }

        Element& assign(Element& x, const Element& y) const {
            for(size_t i=0;i<_rns->_size;i++)
                x._ptr[i*x._stride] = y._ptr[i*y._stride];
            return x;
        }
        std::ostream& write(std::ostream& os, const Element& y) const {
            os<<"[ "<< (long) (y._ptr)[0];
            for(size_t i=1;i<_rns->_size;i++)
                os<<" , "<< (long) (y._ptr)[i*y._stride];
            return os<<" ]";
        }


        std::ostream& write(std::ostream& os) const {
            os<<"RNSInteger with M:=[ "<< (long) _rns->_basis[0];
            for(size_t i=1;i<_rns->_size;i++)
                os<<" , "<< (long) _rns->_basis[i];
            return os<<" ]"<<std::endl;
        }



    }; // end of class Unparametric<rns_double>


}  // end of namespace FFPACK

namespace FFLAS {

    // specialization for the fflas alloc function
    template<>
    inline FFPACK::rns_double_elt_ptr
    fflas_new(const FFPACK::RNSInteger<FFPACK::rns_double> &F, const size_t m, const Alignment align){
        double *ptr=FFLAS::fflas_new<double>(m*F.size(), align);
        return FFPACK::rns_double_elt_ptr(ptr,m);
    }

    template<>
    inline FFPACK::rns_double_elt_ptr
    fflas_new(const FFPACK::RNSInteger<FFPACK::rns_double> &F, const size_t m, const size_t n,  const Alignment align){
        return fflas_new(F, m*n, align);
    }

    // function to convert from integer to RNS (note: this is not the finit function from FFLAS, extra k)
    template<typename RNS>
    void finit_rns(const FFPACK::RNSInteger<RNS> &F, const size_t m, const size_t n, size_t k,
                   const Givaro::Integer *B, const size_t ldb, typename FFPACK::RNSInteger<RNS>::Element_ptr A)
    {
        F.rns().init(m,n,A._ptr,A._stride, B,ldb,k);
    }
    // function to convert from RNS to integer (note: this is not the fconvert function from FFLAS, extra alpha)
    template<typename RNS>
    void fconvert_rns(const FFPACK::RNSInteger<RNS> &F, const size_t m, const size_t n,
                      Givaro::Integer alpha, Givaro::Integer *B, const size_t ldb, typename FFPACK::RNSInteger<RNS>::ConstElement_ptr A)
    {
        F.rns().convert(m,n,alpha,B,ldb,A._ptr,A._stride);
    }


} // end of namespace FFLAS

#endif // __FFPACK_unparametric_rns_integer_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
