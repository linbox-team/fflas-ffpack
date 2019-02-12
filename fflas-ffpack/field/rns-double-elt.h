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

/*! @file field/rns-double-elt.h
 * @ingroup field
 * @brief  rns elt structure with double support
 */


#ifndef __FFLASFFPACK_field_rns_double_elt_INL
#define __FFLASFFPACK_field_rns_double_elt_INL

#include "fflas-ffpack/utils/fflas_memory.h"
#include "fflas-ffpack/utils/cast.h"

namespace FFPACK {

    // forward declaration
    struct rns_double_elt_ptr;
    struct rns_double_elt_cstptr;

    // element of the rns structure (allow virtualization of element from an array of double)
    struct rns_double_elt {
        double    *_ptr;
        size_t  _stride;
        bool     _alloc; // specify wether Element owns its memory; alloc is true only through F.init() and _ptr==NULL (this is to handle Element allocated within a matrix)
        rns_double_elt(): _ptr(NULL), _stride(0), _alloc(false) {}
        ~rns_double_elt(){ if (_alloc) FFLAS::fflas_delete(_ptr);}
        rns_double_elt(double* p, size_t r, size_t a=false) : _ptr(p), _stride(r), _alloc(a) {}
        inline  rns_double_elt_ptr    operator&() ;
        inline  rns_double_elt_cstptr operator&()const ;
        rns_double_elt(const rns_double_elt& x) : _ptr(x._ptr),_stride(x._stride),_alloc(false) {}
    };

    // pointer to element of the rns structure (allow virtualization of element from an array of double)
    struct rns_double_elt_ptr : public rns_double_elt {
        rns_double_elt other;
        rns_double_elt_ptr(){}
        rns_double_elt_ptr(double* p, size_t r)            : rns_double_elt(p,r,false){}
        rns_double_elt_ptr(const rns_double_elt_ptr &x)    : rns_double_elt(x._ptr,x._stride,false){}
        rns_double_elt_ptr(const rns_double_elt_cstptr &x);
        rns_double_elt_ptr(rns_double_elt_ptr &&)=default;
        //inline  operator rns_double_elt_cstptr();
        inline rns_double_elt_ptr* operator&(){return this;}
        inline rns_double_elt&     operator*()  {return static_cast<rns_double_elt&>(*this);}
        inline rns_double_elt     operator[](size_t i) const {return rns_double_elt(_ptr+i,_stride);} // BUGGY
        inline rns_double_elt&     operator[](size_t i) {other=rns_double_elt(_ptr+i,_stride);return other;} // BUGGY
        inline rns_double_elt_ptr  operator++() {return rns_double_elt_ptr(_ptr++,_stride);}
        inline rns_double_elt_ptr  operator--() {return rns_double_elt_ptr(_ptr--,_stride);}
        inline rns_double_elt_ptr  operator+(size_t inc) {return rns_double_elt_ptr(_ptr+inc,_stride);}
        inline rns_double_elt_ptr  operator-(size_t inc) {return rns_double_elt_ptr(_ptr-inc,_stride);}
        inline rns_double_elt_ptr& operator+=(size_t inc) {_ptr+=inc;return *this;}
        inline rns_double_elt_ptr& operator-=(size_t inc) {_ptr-=inc;return *this;}
        inline rns_double_elt_ptr& operator=(const rns_double_elt_ptr& x);
        bool operator< (const rns_double_elt_ptr& x) {return _ptr < x._ptr;}
        bool operator!= (const rns_double_elt_ptr& x) {return _ptr != x._ptr;}
    };
    struct rns_double_elt_cstptr : public rns_double_elt {
        rns_double_elt other;
        rns_double_elt_cstptr(){}
        rns_double_elt_cstptr(double* p, size_t r)            : rns_double_elt(p,r,false){}
        rns_double_elt_cstptr(const rns_double_elt_ptr& x)    : rns_double_elt(x._ptr,x._stride,false){}
        rns_double_elt_cstptr(const rns_double_elt_cstptr& x) : rns_double_elt(x._ptr,x._stride,false){}
        rns_double_elt_cstptr(rns_double_elt_cstptr &&)=default;
        inline rns_double_elt_cstptr* operator&(){return this;}
        inline rns_double_elt&     operator*() const  {
            return *const_cast<rns_double_elt*>(static_cast<const rns_double_elt*>(this));
        }
        inline rns_double_elt      operator[](size_t i)const {return rns_double_elt(_ptr+i,_stride);}
        inline rns_double_elt&     operator[](size_t i) {other=rns_double_elt(_ptr+i,_stride);return other;} // BUGGY

        //inline rns_double_elt&     operator[](size_t i)const {return *((*this)+i);}// BUGGY
        inline rns_double_elt_cstptr  operator++() {return rns_double_elt_cstptr(_ptr++,_stride);}
        inline rns_double_elt_cstptr  operator--() {return rns_double_elt_cstptr(_ptr--,_stride);}
        inline rns_double_elt_cstptr  operator+(size_t inc)const {return rns_double_elt_cstptr(_ptr+inc,_stride);}
        inline rns_double_elt_cstptr  operator-(size_t inc)const {return rns_double_elt_cstptr(_ptr-inc,_stride);}
        inline rns_double_elt_cstptr& operator+=(size_t inc) {_ptr+=inc;return *this;}
        inline rns_double_elt_cstptr& operator-=(size_t inc) {_ptr-=inc;return *this;}
        inline rns_double_elt_cstptr& operator=(const rns_double_elt_cstptr& x);
        bool operator< (const rns_double_elt_cstptr& x) {return _ptr < x._ptr;}
        bool operator!= (const rns_double_elt_cstptr& x) {return _ptr != x._ptr;}
    };

    inline rns_double_elt_ptr& rns_double_elt_ptr::operator=(const rns_double_elt_ptr& x)  {
        if (this != &x){
            if (_alloc) FFLAS::fflas_delete(_ptr);
            _ptr= x._ptr;
            _stride=x._stride;
            _alloc=false;
        }
        return *this;
    }
    inline rns_double_elt_cstptr& rns_double_elt_cstptr::operator=(const rns_double_elt_cstptr& x)  {
        if (this != &x){
            if (_alloc) FFLAS::fflas_delete(_ptr);
            _ptr= x._ptr;
            _stride=x._stride;
            _alloc=false;
        }
        return *this;
    }

    inline rns_double_elt_ptr::rns_double_elt_ptr(const rns_double_elt_cstptr &x)
    : rns_double_elt(x._ptr,x._stride,false){}
    //inline rns_double_elt_ptr::operator rns_double_elt_cstptr(){return rns_double_elt_cstptr(_ptr,_stride);}
    inline rns_double_elt_ptr    rns_double_elt::operator&()       {return 	rns_double_elt_ptr(_ptr,_stride);}
    inline rns_double_elt_cstptr rns_double_elt::operator&() const {return 	rns_double_elt_cstptr(_ptr,_stride);}


    template<>
    inline rns_double_elt_ptr fflas_const_cast (rns_double_elt_cstptr x){return x;}
    template<>
    inline rns_double_elt_cstptr fflas_const_cast (rns_double_elt_ptr x){return x;}


} // end namespace FFPACK:

#endif // __FFLASFFPACK_field_rns_double_elt_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
