/*
 * Copyright (C) 2014 FFLAS-FFPACK group
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

#ifndef __FFLASFFPACK_fscal_mp_INL
#define __FFLASFFPACK_fscal_mp_INL

#include "fflas-ffpack/field/rns-integer.h"
#include "fflas_fscal.h"
#include "fflas_fgemm.inl"
namespace FFLAS {

    /*
     *  specialization for the field RNSInteger<rns_double>
     */

    // level 1 : fscalin
    template<>
    inline void fscalin(const FFPACK::RNSInteger<FFPACK::rns_double> &F,  const size_t n,
                        const FFPACK::rns_double::Element alpha,
                        FFPACK::rns_double::Element_ptr A, const size_t inc)
    {
        for (size_t i=0;i<F.size();i++)
            fscalin(F.rns()._field_rns[i], n, alpha._ptr[i*alpha._stride], A._ptr+i*A._stride,inc);
    }
    // level 1 : fscal
    template<>
    inline void fscal(const FFPACK::RNSInteger<FFPACK::rns_double> &F,  const size_t n,
                      const FFPACK::rns_double::Element alpha,
                      FFPACK::rns_double::ConstElement_ptr A, const size_t Ainc,
                      FFPACK::rns_double::Element_ptr B, const size_t Binc)
    {
        for (size_t i=0;i<F.size();i++)
            fscal(F.rns()._field_rns[i], n, alpha._ptr[i*alpha._stride], A._ptr+i*A._stride,Ainc, B._ptr+i*B._stride,Binc);
    }
    // level 2 : fscalin
    template<>
    inline void fscalin(const FFPACK::RNSInteger<FFPACK::rns_double> &F,  const size_t m, const size_t n,
                        const FFPACK::rns_double::Element alpha,
                        FFPACK::rns_double::Element_ptr A, const size_t lda) {
        for (size_t i=0;i<F.size();i++)
            fscalin(F.rns()._field_rns[i], m, n, alpha._ptr[i*alpha._stride], A._ptr+i*A._stride,lda);
    }
    // level 2 : fscal
    template<>
    inline void fscal(const FFPACK::RNSInteger<FFPACK::rns_double> &F, const size_t m, const size_t n,
                      const FFPACK::rns_double::Element alpha,
                      FFPACK::rns_double::ConstElement_ptr A, const size_t lda,
                      FFPACK::rns_double::Element_ptr B, const size_t ldb) {
        for (size_t i=0;i<F.size();i++)
            fscal(F.rns()._field_rns[i], m, n, alpha._ptr[i*alpha._stride], A._ptr+i*A._stride, lda, B._ptr+i*B._stride, ldb);

    }
}

#include "fflas-ffpack/fflas/fflas_freduce_mp.inl"

namespace FFLAS {
    /*
     *  specialization for the field RNSIntegerMod<rns_double>
     */

    // level 1 : fscalin
    template<>
    inline void fscalin(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F,  const size_t n,
                        const typename FFPACK::RNSIntegerMod<FFPACK::rns_double>::Element alpha,
                        typename FFPACK::RNSIntegerMod<FFPACK::rns_double>::Element_ptr A, const size_t inc)
    {
        fscalin(F.delayed(),n,alpha,A,inc);
        freduce (F, n, A, inc);
    }
    // level 1 : fscal
    template<>
    inline void fscal(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F,  const size_t n,
                      const FFPACK::rns_double::Element alpha,
                      FFPACK::rns_double::ConstElement_ptr A, const size_t Ainc,
                      FFPACK::rns_double::Element_ptr B, const size_t Binc)
    {
        fscal(F.delayed(),n,alpha,A,Ainc,B,Binc);
        freduce (F, n, B, Binc);
    }
    // level 2 : fscalin
    template<>
    inline void fscalin(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F,  const size_t m, const size_t n,
                        const FFPACK::rns_double::Element alpha,
                        FFPACK::rns_double::Element_ptr A, const size_t lda)
    {
        fscalin(F.delayed(),m,n,alpha,A,lda);
        freduce (F, m, n, A, lda);
    }
    // level 2 : fscal
    template<>
    inline void fscal(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F, const size_t m, const size_t n,
                      const FFPACK::rns_double::Element alpha,
                      FFPACK::rns_double::ConstElement_ptr A, const size_t lda,
                      FFPACK::rns_double::Element_ptr B, const size_t ldb)
    {
        fscal(F.delayed(),m,n,alpha,A,lda,B,ldb);
        freduce (F, m, n, B, ldb);
    }

} //end of namespace FFLAS


#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
