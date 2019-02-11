/* fflas/fflas_freduce_mp.inl
 * Copyright (C) 2014 FFLAS FFPACK group
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

#ifndef __FFLASFFPACK_fflas_freduce_mp_INL
#define __FFLASFFPACK_fflas_freduce_mp_INL

#include "fflas-ffpack/field/rns-integer-mod.h"

namespace FFLAS {

    // specialization of the level1 freduce function for the field RNSInteger<rns_double>
    template<>
    inline void freduce (const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F,
                         const size_t n, FFPACK::RNSIntegerMod<FFPACK::rns_double>::Element_ptr A, size_t inc)
    {
        if (n==0) return;
        //cout<<"freduce: "<<n<<" with "<<inc<<endl;
        if (inc==1)
            F.reduce_modp(n,A);
        else
            F.reduce_modp(n,1,A,inc);
    }
    // specialization of the level2 freduce function for the field RNSInteger<rns_double>
    template<>
    inline void freduce (const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F,
                         const size_t m, const size_t n, FFPACK::rns_double::Element_ptr A, size_t lda)
    {
        if (n==0||m==0) return;
        //cout<<"freduce: "<<m<<" x "<<n<<" "<<lda<<endl;
        if (lda == n)
            F.reduce_modp(m*n,A);
        else
            F.reduce_modp(m,n,A,lda);
    }



} // end of namespace FFLAS

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
