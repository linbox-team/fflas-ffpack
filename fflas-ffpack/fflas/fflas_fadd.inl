/*
 * Copyright (C) 2014 FFLAS-FFPACK group
 *
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

#ifndef __FFLASFFPACK_fadd_INL
#define __FFLASFFPACK_fadd_INL

#include "fflas-ffpack/fflas/fflas_simd.h"

namespace FFLAS { namespace vectorised {


    template<class SimdT, class Element, bool positive>
    inline typename std::enable_if<is_simd<SimdT>::value, void>::type
    VEC_ADD(SimdT & C, SimdT & A, SimdT & B, SimdT & Q, SimdT & T, SimdT & P, SimdT & NEGP, SimdT & MIN, SimdT & MAX)
    {
        using simd = Simd<Element>;
        C = simd::add(A, B);
        Q = simd::vand(simd::greater(C, MAX),NEGP);
        if (!positive) {
            T = simd::vand(simd::lesser(C, MIN),P);
            Q = simd::vor(Q, T);
        }
        C = simd::add(C, Q);
    }

    template<bool positive, class Element, class T1, class T2>
    inline typename std::enable_if<FFLAS::support_simd_add<Element>::value, void>::type
    addp(Element * T, const Element * TA, const Element * TB,  size_t n,  Element p,  T1 min_,  T2 max_)
    {
        Element min= (Element)min_, max= (Element)max_;
        using simd = Simd<Element>;
        using vect_t = typename simd::vect_t;

        size_t i = 0;

        if (n < simd::vect_size)
        {
            for (; i < n ; i++)
            {
                T[i] = TA[i] + TB[i];
                T[i] -= (T[i] > max) ? p : 0;
                if (!positive)
                {
                    T[i] += (T[i] < min) ? p : 0;
                }
            }
            return;

        }

        vect_t A,B,C,Q,P,NEGP,TMP,MIN,MAX;
        P   = simd::set1(p);
        NEGP= simd::set1(-p);
        MIN = simd::set1(min);
        MAX = simd::set1(max);
        long st = long(T)%simd::alignment;
        if (st)
        { // the array T is not 32 byte aligned (process few elements s.t. (T+i) is 32 bytes aligned)
            for (size_t j=static_cast<size_t>(st) ; j < simd::alignment ; j += sizeof(Element), i++)
            {
                T[i] = TA[i] + TB[i];
                T[i] -= (T[i] > max) ? p : 0;
                if (!positive)
                    T[i] += (T[i] < min) ? p : 0;
            }
        }
        FFLASFFPACK_check((long(T+i) % simd::alignment == 0));
        if ( (long(TA+i)%simd::alignment==0) && (long(TB+i)%simd::alignment==0))
        {
            // perform the loop using 256 bits SIMD
            for (; i <= n - simd::vect_size ; i += simd::vect_size)
            {
                // C = simd::load(T+i);
                A = simd::load(TA+i);
                B = simd::load(TB+i);
                VEC_ADD<vect_t,Element,positive>(C, A, B, Q, TMP, P, NEGP, MIN, MAX);
                simd::store(T+i, C);
            }
        }
        // perform the last elt from T without SIMD
        for (; i < n ; i++)
        {
            T[i] = TA[i] + TB[i];
            T[i] -= (T[i] > max) ? p : 0;
            if (!positive)
                T[i] += (T[i] < min) ? p : 0;
        }
    }

    template<class SimdT, class Element,bool positive>
    inline typename std::enable_if<is_simd<SimdT>::value, void>::type
    VEC_SUB(SimdT & C, SimdT & A, SimdT & B, SimdT & Q, SimdT & T, SimdT & P, SimdT & NEGP, SimdT & MIN, SimdT & MAX)
    {
        using simd = Simd<Element>;
        C = simd::sub(A, B);
        T = simd::vand(simd::lesser(C, MIN),P);
        if (!positive) {
            Q = simd::vand(simd::greater(C, MAX),NEGP);
            T = simd::vor(Q, T);
        }
        C = simd::add(C, T);
    }

    template<bool positive, class Element, class T1, class T2>
    inline typename std::enable_if<FFLAS::support_simd_add<Element>::value, void>::type
    subp(Element * T, const Element * TA, const Element * TB, const size_t n, const Element p, const T1 min_, const T2 max_)
    {
        Element min = (Element)min_, max = (Element)max_;
        using simd = Simd<Element>;
        using vect_t = typename simd::vect_t;

        size_t i = 0;

        if (n < simd::vect_size)
        {
            for (; i < n ; i++)
            {
                T[i] = TA[i] - TB[i];
                if (!positive)
                    T[i] -= (T[i] > max) ? p : 0;
                T[i] += (T[i] < min) ? p : 0;
            }
            return;

        }
        vect_t A,B,C,Q,P,NEGP,TMP,MIN,MAX;
        P   = simd::set1(p);
        NEGP= simd::set1(-p);
        MIN = simd::set1(min);
        MAX = simd::set1(max);
        long st = long(T) % simd::alignment;
        if (st)
        { // the array T is not 32 byte aligned (process few elements s.t. (T+i) is 32 bytes aligned)
            for (size_t j = static_cast<size_t>(st) ; j < simd::alignment ; j += sizeof(Element), i++)
            {
                T[i] = TA[i] - TB[i];
                if (!positive)
                    T[i] -= (T[i] > max) ? p : 0;
                T[i] += (T[i] < min) ? p : 0;
            }
        }
        FFLASFFPACK_check((long(T+i) % simd::alignment == 0));
        if ( (long(TA+i) % simd::alignment == 0) && (long(TB+i) % simd::alignment == 0))
        {
            // perform the loop using 256 bits SIMD
            for (; i <= n - simd::vect_size ; i += simd::vect_size)
            {
                // C = simd::load(T+i);
                A = simd::load(TA+i);
                B = simd::load(TB+i);
                VEC_SUB<vect_t,Element,positive>(C, A, B, Q, TMP, P, NEGP, MIN, MAX);
                simd::store(T+i, C);
            }
        }

        // perform the last elt from T without SIMD
        for (; i < n ; i++)
        {
            T[i] = TA[i] - TB[i];
            if (!positive)
                T[i] -= (T[i] > max) ? p : 0;
            T[i] += (T[i] < min) ? p : 0;
        }
    }

    template<class Element>
    inline typename std::enable_if<FFLAS::support_simd_add<Element>::value, void>::type
    add(Element * T, const Element * TA, const Element * TB,  size_t n)
    {
        using simd = Simd<Element>;
        using vect_t = typename simd::vect_t;

        size_t i = 0;

        if (n < simd::vect_size)
        {
            for (; i < n ; i++)
            {
                T[i] = TA[i] + TB[i];
            }
            return;
        }

        vect_t A,B,C;
        long st = long(T)%simd::alignment;
        if (st)
        { // the array T is not 32 byte aligned (process few elements s.t. (T+i) is 32 bytes aligned)
            for (size_t j=static_cast<size_t>(st) ; j < simd::alignment ; j += sizeof(Element), i++)
            {
                T[i] = TA[i] + TB[i];
            }
        }
        FFLASFFPACK_check((long(T+i) % simd::alignment == 0));
        if ( (long(TA+i)%simd::alignment==0) && (long(TB+i)%simd::alignment==0))
        {
            // perform the loop using 256 bits SIMD
            for (; i <= n - simd::vect_size ; i += simd::vect_size)
            {
                // C = simd::load(T+i);
                A = simd::load(TA+i);
                B = simd::load(TB+i);
                C = simd::add(A, B);
                simd::store(T+i, C);
            }
        }
        // perform the last elt from T without SIMD
        for (; i < n ; i++)
        {
            T[i] = TA[i] + TB[i];
        }
    }
    template<class Element>
    inline typename std::enable_if<FFLAS::support_simd_add<Element>::value, void>::type
    sub(Element * T, const Element * TA, const Element * TB,  size_t n)
    {
        using simd = Simd<Element>;
        using vect_t = typename simd::vect_t;

        size_t i = 0;

        // CP: is this necessary?
        if (n < simd::vect_size)
        {
            for (; i < n ; i++)
            {
                T[i] = TA[i] - TB[i];
            }
            return;
        }

        vect_t A,B,C;
        long st = long(T)%simd::alignment;
        if (st)
        { // the array T is not 32 byte aligned (process few elements s.t. (T+i) is 32 bytes aligned)
            for (size_t j=static_cast<size_t>(st) ; j < simd::alignment ; j += sizeof(Element), i++)
            {
                T[i] = TA[i] - TB[i];
            }
        }
        FFLASFFPACK_check((long(T+i) % simd::alignment == 0));
        if ( (long(TA+i)%simd::alignment==0) && (long(TB+i)%simd::alignment==0))
        {
            // perform the loop using 256 bits SIMD
            for (; i <= n - simd::vect_size ; i += simd::vect_size)
            {
                // C = simd::load(T+i);
                A = simd::load(TA+i);
                B = simd::load(TB+i);
                C = simd::sub(A, B);
                simd::store(T+i, C);
            }
        }
        // perform the last elt from T without SIMD
        for (; i < n ; i++)
        {
            T[i] = TA[i] - TB[i];
        }
    }

} // vectorised
} //  FFLAS

namespace FFLAS { namespace details {

    /**** Specialised ****/

    template <class Field, bool ADD>
    typename std::enable_if<FFLAS::support_simd_add<typename Field::Element>::value, void>::type
    fadd (const Field & F,  const size_t N,
          typename Field::ConstElement_ptr A, const size_t inca,
          typename Field::ConstElement_ptr B, const size_t incb,
          typename Field::Element_ptr C, const size_t incc
          , FieldCategories::ModularTag
         )
    {
        if (inca == 1 && incb == 1 && incc == 1) {
            typename Field::Element p = (typename Field::Element) F.characteristic();
            if (ADD)
                FFLAS::vectorised::addp<!FieldTraits<Field>::balanced>(C,A,B,N,p,F.minElement(),F.maxElement());
            else
                FFLAS::vectorised::subp<!FieldTraits<Field>::balanced>(C,A,B,N,p,F.minElement(),F.maxElement());
        }
        else {
            typename Field::ConstElement_ptr Ai=A, Bi=B;
            typename Field::Element_ptr Ci=C;
            for (size_t i=0; i<N; i++, Ai+=inca, Bi+=incb, Ci+=incc)
                if (ADD)
                    F.add (*Ci, *Ai, *Bi);
                else
                    F.sub (*Ci, *Ai, *Bi);
        }
    }

    template <class Field, bool ADD>
    typename std::enable_if<!FFLAS::support_simd_add<typename Field::Element>::value, void>::type
    fadd (const Field & F,  const size_t N,
          typename Field::ConstElement_ptr A, const size_t inca,
          typename Field::ConstElement_ptr B, const size_t incb,
          typename Field::Element_ptr C, const size_t incc
          , FieldCategories::ModularTag
         )
    {
        if (inca == 1 && incb == 1 && incc == 1) {
            for (size_t i=0; i<N; i++)
                if (ADD)
                    F.add (C[i], A[i], B[i]);
                else
                    F.sub (C[i], A[i], B[i]);
        }
        else {
            typename Field::ConstElement_ptr Ai=A, Bi=B;
            typename Field::Element_ptr Ci=C;
            for (size_t i=0; i<N; i++, Ai+=inca, Bi+=incb, Ci+=incc)
                if (ADD)
                    F.add (*Ci, *Ai, *Bi);
                else
                    F.sub (*Ci, *Ai, *Bi);
        }
    }



    template <class Field, bool ADD>
    void
    fadd (const Field & F,  const size_t N,
          typename Field::ConstElement_ptr A, const size_t inca,
          typename Field::ConstElement_ptr B, const size_t incb,
          typename Field::Element_ptr C, const size_t incc
          , FieldCategories::GenericTag
         )
    {
        if (inca == 1 && incb == 1 && incc == 1) {
            for (size_t i=0; i<N; i++) {
                if (ADD)
                    F.add (C[i], A[i], B[i]);
                else
                    F.sub (C[i], A[i], B[i]);
            }
        }
        else {
            typename Field::ConstElement_ptr Ai=A, Bi=B;
            typename Field::Element_ptr Ci=C;
            for (size_t i=0; i<N; i++, Ai+=inca, Bi+=incb, Ci+=incc)
                if (ADD)
                    F.add (*Ci, *Ai, *Bi);
                else
                    F.add (*Ci, *Ai, *Bi);
        }
    }

    template <class Field, bool ADD>
    inline typename std::enable_if<!FFLAS::support_simd_add<typename Field::Element>::value, void>::type
    fadd (const Field & F,  const size_t N,
          typename Field::ConstElement_ptr A, const size_t inca,
          typename Field::ConstElement_ptr B, const size_t incb,
          typename Field::Element_ptr C, const size_t incc
          , FieldCategories::UnparametricTag
         )
    {
        typename Field::ConstElement_ptr Ai=A, Bi=B;
        typename Field::Element_ptr Ci=C;
        for (size_t i=0; i<N; i++, Ai+=inca,Bi+=incb,Ci+=incc)
            if (ADD)
                *Ci = *Ai + *Bi;
            else
                *Ci = *Ai - *Bi;
    }

    template <class Field, bool ADD>
    inline typename std::enable_if<FFLAS::support_simd_add<typename Field::Element>::value, void>::type
    fadd (const Field & F,  const size_t N,
          typename Field::ConstElement_ptr A, const size_t inca,
          typename Field::ConstElement_ptr B, const size_t incb,
          typename Field::Element_ptr C, const size_t incc
          , FieldCategories::UnparametricTag
         )
    {
        if (inca == 1 && incb == 1 && incc == 1) {
            if (ADD)
                FFLAS::vectorised::add(C,A,B,N);
            else
                FFLAS::vectorised::sub(C,A,B,N);
        } else {
            typename Field::ConstElement_ptr Ai=A, Bi=B;
            typename Field::Element_ptr Ci=C;
            for (size_t i=0; i<N; i++, Ai+=inca,Bi+=incb,Ci+=incc)
                if (ADD)
                    *Ci = *Ai + *Bi;
                else
                    *Ci = *Ai - *Bi;
        }
    }



} // details
} // FFLAS


#endif // __FFLASFFPACK_fscal_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
