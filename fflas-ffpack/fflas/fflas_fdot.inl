/* fflas_fdot.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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

#ifndef __FFLASFFPACK_fdot_INL
#define __FFLASFFPACK_fdot_INL

#include "fflas-ffpack/fflas/fflas_helpers.inl"
// Default implementation
// Specializations should be written
// to increase efficiency


namespace FFLAS {

    template<class Field>
    inline typename Field::Element
    fdot( const Field& F, const size_t N,
          typename Field::ConstElement_ptr x, const size_t incx,
          typename Field::ConstElement_ptr y, const size_t incy ,
          ModeCategories::DefaultTag& MT)
    {

        typename Field::Element d;
        typename Field::ConstElement_ptr xi = x;
        typename Field::ConstElement_ptr yi = y;
        F.init(d);
        F.assign(d, F.zero);
        for ( ; xi < x+N*incx; xi+=incx, yi+=incy )
            F.axpyin( d, *xi, *yi );
        return d;
    }

    template<class Field>
    inline typename Field::Element
    fdot( const Field& F, const size_t N,
          typename Field::ConstElement_ptr x, const size_t incx,
          typename Field::ConstElement_ptr y, const size_t incy ,
          ModeCategories::DelayedTag& MT)
    {
        typedef typename associatedDelayedField<const Field>::field DelayedField;
        typedef typename associatedDelayedField<const Field>::type DelayedField_t;
        typedef typename DelayedField::Element DFElt;
        DelayedField_t delayedF;

        const DFElt MaxStorableValue = limits<typename DelayedField::Element>::max();
        const DFElt AbsMax = std::max(-F.minElement(), F.maxElement());
        const DFElt r = (MaxStorableValue / AbsMax) / AbsMax;
        size_t delayedDim = FFLAS::Protected::min_types<DFElt>(r);
        if (!delayedDim){
            ModeCategories::DefaultTag DT;
            return fdot(F, N, x, incx, y, incy, DT);
        }
        typename Field::Element d;
        F.init (d);
        F.assign (d, F.zero);
        ModeCategories::DefaultTag DM;
        typename Field::ConstElement_ptr xi = x, yi = y;
        size_t i=delayedDim;
        typename Field::Element dp;
        for (; i<N; i+= delayedDim, xi += incx*delayedDim, yi += incy*delayedDim){
            F.init(dp, fdot(delayedF, delayedDim, xi, incx, yi, incy, DM));
            F.addin(d, dp);
        }
        F.init(dp, fdot(delayedF, N+delayedDim-i, xi, incx, yi, incy, DM));
        F.addin (d, dp);
        return d;
    }

    template<>
    inline Givaro::DoubleDomain::Element
    fdot( const Givaro::DoubleDomain& , const size_t N,
          Givaro::DoubleDomain::ConstElement_ptr x, const size_t incx,
          Givaro::DoubleDomain::ConstElement_ptr y, const size_t incy,
          ModeCategories::DefaultTag& MT)
    {

#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        return cblas_ddot( (int)N, x, (int)incx, y, (int)incy );
    }

    template<>
    inline Givaro::FloatDomain::Element
    fdot( const Givaro::FloatDomain& , const size_t N,
          Givaro::FloatDomain::ConstElement_ptr x, const size_t incx,
          Givaro::FloatDomain::ConstElement_ptr y, const size_t incy,
          ModeCategories::DefaultTag& MT)
    {

#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        return cblas_sdot( (int)N, x, (int)incx, y, (int)incy );
    }

    template<class Field, class T>
    inline typename Field::Element
    fdot( const Field& F, const size_t N,
          typename Field::ConstElement_ptr x, const size_t incx,
          typename Field::ConstElement_ptr y, const size_t incy,
          ModeCategories::ConvertTo<T>& MT)
    {
        typename ModeCategories::DefaultTag mt;
        return fdot (F, N, x, incx, y, incy, mt);
    }

    template<class Field>
    inline typename Field::Element
    fdot( const Field& F, const size_t N,
          typename Field::ConstElement_ptr x, const size_t incx,
          typename Field::ConstElement_ptr y, const size_t incy,
          ModeCategories::DefaultBoundedTag& dbt)
    {
        // Nothing special in fdot for Bounded
        typename ModeCategories::DefaultTag dt;
        return fdot (F, N, x, incx, y, incy, dt);
    }

    template<class Field>
    inline typename Field::Element
    fdot( const Field& F, const size_t N,
          typename Field::ConstElement_ptr x, const size_t incx,
          typename Field::ConstElement_ptr y, const size_t incy,
          const ParSeqHelper::Sequential seq)
    {
        typename ModeTraits<Field>::value mt;
        return fdot (F, N, x, incx, y, incy, mt);
    }

    template<class Field>
    inline typename Field::Element
    fdot( const Field& F, const size_t N,
          typename Field::ConstElement_ptr x, const size_t incx,
          typename Field::ConstElement_ptr y, const size_t incy)
    {
        return fdot (F, N, x, incx, y, incy, FFLAS::ParSeqHelper::Sequential());
    }

} // FFLAS

#endif // __FFLASFFPACK_fdot_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
