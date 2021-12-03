/* fflas/fflas_faxpy.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 * Pierre Karpman <pierre.karpman@univ-grenoble-alpes.fr>
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

#ifndef __FFLASFFPACK_faxpy_INL
#define __FFLASFFPACK_faxpy_INL

/* most of this is a shameless copy from fscal... */

namespace FFLAS { namespace vectorised { namespace unswitch {

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

    // TODO check if returning in register is better??
    template<class Field, class SimdT>
    inline void
    VEC_FMA(typename SimdT::vect_t & A, typename SimdT::vect_t & X, typename SimdT::vect_t & Y, HelperModSimd<Field,SimdT> & H)
    {
        Y = SimdT::fmadd(Y,A,X);
        VEC_MOD(Y, H);
    }

    template<class Field>
    inline typename std::enable_if<FFLAS::support_simd_mod<typename Field::Element>::value, void>::type
    axpyp(const Field &F, const typename Field::Element a, typename Field::ConstElement_ptr X, typename Field::Element_ptr Y, const size_t n, HelperMod<Field> &G)
    {
        typedef typename Field::Element Element;
        using simd = Simd<Element>;
        using vect_t = typename simd::vect_t;
        HelperModSimd<Field,simd> H(F,G);
        vect_t A, U, V;
        A = simd::set1(a);

        size_t i = 0;

        if (n < simd::vect_size)
        {
            for (; i < n ; i++)
            {
                Y[i] = reduce(a*X[i]+Y[i], H);
            }
            return;
        }

        long st = long(Y) % simd::alignment;
        if (st)
        { // the array Y is not aligned (process few elements until (Y+i) is aligned)
            for (size_t j = static_cast<size_t>(st) ; j < simd::alignment ; j+=sizeof(Element), i++)
            {
                Y[i] = reduce(a*X[i]+Y[i], H);
            }
        }
        FFLASFFPACK_check((long(Y+i) % simd::alignment == 0));
        if (long(Y+i)%simd::alignment==0)
        {
            // perform the loop using SIMD
            for (;i <= n - simd::vect_size ; i += simd::vect_size)
            {
                U = simd::loadu(X+i); // could be transformed into a load at the cost of an additional top check. Worth it?
                V = simd::load(Y+i);
                VEC_FMA(A, U, V, H);
                simd::store(Y+i,V);
            }
        }
        // perform the last elt from T without SIMD
        for (; i < n ; i++)
        {
            Y[i] = reduce(a*X[i]+Y[i], H);
        }
    }

#endif // __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

    template<class Field>
    inline typename std::enable_if<!FFLAS::support_simd_mod<typename Field::Element>::value && FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    axpyp(const Field &F, const typename Field::Element a, typename Field::ConstElement_ptr X, typename Field::Element_ptr Y, const size_t n, HelperMod<Field> &H)
    {
        size_t i = 0;

        for (; i < n ; i++)
        {
            Y[i] = reduce(a*X[i] + Y[i], H);
        }
    }


    template<class Field>
    inline typename std::enable_if<FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    axpyp(const Field &F, const typename Field::Element a, typename Field::ConstElement_ptr X, typename Field::Element_ptr Y, const size_t n, const size_t incX, const size_t incY,
              HelperMod<Field> &H)
    {
        typename Field::ConstElement_ptr Xi = X;
        typename Field::Element_ptr Yi=Y;
        for (; Xi < X+n*incX; Xi+=incX, Yi+=incY )
            *Yi = reduce(a*(*Xi) + (*Yi), H);
    }

} // unswitch
} // vectorised
} // FFLAS

namespace FFLAS { namespace vectorised {

    // simd_mod => fast_mod, so this declaration's one in two
    template<class Field>
    inline typename std::enable_if<FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    axpyp(const Field &F, const typename Field::Element a, typename Field::ConstElement_ptr X, typename Field::Element_ptr Y, const size_t n)
    {
        HelperMod<Field> H(F);

        unswitch::axpyp(F,a,X,Y,n,H);
    }

    template<class Field>
    inline typename std::enable_if<FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    axpyp(const Field &F, const typename Field::Element a, typename Field::ConstElement_ptr X, typename Field::Element_ptr Y, const size_t n, const size_t incX, const size_t incY)
    {
        HelperMod<Field> H(F);

        unswitch::axpyp(F,a,X,Y,n,incX,incY,H);
    }

} // vectorised
} // FFLAS

namespace FFLAS { namespace details {

    /*** Entry-points for specialised code ***/
    // simd_mod => fast_mod;
    // will default to the best supported option in unswitch

    template<class Field>
    inline typename std::enable_if<FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    faxpy( const Field & F , const size_t N,
           const typename Field::Element a,
           typename Field::ConstElement_ptr X, const size_t incX,
           typename Field::Element_ptr Y, const size_t incY,
           FieldCategories::ModularTag)
    {
        if((incX == 1) && (incY == 1))
        {
            vectorised::axpyp(F,a,X,Y,N);
        }
        else
        {
            if (N < FFLASFFPACK_COPY_REDUCE) // defined in fflas_freduce.inl
            {
                vectorised::axpyp(F,a,X,Y,N,incX,incY);
            }
            else
            {
                typename Field::Element_ptr Xc;
                typename Field::Element_ptr Yc;
                if (incX != 1)
                {
                    Xc = fflas_new (F,N);
                    fassign(F,N,X,incX,Xc,1);
                }
                else
                {
                    Xc = const_cast<typename Field::Element_ptr>(X); // Oh the horror
                }
                if (incY != 1)
                {
                    Yc = fflas_new (F,N);
                    fassign(F,N,Y,incY,Yc,1);
                }
                else
                {
                    Yc = Y;
                }

                vectorised::axpyp(F,a,Xc,Yc,N);

                if (incY != 1)
                {
                    fassign(F,N,Yc,1,Y,incY);
                    fflas_delete(Yc);
                }
                if (incX != 1)
                {
                    fflas_delete(Xc);
                }
            }
        }
    }

    template<class Field, class FC>
    inline void
    faxpy(const Field& F, const size_t N, const typename Field::Element a,
          typename Field::ConstElement_ptr X, const size_t incX,
          typename Field::Element_ptr Y, const size_t incY, FC)
    {
        if (F.isZero(a))
            return ;

        if (F.isOne(a))
            return faddin(F,N,X,incX,Y,incY);
        //return fassign(F,N,X,incX,Y,incY);

        if (F.isMOne(a))
            return fsubin(F,N,X,incX,Y,incY);
        //return fneg(F,N,X,incX,Y,incY);

        typename Field::ConstElement_ptr Xi = X;
        typename Field::Element_ptr Yi = Y;
        for (; Xi < X+N*incX; Xi+=incX, Yi+=incY )
            F.axpyin( *Yi, a, *Xi );
    }

} // details
} // FFLAS

namespace FFLAS {

    /***************************/
    /*         LEVEL 1         */
    /***************************/

    template<class Field>
    inline void
    faxpy(const Field& F, const size_t n, const typename Field::Element a,
          typename Field::ConstElement_ptr X, const size_t incX,
          typename Field::Element_ptr Y, const size_t incY)
    {
        details::faxpy(F,n,a,X,incX,Y,incY,typename FieldTraits<Field>::category());
    }

    template<>
    inline void
    faxpy( const Givaro::DoubleDomain& , const size_t N,
           const Givaro::DoubleDomain::Element a,
           Givaro::DoubleDomain::ConstElement_ptr x, const size_t incx,
           Givaro::DoubleDomain::Element_ptr y, const size_t incy )
    {

#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_daxpy( (int)N, a, x, (int)incx, y, (int)incy);
    }

    template<>
    inline void
    faxpy( const Givaro::FloatDomain& , const size_t N,
           const Givaro::FloatDomain::Element a,
           Givaro::FloatDomain::ConstElement_ptr x, const size_t incx,
           Givaro::FloatDomain::Element_ptr y, const size_t incy )
    {

#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_saxpy( (int)N, a, x, (int)incx, y, (int)incy);
    }


    /***************************/
    /*         LEVEL 2         */
    /***************************/

    template<class Field>
    inline void
    faxpy( const Field& F, const size_t m, const size_t n,
           const typename Field::Element a,
           typename Field::ConstElement_ptr X, const size_t ldX,
           typename Field::Element_ptr Y, const size_t ldY )
    {

        if (F.isZero(a))
            return ;

        if (F.isOne(a))
            return faddin(F,m,n,X,ldX,Y,ldY);
        //return fassign(F,m,n,X,ldX,Y,ldY);

        if (F.isMOne(a))
            return fsubin(F,m,n,X,ldX,Y,ldY);
        //return fneg(F,m,n,X,ldX,Y,ldY);

        if (n == ldX && n == ldY)
            return faxpy(F,m*n,a,X,1,Y,1);

        typename Field::ConstElement_ptr Xi = X;
        typename Field::Element_ptr Yi=Y;
        for (; Xi < X+m*ldX; Xi+=ldX, Yi+=ldY )
            faxpy(F,n,a,Xi,1,Yi,1);
    }

} // FFLAS

#endif // __FFLASFFPACK_faxpy_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
