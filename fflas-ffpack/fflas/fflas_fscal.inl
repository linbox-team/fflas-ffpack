/* fflas/fflas_fscal.inl
 * Copyright (C) 2014 FFLAS-FFPACK group
 *
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

#ifndef __FFLASFFPACK_fscal_INL
#define __FFLASFFPACK_fscal_INL

namespace FFLAS { namespace vectorised { namespace unswitch {

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

    template<class Field, class SimdT>
    inline void
    VEC_SCAL(typename SimdT::vect_t & C, typename SimdT::vect_t & ALPHA, HelperModSimd<Field,SimdT> & H)
    {
        C = SimdT::mul(C,ALPHA);
        VEC_MOD(C, H);
    }

    template<class Field>
    inline typename std::enable_if<FFLAS::support_simd_mod<typename Field::Element>::value, void>::type
    scalp(const Field &F, typename Field::Element_ptr T, const typename Field::Element alpha, typename Field::ConstElement_ptr U, const size_t n, HelperMod<Field> &G)
    {
        typedef typename Field::Element Element;
        using simd = Simd<Element>;
        using vect_t = typename simd::vect_t;
        HelperModSimd<Field,simd> H(F,G);
        vect_t C,ALPHA;
        ALPHA = simd::set1(alpha);

        size_t i = 0;

        if (n < simd::vect_size)
        {
            for (; i < n ; i++)
            {
                T[i] = reduce(alpha*U[i], H);
            }
            return;

        }

        long st = long(T) % simd::alignment;
        if (st)
        { // the array T is not aligned (process few elements s.t. (T+i) is aligned)
            for (size_t j = static_cast<size_t>(st) ; j < simd::alignment ; j+=sizeof(Element), i++)
            {
                T[i] = reduce(alpha*U[i], H);
            }
        }
        FFLASFFPACK_check((long(T+i) % simd::alignment == 0));
        if ((long(U+i)%simd::alignment==0))
        {
            // perform the loop using SIMD
            for (;i <= n - simd::vect_size ; i += simd::vect_size)
            {
                C = simd::load(U+i);
                VEC_SCAL(C, ALPHA, H);
                simd::store(T+i,C);
            }
        }
        // perform the last elt from T without SIMD
        for (; i < n ; i++)
        {
            T[i] = reduce(alpha*U[i], H);
        }
    }

#endif // __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

    template<class Field>
    inline typename std::enable_if<!FFLAS::support_simd_mod<typename Field::Element>::value && FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    scalp(const Field &F, typename Field::Element_ptr T, const typename Field::Element alpha, typename Field::ConstElement_ptr U, const size_t n, HelperMod<Field> &H)
    {
        size_t i = 0;

        for (; i < n ; i++)
        {
            T[i] = reduce(alpha*U[i], H);
        }
    }

    template<class Field>
    inline typename std::enable_if<FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    scalp(const Field &F, typename Field::Element_ptr T, const typename Field::Element alpha, typename Field::ConstElement_ptr U, const size_t n, const size_t & incX, HelperMod<Field> &H)
    {
        size_t i = 0;
        typename Field::ConstElement_ptr Xi = U;
        for (; Xi < U+n*incX ; Xi+=incX,i+=incX)
        {
            T[i]=reduce(alpha*(*Xi),H);
        }
    }

    template<class Field>
    inline typename std::enable_if<FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    scalp(const Field &F, typename Field::Element_ptr T, const typename Field::Element alpha, typename Field::ConstElement_ptr U, const size_t n,
          const size_t & incX, const size_t & incY, HelperMod<Field> &H)
    {
        typename Field::ConstElement_ptr Xi = U;
        typename Field::Element_ptr Yi = T;
        for (; Xi < U+n*incX ; Xi+=incX,Yi+=incY)
        {
            *Yi = reduce(alpha*(*Xi),H);
        }
    }

} // unswitch
} // vectorised
} // FFLAS

namespace FFLAS { namespace vectorised {

    // simd_mod => fast_mod, so this declaration's one in two
    template<class Field>
    inline typename std::enable_if<FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    scalp(const Field &F, typename Field::Element_ptr T, const typename Field::Element alpha, typename Field::ConstElement_ptr U, const size_t n)
    {
        HelperMod<Field> H(F);

        unswitch::scalp(F,T,alpha,U,n,H);
    }

    template<class Field>
    inline typename std::enable_if<FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    scalp(const Field &F, typename Field::Element_ptr T, const typename Field::Element alpha, typename Field::ConstElement_ptr U, const size_t n, const size_t & incX)
    {
        HelperMod<Field> H(F);

        unswitch::scalp(F,T,alpha,U,n,incX,H);
    }

    template<class Field>
    inline typename std::enable_if<FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    scalp(const Field &F, typename Field::Element_ptr T, const typename Field::Element alpha, typename Field::ConstElement_ptr U, const size_t n, const size_t & incX, const size_t & incY)
    {
        HelperMod<Field> H(F);

        unswitch::scalp(F,T,alpha,U,n,incX,incY,H);
    }


} // vectorised
} // FFLAS

namespace FFLAS { namespace details {

    /*** Entry-points for specialised code ***/
    // simd_mod => fast_mod;
    // will default to the best supported option in unswitch

    template<class Field>
    inline typename std::enable_if<FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    fscalin( const Field & F , const size_t N,
             const typename Field::Element a,
             typename Field::Element_ptr  X, const size_t incX, FieldCategories::ModularTag)
    {
        if(incX == 1)
        {
            vectorised::scalp(F,X,a,X,N);
        }
        else
        {
            if (N < FFLASFFPACK_COPY_REDUCE) // defined in fflas_freduce.inl
            {
                vectorised::scalp(F,X,a,X,N,incX);
            }
            else
            {
                typename Field::Element_ptr Xc = fflas_new (F,N) ;
                fassign (F,N,X,incX,Xc,1);
                vectorised::scalp(F,Xc,a,Xc,N);
                fassign (F,N,Xc,1,X,incX);
                fflas_delete (Xc);
            }
        }
    }

    template<class Field>
    inline typename std::enable_if<FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    fscal( const Field & F , const size_t N,
           const typename Field::Element a,
           typename  Field::ConstElement_ptr X, const size_t incX,
           typename Field::Element_ptr Y, const size_t incY, FieldCategories::ModularTag)
    {
        if (F.isOne(a)) {
            fassign(F,N,X,incX,Y,incY);
            return ;
        }

        if (F.areEqual(a,F.mOne)){
            fneg(F,N,X,incX,Y,incY);
            return;
        }

        if (F.isZero(a)){
            fzero(F,N,Y,incY);
            return;
        }

        if(incX == 1 && incY==1) {
            vectorised::scalp(F,Y,a,X,N);
        }
        else
        {
            if (N < FFLASFFPACK_COPY_REDUCE)
            {
                vectorised::scalp(F,Y,a,X,N,incX,incY);
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

                vectorised::scalp(F,Yc,a,Xc,N);

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
    fscalin (const Field& F, const size_t n, const typename Field::Element a,
             typename Field::Element_ptr X, const size_t incX, FC)
    {
        if (F.isOne(a))
            return ;

        if (F.isMOne(a)){
            fnegin(F,n,X,incX);
            return;
        }

        if (F.isZero(a)){
            fzero(F,n,X,incX);
            return;
        }

        typename Field::Element_ptr Xi = X ;

        if ( incX == 1)
            for (size_t i = 0 ; i < n ; ++i)
                F.mulin( X[i], a);
        else

            for (; Xi < X+n*incX; Xi+=incX )
                F.mulin( *Xi, a);
    }


    template<class Field, class FC>
    inline void
    fscal( const Field& F, const size_t N,
           const typename Field::Element a,
           typename Field::ConstElement_ptr X, const size_t incX,
           typename Field::Element_ptr Y, const size_t incY, FC)
    {
        if (F.isOne(a)) {
            fassign(F,N,X,incX,Y,incY);
            return ;
        }

        if (F.areEqual(a,F.mOne)){
            fneg(F,N,X,incX,Y,incY);
            return;
        }

        if (F.isZero(a)){
            fzero(F,N,Y,incY);
            return;
        }

        typename Field::ConstElement_ptr Xi = X;
        typename Field::Element_ptr Yi = Y;
        if (incX == 1 && incY == 1)
            for (size_t i = 0 ; i < N ; ++i)
                F.mul( Y[i], a, X[i] );
        else
            for (; Xi < X+N*incX; Xi+=incX, Yi+=incY )
                F.mul( *Yi, a, *Xi );
    }
} // details
} // FFLAS

namespace FFLAS {

    /***************************/
    /*         LEVEL 1         */
    /***************************/

    template<class Field>
    inline void
    fscalin (const Field& F, const size_t n, const typename Field::Element a,
             typename Field::Element_ptr X, const size_t incX)
    {
        details::fscalin(F,n,a,X,incX,typename FieldTraits<Field>::category());
    }

    template<class Field>
    inline void
    fscal( const Field& F, const size_t N,
           const typename Field::Element a,
           typename Field::ConstElement_ptr X, const size_t incX,
           typename Field::Element_ptr Y, const size_t incY)
    {
        details::fscal(F,N,a,X,incX,Y,incY,typename FieldTraits<Field>::category());
    }


    template<>
    inline void
    fscal( const Givaro::DoubleDomain& , const size_t N,
           const Givaro::DoubleDomain::Element a,
           Givaro::DoubleDomain::ConstElement_ptr x, const size_t incx,
           Givaro::DoubleDomain::Element_ptr y, const size_t incy )
    {
#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_dcopy( (int)N, x, (int)incy, y, (int)incy);
        cblas_dscal( (int)N, a, y, (int)incy);
    }

    template<>
    inline void
    fscal( const Givaro::FloatDomain& , const size_t N,
           const Givaro::FloatDomain::Element a,
           Givaro::FloatDomain::ConstElement_ptr x, const size_t incx,
           Givaro::FloatDomain::Element_ptr y, const size_t incy )
    {
#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_scopy( (int)N, x, (int)incy, y, (int)incy);
        cblas_sscal( (int)N, a, y, (int)incy);
    }

    template<>
    inline void
    fscalin( const Givaro::DoubleDomain& , const size_t N,
             const Givaro::DoubleDomain::Element a,
             Givaro::DoubleDomain::Element_ptr y, const size_t incy )
    {

#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_dscal( (int)N, a, y, (int)incy);
    }

    template<>
    inline void
    fscalin( const Givaro::FloatDomain& , const size_t N,
             const Givaro::FloatDomain::Element a,
             Givaro::FloatDomain::Element_ptr y, const size_t incy )
    {

#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_sscal( (int)N, a, y, (int)incy);
    }

    /***************************/
    /*         LEVEL 2         */
    /***************************/

    template<class Field>
    inline void
    fscalin (const Field& F, const size_t m , const size_t n,
             const typename Field::Element a,
             typename Field::Element_ptr A, const size_t lda)
    {
        if (F.isOne(a)) {
            return ;
        }
        else if (F.isZero(a)) {
            fzero(F,m,n,A,lda);
        }
        else if (F.isMOne(a)) {
            fnegin(F,m,n,A,lda);
        }
        else {
            if (lda == n) {
                fscalin(F,n*m,a,A,1);
            }
            else {
                for (size_t i = 0 ; i < m ; ++i)
                {
                    fscalin(F,n,a,A+i*lda,1);
                }
            }

            return;
        }
    }

    template<class Field>
    inline void
    fscal (const Field& F, const size_t m , const size_t n,
           const typename Field::Element a,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::Element_ptr B, const size_t ldb)
    {
        if (F.isOne(a)) {
            fassign(F,m,n,A,lda,B,ldb) ;
        }
        else if (F.isZero(a)) {
            fzero(F,m,n,B,ldb);
        }
        else if (F.isMOne(a)) {
            fneg(F,m,n,A,lda,B,ldb);
        }
        else {
            if (n == lda && m == lda)
                fscal(F,m*n,a,A,lda,B,ldb);
            else {
                for (size_t i = 0; i < m ; ++i)
                    fscal(F,n,a,A+i*lda,1,B+i*ldb,1);
            }
        }

        return;
    }

} // FFLAS

#endif // __FFLASFFPACK_fscal_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
