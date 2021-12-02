/* fflas/fflas_freduce.inl
 * Copyright (C) 2014 Pascal Giorgi
 *
 * Written by Pascal Giorgi <Pascal.Giorgi@lirmm.fr>
 * Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 * Pierre Karpman <pierre.karpman@univ-grenoble-alpes.fr>
 *
 * Part of this code is taken from http://libdivide.com/
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

#ifndef __FFLASFFPACK_fflas_freduce_INL
#define __FFLASFFPACK_fflas_freduce_INL

#include <givaro/udl.h>

#include "fflas-ffpack/fflas/fflas_fassign.h"

#define FFLASFFPACK_COPY_REDUCE 32 /*  TODO TO BENCMARK LATER */

namespace FFLAS { namespace vectorised { /*  for casts (?) */

    template<class T>
    inline typename std::enable_if< ! std::is_integral<T>::value, T>::type
    reduce(T A, T B)
    {
        return fmod(A,B);
    }

    template<class T>
    inline typename std::enable_if< std::is_integral<T>::value, T>::type
    reduce(T  A, T B)
    {
        return A % B; // B > 0
    }

    template<>
    inline Givaro::Integer reduce(Givaro::Integer  A, Givaro::Integer B) // @bug B is not integer, but uint64_t usually
    {
        return A % B; // B > 0
    }

    inline float reduce(float A, float B, float invB, float min, float max)
    {
        float Q = A * invB;
        Q = floorf(Q);
        A = A - Q*B;
        A = A < min ? A + B : A;
        A = A > max ? A - B : A;
        return A;
    }

    inline double reduce(double A, double B, double invB, double min, double max)
    {
        //std::cerr<<"fmod"<<std::endl;
        double Q = A * invB;
        Q = floor(Q);
        A = A - Q*B;
        A = A < min ? A + B : A;
        A = A > max ? A - B : A;
        return A;
    }

    inline int64_t reduce(int64_t A, int64_t p, double invp, double min, double max, int64_t pow50rem)
    {
        // assert(p < 1LL << 33);
        // nothing so special with 50; could be something else

        int64_t Aq50 = A >> 50;                         // Aq50 < 2**14
        int64_t Ar50 = A & 0x3FFFFFFFFFFFFLL;           // Ar50 < 2**50

        int64_t Aeq  = Aq50 * pow50rem + Ar50;          // Aeq < 2**47 + 2**50 < 2**51; Aeq ~ A mod p

        return static_cast<int64_t>(reduce(static_cast<double>(Aeq),static_cast<double>(p), invp, min, max));
    }


} // vectorised
} // FFLAS

namespace FFLAS { namespace vectorised {

    template<class Field, class ElementTraits = typename ElementTraits<typename Field::Element>::value>
    struct HelperMod  ;


    template<class Field>
    struct HelperMod<Field, ElementCategories::MachineIntTag> {
        typename Field::Element p;
        double invp;
        double min;
        double max;
        int64_t pow50rem;

        HelperMod()
        {
//             std::cout << "empty cstor called" << std::endl;
        } ;

        HelperMod( const Field & F)
        {
//             std::cout << "field cstor called" << std::endl;
            p =  (typename Field::Element) F.characteristic();
            pow50rem = (1LL << 50) % p;
            invp = 1/static_cast<double>(p);
            min = static_cast<double>(F.minElement());
            max = static_cast<double>(F.maxElement());
        }

    } ;

    template<class Field>
    struct HelperMod<Field, FFLAS::ElementCategories::MachineFloatTag> {
        typename Field::Element p;
        typename Field::Element invp;
        typename Field::Element min ;
        typename Field::Element max ;

        HelperMod() {} ;

        HelperMod( const Field & F)
        {
            p = (typename Field::Element) F.characteristic();
            invp = (typename Field::Element)1/p;
            min = F.minElement();
            max = F.maxElement();
        }

    } ;

    template<class Field>
    struct HelperMod<Field, FFLAS::ElementCategories::ArbitraryPrecIntTag> {
        typename Field::Element p;
        // typename Field::Element invp;
        // typename Field::Elmeent min ;
        // typename Field::Elmeent max ;

        HelperMod() {} ;

        HelperMod( const Field & F)
        {
            p = (typename Field::Element) F.characteristic();
            // invp = (typename Field::Element)1/p;
            // min = F.minElement();
            // max = F.maxElement();
        }

    } ;

    template<class Field>
    struct HelperMod<Field, FFLAS::ElementCategories::FixedPrecIntTag> {
        typename Field::Element p;
        // typename Field::Element invp;
        // typename Field::Elmeent min ;
        // typename Field::Elmeent max ;

        HelperMod() {} ;

        HelperMod( const Field & F)
        {
            p = (typename Field::Element) F.characteristic();
            // invp = (typename Field::Element)1/p;
            // min = F.minElement();
            // max = F.maxElement();
        }

    } ;


#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
    template<class Field, class SimdT, class ElementTraits = typename ElementTraits<typename Field::Element>::value>
    struct HelperModSimd  ;

    template<class Field, class SimdT>
    struct HelperModSimd<Field, SimdT, ElementCategories::MachineIntTag> : public HelperMod<Field> {
        typedef typename SimdT::vect_t vect_t ;
#ifndef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
        // with AVX but not AVX2, integral vectors are on 128 bits only
        Simd128<double>::vect_t P ;
        Simd128<double>::vect_t MIN ;
        Simd128<double>::vect_t MAX ;
        Simd128<double>::vect_t NEGP ;
        Simd128<double>::vect_t Q ;
        Simd128<double>::vect_t T ;
        Simd128<double>::vect_t INVP;
#else
        Simd<double>::vect_t P ;
        Simd<double>::vect_t MIN ;
        Simd<double>::vect_t MAX ;
        Simd<double>::vect_t NEGP ;
        Simd<double>::vect_t Q ;
        Simd<double>::vect_t T ;
        Simd<double>::vect_t INVP;
#endif
        vect_t POW50REM;

        HelperModSimd ( const Field & F) :
            HelperMod<Field>(F)
        {
//             std::cout << "HelperMod constructed " << this->shift << std::endl;
#ifndef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
            P       = Simd128<double>::set1(static_cast<double>(this->p));
            NEGP    = Simd128<double>::set1(-static_cast<double>(this->p));
            MIN     = Simd128<double>::set1(this->min);
            MAX     = Simd128<double>::set1(this->max);
            INVP    = Simd128<double>::set1(this->invp);
#else
            P       = Simd<double>::set1(static_cast<double>(this->p));
            NEGP    = Simd<double>::set1(-static_cast<double>(this->p));
            MIN     = Simd<double>::set1(this->min);
            MAX     = Simd<double>::set1(this->max);
            INVP    = Simd<double>::set1(this->invp);
#endif
            POW50REM= SimdT::set1(this->pow50rem);
        }

        HelperModSimd( const Field & F, const HelperMod<Field> & G)
        {
            this->p         = G.p;
            this->invp      = G.invp;
            this->min       = G.min;
            this->max       = G.max;
            this->pow50rem  = G.pow50rem;
#ifndef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
            P               = Simd128<double>::set1(static_cast<double>(this->p));
            NEGP            = Simd128<double>::set1(-static_cast<double>(this->p));
            MIN             = Simd128<double>::set1(this->min);
            MAX             = Simd128<double>::set1(this->max);
            INVP            = Simd128<double>::set1(this->invp);
#else
            P               = Simd<double>::set1(static_cast<double>(this->p));
            NEGP            = Simd<double>::set1(-static_cast<double>(this->p));
            MIN             = Simd<double>::set1(this->min);
            MAX             = Simd<double>::set1(this->max);
            INVP            = Simd<double>::set1(this->invp);
#endif
            POW50REM        = SimdT::set1(this->pow50rem);
        }

    } ;

    template<class Field, class SimdT>
    struct HelperModSimd<Field, SimdT, ElementCategories::MachineFloatTag>  : public HelperMod<Field> {
        typedef typename SimdT::vect_t vect_t ;
        vect_t INVP;
        vect_t MIN ;
        vect_t MAX ;
        vect_t NEGP ;
        vect_t P ;
        vect_t Q ;
        vect_t T ;

        HelperModSimd( const Field & F) :
            HelperMod<Field>(F)
        {
            P   = SimdT::set1(this->p);
            NEGP= SimdT::set1(-(this->p));
            MIN = SimdT::set1(this->min);
            MAX = SimdT::set1(this->max);
            INVP= SimdT::set1(this->invp);
        }

        HelperModSimd( const Field & F, const HelperMod<Field> & G)
        {
            this->p     = G.p;
            this->invp  = G.invp;
            this->min   = G.min;
            this->max   = G.max;
            P   = SimdT::set1(this->p);
            NEGP= SimdT::set1(-this->p);
            MIN = SimdT::set1(this->min);
            MAX = SimdT::set1(this->max);
            INVP= SimdT::set1(this->invp);
        }
    } ;
#endif // __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS


    //TODO why should this be in this ifdef??
#ifdef __x86_64__
    template<class Field>
    inline typename std::enable_if< std::is_same<typename Field::Element,int64_t>::value , int64_t>::type
    reduce (typename Field::Element A, HelperMod<Field,ElementCategories::MachineIntTag> & H)
    {
        return reduce(A, H.p, H.invp, H.min, H.max, H.pow50rem);
    }
#endif // __x86_64__


    template<class Field>
#ifdef __x86_64__
    inline typename std::enable_if< ! std::is_same<typename Field::Element,int64_t>::value , typename Field::Element>::type
#else
    inline typename Field::Element
#endif // __x86_64__
    reduce (typename Field::Element A, HelperMod<Field,ElementCategories::MachineIntTag> & H)
    {
        bool positive = !FieldTraits<Field>::balanced;
        typename Field::Element r = reduce(A,H.p);

        if(!positive)
        {
            r = r > H.max ? r - H.p : r;
        }
        r = r < H.min ? r + H.p : r;

        return r;
    }

    template<class Field>
    inline
    typename Field::Element reduce (typename Field::Element A, HelperMod<Field,ElementCategories::MachineFloatTag> & H)
    {
        return reduce(A, H.p, H.invp, H.min, H.max);
    }

    template<class Field>
    inline typename Field::Element reduce (typename Field::Element A, HelperMod<Field,ElementCategories::ArbitraryPrecIntTag> & H)
    {
        bool positive = !FieldTraits<Field>::balanced;
        typename Field::Element r = reduce(A,H.p);

        if(!positive)
        {
            r = r > H.max ? r - H.p : r;
        }
        r = r < H.min ? r + H.p : r;

        return r;
    }

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

    template<class Field, class SimdT>
    inline void
    VEC_MOD(typename SimdT::vect_t & C, HelperModSimd<Field,SimdT,ElementCategories::MachineFloatTag> & H)
    {
        C = SimdT::mod( C, H.P, H.INVP, H.NEGP, H.MIN, H.MAX, H.Q, H.T );
    }

    template<class Field, class SimdT>
    inline void
    VEC_MOD(typename SimdT::vect_t & C, HelperModSimd<Field,SimdT,ElementCategories::MachineIntTag> & H)
    {
        C = SimdT::mod(C, H.P, H.INVP, H.NEGP, H.POW50REM, H.MIN, H.MAX, H.Q, H.T);
    }

#endif // __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

} // vectorised
} // FFLAS

namespace FFLAS  { namespace vectorised { namespace unswitch  {

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
    template<class Field>
    inline typename std::enable_if<FFLAS::support_simd_mod<typename Field::Element>::value, void>::type
    modp(const Field &F, typename Field::ConstElement_ptr U, const size_t & n,
         typename Field::Element_ptr T
         , HelperMod<Field> & G
        )
    {

        // std::cerr<<"modp vectorized"<<std::endl;
        typedef typename Field::Element Element;
        using simd = Simd<Element>;
        using vect_t = typename simd::vect_t;
        HelperModSimd<Field,simd> H(F,G);

        size_t i = 0;
        if (n < simd::vect_size)
        {
            // std::cerr<< n<< " < "<<simd::vect_size<<std::endl;
            for (; i < n ; i++)
            {
                T[i]=reduce(U[i],H);
            }
            return;
        }

        long st = long(T) % simd::alignment;

        // the array T is not aligned (process few elements s.t. (T+i) is 32 bytes aligned)
        if (st)
        {
            // std::cerr<< st << " not aligned on  "<<simd::alignment<<std::endl;

            for (size_t j = static_cast<size_t>(st) ; j < simd::alignment ; j += sizeof(Element), i++)
            {
                T[i] = reduce(U[i],H);
            }
        }

        FFLASFFPACK_check((long(T+i) % simd::alignment == 0));

        vect_t C ;

        if((long(U+i) % simd::alignment == 0))
        {
            // perform the loop using SIMD
            for (; i<= n - simd::vect_size ; i += simd::vect_size)
            {
                C = simd::load(U + i);

                VEC_MOD<Field,simd>(C,H);
                simd::store(T+i, C);
            }
        }

        // perform the last elt from T without SIMD
        // std::cerr<< n-i<< " unaligned elements left "<<std::endl;
        for (;i<n;i++)
        {

            T[i] = reduce(U[i],H);
        }
    }
#endif

    // not vectorised but allows better code than % or fmod via helper
    template<class Field>
    inline typename std::enable_if<!FFLAS::support_simd_mod<typename Field::Element>::value && FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    modp(const Field &F, typename Field::ConstElement_ptr U, const size_t & n,
         typename Field::Element_ptr T
         , HelperMod<Field> & H
        )
    {
        size_t i = 0;
        for (; i < n ; i++)
        {
            T[i]=reduce(U[i],H);
        }
    }

    // Using fast helper-based reduce for non-local input
    template<class Field>
    inline typename std::enable_if<FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    modp(const Field &F, typename Field::ConstElement_ptr U, const size_t & n, const size_t & incX,
         typename Field::Element_ptr T
         , HelperMod<Field> & H
        )
    {
        size_t i = 0;
        typename Field::ConstElement_ptr Xi = U;
        for (; Xi < U+n*incX ; Xi+=incX,i+=incX)
        {
            T[i]=reduce(*Xi,H);
        }
    }

} // unswitch
} // vectorised
} // FFLAS

namespace FFLAS { namespace vectorised {


    // simd_mod => fast_mod, so this declaration's one in two
    template<class Field>
    inline typename std::enable_if<FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    modp(const Field &F, typename Field::ConstElement_ptr U, const size_t & n,
         typename Field::Element_ptr T)
    {
        HelperMod<Field> H(F);

        unswitch::modp(F,U,n,T,H);
    }

    template<class Field>
    inline typename std::enable_if<FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    modp(const Field &F, typename Field::ConstElement_ptr U, const size_t & n, const size_t & incX,
         typename Field::Element_ptr T)
    {
        HelperMod<Field> H(F);

        unswitch::modp(F,U,n,incX,T,H);
    }

} // vectorised
} // FFLAS


namespace FFLAS { namespace details {


    /*** Entry-points for specialised code ***/
    // simd_mod => fast_mod;
    // will default to the best supported option in unswitch

    template<class Field>
    inline typename std::enable_if<FFLAS::support_fast_mod<typename Field::Element>::value  , void>::type
    freduce (const Field & F, const size_t m,
             typename Field::Element_ptr A, const size_t incX, FieldCategories::ModularTag)
    {
        if(incX == 1)
        {
            vectorised::modp(F,A,m,A);
        }
        else
        { /*  faster with copy, use incX=1, copy back ? */
            if (m < FFLASFFPACK_COPY_REDUCE)
            {
                vectorised::modp(F,A,m,incX,A);
            }
            else
            {
                typename Field::Element_ptr Ac = fflas_new (F,m) ;
                fassign (F,m,A,incX,Ac,1);
                vectorised::modp(F,Ac,m,Ac);
                fassign (F,m,Ac,1,A,incX);
                fflas_delete (Ac);
            }
        }
    }

    template<class Field>
    inline typename std::enable_if< FFLAS::support_fast_mod<typename Field::Element>::value, void>::type
    freduce (const Field & F, const size_t m,
             typename Field::ConstElement_ptr  B, const size_t incY,
             typename Field::Element_ptr A, const size_t incX,
             FieldCategories::ModularTag)
    {
        if(incX == 1 && incY == 1)
        {
            vectorised::modp(F,B,m,A);
        }
        else
        {
            if (m < FFLASFFPACK_COPY_REDUCE)
            {
                // TODO also write a ``fast'' modp that handles an incY
                typename Field::Element_ptr Xi = A ;
                typename Field::ConstElement_ptr Yi = B ;
                for (; Xi < A+m*incX; Xi+=incX, Yi += incY )
                    F.reduce (*Xi , *Yi);
            }
            else
            {
                typename Field::Element_ptr Ac;
                typename Field::Element_ptr Bc;
                if (incX != 1)
                {
                    Ac = fflas_new (F,m);
                    fassign(F,m,A,incX,Ac,1);
                }
                else
                {
                    Ac = A;
                }
                if (incY != 1)
                {
                    Bc = fflas_new (F,m);
                    fassign(F,m,B,incY,Bc,1);
                }
                else
                {
                    Bc = const_cast<typename Field::Element_ptr>(B); // Oh the horror
                }

                vectorised::modp(F,Bc,m,Ac);

                if (incX != 1)
                {
                    fassign(F,m,Ac,1,A,incX);
                    fflas_delete(Ac);
                }
                if (incY != 1)
                {
                    fflas_delete(Bc);
                }
            }
        }
    }

    /*** Non-specialised code ***/

    template<class Field, class FC>
    inline void
    freduce (const Field & F, const size_t m,
             typename Field::Element_ptr A, const size_t incX,
             FC)
    {
        typename Field::Element_ptr Xi = A ;
        for (; Xi < A+m*incX; Xi+=incX )
            F.reduce (*Xi);
    }

    template<class Field, class FC>
    inline void
    freduce (const Field & F, const size_t m,
             typename Field::ConstElement_ptr  B, const size_t incY,
             typename Field::Element_ptr A, const size_t incX,
             FC)
    {
        typename Field::Element_ptr Xi = A ;
        typename Field::ConstElement_ptr Yi = B ;
        for (; Xi < A+m*incX; Xi+=incX, Yi += incY )
            F.reduce (*Xi , *Yi);
    }

} // details
} // FFLAS


#endif // __FFLASFFPACK_fflas_freduce_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
