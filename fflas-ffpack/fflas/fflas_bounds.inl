/* fflas/fflas_bounds.inl
 * Copyright (C) 2008 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *            Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

#ifndef __FFLASFFPACK_fflas_bounds_INL
#define __FFLASFFPACK_fflas_bounds_INL

#define FFLAS_INT_TYPE uint64_t

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/utils/flimits.h"

#include <givaro/udl.h>
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>

namespace FFLAS { namespace Protected {

    template <class Field>
    inline double computeFactorClassic (const Field& F)
    {
        //FFLAS_INT_TYPE p=0;
        Givaro::Integer p=0;
        F.characteristic(p);
        return (double) (p-1);
    }

    /*************************************************************************************
     * Specializations for ModularPositive and ModularBalanced over double and float
     *************************************************************************************/
    template <>
    inline double computeFactorClassic (const Givaro::ModularBalanced<double>& F)
    {
        //FFLAS_INT_TYPE p;
        Givaro::Integer p;
        F.characteristic(p);
        return double((p-1) >> 1);
    }

    //BB: ajout, pourquoi pas ?
    template <>
    inline double computeFactorClassic (const Givaro::ModularBalanced<float>& F)
    {
        //FFLAS_INT_TYPE p;
        Givaro::Integer p;
        F.characteristic(p);
        return double((p-1) >> 1);
    }

    template <class Field>
    inline size_t DotProdBoundClassic (const Field& F,
                                       const typename Field::Element& beta
                                      )
    {

        //FFLAS_INT_TYPE p=0;
        Givaro::Integer p=0;
        F.characteristic(p);

        //unsigned long mantissa = Protected::Mantissa<typename Field::Element>();

        if (p == 0)
            return std::numeric_limits<size_t>::max();

        double kmax;
        {

            double c = computeFactorClassic(F);

            double cplt=0;
            if (!F.isZero (beta)){
                if (F.isOne (beta) || F.areEqual (beta, F.mOne)) cplt = c;
                else{
                    double be;
                    F.convert(be, beta);
                    cplt = fabs(be)*c;
                }
            }
            kmax = floor ( (double (double(limits<typename Field::Element>::max()) + 1 - cplt)) / (c*c));
            if (kmax  <= 1) return 1;
        }

        //kmax--; // we computed a strict upper bound
        return  (size_t) std::min ((uint64_t)kmax, 1_ui64 << 31);
    }

} // FFLAS
} // Protected

namespace FFLAS {

    inline Givaro::Integer
    InfNorm (const size_t M, const size_t N, const Givaro::Integer* A, const size_t lda){
        Givaro::Integer max = 0;
        for (size_t i=0; i<M; ++i)
            for (size_t j=0; j<N; ++j) {
                const Givaro::Integer & x(A[i*lda+j]);
                if (Givaro::absCompare(x,max)>0) max = x;
            }
        return abs(max);
    }

    namespace Protected {


        /**
         * TRSMBound
         *
         * \brief  computes the maximal size for delaying the modular reduction
         *         in a triangular system resolution
         *
         * This is the default version over an arbitrary field.
         * It is currently never used (the recursive algorithm is run until n=1 in this case)
         *
         * \param F Finite Field/Ring of the computation
         *
         */
        template <class Field>
        inline size_t TRSMBound (const Field&)
        {
            return 1;
        }

        // /**
        //  * Specialization for positive modular representation over double
        //  * Computes nmax s.t. (p-1)/2*(p^{nmax-1} + (p-2)^{nmax-1}) < 2^53
        //  * See [Dumas Giorgi Pernet 06, arXiv:cs/0601133]
        //  */
        // template<>
        // inline size_t TRSMBound (const Givaro::Modular<double>& F)
        // {

        // 	FFLAS_INT_TYPE pi;
        // 	F.characteristic(pi);
        // 	unsigned long p = pi;
        // 	unsigned long long p1(1), p2(1);
        // 	size_t nmax = 0;
        // 	unsigned long long max = ( (1 << (DBL_MANT_DIG + 1) ) / ((unsigned long long)(p - 1)));
        // 	while ( (p1 + p2) < max ){
        // 		p1*=p;
        // 		p2*=p-2;
        // 		nmax++;
        // 	}
        // 	return nmax;
        // }


        /**
         * Specialization for positive modular representation over float.
         * Computes nmax s.t. (p-1)/2*(p^{nmax-1} + (p-2)^{nmax-1}) < 2^24
         * @pbi
         * See [Dumas Giorgi Pernet 06, arXiv:cs/0601133]
         */
        template<class Element>
        inline size_t TRSMBound (const Givaro::Modular<Element>& F)
        {

            FFLAS_INT_TYPE pi;
            F.characteristic(pi);
            double p = pi;
            double p1 = 1.0, p2 = 1.0;
            double pm1 = (p - 1) / 2;
            size_t nmax = 0;
            unsigned long long max = limits<Element>::max();
            while ( (p1 + p2)*pm1 <= max ){
                p1*=p;
                p2*=p-2;
                nmax++;
            }
            return std::max((size_t)1,nmax);
        }

        /**
         * Specialization for balanced modular representation over double.
         * Computes nmax s.t. (p-1)/2*(((p+1)/2)^{nmax-1}) < 2^53
         * @bib
         * - Dumas Giorgi Pernet 06, arXiv:cs/0601133
         */
        template<class Element>
        inline size_t TRSMBound (const Givaro::ModularBalanced<Element>& F)
        {

            FFLAS_INT_TYPE pi;
            F.characteristic (pi);
            double pp1 = (pi + 1) / 2;
            double pm1 = (pi - 1) / 2;
            double p1 = 1.0;
            size_t nmax = 0;
            double max = limits<Element>::max();
            while (pm1*p1 <= max){
                p1 *= pp1;
                nmax++;
            }
            return std::max((size_t) 1,nmax);
        }

        // /**
        //  * Specialization for balanced modular representation over float
        //  * Computes nmax s.t. (p-1)/2*(((p+1)/2)^{nmax-1}) < 2^24
        //  * See [Dumas Giorgi Pernet 06, arXiv:cs/0601133]
        //  */
        // template<>
        // inline size_t TRSMBound (const Givaro::ModularBalanced<float>& F)
        // {

        // 	FFLAS_INT_TYPE pi;
        // 	F.characteristic (pi);
        // 	unsigned long p = (pi + 1) / 2;
        // 	unsigned long long p1(1);
        // 	size_t nmax = 0;
        // 	unsigned long long max = (1 << (FLT_MANT_DIG + 1)) ;
        // 	while ((pi-1)*p1 < max){
        // 		p1 *= p;
        // 		nmax++;
        // 	}
        // 	return nmax;

        // }
    } // Protected
} // FFLAS

#endif // __FFLASFFPACK_fflas_bounds_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
