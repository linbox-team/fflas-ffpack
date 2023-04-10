/* fflas/fflas_fgemm.inl
 * Copyright (C) 2005 Clement Pernet
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by Clement Pernet  < Clement.Pernet@imag.fr >
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

#ifndef __FFLASFFPACK_fgemm_INL
#define __FFLASFFPACK_fgemm_INL

#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include "fflas-ffpack/utils/debug.h"

namespace FFLAS { namespace Protected{

    template <class NewField, class Field, class FieldMode>
    inline typename Field::Element_ptr
    fgemm_convert (const Field& F,
                   const FFLAS_TRANSPOSE ta,
                   const FFLAS_TRANSPOSE tb,
                   const size_t m, const size_t n, const size_t k,
                   const typename Field::Element alpha,
                   typename Field::ConstElement_ptr A,const size_t lda,
                   typename Field::ConstElement_ptr B,const size_t ldb,
                   const typename Field::Element beta,
                   typename Field::Element_ptr C, const size_t ldc,
                   MMHelper<Field, MMHelperAlgo::Winograd, FieldMode> & H)
    {
        // CP: lda, ldb, ldc can be zero (if m,n or k is 0) and since  this may have not
        // been checked by the caller at this point.
        // FFLASFFPACK_check(lda);
        // FFLASFFPACK_check(ldb);
        // FFLASFFPACK_check(ldc);
        typedef typename NewField::Element FloatElement;
        NewField G((FloatElement) F.characteristic());
        FloatElement tmp,alphaf, betaf;
        // This conversion is quite tricky, but convert and init are required
        // in sequence e.g. for when F is a ModularBalanced field and alpha == -1
        F.convert (tmp, beta);
        G.init(betaf, tmp);
        F.convert (tmp, alpha);
        G.init(alphaf, tmp);

        FloatElement* Af = FFLAS::fflas_new(G, m, k);
        FloatElement* Bf = FFLAS::fflas_new(G, k, n);
        FloatElement* Cf = FFLAS::fflas_new(G, m, n);

        size_t ma, ka, kb, nb; //mb, na
        if (ta == FflasTrans) { ma = k; ka = m; }
        else { ma = m; ka = k; }
        if (tb == FflasTrans) { kb = n; nb = k; }
        else {  kb = k; nb = n; }
        size_t ldaf = ka, ldbf = nb, ldcf= n;

        fconvert(F, ma, ka, Af, ka, A, lda);
        freduce(G, ma, ka, Af, ka);
        fconvert(F, kb, nb, Bf, nb, B, ldb);
        freduce(G, kb, nb, Bf, nb);

        if (!F.isZero(beta)){
            fconvert(F, m, n, Cf, n, C, ldc);
            freduce (G, m, n, Cf, n);
        }
        MMHelper<NewField, MMHelperAlgo::Winograd> HG(G,H.recLevel, ParSeqHelper::Sequential());
        fgemm (G, ta, tb, m, n, k, alphaf, Af, ldaf, Bf, ldbf, betaf, Cf, ldcf, HG);

        finit (F, m, n, Cf, n, C, ldc);

        fflas_delete (Af);
        fflas_delete (Bf);
        fflas_delete (Cf);
        return C;
    }
}//Protected
}//FFLAS

namespace FFLAS{ namespace Protected{
    template <class Field, class Element, class AlgoT, class ParSeqTrait>
    inline bool NeedPreAddReduction (Element& Outmin, Element& Outmax,
                                     Element& Op1min, Element& Op1max,
                                     Element& Op2min, Element& Op2max,
                                     MMHelper<Field, AlgoT, ModeCategories::LazyTag, ParSeqTrait >& WH)
    {
        Outmin = Op1min + Op2min;
        Outmax = Op1max + Op2max;
        if (WH.MaxStorableValue - Op1max < Op2max ||
            WH.MaxStorableValue + Op1min < -Op2min){
            // Reducing both Op1 and Op2
            Op1min = Op2min = WH.FieldMin;
            Op1max = Op2max = WH.FieldMax;
            Outmin = 2*WH.FieldMin;
            Outmax = 2*WH.FieldMax;
            return true;
        } else return false;
    }

    template <class Field, class Element, class AlgoT, class ModeT, class ParSeqTrait>
    inline bool NeedPreAddReduction (Element& Outmin, Element& Outmax,
                                     Element& Op1min, Element& Op1max,
                                     Element& Op2min, Element& Op2max,
                                     MMHelper<Field, AlgoT, ModeT, ParSeqTrait >& WH)
    {
        Outmin = WH.FieldMin;
        Outmax = WH.FieldMax;
        return false;
    }

    template <class Field, class Element, class AlgoT, class ParSeqTrait>
    inline bool NeedPreSubReduction (Element& Outmin, Element& Outmax,
                                     Element& Op1min, Element& Op1max,
                                     Element& Op2min, Element& Op2max,
                                     MMHelper<Field, AlgoT, ModeCategories::LazyTag, ParSeqTrait >& WH)
    {
        Outmin = Op1min - Op2max;
        Outmax = Op1max - Op2min;
        if (WH.MaxStorableValue - Op1max < -Op2min ||
            WH.MaxStorableValue - Op2max < -Op1min){
            // Reducing both Op1 and Op2
            Op1min = Op2min = WH.FieldMin;
            Op1max = Op2max = WH.FieldMax;
            Outmin = WH.FieldMin-WH.FieldMax;
            Outmax = -Outmin;
            return true;
        } else return false;
    }

    template <class Field, class Element, class AlgoT, class ModeT, class ParSeqTrait>
    inline bool NeedPreSubReduction (Element& Outmin, Element& Outmax,
                                     Element& Op1min, Element& Op1max,
                                     Element& Op2min, Element& Op2max,
                                     MMHelper<Field, AlgoT, ModeT, ParSeqTrait >& WH)
    {
        // Necessary? -> CP: Yes, for generic Mode of op
        Outmin = WH.FieldMin;
        Outmax = WH.FieldMax;
        return false;
    }

    //Probable bug here due to overflow of int64_t
    template<class Field, class Element, class AlgoT, class ParSeqTrait>
    inline bool NeedDoublePreAddReduction (Element& Outmin, Element& Outmax,
                                           Element& Op1min, Element& Op1max,
                                           Element& Op2min, Element& Op2max, Element beta,
                                           MMHelper<Field, AlgoT, ModeCategories::LazyTag, ParSeqTrait >& WH)
    {
        // Testing if P5 need to be reduced
        Outmin =  std::min(beta*Op2min,beta*Op2max);
        Outmax =  std::max(beta*Op2min,beta*Op2max);
        if (Op1max > WH.MaxStorableValue-Outmax ||
            -Op1min > WH.MaxStorableValue+Outmin){
            Outmin += WH.FieldMin;
            Outmax += WH.FieldMax;
            return true;
        } else{
            Outmin += Op1min;
            Outmax += Op1max;
            return false;
        }
    }

    template<class Field, class Element, class AlgoT, class ModeT, class ParSeqTrait>
    inline bool NeedDoublePreAddReduction (Element& Outmin, Element& Outmax,
                                           Element& Op1min, Element& Op1max,
                                           Element& Op2min, Element& Op2max, Element beta,
                                           MMHelper<Field, AlgoT, ModeT, ParSeqTrait>& WH)
    {
        Outmin = WH.FieldMin;
        Outmax = WH.FieldMax;
        return false;
    }

    template <class Field, class AlgoT, class ParSeqTrait>
    inline void ScalAndReduce (const Field& F, const size_t N,
                               const typename Field::Element alpha,
                               typename Field::Element_ptr X, const size_t incX,
                               const MMHelper<Field, AlgoT, ModeCategories::LazyTag, ParSeqTrait >& H)
    {
        if (!F.isOne(alpha) && !F.isMOne(alpha)){
            typename MMHelper<Field, AlgoT, ModeCategories::LazyTag, ParSeqTrait >::DFElt al;
            F.convert(al, alpha);
            if (al < 0) al = -al;
            if (std::max(-H.Outmin, H.Outmax) > H.MaxStorableValue/al){
                freduce (F, N, X, incX);
                fscalin (F, N, alpha, X, incX);
            } else {
                fscalin (H.delayedField, N, alpha, X, incX);
                freduce (F, N, X, incX);
            }
        } else
            freduce (F, N, X, incX);
    }

    template <class Field, class AlgoT, class ParSeqTrait>
    inline void ScalAndReduce (const Field& F, const size_t M, const size_t N,
                               const typename Field::Element alpha,
                               typename Field::Element_ptr A, const size_t lda,
                               const MMHelper<Field, AlgoT, ModeCategories::LazyTag, ParSeqTrait >& H)
    {
        if (!F.isOne(alpha) && !F.isMOne(alpha)){
            typename MMHelper<Field, AlgoT, ModeCategories::LazyTag, ParSeqTrait >::DFElt al;
            F.convert(al, alpha);
            if (al<0) al = -al;
            if (std::max(-H.Outmin, H.Outmax) > H.MaxStorableValue/al){
                freduce (F, M, N, A, lda);
                fscalin (F, M, N, alpha, A, lda);
            } else {
                fscalin (H.delayedField, M, N, alpha, (typename MMHelper<Field, AlgoT, ModeCategories::LazyTag, ParSeqTrait >::DFElt*)A, lda);
                freduce (F, M, N, A, lda);
            }
        } else
            freduce (F, M, N, A, lda);
    }

} // Protected
} // FFLAS

namespace FFLAS {

    template<class Field>
    inline  typename Field::Element_ptr
    fgemm (const Field& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m, const size_t n, const size_t k,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr B, const size_t ldb,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           MMHelper<Field, MMHelperAlgo::Winograd, ModeCategories::ConvertTo<ElementCategories::MachineFloatTag>, ParSeqHelper::Sequential> & H)
    {
        if (!std::is_same<Field,Givaro::Modular<float> >::value){
            if (F.cardinality() == 2)
                return Protected::fgemm_convert<Givaro::Modular<float>,Field>(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H);
            else if (!std::is_same<Field,Givaro::ModularBalanced<float> >::value){
                if (F.cardinality() < DOUBLE_TO_FLOAT_CROSSOVER)
                    return Protected::fgemm_convert<Givaro::ModularBalanced<float>,Field>(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H);
                else if (!std::is_same<Field,Givaro::ModularBalanced<double> >::value &&  16*F.cardinality() < Givaro::ModularBalanced<double>::maxCardinality())
                    return Protected::fgemm_convert<Givaro::ModularBalanced<double>,Field>(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H);
            }
        }
        // else if (Protected::AreEqual<typename Field::Element,int64_t>::value) {
        // 	    // Stays over int64_t
        // 	MMHelper<Field, MMHelperAlgo::Winograd, ModeCategories::DelayedTag, ParSeqHelper::Sequential> HG(H);
        // 	H.Outmin=HG.Outmin;
        // 	H.Outmax=HG.Outmax;
        // 	return fgemm(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,HG);

        //	}
        else {
            // Fall back case: used
            FFPACK::failure()(__func__,__LINE__,"Invalid ConvertTo Mode for this field");
        }
        return C;
    }
}// FFLAS


// fgemm
namespace FFLAS {

    template<typename Field>
    inline typename Field::Element_ptr
    fgemm( const Field& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr B, const size_t ldb,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           const ParSeqHelper::Sequential seq)
    {
        MMHelper<Field, MMHelperAlgo::Auto, typename FFLAS::ModeTraits<Field>::value, ParSeqHelper::Sequential > HW (F, m, k, n, seq);
        return 	fgemm (F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, HW);
    }

    template<typename Field,class Cut,class Param>
    inline typename Field::Element_ptr
    fgemm( const Field& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr B, const size_t ldb,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           const ParSeqHelper::Parallel<Cut,Param> par)
    {
        MMHelper<Field, MMHelperAlgo::Auto, typename FFLAS::ModeTraits<Field>::value, ParSeqHelper::Parallel<Cut,Param> > HW (F, m, k, n, par);
        return 	fgemm (F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, HW);
    }

    template<typename Field>
    inline typename Field::Element_ptr
    fgemm( const Field& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr B, const size_t ldb,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc)
    {
        if (!m || !n) {return C;}

        if (!k || F.isZero (alpha)){
            fscalin(F, m, n, beta, C, ldc);
            return C;
        }
        Checker_fgemm<Field> checker(F,m,n,k,beta,C,ldc);
        fgemm(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,FFLAS::ParSeqHelper::Sequential());
        checker.check(ta,tb,alpha,A,lda,B,ldb,C);
        return C;
    }

    template<typename Field, class ModeT, class ParSeq>
    inline typename Field::Element_ptr
    fgemm( const Field& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr B, const size_t ldb,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           MMHelper<Field, MMHelperAlgo::Auto, ModeT, ParSeq> & H)
    {
        MMHelper<Field, typename AlgoChooser<ModeT, ParSeq>::value, ModeT, ParSeq> HW (H);
        return fgemm(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,HW);
    }

    template<class Field>
    inline typename Field::Element_ptr
    fgemm( const Field& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr B, const size_t ldb,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           MMHelper<Field, MMHelperAlgo::Winograd, ModeCategories::DelayedTag, ParSeqHelper::Sequential> & H)
    {
        if (!m || !n) {return C;}

        if (!k || F.isZero (alpha)){
            fscalin(F, m, n, beta, C, ldc);
            return C;
        }

#ifndef NDEBUG
        /*  check if alpha is invertible.
         *  XXX do it in F.isInvertible(Element&) ?
         *  XXX do it in return status of F.inv(Element&,Element&)
         */
        typename Field::Element e ;
        F.assign(e,beta);
        F.divin(e,alpha);
        F.mulin(e,alpha);
        FFLASFFPACK_check(F.areEqual(e,beta));
#endif

#if 0
        // detect fgemv
        if (n == 1 and ...) {}
        // detect fger
        if (k==1 and ...) {}
#endif
        if (!std::is_same<Field,Givaro::Modular<float> >::value){
            if (F.cardinality() == 2)
                return Protected::fgemm_convert<Givaro::Modular<float>,Field>(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H);
            else if (!std::is_same<Field,Givaro::ModularBalanced<float> >::value){
                if (F.characteristic() < DOUBLE_TO_FLOAT_CROSSOVER)
                    return Protected::fgemm_convert<Givaro::ModularBalanced<float>,Field>(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H);
                else if (!std::is_same<Field,Givaro::ModularBalanced<double> >::value &&
                         !std::is_same<Field,Givaro::Modular<double> >::value &&
                         16*F.cardinality() < Givaro::ModularBalanced<double>::maxCardinality())
                    return Protected::fgemm_convert<Givaro::ModularBalanced<double>,Field>(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H);
            }
        }
        typename Field::Element alpha_,beta_;
        if ( !F.isOne(alpha) && !F.isMOne(alpha)){
            F.assign (alpha_, F.one);
            F.div (beta_, beta, alpha);
        } else {
            F.assign (alpha_,alpha);
            F.assign (beta_,beta);
        }
        MMHelper<Field, MMHelperAlgo::Winograd, ModeCategories::LazyTag>  HD(H);
        // std::cerr<<"\n Delayed -> Lazy alpha_ = "<<alpha_<<std::endl;
        // std::cerr<<" A = "<<*A<<"\n B = "<<*B<<"\n C = "<<*C<<"\n alpha, beta ="<<alpha<<" "<<beta<<std::endl;

        fgemm (F, ta, tb, m, n, k, alpha_, A, lda, B, ldb, beta_, C, ldc, HD);

        Protected::ScalAndReduce (F, m, n, alpha, C, ldc, HD);
        
        H.initOut();

        return C;
    }
} // FFLAS

// #include "fflas_fgemm/matmul_algos.inl"
#include "fflas_fgemm/fgemm_classical.inl"
#include "fflas_fgemm/fgemm_winograd.inl"
// #include "fflas_fgemm/gemm_bini.inl"

// fsquare
namespace FFLAS {
    template < class Field >
    inline typename Field::Element_ptr
    fsquare (const Field& F,
             const FFLAS_TRANSPOSE ta,
             const size_t n, const typename Field::Element alpha,
             typename Field::ConstElement_ptr A, const size_t lda,
             const typename Field::Element beta,
             typename Field::Element_ptr C, const size_t ldc)
    {

        double alphad, betad;
        F.convert (alphad, alpha);
        if (F.isMOne (beta))
            betad = -1.0;
        else
            F.convert (betad, beta);

        //! @bug why double ?
        // Double  matrices initialisation
        Givaro::DoubleDomain::Element_ptr Ad = fflas_new (Givaro::DoubleDomain(),n,n);
        Givaro::DoubleDomain::Element_ptr Cd = fflas_new (Givaro::DoubleDomain(),n,n);
        // Conversion finite Field = >  double
        fconvert (F, n, n, Ad, n, A, lda);
        if (!F.isZero(beta)) fconvert(F, n, n, Cd, n, C, ldc);

        // Call to the blas Multiplication
        FFLASFFPACK_check(n);
#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_dgemm (CblasRowMajor, (CBLAS_TRANSPOSE)ta,
                     (CBLAS_TRANSPOSE)ta, (int)n, (int)n, (int)n,
                     (Givaro::DoubleDomain::Element) alphad, Ad, (int)n, Ad, (int)n,
                     (Givaro::DoubleDomain::Element) betad, Cd, (int)n);
        // Conversion double = >  Finite Field
        fflas_delete (Ad);
        finit (F,n,n, Cd, n, C, ldc);
        fflas_delete (Cd);
        return C;
    }

    namespace Protected {

        // F is Modular(Balanced)<float/double>
        template < class Field >
        inline typename Field::Element_ptr
        fsquareCommon (const Field& F,
                       const FFLAS_TRANSPOSE ta,
                       const size_t n, const typename Field::Element alpha,
                       typename Field::ConstElement_ptr A, const size_t lda,
                       const typename Field::Element beta,
                       typename Field::Element_ptr C, const size_t ldc)
        {
            if (C==A) {
                typename Field::Element_ptr Ad = fflas_new (F, n, n);
                fassign(F,n,n,A,lda,Ad,n);
                fgemm (F, ta, ta, n, n, n, alpha, Ad, n, Ad, n, beta, C, ldc);
                fflas_delete (Ad);
            }
            else
                fgemm (F, ta, ta, n, n, n, alpha, A, lda, A, lda, beta, C, ldc);
            return C;

        }

    } // Protected

    template <>
    inline double* fsquare (const  Givaro::ModularBalanced<double> & F,
                            const FFLAS_TRANSPOSE ta,
                            const size_t n, const double alpha,
                            const double* A, const size_t lda,
                            const double beta,
                            double* C, const size_t ldc)
    {
        return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
    }

    template <>
    inline float * fsquare (const  Givaro::ModularBalanced<float> & F,
                            const FFLAS_TRANSPOSE ta,
                            const size_t n, const float alpha,
                            const float* A, const size_t lda,
                            const float beta,
                            float* C, const size_t ldc)
    {
        return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
    }

    template <>
    inline double* fsquare (const  Givaro::Modular<double> & F,
                            const FFLAS_TRANSPOSE ta,
                            const size_t n, const double alpha,
                            const double* A, const size_t lda,
                            const double beta,
                            double* C, const size_t ldc)
    {
        return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
    }

    template <>
    inline float * fsquare (const  Givaro::Modular<float> & F,
                            const FFLAS_TRANSPOSE ta,
                            const size_t n, const float alpha,
                            const float* A, const size_t lda,
                            const float beta,
                            float* C, const size_t ldc)
    {
        return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
    }

} // FFLAS

#endif // __FFLASFFPACK_fgemm_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
