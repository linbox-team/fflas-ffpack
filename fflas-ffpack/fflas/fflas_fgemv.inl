/* fflas/fflas_fgemv.inl
 * Copyright (C) 2005 Clement Pernet
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

#ifndef __FFLASFFPACK_fgemv_INL
#define __FFLASFFPACK_fgemv_INL

#include <givaro/zring.h> // DoubleDomain

#if defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS) and defined(__x86_64__)
#include "fflas-ffpack/fflas/fflas_igemm/igemm.h"
#endif

namespace FFLAS{ namespace Protected {
    template <typename FloatElement, class Field>
    inline typename Field::Element_ptr
    fgemv_convert (const Field& F,
                   const FFLAS_TRANSPOSE ta,
                   const size_t M, const size_t N,
                   const typename Field::Element alpha,
                   typename Field::ConstElement_ptr A,const size_t lda,
                   typename Field::ConstElement_ptr X,const size_t incX,
                   const typename Field::Element beta,
                   typename Field::Element_ptr Y, const size_t incY)
    {
        FFLASFFPACK_check(lda);

        Givaro::ModularBalanced<FloatElement> G((FloatElement) F.characteristic());
        FloatElement tmp,alphaf, betaf;
        F.convert (tmp, beta);
        G.init(betaf,tmp);
        F.convert (tmp, alpha);
        G.init(alphaf, tmp);
        size_t ma, na;
        if (ta == FflasTrans) { ma = N; na = M; }
        else { ma = M; na = N; }
        // sizet ldaf = na;
        FloatElement* Af = FFLAS::fflas_new<FloatElement>(M*N);
        FloatElement* Xf = FFLAS::fflas_new<FloatElement>(na);
        FloatElement* Yf = FFLAS::fflas_new<FloatElement>(ma);

        fconvert(F, M, N, Af, N, A, lda);
        freduce (G, M, N, Af, N);
        fconvert(F, na, Xf, 1, X, incX);
        freduce (G, na, Xf, 1);

        if (!F.isZero(beta)){
            fconvert (F, ma, Yf, 1, Y, incY);
            freduce (G, ma, Yf, 1);
        }

        fgemv (G, ta, M, N, alphaf, Af, N, Xf, 1, betaf, Yf, 1);

        finit(F, ma, Yf, 1, Y, incY);
        fflas_delete (Af);
        fflas_delete (Xf);
        fflas_delete (Yf);
        return Y;
    }
}// Protected
}// FFLAS

namespace FFLAS {
    template<class Field>
    inline  typename Field::Element_ptr
    fgemv (const Field& F,
           const FFLAS_TRANSPOSE ta,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr X, const size_t incX,
           const typename Field::Element beta,
           typename Field::Element_ptr Y, const size_t incY,
           MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> > & H)
    {
        if (F.cardinality() < DOUBLE_TO_FLOAT_CROSSOVER && F.cardinality() > 2)
            return Protected::fgemv_convert<float,Field>(F,ta,M,N,alpha,A,lda,X, incX, beta,Y,incY);
        else if (16*F.cardinality() < Givaro::ModularBalanced<double>::maxCardinality() && F.cardinality() > 2)
            return Protected::fgemv_convert<double,Field>(F,ta,M,N,alpha,A,lda,X, incX, beta,Y,incY);
        else {
            FFPACK::failure()(__func__,__LINE__,"Invalid ConvertTo Mode for this field");
        }
        return Y;
    }
}// FFLAS

namespace FFLAS {

    //---------------------------------------------------------------------
    // fgemv: GEneral Matrix Vector Multiplication
    // Computes  Y <- alpha.op(A).X + beta.Y
    // A is M*N,
    //---------------------------------------------------------------------

    template<class Field>
    inline typename Field::Element_ptr
    fgemv (const Field& F, const FFLAS_TRANSPOSE ta,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr X, const size_t incX,
           const typename Field::Element beta,
           typename Field::Element_ptr Y, const size_t incY,
           MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DelayedTag> & H)
    {

        if (!M) {return Y;}
        size_t Ydim = (ta == FflasNoTrans)?M:N;
        size_t Xdim = (ta == FflasNoTrans)?N:M;
        if (!Xdim || F.isZero (alpha)){
            if (F.isZero (beta))
                fzero (F, Ydim, Y, incY);
            else
                fscalin(F, Ydim, beta, Y, incY);
            return Y;
        }

        typename Field::Element alpha_,beta_;
        F.assign (alpha_,alpha);
        F.assign (beta_,beta);
        if (Protected::AreEqual<Field, Givaro::Modular<double> >::value ||
            Protected::AreEqual<Field, Givaro::ModularBalanced<double> >::value){
            //Givaro::Modular<double> need to switch to float if p too small
            if (F.characteristic() < DOUBLE_TO_FLOAT_CROSSOVER && F.cardinality() > 2)
                return Protected::fgemv_convert<float,Field>(F,ta,M,N,alpha,A,lda,X,incX,beta,Y,incY);
        }

        if (Protected::AreEqual<Field, Givaro::Modular<int64_t> >::value ||
            Protected::AreEqual<Field, Givaro::ModularBalanced<int64_t> >::value){
            if (16*F.cardinality() < Givaro::ModularBalanced<double>::maxCardinality() && F.cardinality()>2)
                return Protected::fgemv_convert<double,Field>(F,ta,M,N,alpha,A,lda,X, incX,beta,Y,incY);
            else{
                // Stay over int64_t
                MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::LazyTag, ParSeqHelper::Sequential> HG(H);
                HG.recLevel = 0;
                if (ta == FflasNoTrans)
                    fgemm(F,FflasNoTrans,FflasNoTrans,M,1,N,alpha,A,lda,X,incX,beta,Y,incY,HG);
                else
                    fgemm(F,FflasTrans,FflasNoTrans,N,1,M,alpha,A,lda,X,incX,beta,Y,incY,HG);
                freduce(F,(ta==FflasNoTrans)?M:N, Y,incY);
                H.initOut();
                return Y;
            }
        }
        if ( !F.isOne(alpha) && !F.isMOne(alpha)){
            F.assign (alpha_, F.one);
            F.div (beta_, beta, alpha);
        }
        MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::LazyTag> HD(F,0);

        fgemv (F, ta, M, N, alpha_,
               FFPACK::fflas_const_cast<typename Field::Element_ptr>(A), lda,
               FFPACK::fflas_const_cast<typename Field::Element_ptr>(X), incX,
               beta_, Y, incY, HD);

        Protected::ScalAndReduce (F, Ydim, alpha, Y, incY, HD);
        H.initOut();

        return Y;
    }



}

namespace FFLAS{
    template<class Field>
    inline typename Field::Element_ptr
    fgemv (const Field& F, const FFLAS_TRANSPOSE ta,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr X, const size_t incX,
           const typename Field::Element beta,
           typename Field::Element_ptr Y, const size_t incY,
           MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DefaultTag> & H)

    {
        size_t Ydim = (ta==FflasNoTrans)?M:N;

        if (F.isZero (beta))
            fzero (F, Ydim, Y, incY);
        else {
            typename Field::Element betadivalpha;
            FFLASFFPACK_check(!F.isZero(alpha));
            F.div (betadivalpha, beta, alpha);
            fscalin (F, Ydim, betadivalpha, Y, incY);
        }
        if (ta == FflasNoTrans)
            for (size_t i = 0; i < Ydim; ++i)
                F.addin (Y[i*incY], fdot(F, N, A+i*lda, 1, X, incX));
        else
            for (size_t i = 0; i < Ydim; ++i)
                F.addin (Y[i*incY], fdot(F, M, A+i, lda, X, incX));
        fscalin (F, Ydim, alpha, Y, incY);

        return Y;
    }
}

namespace FFLAS{
    template<class Field>
    inline typename Field::Element_ptr
    fgemv (const Field& F, const FFLAS_TRANSPOSE ta,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr X, const size_t incX,
           const typename Field::Element beta,
           typename Field::Element_ptr Y, const size_t incY,
           MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::LazyTag> & H)
    {
        typedef MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::LazyTag> HelperType;
        typedef typename HelperType::DelayedField::Element DFElt;
        typedef typename HelperType::DelayedField::Element_ptr DFElt_ptr;
        typedef typename HelperType::DelayedField::ConstElement_ptr DFCElt_ptr;
        DFElt alphadf=alpha, betadf=beta;
        size_t Ydim = (ta==FflasNoTrans)?M:N;
        size_t Xdim = (ta==FflasNoTrans)?N:M;
        if (F.isMOne (alpha)) alphadf =  -F.one;
        else {
            alphadf = F.one;
            if (! F.isOne( alpha)) {
                // Compute y = A*x + beta/alpha.y, then y *= alpha
                FFLASFFPACK_check(!F.isZero(alpha));
                typename Field::Element betadalpha;
                F.init(betadalpha);
                F.div (betadalpha, beta, alpha);
                betadf=betadalpha;
            }
        }
        if (F.isMOne(betadf)) betadf = -F.one;

        size_t kmax = H.MaxDelayedDim (betadf);

        if (kmax <=  Xdim/2 ){
            // Might as well reduce inputs
            if (H.Amin < H.FieldMin || H.Amax>H.FieldMax){
                H.initA();
                freduce_constoverride (F, M, N, A, lda);
            }
            if (H.Bmin < H.FieldMin || H.Bmax>H.FieldMax){
                H.initB();
                freduce_constoverride (F, Xdim, X, incX);
            }
            if (H.Cmin < H.FieldMin || H.Cmax>H.FieldMax){
                H.initC();
                freduce (F, Ydim, Y, incY);
            }
            kmax = H.MaxDelayedDim (betadf);
        }

        if (!kmax){
            MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DefaultTag> HG(H);
            H.initOut();
            return fgemv (F, ta, M, N, alpha, A, lda, X, incX, beta, Y, incY, HG);
        }
        size_t k2 = std::min (Xdim, kmax);
        size_t nblock = Xdim / kmax;
        size_t remblock = Xdim % kmax;
        if (!remblock) {
            remblock = kmax;
            --nblock;
        }
        size_t shiftA, M1, N1, Mi, Ni;
        if (ta == FflasTrans) {
            shiftA = k2*lda;
            M1 = remblock;
            Mi = k2;
            Ni = N1 = N;
        }else {
            shiftA = k2;
            Mi = M1 = M;
            N1 = remblock;
            Ni = k2;
        }
        MMHelper<typename associatedDelayedField<const Field>::field, MMHelperAlgo::Classic, ModeCategories::DefaultBoundedTag> Hfp(H);


        fgemv (H.delayedField, ta, M1, N1, alphadf, (DFCElt_ptr)A+nblock*shiftA, lda,
               (DFCElt_ptr)X+nblock*k2*incX, incX, betadf, (DFElt_ptr)Y, incY, Hfp);

        for (size_t i = 0; i < nblock; ++i) {
            freduce (F, Ydim ,Y, incY);
            Hfp.initC();
            fgemv (H.delayedField, ta, Mi, Ni, alphadf, (DFCElt_ptr)A+i*shiftA, lda,
                   (DFCElt_ptr)X+i*k2*incX, incX, F.one, (DFElt_ptr)Y, incY, Hfp);
        }

        if (!F.isOne(alpha) && !F.isMOne(alpha)){
            DFElt al; F.convert(al, alpha);
            if (al<0) al = -al;
            if (std::max(-Hfp.Outmin, Hfp.Outmax) > Hfp.MaxStorableValue/al){
                freduce (F, Ydim, Y, incY);
                Hfp.initOut();
            }
            fscalin (H.delayedField, Ydim, alpha, (DFElt_ptr)Y, incY);
            if (alpha>0){
                H.Outmin = al*Hfp.Outmin;
                H.Outmax = al*Hfp.Outmax;
            } else {
                H.Outmin = -al*Hfp.Outmax;
                H.Outmax = -al*Hfp.Outmin;
            }
        }else {
            H.Outmin = Hfp.Outmin;
            H.Outmax = Hfp.Outmax;
        }
        return Y;
    }
}

namespace FFLAS{
    template<class Field>
    inline typename Field::Element_ptr
    fgemv (const Field& F, const FFLAS_TRANSPOSE ta,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           typename Field::ConstElement_ptr X, const size_t incX,
           const typename Field::Element beta,
           typename Field::Element_ptr Y, const size_t incY)
    {
        if (!M) {return Y;}
        size_t Ydim = (ta == FflasNoTrans)?M:N;
        size_t Xdim = (ta == FflasNoTrans)?N:M;
        if (!Xdim || F.isZero (alpha)){
            if (F.isZero (beta))
                fzero (F, Ydim, Y, incY);
            else
                fscalin(F, Ydim, beta, Y, incY);
            return Y;
        }
        MMHelper<Field, MMHelperAlgo::Classic > HW (F, 0);
        return 	fgemv (F, ta, M, N, alpha,
                       FFPACK::fflas_const_cast<typename Field::Element_ptr>(A), lda,
                       FFPACK::fflas_const_cast<typename Field::Element_ptr>(X), incX,
                       beta, Y, incY, HW);
    }
}


namespace FFLAS{
    inline Givaro::ZRing<int64_t>::Element_ptr
    fgemv (const Givaro::ZRing<int64_t>& F, const FFLAS_TRANSPOSE ta,
           const size_t M, const size_t N,
           const int64_t alpha,
           const int64_t* A, const size_t lda,
           const int64_t* X, const size_t incX,
           const int64_t beta,
           int64_t* Y, const size_t incY,
           MMHelper<Givaro::ZRing<int64_t>, MMHelperAlgo::Classic, ModeCategories::DefaultTag> & H)
    {
        FFLASFFPACK_check(lda);

#if defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS) and defined(__x86_64__)
        if (ta == FflasNoTrans)
            igemm_ (FflasRowMajor, ta, FflasNoTrans,M,1,N,alpha,A,lda,X,incX,beta,Y,incY);
        else
            igemm_ (FflasRowMajor, ta, FflasNoTrans,N,1,M,alpha,A,lda,X,incX,beta,Y,incY);
#else
        if (ta == FflasNoTrans){
            int64_t* Yi=Y;

            for (size_t i=0;i<M;i++, Yi+=incY){
                *Yi *= beta * *Yi;
                const int64_t* Xj=X;
                for (size_t j=0; j < N; j++, Xj += incX)
                    *Yi += alpha*A[i*lda+j] * *Xj;
            }
        } else {
            int64_t* Yi=Y;

            for (size_t i=0;i<N;i++, Yi+=incY){
                *Yi *= beta * *Yi;
                const int64_t* Xj=X;
                for (size_t j=0; j < M; j++, Xj += incX)
                    *Yi += alpha*A[i+j*lda] * *Xj;
            }
        }
#endif
        return Y;
    }
    inline Givaro::DoubleDomain::Element_ptr
    fgemv (const Givaro::DoubleDomain& F, const FFLAS_TRANSPOSE ta,
           const size_t M, const size_t N,
           const Givaro::DoubleDomain::Element alpha,
           const Givaro::DoubleDomain::ConstElement_ptr A, const size_t lda,
           const Givaro::DoubleDomain::ConstElement_ptr X, const size_t incX,
           const Givaro::DoubleDomain::Element beta,
           Givaro::DoubleDomain::Element_ptr Y, const size_t incY,
           MMHelper<Givaro::DoubleDomain, MMHelperAlgo::Classic, ModeCategories::DefaultTag> & H)
    {
        FFLASFFPACK_check(lda);

#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_dgemv (CblasRowMajor, (CBLAS_TRANSPOSE) ta,
                     (int)M, (int)N, (Givaro::DoubleDomain::Element) alpha,
                     A, (int)lda, X, (int)incX, (Givaro::DoubleDomain::Element) beta, Y, (int)incY);
        return Y;
    }

    template <class Field>
    inline typename Field::Element_ptr
    fgemv (const Field& F, const FFLAS_TRANSPOSE ta,
           const size_t M, const size_t N,
           const typename Field::Element alpha,
           const typename Field::ConstElement_ptr A, const size_t lda,
           const typename Field::ConstElement_ptr X, const size_t incX,
           const typename Field::Element beta,
           typename Field::Element_ptr Y, const size_t incY,
           MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DefaultBoundedTag> & H)
    {
        H.setOutBounds((ta ==FflasNoTrans)?N:M, alpha, beta);
        MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DefaultTag> Hb(F,0);

        return fgemv(F, ta, M, N, alpha, A, lda, X, incX, beta, Y, incY, Hb);
    }

    inline Givaro::FloatDomain::Element_ptr
    fgemv (const Givaro::FloatDomain& F, const FFLAS_TRANSPOSE ta,
           const size_t M, const size_t N,
           const Givaro::FloatDomain::Element alpha,
           const Givaro::FloatDomain::ConstElement_ptr A, const size_t lda,
           const Givaro::FloatDomain::ConstElement_ptr X, const size_t incX,
           const Givaro::FloatDomain::Element beta,
           Givaro::FloatDomain::Element_ptr Y, const size_t incY,
           MMHelper<Givaro::FloatDomain, MMHelperAlgo::Classic, ModeCategories::DefaultTag> & H)
    {
        FFLASFFPACK_check(lda);

#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_sgemv (CblasRowMajor, (CBLAS_TRANSPOSE) ta,
                     (int)M, (int)N, (Givaro::FloatDomain::Element) alpha,
                     A, (int)lda, X, (int)incX, (Givaro::FloatDomain::Element) beta, Y, (int)incY);
        return Y;
    }

    template<class Field, class Cut, class Param>
    typename Field::Element_ptr
    fgemv(const Field& F,
           const FFLAS_TRANSPOSE ta,
           const size_t m,
           const size_t n,
           const typename Field::Element alpha,
           const typename Field::ConstElement_ptr A, const size_t lda,
           const typename Field::ConstElement_ptr X, const size_t incX,
           const typename Field::Element beta,
           typename Field::Element_ptr Y, const size_t incY,
           ParSeqHelper::Parallel<Cut,Param>& parH){
        MMHelper<Field, MMHelperAlgo::Auto, typename FFLAS::ModeTraits<Field>::value, ParSeqHelper::Parallel<Cut,Param> > pH (F,m,n,1,parH);
        return fgemv(F, ta, m, n, alpha, A, lda, X, incX, beta, Y, incY, pH);
    }

    template<class Field>
    typename Field::Element_ptr
    fgemv(const Field& F,
           const FFLAS_TRANSPOSE ta,
           const size_t m,
           const size_t n,
           const typename Field::Element alpha,
           const typename Field::ConstElement_ptr A, const size_t lda,
           const typename Field::ConstElement_ptr X, const size_t incX,
           const typename Field::Element beta,
           typename Field::Element_ptr Y, const size_t incY,
           ParSeqHelper::Sequential& seqH ){
        MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DefaultTag> pH(F,m,n,1,seqH);
        return fgemv(F, ta, m, n, alpha, A, lda, X, incX, beta, Y, incY, pH);
    }
}

#endif //  __FFLASFFPACK_fgemv_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
