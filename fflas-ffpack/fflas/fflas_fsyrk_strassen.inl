/*
 * Copyright (C) 2019 the FFLAS-FFPACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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

#ifndef __FFLASFFPACK_fflas_fsyrk_strassen_INL
#define __FFLASFFPACK_fflas_fsyrk_strassen_INL

#include <givaro/givintsqrootmod.h>

namespace FFLAS{ namespace Protected{
    template <class Field, class Element, class AlgoT, class ParSeqTrait>
    inline bool NeedPreScalReduction (Element& Outmin, Element& Outmax,
                                      Element& Op1min, Element& Op1max,
                                      const Element& x,
                                      MMHelper<Field, AlgoT, ModeCategories::LazyTag, ParSeqTrait >& WH)
    {
        if (x<0){
            if ( std::max(-Op1min, Op1max) > WH.MaxStorableValue / (-x)){
                Outmin = x * WH.FieldMax;
                Outmax = x * WH.FieldMin;
                Op1min = WH.FieldMin;
                Op1max = WH.FieldMax;
                return true;
            }else{
                Outmin = x * Op1max;
                Outmax = x * Op1min;
                return false;
            }
        } else {
            if ( std::max(-Op1min, Op1max) > WH.MaxStorableValue / x){
                Outmin = x * WH.FieldMin;
                Outmax = x * WH.FieldMax;
                Op1min = WH.FieldMin;
                Op1max = WH.FieldMax;
                return true;
            }else{
                Outmin = x * Op1min;
                Outmax = x * Op1max;
                return false;
            }
        }
    }

    template <class Field, class Element, class AlgoT, class ModeT, class ParSeqTrait>
    inline bool NeedPreScalReduction (Element& Outmin, Element& Outmax,
                                      Element& Op1min, Element& Op1max,
                                      const Element& x,
                                      MMHelper<Field, AlgoT, ModeT, ParSeqTrait >& WH)
    {
        Outmin = WH.FieldMin;
        Outmax = WH.FieldMax;
        return false;
    }
        
    template <class Field, class Element, class AlgoT, class ParSeqTrait>
    inline bool NeedPreAxpyReduction (Element& Outmin, Element& Outmax,
                                      Element& Op1min, Element& Op1max,
                                      Element& Op2min, Element& Op2max,
                                      const Element& x,
                                      MMHelper<Field, AlgoT, ModeCategories::LazyTag, ParSeqTrait >& WH)
    // Out <- Op1 + x.Op2
    // x is assumed to be reduced
    {
        if (x<0){
            if (Op2min < (WH.MaxStorableValue - Op1max) / x ||
                Op2max > (-WH.MaxStorableValue - Op1min) / x ){
                Outmin = WH.FieldMin + x * WH.FieldMax ;
                Outmax = WH.FieldMax + x * WH.FieldMin;
                Op1min = WH.FieldMin;
                Op1max = WH.FieldMax;
                Op2min = WH.FieldMin;
                Op2max = WH.FieldMax;
                return true;
            }else{
                Outmin = Op1min + x * Op2max;
                Outmax = Op1max + x * Op2min;
                return false;
            }
        } else {
            if (Op2max > (WH.MaxStorableValue - Op1max) / x ||
                Op2min < (-WH.MaxStorableValue - Op1min) / x ){
                Outmin = (x+1) * WH.FieldMin ;
                Outmax = (x+1) * WH.FieldMax;
                Op1min = WH.FieldMin;
                Op1max = WH.FieldMax;
                Op2min = WH.FieldMin;
                Op2max = WH.FieldMax;
                return true;
            }else{
                Outmin = Op1min + x * Op2min;
                Outmax = Op1max + x * Op2max;
                return false;
            }
        }
    }

    template <class Field, class Element, class AlgoT, class ModeT, class ParSeqTrait>
    inline bool NeedPreAxpyReduction (Element& Outmin, Element& Outmax,
                                      Element& Op1min, Element& Op1max,
                                      Element& Op2min, Element& Op2max,
                                      const Element& x,
                                      MMHelper<Field, AlgoT, ModeT, ParSeqTrait >& WH)
    {
        Outmin = WH.FieldMin;
        Outmax = WH.FieldMax;
        return false;
    }
    }
}
namespace FFLAS {

    template<class Field, class FieldTrait>
    inline void
    computeS1S2 (const Field& F,
                 const FFLAS_TRANSPOSE trans,
                 const size_t N,
                 const size_t K,
                 const typename Field::Element x,
                 const typename Field::Element y,
                 typename Field::Element_ptr A, const size_t lda,
                 typename Field::Element_ptr S, const size_t lds,
                 typename Field::Element_ptr T, const size_t ldt,
                 MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait>& WH)  {
            // WH stores and maintain bounds on A11 in Amin,Amax, A21 in Bmin, Bmax and A22, in Cmin,Cmax
            // Bounds on S and T are stored in Outmin and Outmax

            // Computes (when trans = NoTrans)
            // S = (A21-A11) x Y in S
            // T =  A22 - A21 x Y  in T
            // where Y = [ x.I  y.I]
            //           [ -y.I x.I]
        typedef MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > MMH_t;
        const typename MMH_t::DelayedField & DF = WH.delayedField;
        typedef typename  MMH_t::DelayedField::Element DFElt;

        size_t N2 = N>>1;
        size_t K2 = K>>1;
        size_t K4 = K2>>1;
        size_t Axxrows, Axxcols;
        typename Field::Element_ptr A11 = A, A21, A22, A11r, A21r, A22r;
        typename Field::Element_ptr Sr, Tr;
        if (trans==FflasNoTrans){
            A21 = A + N2*lda;
            A22 = A21 + K2;
            A22r = A22 + K4;
            A11r = A11 + K4;
            A21r = A21 + K4;
            Sr = S + K4;
            Tr = T+K4;
            Axxrows=N2;
            Axxcols=K2;
        } else { // Trans
            A21 = A + N2;
            A22 = A21 + K2*lda;
            A22r = A22 + K4*lda;
            A11r = A11 + K4*lda;
            A21r = A21 + K4*lda;
            Sr = S + K4*lds;
            Tr = T + K4*ldt;
            Axxrows=K2;
            Axxcols=N2;
        }

        typename Field::Element negx, negy;
        F.init(negx);
        F.neg(negx, x);
        F.init(negy);
        F.neg(negy, y); 


        // Operations:
        // S1 <- A21 x x.In
        // S2 <- S1 + A21 x y [ 0    In/2 ]
        //                    [ -In/2  0   ]
        // T <- A22 - S2
        // S3 <- S2 - A11 x x.In
        // S  <- S3 - A11 x y. [ 0    In/2 ]
        //                     [ -In/2  0   ]

        DFElt S1min, S1max, S2min, S2max, S2minr, S2maxr, S3min, S3max, Smin, Smax, Sminr, Smaxr, Tmin, Tmax;
            // Should we reduce Axx and Sx beforehand?
        bool reduceA21 = false;
        if (Protected::NeedPreScalReduction (S1min, S1max, WH.Bmin, WH.Bmax, x, WH)){
            reduceA21 = true;
        }

        bool reduceS1 = false;
        if (y){
            bool redS1 = Protected::NeedPreAxpyReduction (S2min, S2max, S1min, S1max, WH.Bmin, WH.Bmax, y, WH);
            bool redS1r = Protected::NeedPreAxpyReduction (S2minr, S2maxr, S1min, S1max, WH.Bmin, WH.Bmax, negy, WH);
            if (redS1 || redS1r){ // Need to avoid lazy bool eval in order to always evaluate S2minr S2maxr
                reduceA21 = true;
                reduceS1 = true;
            }
            if (S2minr < S2min) S2min = S2minr;
            if (S2maxr > S2max) S2max = S2maxr;
        } else {
            S2min = S1min;
            S2max = S1max;
        }

        bool reduceS2 = false;
        bool reduceA11 = false;
        if (Protected::NeedPreAxpyReduction (S3min, S3max, S2min, S2max, WH.Amin, WH.Amax, negx, WH)){
            reduceS2 = true;
            reduceA11 = true;
        }

        bool reduceS3 = false;
        if (y){
            bool redS3 = Protected::NeedPreAxpyReduction (Smin, Smax, S3min, S3max, WH.Amin, WH.Amax, negy, WH);
            bool redS3r = Protected::NeedPreAxpyReduction (Sminr, Smaxr, S3min, S3max, WH.Amin, WH.Amax, y, WH);
            if (redS3 || redS3r) {
                reduceS3 = true;
                reduceA11=true;
            }
            if (Smin > Sminr) Smin = Sminr;
            if (Smax < Smaxr) Smax = Smaxr;
        } else {
            Smin = S3min;
            Smax = S3max;
        }
        bool reduceA22 = false;
        if (Protected::NeedPreSubReduction (Tmin, Tmax, WH.Cmin, WH.Cmax, S2min, S2max, WH)) {
            reduceA22 = true;
                // TODO: shouldn't we also reduce S2?
        }
        WH.initOut();
        
            // Note that the field must ensure that an axpy (a.b+c) can be stored in a storage_t.
            // This is more constraining than what Modular<int64_t,uint64_t>::MaxCardinality ensures.
            // A specific test is add in utils/test-utils.h to avoid generating such field in the testsuite.
            // TODO: add asserts here and there to warn user against this issue.

        if (reduceA21) freduce (F, Axxrows, Axxcols, A21, lda);
        if (reduceA11) freduce (F, Axxrows, Axxcols, A, lda);
        if (reduceA22) freduce (F, Axxrows, Axxcols, A22, lda);

        if (trans==FflasNoTrans){
                // S <- A21 Y
            for (size_t i=0; i<N2; ++i, A11+=lda, A11r+=lda, A21+=lda, A21r+=lda, A22+=lda, S+=lds, Sr+=lds, T+=ldt){
                if (reduceS1){
                    fscal (F, K2, x, A21, 1, S, 1);
                } else
                    fscal (DF, K2, x, A21, 1, S, 1);

                if (!F.isZero(y)){
                        faxpy (DF, K4, negy, A21r, 1, S, 1);
                        faxpy (DF, K4, y, A21, 1, Sr, 1);
                }
                if (reduceS2)
                    freduce (F, K2, S, 1);
                    // T <- A22 -S
                if (reduceA22 && reduceS2){
                    fsub (F, K2, A22, 1, S, 1, T, 1);
                } else {
                        fsub (DF, K2, A22, 1, S, 1, T, 1);
                        freduce (F, K2, T, 1);
                }
                    // S <- S - A11 Y
                faxpy (DF, K2, negx, A11, 1, S, 1);
                if (reduceS3 || F.isZero(y))
                    freduce (F, K2, S, 1);
                if (!F.isZero(y)){
                    faxpy (DF, K4, y, A11r, 1, S, 1);
                    faxpy (DF, K4, negy, A11, 1, Sr, 1);
                    freduce (F, K2, S, 1);
                }
            }
        } else { // FflasTrans         
            for (size_t i=0; i<K4; ++i, A11+=lda, A11r+=lda, A21+=lda, A21r+=lda, A22+=lda, A22r+=lda, S+=lds, Sr+=lds, T+=ldt, Tr+=ldt){
                // S <- Y A21
                
                if (reduceS1){
                    fscal (F, N2, x, A21, 1, S, 1);
                    fscal (F, N2, x, A21r, 1, Sr, 1);
                } else {
                    fscal (DF, N2, x, A21, 1, S, 1);
                    fscal (DF, N2, x, A21r, 1, Sr, 1);
                }
                if (!F.isZero(y)){
                    faxpy (DF, N2, negy, A21r, 1, S, 1);
                    faxpy (DF, N2, y, A21, 1, Sr, 1);
                }
                if (reduceS2){
                    freduce(F, N2, S, 1);
                    freduce(F, N2, Sr, 1);
                }
                    // T <- A22 -S
                if (reduceA22 && reduceS2){
                    fsub (F, N2, A22, 1, S, 1, T, 1);
                    fsub (F, N2, A22r, 1, Sr, 1, Tr, 1);
                } else {
                    fsub (DF, N2, A22, 1, S, 1, T, 1);
                    freduce (F, N2, T, 1);
                    fsub (DF, N2, A22r, 1, Sr, 1, Tr, 1);
                    freduce (F, N2, Tr, 1);
                }
                    // S <- S - Y A11
                faxpy (DF, N2, negx, A11, 1, S, 1);
                if (reduceS3 || F.isZero(y)) freduce (F, N2, S, 1);
                faxpy (DF, N2, negx, A11r, 1, Sr, 1);
                if (reduceS3 || F.isZero(y)) freduce (F, N2, Sr, 1);
                
                if (!F.isZero(y)){
                    faxpy (DF, N2, y, A11r, 1, S, 1);
                    freduce (F,N2,S, 1);
                    faxpy (DF, N2, negy, A11, 1, Sr, 1);
                    freduce (F,N2,Sr, 1);
                }
            }
        }
    }

    template<class Field>
    inline typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t N,
           const size_t K,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           MMHelper<Field, MMHelperAlgo::Winograd,  ModeCategories::DelayedTag, ParSeqHelper::Sequential> & H){
        if (!N) return C;
        if (!K || F.isZero(alpha)){
            fscalin (F, N, N, beta, C, ldc); // TODO UpLo
            return C;
        }
            // TODO convert from XXX to float/double or double to float as in fgemm
            //...
        typename Field::Element alpha_,beta_;
        if ( !F.isOne(alpha) && !F.isMOne(alpha)){
            F.assign (alpha_, F.one);
            F.div (beta_, beta, alpha);
        } else {
            F.assign (alpha_,alpha);
            F.assign (beta_,beta);
        }
        MMHelper<Field, MMHelperAlgo::Winograd, ModeCategories::LazyTag>  HD(H);

        fsyrk(F, UpLo, trans, N, K, alpha_, A, lda, beta_, C, ldc, HD);

        Protected::ScalAndReduce (F, N, N, alpha, C, ldc, HD);

        H.initOut();
        
        return C;
    }

    template<class Field, class Mode>
    inline typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t N,
           const size_t K,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           MMHelper<Field, MMHelperAlgo::Winograd, Mode> & H){

        if (!H.recLevel){
            MMHelper<Field, MMHelperAlgo::Classic, Mode>  H2 (H);
            fsyrk (F, UpLo, trans, N, K, alpha, A, lda, beta, C, ldc, H2);
            H.Outmin = H2.Outmin;
            H.Outmax = H2.Outmax;
            return C;
        }
        size_t q = 1 << (H.recLevel+1);
        size_t Ns = (N/q)*q;
        size_t Ks = (K/q)*q;
        
            // find a, b such that a^2 + b^2 = -1 mod p
        Givaro::Integer a,b;
        Givaro::IntSqrtModDom<> ISM;
        Givaro::ZRing<Givaro::Integer> Z;
        Z.init(a);
        Z.init(b);
        ISM.sumofsquaresmodprime (a, b, -1, F.characteristic());
        typename Field::Element y1, y2;
        F.init (y1, a);
        F.init (y2, b);

        typename Field::ConstElement_ptr A12, A21;
        if (trans==FflasNoTrans){
            A12 = A + Ks; A21 = A + Ns*lda;
        } else {
            A12 = A + Ks*lda; A21 = A + Ns;
        }
            // C11 = A11 x A11^T
        fsyrk_strassen (F, UpLo, trans, Ns, Ks, y1, y2, alpha,
                        FFPACK::fflas_const_cast<typename Field::Element_ptr>(A), lda, beta, C, ldc, H);

            // C11 += A12 x A12 ^T
        MMHelper<Field, MMHelperAlgo::Classic, Mode>  H2 (H);
        H2.Cmin = H.Outmin;
        H2.Cmax = H.Outmax;
        fsyrk (F, UpLo, trans, Ns, K-Ks, alpha, A12, lda, F.one, C, ldc, H2);

            // C22 = [A21 A22] x [A21 A22]^T
        fsyrk (F, UpLo, trans, N-Ns, K, alpha, A21, lda, beta, C+Ns*(ldc+1), ldc);

            // C21 = A21 x A11^T
        if (UpLo == FflasLower)
            fgemm (F, trans, (trans == FflasNoTrans)? FflasTrans : FflasNoTrans, N-Ns, Ns, K,
                   alpha, A21, lda, A, lda, beta, C+Ns*ldc, ldc);
        else
            fgemm (F, trans, (trans == FflasNoTrans)? FflasTrans : FflasNoTrans, Ns, N-Ns, K,
                   alpha, A, lda, A21, lda, beta, C+Ns, ldc);
        return C;
    }

        // Assumes that 2^(reclevel+1) divides N and K
    template<class Field, class FieldTrait>
    inline typename Field::Element_ptr
    fsyrk_strassen (const Field& F,
                    const FFLAS_UPLO uplo,
                    const FFLAS_TRANSPOSE trans,
                    const size_t N,
                    const size_t K,
                    const typename Field::Element y1,
                    const typename Field::Element y2,
                    const typename Field::Element alpha,
                    typename Field::Element_ptr A, const size_t lda,
                    const typename Field::Element beta,
                    typename Field::Element_ptr C, const size_t ldc,
                    MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait> & WH
                    ) {
            // Comments are written for the NoTrans, Lower version
        if (WH.recLevel == 0){
            MMHelper<Field, MMHelperAlgo::Classic, FieldTrait>  CH(WH);
            fsyrk (F, uplo, trans, N, K, alpha, A, lda, beta, C, ldc,CH);
            WH.Outmin = CH.Outmin;
            WH.Outmax = CH.Outmax;
            return C;
        }

        typedef MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > MMH_t;
        typedef typename  MMH_t::DelayedField::Element DFElt;

        const typename MMH_t::DelayedField & DF = WH.delayedField;

        size_t N2 = N>>1;
        size_t K2 = K>>1;
        size_t Arows, Acols;
        typename Field::Element_ptr A11 = A, A12, A21, A22;

        if (trans == FflasNoTrans){
            A12 = A + K2; A21 = A + N2*lda; A22 = A21 + K2; Arows = N2; Acols = K2;
        } else {
            A12 = A + K2*lda; A21 = A + N2; A22 = A12 + N2; Arows = K2; Acols = N2;
        }
        typename Field::Element_ptr C11 = C, C12, C21, C22;
        if (uplo == FflasLower){
            C12 = C11 + N2; C21  = C11 + N2*ldc; C22 = C21 + N2;
        } else {
            C12 = C11 + N2*ldc; C21  = C11 + N2; C22 = C12 + N2;
        }
        typename Field::Element_ptr S1, S2, S4;
        size_t lds;
        FFLAS_TRANSPOSE OppTrans = (trans == FflasNoTrans)? FflasTrans : FflasNoTrans;

        typename Field::Element negalpha;
        F.init(negalpha);
        F.neg(negalpha,alpha);

        if (F.isZero(beta)){ // no accumulation, schedule without extra temp

            if (K>N){ // Not possible to store S2 in C12 nor S1 in C21
                S1 = fflas_new(F, N2, K2);
                S2 = fflas_new(F, N2, K2);
                S4 = S1;
                if (trans==FflasNoTrans) lds = K2;
                else  lds = N2;
            } else{
                S1 = C21;
                S2 = C12;
                S4 = C11;
                lds = ldc;
            }
                // S1 = (A21-A11) x Y in S1 (C21)
                // S2 =  A22 - A21 x Y  in S2 (C12)
            MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait>  PH(WH);
            PH.Bmin = WH.Amin; PH.Bmax = WH.Amax; PH.Cmin = WH.Amin; PH.Cmax = WH.Amax;

            computeS1S2 (F, trans, N, K, y1, y2, A,lda, S1, lds, S2, lds, PH);

                //  P4^T =  S2 x S1^T in  C22
            MMH_t H4 (F, -1, PH.Outmin, PH.Outmax, PH.Outmin, PH.Outmax, 0,0);
            if (uplo == FflasLower)
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, S2, lds, S1, lds, F.zero, C22, ldc, H4);
            else
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, S1, lds, S2, lds, F.zero, C22, ldc, H4);

            if (K>N){ fflas_delete(S2); }

                // S3 = S1 - A22 in S1
            fsubin (DF, Arows, Acols, A22, lda, S1, lds);

                // P5 = S3 x S3^T in C12
            MMH_t H5 (F, WH.recLevel-1, PH.Outmin-PH.Cmax, PH.Outmax-PH.Cmin,
                      PH.Outmin-PH.Cmax, PH.Outmax-PH.Cmin, 0,0);
            fsyrk_strassen (F, uplo, trans, N2, K2, y1, y2, alpha, S1, lds, F.zero, C12, ldc, H5);

                // S4 = S3 + A12 in S4
            fadd (DF, Arows, Acols, S1, lds, A12, lda, S4, lds);

                // P3 = A22 x S4^T in C21
            MMH_t H3 (F, WH.recLevel-1, PH.Cmin, PH.Cmax, PH.Outmin+WH.Bmin-PH.Cmax, PH.Outmax+WH.Bmax-PH.Cmin, 0,0);
            if (uplo == FflasLower)
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, A22, lda, S4, lds, F.zero, C21, ldc, H3);
            else{
                H3.Amin = PH.Outmin+WH.Amin-PH.Cmax;
                H3.Amax = PH.Outmax+WH.Amax-PH.Cmin;
                H3.Bmin = PH.Cmin;
                H3.Bmax = PH.Cmax;
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, S4, lds, A22, lda, F.zero, C21, ldc, H3);
            }

            if (K>N){ fflas_delete(S1); }

                // P1 = A11 x A11^T in C11
            MMH_t H1 (F, WH.recLevel-1, PH.Amin, PH.Amax, PH.Amin, PH.Amax, 0,0);
            fsyrk_strassen (F, uplo, trans, N2, K2, y1, y2, alpha, A11, lda, F.zero, C11, ldc, H1);

                // U1 = P1 + P5 in C12
            DFElt U1Min, U1Max;
            if (Protected::NeedPreAddReduction (U1Min, U1Max, H1.Outmin, H1.Outmax, H5.Outmin, H5.Outmax, WH)){
                freduce(F,uplo,N2,C12,ldc);
                freduce(F,uplo,N2,C11,ldc);
            }
            faddin (DF, uplo, N2, C11, ldc, C12, ldc);

                // U2 = U1 + P4 in C12
            DFElt U2Min, U2Max;
            if (Protected::NeedPreAddReduction (U2Min, U2Max, U1Min, U1Max, H4.Outmin, H4.Outmax, WH)){
                freduce(F,uplo,N2,C12,ldc);
                freduce(F,N2,N2,C22,ldc);
            }
                // make U1 explicit: Up(U1)=Low(U1)^T
            if (uplo == FflasLower)
                for (size_t i=0; i<N2; i++)
                    fassign(F, i, C12+i*ldc, 1, C12+i, ldc);
            else
                for (size_t i=0; i<N2; i++)
                    fassign(F, i, C12+i, ldc, C12+i*ldc, 1);

            for (size_t i=0; i<N2; ++i){
                faddin (DF,  N2, C22+i*ldc, 1, C12+i, ldc);
            }

                // U4 = U2 + P3 in C21
            DFElt U4Min, U4Max;
            if (Protected::NeedPreAddReduction (U4Min, U4Max, U2Min, U2Max, H3.Outmin, H3.Outmax, WH)){
                freduce(F,N2,N2,C21,ldc);
                freduce(F,N2,N2,C12,ldc);
            }
            faddin (DF, N2, N2, C12, ldc, C21, ldc);

                // U5 = U2 + P4^T in C22
            DFElt U5Min, U5Max;
            if (Protected::NeedPreAddReduction (U5Min, U5Max, U2Min, U2Max, H4.Outmin, H4.Outmax, WH)){
                freduce(F,uplo,N2,C22,ldc);
                freduce(F,uplo,N2,C12,ldc);
            }
            faddin (DF, uplo, N2, C12, ldc, C22, ldc);

                // P2 = A12 x A12^T in C12
            MMH_t H2 (F, WH.recLevel-1, WH.Amin, WH.Amax, WH.Bmin, WH.Bmax, 0,0);
            fsyrk_strassen (F, uplo, trans, N2, K2, y1, y2, alpha, A12, lda, F.zero , C12, ldc, H2);

                // U3 = P1 + P2 in C11
            DFElt U3Min, U3Max;
            if (Protected::NeedPreAddReduction (U3Min, U3Max, H1.Outmin, H1.Outmax, H2.Outmin, H2.Outmax, WH)){
                freduce(F,uplo,N2,C11,ldc);
                freduce(F,uplo,N2,C12,ldc);
            }
            faddin (DF, uplo, N2, C12, ldc, C11, ldc);
                // Updating WH with Outmin, Outmax of the result
            WH.Outmin = min3 (U3Min, U4Min, U5Min);
            WH.Outmax = max3 (U3Max, U4Max, U5Max);
        } else { // with accumulation, schedule with 1 temp

            typename Field::Element negbeta;
            F.init (negbeta);
            F.neg (negbeta, beta);
            typename Field::Element_ptr T = fflas_new (F, N2, std::max(N2,K2));
            size_t ldt = (trans == FflasNoTrans)?std::max(N2,K2) : N2;

            if (K>N){ // Not possible to store S2 in C12
                S2 = fflas_new(F, N2, K2);
                if (trans == FflasNoTrans) lds = K2;
                else lds = N2;
            } else{
                S2 = C12;
                lds = ldc;
            }
                // S1 = (A21-A11) x Y^T in T1
                // S2 = A22 - A21 x Y^T in C12
            MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait>  PH(WH);
            PH.Bmin = WH.Amin; PH.Bmax = WH.Amax; PH.Cmin = WH.Amin; PH.Cmax = WH.Amax;

            computeS1S2 (F, trans, N, K, y1, y2, A, lda, T, ldt, S2, lds, PH);

                // Up(C11) = Low(C22) (saving C22)
            if (uplo == FflasLower)
                for (size_t i=1; i<N2; ++i)
                    fassign (F, N2-i, C22 + (N2-i)*ldc, 1, C11 + 1 + (i-1)*(ldc+1), 1);
            else 
                for (size_t i=1; i<N2; ++i)
                    fassign (F, N2-i, C22 + 1 + (i-1)*(ldc+1), 1,  C11 + (N2-i)*ldc, 1);

                // temp for storing the diagonal of C22
            typename Field::Element_ptr DC22 = fflas_new(F,N2);
            for (size_t i=0; i<N2; ++i)
                F.assign (DC22[i], C22[i*(ldc+1)]);

                // P4^T = S2 x S1^T in C22
            MMH_t H4 (F,WH.recLevel-1, PH.Outmin, PH.Outmax, PH.Outmin, PH.Outmax, 0,0);
            if (uplo == FflasLower)
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, S2, lds, T, ldt, F.zero, C22, ldc, H4);
            else
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, T, ldt, S2, lds, F.zero, C22, ldc, H4);

            if (K>N){ fflas_delete(S2); }

                // S3 = S1 - A22 in T
            fsubin (DF, Arows, Acols, A22, lda, T, ldt);

                // P5 = S3 x S3^T in C12
            MMH_t H5 (F, WH.recLevel-1, PH.Outmin-PH.Cmax, PH.Outmax-PH.Cmin, PH.Outmin-PH.Cmax, PH.Outmax-PH.Cmin, 0,0);
            fsyrk_strassen (F, uplo, trans, N2, K2, y1, y2, alpha, T, ldt, F.zero, C12, ldc, H5);

                // S4 = S3 + A12 in T1
            faddin (DF, Arows, Acols, A12, lda, T, ldt);

                //  P3 = A22 x S4^T + beta C21 in C21
            MMH_t H3 (F, WH.recLevel-1, PH.Cmin, PH.Cmax, PH.Outmin+WH.Bmin-PH.Cmax, PH.Outmax+WH.Bmax-PH.Cmin, WH.Cmin, WH.Cmax);
            if (uplo == FflasLower)
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, A22, lda, T, ldt, beta, C21, ldc, H3);
            else{
                H3.Amin = PH.Outmin+WH.Amin-WH.Amax; 
                H3.Amax = PH.Outmax+WH.Amax-WH.Amin; 
                H3.Bmin = PH.Cmin;
                H3.Bmax = PH.Cmax;
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, T, ldt, A22, lda, beta, C21, ldc, H3);
            }

                // P1 = A11 x A11^T in T1
            MMH_t H1 (F, WH.recLevel-1, PH.Amin, PH.Amax, PH.Amin, PH.Amax, 0, 0);
            fsyrk_strassen (F, uplo, trans, N2, K2, y1, y2, alpha, A11, lda, F.zero, T, ldt, H1);

                // U1 = P5 + P1  in C12 // Still symmetric
            DFElt U1Min, U1Max;
            if (Protected::NeedPreAddReduction (U1Min, U1Max, H5.Outmin, H5.Outmax, H1.Outmin, H1.Outmax, WH)){
                freduce(F,uplo,N2,C12,ldc);
                freduce(F,uplo,N2,T,ldt);
           }
            faddin (DF, uplo, N2, T, ldt, C12, ldc);

                // U2 = U1 + P4 in C12
            DFElt U2Min, U2Max;
            if (Protected::NeedPreAddReduction (U2Min, U2Max, U1Min, U1Max, H4.Outmin, H4.Outmax, WH)){
                freduce(F,uplo,N2,C12,ldc);
                freduce(F,N2,N2,C22,ldc);
            }
                // Make U1 explicit (copy the N^2/2 missing part)
            if (uplo == FflasLower)
                for (size_t i=0; i<N2; ++i)
                    fassign (DF, i, C12 + i*ldc, 1, C12 + i, ldc);
            else
                for (size_t i=0; i<N2; ++i)
                    fassign (DF, i, C12 + i, ldc, C12 + i*ldc, 1);                

            for (size_t i=0; i<N2; i++)
                faddin (DF, N2, C22 + i*ldc, 1, C12 + i, ldc);

                // U4 = U2 + P3 in C21
            DFElt U4Min, U4Max;
            if (Protected::NeedPreAddReduction (U4Min, U4Max, U2Min, U2Max, H3.Outmin, H3.Outmax, WH)){
                freduce(F,N2,N2,C21,ldc);
                freduce(F,N2,N2,C12,ldc);
            }
            faddin (DF, N2, N2, C12, ldc, C21, ldc);

                // U5 = U2 + P4^T  in C22 (only the lower triang part)
            DFElt U5Min, U5Max;
            if (Protected::NeedPreAddReduction (U5Min, U5Max, U2Min, U2Max, H4.Outmin, H4.Outmax, WH)){
                freduce(F,uplo,N2,C22,ldc);
                freduce(F,uplo,N2,C12,ldc);
            }
            faddin (DF, uplo, N2, C12, ldc, C22, ldc);
            freduce(F,N2,N2,C22,ldc);

                // U5' = U5 +  beta Up(C11)^T in C22
                // TODO use delayed field and a needPreAXPYReduction
                // Suspicious BUG: if C22 is not reduced, possible overflow???
            if (uplo == FflasLower)
                for (size_t i=1; i<N2; i++){ // TODO factorize out in a triple add
                    faxpy (F, N2-i, beta, C11 + 1+(i-1)*(ldc+1), 1, C22 + (N2-i)*ldc, 1);
                    F.axpyin (C22[(i-1)*(ldc+1)], beta, DC22[i-1]);
                }
            else
                for (size_t i=1; i<N2; i++){ // TODO factorize out in a triple add
                    faxpy (F, N2-i, beta, C11 + (N2-i)*ldc, 1, C22 + 1 + (i-1)*(ldc+1), 1);
                    F.axpyin (C22[(i-1)*(ldc+1)], beta, DC22[i-1]);
                }
            F.axpyin (C22[(N2-1)*(ldc+1)], beta, DC22[N2-1]);
            U5Min = WH.FieldMin;
            U5Max = WH.FieldMax;

            fflas_delete(DC22);

                // P2 = A12 x A12^T + beta C11 in C11
            MMH_t H2 (F, WH.recLevel-1, WH.Amin, WH.Amax, WH.Bmin, WH.Bmax, WH.Cmin, WH.Cmax);
            fsyrk_strassen (F, uplo, trans, N2, K2, y1, y2, alpha, A12, lda, beta, C11, ldc, H2);

                // U3 = P1 + P2 in C11
            DFElt U3Min, U3Max;
            if (Protected::NeedPreAddReduction (U3Min, U3Max, H1.Outmin, H1.Outmax, H2.Outmin, H2.Outmax, WH)){
                freduce(F,uplo,N2,C11,ldc);
                freduce(F,uplo,N2,T,ldt);
            }
            faddin (DF, uplo, N2, T, ldt, C11, ldc);

            fflas_delete (T);
            WH.Outmin = min3 (U3Min, U4Min, U5Min);
            WH.Outmax = max3 (U3Max, U4Max, U5Max);
        }
        WH.checkOut(F, uplo, N, N, C, ldc);
        return C;
    }
}


#endif // __FFLASFFPACK_fflas_fsyrk_strassen_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
