/*
 * Copyright (C) 2016 the FFLAS-FFPACK group
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

#ifndef __FFLASFFPACK_fflas_fsyrk_INL
#define __FFLASFFPACK_fflas_fsyrk_INL
namespace FFLAS{namespace Protected{
    template<class NewField, class Field, class FieldMode>
    inline typename Field::Element_ptr
    fsyrk_convert (const Field& F,
                   const FFLAS_UPLO UpLo,
                   const FFLAS_TRANSPOSE trans,
                   const size_t N,
                   const size_t K,
                   const typename Field::Element alpha,
                   typename Field::ConstElement_ptr A, const size_t lda,
                   const typename Field::Element beta,
                   typename Field::Element_ptr C, const size_t ldc,
                   MMHelper<Field, MMHelperAlgo::Classic, FieldMode> & H)
    {
        typedef typename NewField::Element FloatElement;
        NewField G((FloatElement) F.characteristic());
        FloatElement tmp,alphaf, betaf;
        // This conversion is quite tricky, but convert and init are required
        // in sequence e.g. for when F is a ModularBalanced field and alpha == -1
        F.convert (tmp, beta);
        G.init(betaf, tmp);
        F.convert (tmp, alpha);
        G.init(alphaf, tmp);

        FloatElement* Af = FFLAS::fflas_new(G, N, K);
        FloatElement* Cf = FFLAS::fflas_new(G, N, N);

        size_t ma, ka;
        if (trans == FflasTrans) { ma = K; ka = N; }
        else { ma = N; ka = K; }

        fconvert(F, ma, ka, Af, ka, A, lda);
        freduce(G, ma, ka, Af, ka);

        if (!F.isZero(beta)){
            fconvert(F, N, N, Cf, N, C, ldc); // @todo: take advantage of the symmetry
            freduce (G, UpLo, N, N, Cf, N);
        }
        MMHelper<NewField, MMHelperAlgo::Classic> HG(G,H.recLevel, ParSeqHelper::Sequential());
        fsyrk (G, UpLo, trans, N, K, alphaf, Af, ka, betaf, Cf, N, HG);

        finit (F, N, N, Cf, N, C, ldc); // @todo: take advantage of the symmetry
        fflas_delete (Af);
        fflas_delete (Cf);
        return C;
    }
}// Protected
}// FFLAS
namespace FFLAS {
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
           typename Field::Element_ptr C, const size_t ldc)
    {
        if (!N) return C;
        if (!K || F.isZero (alpha)){
            fscalin(F, N, N, beta, C, ldc);
            return C;
        }
        fsyrk(F,UpLo, trans, N, K, alpha, A, lda, beta, C, ldc, FFLAS::ParSeqHelper::Sequential());
        return C;
        
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
           const ParSeqHelper::Sequential seq )
    {
            //std::cerr<<"fsyrk PSH::Seq"<<std::endl;
        MMHelper<Field, MMHelperAlgo::Classic, typename FFLAS::ModeTraits<Field>::value, ParSeqHelper::Sequential > HW (F, N, K, N, seq);
        return fsyrk (F, UpLo, trans, N, K, alpha, A, lda, beta, C, ldc, HW);
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
           MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DefaultTag> & H){
            // The default implementation
        if (UpLo == FflasUpper){
            if (trans == FflasNoTrans){
                for (size_t i=0; i<N; i++)
                    for (size_t j=i; j<N; j++){
                        F.mulin (C[i*ldc+j],beta);
                        F.axpyin (C[i*ldc+j],
                                  alpha,
                                  fdot(F, K, A+i*lda, 1, A+j*lda, 1));
                    }
            } else { // Trans
                for (size_t i=0; i<N; i++)
                    for (size_t j=i; j<N; j++){
                        F.mulin (C[i*ldc+j],beta);
                        F.axpyin (C[i*ldc+j],
                                  alpha,
                                  fdot(F, K, A+i, lda, A+j, lda));
                    }
            }
        } else { // Lower
            if (trans == FflasNoTrans){
                for (size_t i=0; i<N; i++)
                    for (size_t j=0; j<=i; j++){
                        F.mulin (C[i*ldc+j],beta);
                        F.axpyin (C[i*ldc+j],
                                  alpha,
                                  fdot(F, K, A+i*lda, 1, A+j*lda, 1));
                    }
            } else { // Trans
                for (size_t i=0; i<N; i++)
                    for (size_t j=0; j<=i; j++){
                        F.mulin (C[i*ldc+j],beta);
                        F.axpyin (C[i*ldc+j],
                                  alpha,
                                  fdot(F, K, A+i, lda, A+j, lda));
                    }
            }
        }
        return C;
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
           MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::ConvertTo<ElementCategories::MachineFloatTag>, ParSeqHelper::Sequential> & H)
     {
             //std::cerr<<"fsyrk Classic ConvertTo"<<std::endl;
         if (!std::is_same<Field,Givaro::Modular<float> >::value){
            if (F.cardinality() == 2)
                return Protected::fsyrk_convert<Givaro::Modular<float>,Field>(F,UpLo,trans,N,K,alpha,A,lda,beta,C,ldc,H);
            else if (!std::is_same<Field,Givaro::ModularBalanced<float> >::value){
                if (F.cardinality() < DOUBLE_TO_FLOAT_CROSSOVER)
                    return Protected::fsyrk_convert<Givaro::ModularBalanced<float>,Field>(F,UpLo,trans,N,K,alpha,A,lda,beta,C,ldc,H);
                else if (!std::is_same<Field,Givaro::ModularBalanced<double> >::value &&  16*F.cardinality() < Givaro::ModularBalanced<double>::maxCardinality())
                    return Protected::fsyrk_convert<Givaro::ModularBalanced<double>,Field>(F,UpLo, trans,N, K,alpha,A,lda,beta,C,ldc,H);
            }
        } else {
                 // Fall back case
             FFPACK::failure()(__func__,__LINE__,"Invalid ConvertTo Mode for this field");
        }
        return C;

     }


    namespace Protected {
        template <class Field, class AlgoT, class ParSeqTrait>
        inline void ScalAndReduce (const Field& F, const FFLAS_UPLO UpLo, const size_t N,
                                   const typename Field::Element alpha,
                                   typename Field::Element_ptr A, const size_t lda,
                                   const MMHelper<Field, AlgoT, ModeCategories::LazyTag, ParSeqTrait >& H)
        {
            if (!F.isOne(alpha) && !F.isMOne(alpha)){
                typename MMHelper<Field, AlgoT, ModeCategories::LazyTag, ParSeqTrait >::DFElt al;
                F.convert(al, alpha);
                if (al<0) al = -al;
                if (std::max(-H.Outmin, H.Outmax) > H.MaxStorableValue/al){
                    freduce (F, UpLo, N, A, lda);
                    fscalin (F, N, N, alpha, A, lda);
                } else {
                    fscalin (H.delayedField, N, N, alpha, (typename MMHelper<Field, AlgoT, ModeCategories::LazyTag, ParSeqTrait >::DFElt*)A, lda);
                    freduce (F, UpLo, N, A, lda);
                }
            } else
                freduce (F, UpLo, N, A, lda);
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
           MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DelayedTag> & H)
    {
            //std::cerr<<"fsyrk Classic Delayed"<<std::endl;
        if (!N) return C;
        if (!K || F.isZero (alpha)){
            fscalin(F, N, N, beta, C, ldc);
            return C;
        }
        
        typename Field::Element alpha_,beta_;
        if ( !F.isOne(alpha) && !F.isMOne(alpha)){
            F.assign (alpha_, F.one);
            F.div (beta_, beta, alpha);
        } else {
            F.assign (alpha_,alpha);
            F.assign (beta_,beta);
        }
        MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::LazyTag>  HD(H);

        fsyrk (F, UpLo, trans, N, K, alpha_, A, lda, beta_, C, ldc, HD);

        Protected::ScalAndReduce (F, UpLo, N, alpha, C, ldc, HD);

        H.initOut();

        return C;
    }

// F is a field supporting delayed reductions
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
           MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::LazyTag> & H)
    {
        if (!N) return C;
        if (!K || F.isZero (alpha)){
            fscalin(F, N, N, beta, C, ldc);
            if (beta> 0){
                H.Outmin = H.Cmin * beta;
                H.Outmax = H.Cmax * beta;
            }else{
                H.Outmin = H.Cmax * beta;
                H.Outmax = H.Cmin * beta;
            }
            return C;
        }
            //std::cerr<<"fsyrk Classic Lazy"<<std::endl;
        // Input matrices are unreduced: need to figure out the best option between:
        // - reducing them
        // - making possibly more blocks (smaller kmax)
        typedef MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::LazyTag> HelperType;
        typename HelperType::DelayedField::Element alphadf, betadf;
        betadf = beta;
        if (F.isMOne (alpha)) {
            alphadf = -H.delayedField.one;
        } else {
            alphadf = F.one;
            if (! F.isOne( alpha)) {
                // Compute y = A*x + beta/alpha.y
                // and after y *= alpha
                FFLASFFPACK_check(!F.isZero(alpha));
                typename Field::Element betadalpha;
                F.init(betadalpha);
                F.div (betadalpha, beta, alpha);
                betadf = betadalpha;
            }
        }

        if (F.isMOne(betadf)) betadf = -F.one;

        size_t kmax = H.MaxDelayedDim (betadf);
        H.checkA(F, trans, N, K, A, lda);
        if (kmax <=  K/2 || H.Aunfit() ){
            // Might as well reduce inputs
            if (H.Amin < H.FieldMin || H.Amax>H.FieldMax){
                H.initA();
                H.initB();
                freduce_constoverride (F, (trans==FflasNoTrans)?N:K, (trans==FflasNoTrans)?K:N, A, lda);
            }
            if (!F.isZero(beta) && (H.Cmin < H.FieldMin || H.Cmax>H.FieldMax)){
                H.initC();
                freduce (F, UpLo, N, C, ldc);
            }
            kmax = H.MaxDelayedDim (betadf);
        }

        if (!kmax){
            MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DefaultTag> HG(H);
            H.initOut();
            return fsyrk (F, UpLo, trans, N, K, alpha, A, lda, beta, C, ldc, HG);
        }

        size_t k2 = std::min(K,kmax);
        size_t nblock = K / kmax;
        size_t remblock = K % kmax;
        if (!remblock) {
            remblock = kmax;
            --nblock;
        }
        size_t shiftA;
        if (trans == FflasTrans) shiftA = k2*lda;
        else shiftA = k2;

        typedef MMHelper<typename HelperType::DelayedField, MMHelperAlgo::Classic, ModeCategories::DefaultBoundedTag> DelayedHelper_t;
        DelayedHelper_t Hfp(H);
        typedef typename HelperType::DelayedField::Element DFElt;
        typedef typename HelperType::DelayedField::Element_ptr DFElt_ptr;
        typedef typename HelperType::DelayedField::ConstElement_ptr DFCElt_ptr;

        fsyrk (H.delayedField, UpLo, trans, N, remblock, alphadf,
               (DFCElt_ptr)A +nblock*shiftA, lda,
               betadf, (DFElt_ptr)C, ldc, Hfp);

        for (size_t i = 0; i < nblock; ++i) {
            freduce (F, UpLo, N, C, ldc);
            Hfp.initC();
            fsyrk (H.delayedField, UpLo, trans, N, k2, alphadf,
                   (DFCElt_ptr)A +i*shiftA, lda,
                   F.one, (DFElt_ptr)C, ldc, Hfp);
        }

        if (!F.isOne(alpha) && !F.isMOne(alpha)){
            DFElt al; F.convert(al, alpha);
            if (al<0) al = -al;
            // This cast is needed when Outmin base type is int8/16_t,
            // getting -Outmin returns a int, not the same base type.
            if (std::max(static_cast<const decltype(Hfp.Outmin)&>(-Hfp.Outmin), Hfp.Outmax)
                >Hfp.MaxStorableValue/al){
                freduce (F, UpLo, N, C, ldc);
                Hfp.initOut();
            }

            fscalin(H.delayedField, N,N,alpha,(typename DelayedHelper_t::DelayedField_t::Element_ptr)C,ldc);

            if (alpha>0){
                H.Outmin = (const DFElt)(alpha) * Hfp.Outmin;
                H.Outmax = (const DFElt)alpha * Hfp.Outmax;
            } else {
                H.Outmin = (const DFElt)alpha * Hfp.Outmax;
                H.Outmax = (const DFElt)alpha * Hfp.Outmin;
            }
        }else {
            H.Outmin = Hfp.Outmin;
            H.Outmax = Hfp.Outmax;
        }
        H.checkOut(F, UpLo, N, N, C, ldc);
        return C;
    }

    template<class Field, typename Mode>
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
           MMHelper<Field, MMHelperAlgo::DivideAndConquer, Mode> & H) {

        //@TODO: write an optimized iterative basecase
        // if (N==1){ // Base case
        //     F.mulin (*C, beta);
        //     size_t incA = (trans==FFLAS::FflasNoTrans)?1:lda;
        //     F.axpyin (*C, alpha, fdot (F, K, A, incA, A, incA));
        //     return C;
        if (H.recLevel == 0){
            MMHelper<Field, MMHelperAlgo::Classic,Mode> CH (H);
            return fsyrk(F, UpLo, trans, N, K, alpha, A, lda, beta, C, ldc, CH);
        } else {
            size_t N1 = N>>1;
            size_t N2 = N - N1;
            // Comments written for the case UpLo==FflasUpper, trans==FflasNoTrans
            FFLAS_TRANSPOSE oppTrans;
            if (trans==FflasNoTrans) {oppTrans=FflasTrans;}
            else {oppTrans=FflasNoTrans;}

            typename Field::ConstElement_ptr A2 = A + N1*(trans==FflasNoTrans?lda:1);
            typename Field::Element_ptr C12 = C + N1;
            typename Field::Element_ptr C21 = C + N1*ldc;
            typename Field::Element_ptr C22 = C12 + N1*ldc;
            // C11 <- alpha A1 x A1^T + beta C11
            MMHelper<Field, MMHelperAlgo::DivideAndConquer, Mode> CH (F, H.recLevel - 1);
            fsyrk (F, UpLo, trans, N1, K, alpha, A, lda, beta, C, ldc, CH);
            // C22 <- alpha A2 x A2^T + beta C22
            fsyrk (F, UpLo, trans, N2, K, alpha, A2, lda, beta, C22, ldc, CH);
            if (UpLo == FflasUpper) {
                // C12 <- alpha A1 * A2^T + beta C12
                fgemm (F, trans, oppTrans, N1, N2, K, alpha, A, lda, A2, lda, beta, C12, ldc);
            } else {
                // C21 <- alpha A2 * A1^T + beta C21
                fgemm (F, trans, oppTrans, N2, N1, K, alpha, A2, lda, A, lda, beta, C21, ldc);
            }
            return C;
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
           MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DefaultBoundedTag> & H) {

            //std::cerr<<"fsyrk Classic DefaultBounded"<<std::endl;
        MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DefaultTag>  Hd(H);
        fsyrk (F, UpLo, trans, N, K, alpha, A, lda, beta, C, ldc, Hd);
        H.setOutBounds (K,alpha,beta);
       return C;
    }

    inline Givaro::FloatDomain::Element_ptr
    fsyrk (const Givaro::FloatDomain& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t N,
           const size_t K,
           const Givaro::FloatDomain::Element alpha,
           Givaro::FloatDomain::ConstElement_ptr A, const size_t lda,
           const Givaro::FloatDomain::Element beta,
           Givaro::FloatDomain::Element_ptr C, const size_t ldc,
           MMHelper<Givaro::FloatDomain, MMHelperAlgo::Classic, ModeCategories::DefaultTag> &H) {
            //std::cerr<<"fsyrk Classic Default FloatDomain"<<std::endl;
         cblas_ssyrk (CblasRowMajor, (CBLAS_UPLO) UpLo, (CBLAS_TRANSPOSE) trans, N, K, alpha, A, lda, beta, C, ldc);
        return C;
    }

    inline Givaro::DoubleDomain::Element_ptr
    fsyrk (const Givaro::DoubleDomain& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t N,
           const size_t K,
           const Givaro::DoubleDomain::Element alpha,
           Givaro::DoubleDomain::ConstElement_ptr A, const size_t lda,
           const Givaro::DoubleDomain::Element beta,
           Givaro::DoubleDomain::Element_ptr C, const size_t ldc,
           MMHelper<Givaro::DoubleDomain, MMHelperAlgo::Classic, ModeCategories::DefaultTag> &H) {
            //std::cerr<<"fsyrk Classic Default DoubleDomain"<<std::endl;
        cblas_dsyrk (CblasRowMajor, (CBLAS_UPLO) UpLo, (CBLAS_TRANSPOSE) trans, N, K, alpha, A, lda, beta, C, ldc);
        return C;
    }

    template<class Field>
    inline typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t N,
           const size_t K,
           const typename Field::Element alpha,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::ConstElement_ptr D, const size_t incD,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc, const size_t threshold){
        return fsyrk (F, UpLo, trans, N, K, alpha, A, lda, D, incD, beta, C, ldc, ParSeqHelper::Sequential(), threshold);
    }
    template<class Field>
    inline typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t N,
           const size_t K,
           const typename Field::Element alpha,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::ConstElement_ptr D, const size_t incD,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           const ParSeqHelper::Sequential seq,
           const size_t threshold){
        size_t incRow,incCol;
        FFLAS_TRANSPOSE oppTrans;
        if (trans==FflasNoTrans) {incRow=lda;incCol=1;oppTrans=FflasTrans;}
        else {incRow = 1; incCol = lda;oppTrans=FflasNoTrans;}

        if (N <= threshold){
            typename Field::Element_ptr temp = fflas_new (F, N, K);
            size_t ldt;
            if (trans==FFLAS::FflasNoTrans){
                ldt =K;
                fassign(F, N, K, A, lda, temp, ldt);
            } else{
                ldt = N;
                fassign(F, K, N, A, lda, temp, ldt);
            }
            typename Field::Element_ptr Ai = A;
            typename Field::ConstElement_ptr Di = D;
            for (; Ai != A + K*incCol; Ai += incCol, Di+=incD){
                fscalin (F, N, *Di, Ai, incRow);
            }
            FFLAS::fgemm (F, trans, oppTrans, N, N, K, alpha, A, lda, temp, ldt, beta, C, ldc);
            FFLAS::fflas_delete(temp);
            return C;

        }else {
            size_t N1 = N>>1;
            size_t N2 = N - N1;
            // Comments written for the case UpLo==FflasUpper, trans==FflasNoTrans

            typename Field::Element_ptr A2 = A + N1*incRow;
            typename Field::Element_ptr C12 = C + N1;
            typename Field::Element_ptr C21 = C + N1*ldc;
            typename Field::Element_ptr C22 = C12 + N1*ldc;

            typename Field::Element_ptr temp = fflas_new (F, std::max(N2,N1),K);
            size_t ldt, incRowT,incColT;
            if (trans==FflasNoTrans) {ldt=K; incRowT=ldt; incColT=1;}
            else {ldt = N2; incRowT=1; incColT=ldt;}

            // temp <- A2 x D1
            typename Field::Element_ptr Ai = A2, Ti = temp;
            typename Field::ConstElement_ptr Di = D;
            for (; Ai != A2 + K*incCol; Ai += incCol, Ti += incColT, Di+=incD){
                fscal (F, N2, *Di, Ai, incRow, Ti, incRowT);
            }
            if (UpLo == FflasUpper) {
                // C12 <- alpha A1 x temp^T + beta C12
                fgemm (F, trans, oppTrans, N1, N2, K, alpha, A, lda, temp, ldt, beta, C12, ldc);
            } else {
                // C21 <- alpha temp x A11^T + beta C21
                fgemm (F, trans, oppTrans, N2, N1, K, alpha, temp, ldt, A, lda, beta, C21, ldc);
            }
            fflas_delete (temp);

            // C11 <- alpha A1 x D1 x A1^T + beta C11 and A1 <- A1 x D1
            fsyrk (F, UpLo, trans, N1, K, alpha, A, lda, D, incD, beta, C, ldc, threshold);
            // C22 <- alpha A2 x D1 x A2^T + beta C22 and A2 <- A2 x D1
            fsyrk (F, UpLo, trans, N2, K, alpha, A2, lda, D, incD, beta, C22, ldc, threshold);

            return C;
        }
    }
    template<class Field, class Cut, class Param>
    inline typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t N,
           const size_t K,
           const typename Field::Element alpha,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::ConstElement_ptr D, const size_t incD,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           const ParSeqHelper::Parallel<Cut,Param> par,
           const size_t threshold){

        size_t incRow,incCol;
        FFLAS_TRANSPOSE oppTrans;
        if (trans==FflasNoTrans) {incRow=lda;incCol=1;oppTrans=FflasTrans;}
        else {incRow = 1; incCol = lda;oppTrans=FflasNoTrans;}

        if (N <= threshold){
            typename Field::Element_ptr temp = fflas_new (F, N, K);
            size_t ldt;
            if (trans==FFLAS::FflasNoTrans){
                ldt =K;
                fassign(F, N, K, A, lda, temp, ldt);
            } else{
                ldt = N;
                fassign(F, K, N, A, lda, temp, ldt);
            }
            typename Field::Element_ptr Ai = A;
            typename Field::ConstElement_ptr Di = D;
            for (; Ai != A + K*incCol; Ai += incCol, Di+=incD){
                fscalin (F, N, *Di, Ai, incRow);
            }
            FFLAS::fgemm (F, trans, oppTrans, N, N, K, alpha, A, lda, temp, ldt, beta, C, ldc, par);
            FFLAS::fflas_delete(temp);
            return C;

        }else {
            size_t N1 = N>>1;
            size_t N2 = N - N1;
            // Comments written for the case UpLo==FflasUpper, trans==FflasNoTrans

            typename Field::Element_ptr A2 = A + N1*incRow;
            typename Field::Element_ptr C12 = C + N1;
            typename Field::Element_ptr C21 = C + N1*ldc;
            typename Field::Element_ptr C22 = C12 + N1*ldc;

            size_t nt = par.numthreads();
            if (nt == 1)
                return fsyrk(F, UpLo, trans, N, K, alpha, A, lda, D, incD, beta, C, ldc, ParSeqHelper::Sequential(), threshold);
            size_t nt2 = nt >> 1;
            size_t ntr = nt - nt2;
            ParSeqHelper::Parallel<Cut, Param> ps_rec1(nt2);
            ParSeqHelper::Parallel<Cut, Param> ps_rec2(ntr);
            ParSeqHelper::Parallel<CuttingStrategy::Block,StrategyParameter::Threads> ps_fgemm (nt);

            typename Field::Element_ptr temp = fflas_new (F, std::max(N2,N1),K);
            size_t ldt, incRowT,incColT;
            if (trans==FflasNoTrans) {ldt=K; incRowT=ldt; incColT=1;}
            else {ldt = N2; incRowT=1; incColT=ldt;}

            // temp <- A2 x D1
            typename Field::Element_ptr Ai = A2, Ti = temp;
            typename Field::ConstElement_ptr Di = D;
            //TODO parallelize fscal
            //TODO more cache efficient fscal
            for (; Ai != A2 + K*incCol; Ai += incCol, Ti += incColT, Di+=incD){
                fscal (F, N2, *Di, Ai, incRow, Ti, incRowT);
            }
            //TODO parallelize fgemm
            if (UpLo == FflasUpper) {
                // C12 <- alpha A1 x temp^T + beta C12
                fgemm (F, trans, oppTrans, N1, N2, K, alpha, A, lda, temp, ldt, beta, C12, ldc, ps_fgemm);
            } else {
                // C21 <- alpha temp x A11^T + beta C21
                fgemm (F, trans, oppTrans, N2, N1, K, alpha, temp, ldt, A, lda, beta, C21, ldc, ps_fgemm);
            }
            fflas_delete (temp);

            SYNCH_GROUP(
                        // C11 <- alpha A1 x D1 x A1^T + beta C11 and A1 <- A1 x D1
                        TASK(MODE(READ(D[0]) READWRITE(A[0]) WRITE(C[0]) CONSTREFERENCE(A, D, C, F, ps_rec1)),
                             fsyrk (F, UpLo, trans, N1, K, alpha, A, lda, D, incD, beta, C, ldc, ps_rec1, threshold));
                        // C22 <- alpha A2 x D1 x A2^T + beta C22 and A2 <- A2 x D1
                        TASK(MODE(READ(D[0]) READWRITE(A2[0]) WRITE(C22[0]) CONSTREFERENCE(A2, D, C22, F, ps_rec2)),
                             fsyrk (F, UpLo, trans, N2, K, alpha, A2, lda, D, incD, beta, C22, ldc, ps_rec2, threshold));
                       );
            return C;
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
           typename Field::Element_ptr A, const size_t lda,
           typename Field::ConstElement_ptr D, const size_t incD,
           const std::vector<bool>& twoBlocks,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc, const size_t threshold){

        size_t incRow,incCol;
        FFLAS_TRANSPOSE oppTrans;
        if (trans==FflasNoTrans) {incRow=lda;incCol=1;oppTrans=FflasTrans;}
        else {incRow = 1; incCol = lda;oppTrans=FflasNoTrans;}

        if (N <= threshold){
            typename Field::Element_ptr temp = fflas_new (F, N, K);
            size_t ldt,incRowT,incColT;
            if (trans==FFLAS::FflasNoTrans){
                ldt =K;
                incRowT=ldt; incColT=1;
                fassign(F, N, K, A, lda, temp, ldt);
            } else{
                ldt = N;
                incRowT=1; incColT=ldt;
                fassign(F, K, N, A, lda, temp, ldt);
            }
            typename Field::Element_ptr Ai = A;
            typename Field::Element_ptr tempi = temp;
            typename Field::ConstElement_ptr Di = D;
            for (size_t i=0; i<K;i++, Ai += incCol, Di+=incD, tempi+=incColT){
                if (!twoBlocks[i])
                    fscalin (F, N, *Di, Ai, incRow);
                else{
                    fscal (F, N, *Di, Ai+incCol,incRow, Ai,incRow);
                    fscal (F, N, *Di, tempi,incRowT, Ai+incCol,incRow);
                    Ai+=incCol; Di+=incD; tempi+=incColT;i++;
                }
            }
            FFLAS::fgemm (F, trans, oppTrans, N, N, K, alpha, A, lda, temp, ldt, beta, C, ldc);
            FFLAS::fflas_delete(temp);
            return C;

        }else {
            size_t N1 = N>>1;
            if (twoBlocks[N1-1]) N1++; // don't split a 2x2 block
            size_t N2 = N - N1;
            // Comments written for the case UpLo==FflasUpper, trans==FflasNoTrans

            typename Field::Element_ptr A2 = A + N1*incRow;
            typename Field::Element_ptr C12 = C + N1;
            typename Field::Element_ptr C21 = C + N1*ldc;
            typename Field::Element_ptr C22 = C12 + N1*ldc;

            typename Field::Element_ptr temp = fflas_new (F, std::max(N2,N1),K);
            size_t ldt, incRowT,incColT;
            if (trans==FflasNoTrans) {ldt=K; incRowT=ldt; incColT=1;}
            else {ldt = N2; incRowT=1; incColT=ldt;}

            // temp <- A2 x D1
            typename Field::Element_ptr Ai = A2, Ti = temp;
            typename Field::ConstElement_ptr Di = D;
            for (size_t i=0; i<K; Ai += incCol, Ti += incColT, Di+=incD,i++){
                if (!twoBlocks[i])
                    fscal (F, N2, *Di, Ai, incRow, Ti, incRowT);
                else {
                    fscal (F, N2, *Di, Ai, incRow, Ti+incColT, incRowT);
                    fscal (F, N2, *Di, Ai+incCol, incRow, Ti, incRowT);
                    Ti+=incColT; Ai+=incCol; Di+=incD; i++;
                }

            }
            if (UpLo == FflasUpper) {
                // C12 <- alpha A1 x temp^T + beta C12
                fgemm (F, trans, oppTrans, N1, N2, K, alpha, A, lda, temp, ldt, beta, C12, ldc);
            } else {
                // C21 <- alpha temp x A11^T + beta C21
                fgemm (F, trans, oppTrans, N2, N1, K, alpha, temp, ldt, A, lda, beta, C21, ldc);
            }
            fflas_delete (temp);

            // C11 <- alpha A1 x D1 x A1^T + beta C11 and A1 <- A1 x D1
            fsyrk (F, UpLo, trans, N1, K, alpha, A, lda, D, incD, twoBlocks, beta, C, ldc, threshold);
            // C22 <- alpha A2 x D1 x A2^T + beta C22 and A2 <- A2 x D1
            fsyrk (F, UpLo, trans, N2, K, alpha, A2, lda, D, incD, twoBlocks, beta, C22, ldc, threshold);

            return C;
        }
    }
}

#endif //__FFLASFFPACK_fflas_fsyrk_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
