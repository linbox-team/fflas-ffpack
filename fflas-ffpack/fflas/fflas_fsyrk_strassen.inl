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

namespace FFLAS {

    // template<class Field>
    // inline void
    // S1S2quadrant (const Field& F,
    //               const size_t N,
    //               const size_t K,
    //               const typename Field::Element x,
    //               typename Field::ConstElement_ptr A, const size_t lda,
    //               typename Field::Element_ptr S, const size_t lds,
    //               typename Field::Element_ptr T, const size_t ldt){

    //     typename Field::ConstElement_ptr A11=A;
    //     typename Field::ConstElement_ptr A21=A+N*lda;
    //     typename Field::ConstElement_ptr A22=A21+K;

    //     size_t N2 = N>>1;
    //     size_t K2 = K>>1;
    //     for (size_t i=0; i<N2; i++, A11+=lda, A21+=lda, A22+=lda, S+=lds, T+=ldt){
    //             // @todo: vectorize this inner loop
    //         for (size_t j=0; j<K2; j++){
    //             typename Field::Element t;
    //             F.init(t);
    //             F.sub(t,A11[j],A21[j]);
    //             F.mul(S[j],t,x);
    //             F.maxpy(T[j], A21[j], x, A22[j]);
    //         }
    //     }
    // }

    // template<class Field>
    // inline void
    // computeS1S2 (const Field& F,
    //              const size_t N,
    //              const size_t K,
    //              const typename Field::Element x,
    //              const typename Field::Element y,
    //              typename Field::ConstElement_ptr A, const size_t lda,
    //              typename Field::Element_ptr S, const size_t lds,
    //              typename Field::Element_ptr T, const size_t ldt){
    //         // S1 = (A11-A21) x Y^T in S
    //         // S2 = A22 - A21 x Y^T in T
    //         // where Y = [ x.I  y.I]
    //         //           [ -y.I x.I]
    //     size_t N2 = N>>1;
    //     typename Field::Element negy;
    //     F.init(negy);
    //     F.neg(negy, y);
    //     S1S2quadrant (F, N, K, x, A, lda, S, lds, T, ldt);
    //     S1S2quadrant (F, N, K, negy, A+N2, lda, S+N2, lds, T+N2, ldt);
    //     S1S2quadrant (F, N, K, y, A+N2*lda, lda, S+N2*lds, lds, T+N2*ldt, ldt);
    //     S1S2quadrant (F, N, K, x, A+N2*lda+N2, lda, S+N2*lds+N2, lds, T+N2*ldt+N2, ldt);
    // }
    //     ///////////
    
    template<class Field, class FieldTrait>
    inline void
    computeS1S2 (const Field& F,
                 const FFLAS_TRANSPOSE trans,
                 const size_t N,
                 const size_t K,
                 const typename Field::Element x,
                 const typename Field::Element y,
                 typename Field::ConstElement_ptr A, const size_t lda,
                 typename Field::Element_ptr S, const size_t lds,
                 typename Field::Element_ptr T, const size_t ldt,
                 MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait> WH)  {
            // S1 = (A21-A11) x Y in S
            // S2 =  A22 - A21 x Y  in T
            // where Y = [ x.I  y.I]
            //           [ -y.I x.I]
        typedef MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > MMH_t;
            //typedef typename  MMH_t::DelayedField::Element DFElt;
        const typename MMH_t::DelayedField & DF = WH.delayedField;

        size_t N2 = N>>1;
        size_t K2 = K>>1;
        size_t K4 = K2>>1;
        typename Field::ConstElement_ptr A11 = A, A21, A22, A11r, A21r, A22r;
        typename Field::Element_ptr Sr, Tr;
        if (trans==FflasNoTrans){
            A21 = A + N2*lda;
            A22 = A21 + K2;
            A22r = A22 + K4;
            A11r = A11 + K4;
            A21r = A21 + K4;
            Sr = S + K4;
            Tr = T+K4;
        } else { // Trans
            A21 = A + N2;
            A22 = A21 + K2*lda;
            A22r = A22 + K4*lda;
            A11r = A11 + K4*lda;
            A21r = A21 + K4*lda;
            Sr = S + K4*lds;
            Tr = T + K4*ldt;
        }
        typename Field::Element negx, negy, y1, y2;
        F.init(negx);
        F.init(y1);
        F.init(y2);
        F.neg(negx, x);
        F.init(negy);
        F.neg(negy, y); 
        //     // T <-  A21 Y^T
        // fscal (F, N2, K2, x, A21, lda, T, ldt);
        // if (!F.isZero(y)){
        //     faxpy (F, N2, K4, y, A21r, lda, T, ldt);
        //     faxpy (F, N2, K4, negy, A21, lda, Tr, ldt);
        // }
        //     // S <- T + A11 Y^T
        //     // introduce a faxpy with 4 operands
        //       //  S <- A11 Y^T
        // fscal (F, N2, K2, x, A11, lda, S, lds);
        // if (!F.isZero(y)){
        //     faxpy (F, N2, K4, y, A11r, lda, S, lds);
        //     faxpy (F, N2, K4, negy, A11, lda, Sr, lds);
        // }
        //     // S <- S - T
        // fsubin (F, N2, K2, T, ldt, S, lds);

        //     // T <- T - A22
        // fsubin (F, N2, K2, A22, lda, T, ldt);

        if (trans==FflasNoTrans){
            F.assign(y1, y);
            F.assign(y2, negy);
                // TODO: write a distinct loop for the trans==FflasTrans case for better cache efficiency
                // S <- A21 Y
            for (size_t i=0; i<N2; ++i, A11+=lda, A11r+=lda, A21+=lda, A21r+=lda, A22+=lda, S+=lds, Sr+=lds, T+=ldt){
                fscal (DF, K2, x, A21, 1, S, 1);
                if (!F.isZero(y)){
                    faxpy (DF, K4, negy, A21r, 1, S, 1);
                    faxpy (DF, K4, y, A21, 1, Sr, 1);
                }
                    // T <- A22 -S
                fsub (DF, K2, A22, 1, S, 1, T, 1);
                
                    // S <- S - Y A11 
                faxpy (DF, K2, negx, A11, 1, S, 1);
                if (!F.isZero(y)){
                    faxpy (DF, K4, y, A11r, 1, S, 1);
                    faxpy (DF, K4, negy, A11, 1, Sr, 1);
                }
                
                freduce (F,K2,S, 1);
                freduce (F,K2,T, 1);
            }
        } else { // FflasTrans
            F.assign(y2, y);
            F.assign(y1, negy);
            
            for (size_t i=0; i<K4; ++i, A11+=lda, A11r+=lda, A21+=lda, A21r+=lda, A22+=lda, A22r+=lda, S+=lds, Sr+=lds, T+=ldt, Tr+=ldt){
                // S <- Y A21
                fscal (DF, N2, x, A21, 1, S, 1);
                fscal (DF, N2, x, A21r, 1, Sr, 1);
                if (!F.isZero(y)){
                    faxpy (DF, N2, negy, A21r, 1, S, 1);
                    faxpy (DF, N2, y1, A21, 1, Sr, 1);
                }
                    // T <- A22 -S
                fsub (DF, N2, A22, 1, S, 1, T, 1);
                fsub (DF, N2, A22r, 1, S, 1, Tr, 1);
                
                    // S <- S - A11 Y^T
                faxpy (DF, N2, negx, A11, 1, S, 1);
                faxpy (DF, N2, negx, A11r, 1, Sr, 1);
                if (!F.isZero(y)){
                    faxpy (DF, K4, y1, A11r, 1, S, 1);
                    faxpy (DF, K4, y2, A11, 1, Sr, 1);
                }
                
                freduce (F,N2,S, 1);
                freduce (F,N2,Sr, 1);
                freduce (F,N2,T, 1);
                freduce (F,N2,Tr, 1);
            }
            
        }
        // } else { // -1 is not a square in F
//         size_t K4 = K2>>1;
//         typename Field::ConstElement_ptr A11r = A11 + K4;
//             typename Field::ConstElement_ptr A21r = A21 + K4;
//             typename Field::ConstElement_ptr A22r = A22 + K4;
//             typename Field::Element_ptr Sr = S + K4;
//             typename Field::Element_ptr Tr = T + K4;
//             for (size_t i=0; i<N2; i++, A11+=lda, A21+=lda, A22+=lda, A11r+=lda, A21r+=lda, A22r+=lda, S+=lds, Sr+=lds, T+=ldt, Tr+=ldt){
//                     // @todo: vectorize this inner loop
//                 for (size_t j=0; j<K4; j++){
//                     typename Field::Element t, tr;
//                     F.init(t);
//                     F.init(tr);
//                     F.sub(t,A11[j],A21[j]);
//                     F.sub(tr,A11r[j],A21r[j]);
//                     F.mul(S[j],t,x);
//                     F.mul(Sr[j],t,negy);
//                     F.axpyin(S[j],tr,y);
//                     F.axpyin(Sr[j],tr,x);

//                     F.maxpy(T[j], A21[j], x, A22[j]);
//                     F.maxpy(Tr[j], A21[j], negy, A22r[j]);

//                     F.maxpyin(T[j], A21r[j], y);
//                     F.maxpyin(Tr[j], A21r[j], x);
//                 }
//             }
//         }
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
            //std::cerr<<"fsyrk Wino Delayed"<<std::endl;
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
        // std::cerr<<"\n Delayed -> Lazy alpha_ = "<<alpha_<<std::endl;
        // std::cerr<<" A = "<<*A<<"\n B = "<<*B<<"\n C = "<<*C<<"\n alpha, beta ="<<alpha<<" "<<beta<<std::endl;

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

            //std::cerr<<"Wino Mode, n quelconque"<<std::endl;
        size_t q = 1 << (H.recLevel+1); // 
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
        fsyrk_strassen (F, UpLo, trans, Ns, Ks, y1, y2, alpha, A, lda, beta, C, ldc, H);

            // C11 += A12 x A12 ^T
        MMHelper<Field, MMHelperAlgo::Classic, Mode>  H2 (H);
        H2.Cmin = H.Outmin;
        H2.Cmax = H.Outmax;
            //WriteMatrix (std::cerr<<"---------------"<<std::endl<<"C11_strassen = "<<std::endl, F, Ns, Ns, C, ldc);
        fsyrk (F, UpLo, trans, Ns, K-Ks, alpha, A12, lda, F.one, C, ldc, H2);
            //WriteMatrix (std::cerr<<"---------------"<<std::endl<<"C11++_strassen = "<<std::endl, F, Ns, Ns, C, ldc);

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
                    typename Field::ConstElement_ptr A, const size_t lda,
                    const typename Field::Element beta,
                    typename Field::Element_ptr C, const size_t ldc,
                    MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait> & WH
                    ) {
        // std::cerr<<"Entree fsyrk_strassen"
        //          <<" WH = "<<WH
        //          <<std::endl;
            // written for NoTrans, Lower
        if (WH.recLevel == 0){
            MMHelper<Field, MMHelperAlgo::Classic, FieldTrait>  CH(WH);
            fsyrk (F, uplo, trans, N, K, alpha, A, lda, beta, C, ldc,CH);
            WH.Outmin = CH.Outmin;
            WH.Outmax = CH.Outmax;
            return C;
        }

        typedef MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > MMH_t;
            //       typedef typename  MMH_t::DelayedField::Element_ptr DFEptr;
            //       typedef typename  MMH_t::DelayedField::ConstElement_ptr DFCEptr;
        typedef typename  MMH_t::DelayedField::Element DFElt;

        const typename MMH_t::DelayedField & DF = WH.delayedField;

        size_t N2 = N>>1;
        size_t K2 = K>>1;
        typename Field::ConstElement_ptr A11 = A, A12, A21, A22;

        if (trans == FflasNoTrans){
            A12 = A + K2; A21 = A + N2*lda; A22 = A21 + K2;
        } else {
            A12 = A + K2*lda; A21 = A + N2; A22 = A12 + N2;
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
                lds = K2;
            } else{
                S1 = C21;
                S2 = C12;
                S4 = C11;
                lds = ldc;
            }
                // S1 = (A21-A11) x Y in S1 (C21)
                // S2 =  A22 - A21 x Y  in S2 (C12)
            computeS1S2 (F, trans, N, K, y1, y2, A,lda, S1, lds, S2, lds, WH);
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"A21 = "<<std::endl, F, N2, K2, A21, lda);
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"A11 = "<<std::endl, F, N2, K2, A11, lda);
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"S1 = "<<std::endl, F, N2, K2, C21, ldc);
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"S2 = "<<std::endl, F, N2, K2, C12, ldc);
            //std::cerr<<"x = "<<y1<< " y = "<<y2<<std::endl;
                //  P4^T =  S2 x S1^T in  C22
            MMH_t H4 (F, -1, WH.Amin, WH.Amax, WH.Bmin, WH.Bmax, 0,0);
            if (uplo == FflasLower)
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, S2, lds, S1, lds, F.zero, C22, ldc, H4);
            else
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, S1, lds, S2, lds, F.zero, C22, ldc, H4);
                
            // std::cerr<<"alpha = "<<alpha<<std::endl;
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"P4^T = "<<std::endl, F, N2, N2, C22, ldc);

            if (K>N){ fflas_delete(S2); }

                // S3 = S1 - A22 in S1
            fsubin (DF, N2, K2, A22, lda, S1, lds);
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"S3 = "<<std::endl, F, N2, K2, S1, lds);

                // P5 = S3 x S3^T in C12
            MMH_t H5 (F, WH.recLevel-1, 2*WH.Amin, 2*WH.Amax, 2*WH.Bmin, 2*WH.Bmax, 0,0);
            fsyrk_strassen (F, uplo, trans, N2, K2, y1, y2, alpha, S1, lds, F.zero, C12, ldc, H5);
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"P5 = "<<std::endl, F, N2, N2, C12, ldc);

                // S4 = S3 + A12 in S4
            fadd (DF, N2, K2, S1, lds, A12, lda, S4, lds);
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"S4 = "<<std::endl, F, N2, K2, S4, lds);

                // P3 = A22 x S4^T in C21
            MMH_t H3 (F, WH.recLevel-1, WH.Amin, WH.Amax, 2*WH.Bmin-WH.Bmax, 2*WH.Bmax-WH.Bmin, 0,0);
            if (uplo == FflasLower)
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, A22, lda, S4, lds, F.zero, C21, ldc, H3);
            else{
                H3.Amin = 2*WH.Bmin-WH.Bmax; 
                H3.Amax = 2*WH.Bmax-WH.Bmin; 
                H3.Bmin = WH.Amin;
                H3.Bmax = WH.Amax;
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, S4, lds, A22, lda, F.zero, C21, ldc, H3);
            }
// WriteMatrix (std::cerr<<"---------------"<<std::endl<<"P3 = "<<std::endl, F, N2, N2, C21, ldc);

            if (K>N){ fflas_delete(S1); }

                // P1 = A11 x A11^T in C11
            MMH_t H1 (F, WH.recLevel-1, WH.Amin, WH.Amax, WH.Bmin, WH.Bmax, 0,0);
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"A11 = "<<std::endl, F, N2, K2, A11, lda);
            fsyrk_strassen (F, uplo, trans, N2, K2, y1, y2, alpha, A11, lda, F.zero, C11, ldc, H1);
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"P1 = "<<std::endl, F, N2, N2, C11, ldc);

                // U1 = P1 + P5 in C12
            DFElt U1Min, U1Max;
            if (Protected::NeedPreAddReduction (U1Min, U1Max, H1.Outmin, H1.Outmax, H5.Outmin, H5.Outmax, WH)){
                freduce(F,N2,N2,C12,ldc);
                freduce(F,N2,N2,C11,ldc);
            }
           // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"after reduce P1 = "<<std::endl, F, N2, N2, C11, ldc);
            faddin (DF, uplo, N2, C11, ldc, C12, ldc); // TODO triangular addin (to be implemented)
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"U1 = "<<std::endl, F, N2, N2, C12, ldc);

                // make U1 explicit: Up(U1)=Low(U1)^T
            if (uplo == FflasLower)
                for (size_t i=0; i<N2; i++)
                    fassign(F, i, C12+i*ldc, 1, C12+i, ldc);
            else
                for (size_t i=0; i<N2; i++)
                    fassign(F, i, C12+i, ldc, C12+i*ldc, 1);
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"U1 = "<<std::endl, F, N2, N2, C12, ldc);

                // U2 = U1 + P4 in C12
            DFElt U2Min, U2Max;
            if (Protected::NeedPreAddReduction (U2Min, U2Max, U1Min, U1Max, H4.Outmin, H4.Outmax, WH)){
                freduce(F,N2,N2,C12,ldc);
                freduce(F,N2,N2,C22,ldc);
            }
            for (size_t i=0; i<N2; ++i){
                faddin (DF,  N2, C22+i*ldc, 1, C12+i, ldc);
                // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"WIP U2 = "<<std::endl, F, N2, N2, C12, ldc);
            }
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"U2 = "<<std::endl, F, N2, N2, C12, ldc);

                // U4 = U2 + P3 in C21
            DFElt U4Min, U4Max;
            if (Protected::NeedPreAddReduction (U4Min, U4Max, U2Min, U2Max, H3.Outmin, H3.Outmax, WH)){
                freduce(F,N2,N2,C21,ldc);
                freduce(F,N2,N2,C12,ldc);
            }
            faddin (DF, N2, N2, C12, ldc, C21, ldc);
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"U4 = "<<std::endl, F, N2, N2, C21, ldc);

                // U5 = U2 + P4^T in C22
            DFElt U5Min, U5Max;
            if (Protected::NeedPreAddReduction (U5Min, U5Max, U2Min, U2Max, H4.Outmin, H4.Outmax, WH)){
                freduce(F,N2,N2,C22,ldc);
                freduce(F,N2,N2,C12,ldc);
            }
            faddin (DF, uplo, N2, C12, ldc, C22, ldc);
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"U5 = "<<std::endl, F, N2, N2, C22, ldc);

                // P2 = A12 x A12^T in C12
            MMH_t H2 (F, WH.recLevel-1, WH.Amin, WH.Amax, WH.Bmin, WH.Bmax, 0,0);
            fsyrk_strassen (F, uplo, trans, N2, K2, y1, y2, alpha, A12, lda, F.zero , C12, ldc, H2);
            // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"P2 = "<<std::endl, F, N2, N2, C12, ldc);

                // U3 = P1 + P2 in C11
            DFElt U3Min, U3Max;
            if (Protected::NeedPreAddReduction (U3Min, U3Max, H1.Outmin, H1.Outmax, H2.Outmin, H2.Outmax, WH)){
                freduce(F,N2,N2,C11,ldc);
                freduce(F,N2,N2,C12,ldc);
            }
            faddin (DF, uplo, N2, C12, ldc, C11, ldc);
           // WriteMatrix (std::cerr<<"---------------"<<std::endl<<"U3 = "<<std::endl, F, N2, N2, C11, ldc);
                // Updating WH with Outmin, Outmax of the result
            WH.Outmin = min3 (U3Min, U4Min, U5Min);
            WH.Outmax = max3 (U3Max, U4Max, U5Max);
        } else { // with accumulation, schedule with 1 temp

            typename Field::Element negbeta;
            F.init (negbeta);
            F.neg (negbeta, beta);
            typename Field::Element_ptr T = fflas_new (F, N2, std::max(N2,K2));
            size_t ldt = std::max(N2,K2);

            if (K>N){ // Not possible to store S2 in C12
                S2 = fflas_new(F, N2, K2);
                lds = K2;
            } else{
                S2 = C12;
                lds = ldc;
            }
                // S1 = (A21-A11) x Y^T in T1
                // S2 = A22 - A21 x Y^T in C12
            computeS1S2 (F, trans, N, K, y1, y2, A, lda, T, ldt, S2, lds, WH);

                // Up(C11) = Low(C22) (saving C22)
            if (uplo == FflasLower)
                for (size_t i=0; i<N2-1; ++i)
                    fassign (F, N2-i-1, C22 + (N2-i-1)*ldc, 1, C11 + 1 + i*(ldc+1), 1);
            else 
                for (size_t i=0; i<N2-1; ++i)
                    fassign (F, N2-i-1, C22 + 1 + i*(ldc+1), 1,  C11 + (N2-i-1)*ldc, 1);

                // temp for storing the diagonal of C22
            typename Field::Element_ptr DC22 = fflas_new(F,N2);
            for (size_t i=0; i<N2; ++i)
                F.assign (DC22[i], C22[i*(ldc+1)]);

                // P4^T = S2 x S1^T in C22
            MMH_t H4 (F,WH.recLevel-1, WH.Amin, WH.Amax, WH.Bmin, WH.Bmax, 0,0);
            if (uplo == FflasLower)
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, S2, lds, T, ldt, F.zero, C22, ldc, H4);
            else
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, T, ldt, S2, lds, F.zero, C22, ldc, H4);

            if (K>N){ fflas_delete(S2); }

                // S3 = S1 - A22 in T
            fsubin (DF, N2, K2, A22, lda, T, ldt);

                // P5 = S3 x S3^T in C12
            MMH_t H5 (F, WH.recLevel-1, WH.Amin-WH.Amax, WH.Amax-WH.Amin, WH.Bmin-WH.Bmax, WH.Bmax-WH.Bmin, 0,0);
            fsyrk_strassen (F, uplo, trans, N2, K2, y1, y2, alpha, T, ldt, F.zero, C12, ldc, H5);

                // S4 = S3 + A12 in T1
            faddin (DF, N2, K2, A12, lda, T, ldt);

                //  P3 = A22 x S4^T + beta C21 in C21
            MMH_t H3 (F, WH.recLevel-1, WH.Amin, WH.Amax, 2*WH.Bmin-WH.Bmax, 2*WH.Bmax-WH.Bmin, WH.Cmin, WH.Cmax);
            if (uplo == FflasLower)
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, A22, lda, T, ldt, beta, C21, ldc, H3);
            else{
                H3.Amin = 2*WH.Bmin-WH.Bmax; 
                H3.Amax = 2*WH.Bmax-WH.Bmin; 
                H3.Bmin = WH.Amin;
                H3.Bmax = WH.Amax;
                fgemm (F, trans, OppTrans, N2, N2, K2, alpha, T, ldt, A22, lda, beta, C21, ldc, H3);
            }
                // P1 = A11 x A11^T in T1
            MMH_t H1 (F, WH.recLevel-1, WH.Amin, WH.Amax, WH.Bmin, WH.Bmax, 0, 0);
            fsyrk_strassen (F, uplo, trans, N2, K2, y1, y2, alpha, A11, lda, F.zero, T, ldt, H1);
                //WriteMatrix (std::cerr<<"---------------"<<std::endl<<"P1 = "<<std::endl, F, N2, N2, T, ldt);

                // U1 = P5 + P1  in C12 // Still symmetric
            DFElt U1Min, U1Max;
            if (Protected::NeedPreAddReduction (U1Min, U1Max, H5.Outmin, H5.Outmax, H1.Outmin, H1.Outmax, WH)){
                freduce(F,N2,N2,C12,ldc);
                freduce(F,N2,N2,T,ldt);
           }
            faddin (DF, uplo, N2, T, ldt, C12, ldc);

                // Make U1 explicit (copy the N/2 missing part)
            if (uplo == FflasLower)
                for (size_t i=0; i<N2; ++i)
                    fassign (DF, i, C12 + i*ldc, 1, C12 + i, ldc);
            else
                for (size_t i=0; i<N2; ++i)
                    fassign (DF, i, C12 + i, ldc, C12 + i*ldc, 1);                

                // U2 = U1 + P4 in C12
            DFElt U2Min, U2Max;
            if (Protected::NeedPreAddReduction (U2Min, U2Max, U1Min, U1Max, H4.Outmin, H4.Outmax, WH)){
                freduce(F,N2,N2,C12,ldc);
                freduce(F,N2,N2,C22,ldc);
            }
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
                freduce(F,N2,N2,C22,ldc);
                freduce(F,N2,N2,C12,ldc);
            }
            faddin (DF, uplo, N2, C12, ldc, C22, ldc);
            freduce(F,N2,N2,C22,ldc);

                // U5' = U5 +  beta Up(C11)^T in C22
                // TODO use delayed field and a needPreAXPYReduction
                // Suspicious BUG: if C22 is not reduced, possible overflow???
            if (uplo == FflasLower)
                for (size_t i=0; i<N2-1; i++){ // TODO factorize out in a triple add
                    faxpy (F, N2-i-1, beta, C11 + 1+i*(ldc+1), 1, C22 + (N2-i-1)*ldc, 1);
                    F.axpyin (C22[i*(ldc+1)], beta, DC22[i]);
                }
            else
                for (size_t i=0; i<N2-1; i++){ // TODO factorize out in a triple add
                    faxpy (F, N2-i-1, beta, C11 + (N2-i-1)*ldc, 1, C22 + 1 + i*(ldc+1), 1);
                    F.axpyin (C22[i*(ldc+1)], beta, DC22[i]);
                }

            F.axpyin (C22[(N2-1)*(ldc+1)], beta, DC22[N2-1]);

            fflas_delete(DC22);

                // P2 = A12 x A12^T + beta C11 in C11
            MMH_t H2 (F, WH.recLevel-1, WH.Amin, WH.Amax, WH.Bmin, WH.Bmax, WH.Cmin, WH.Cmax);
            fsyrk_strassen (F, uplo, trans, N2, K2, y1, y2, alpha, A12, lda, beta, C11, ldc, H2);
                //WriteMatrix (std::cerr<<"---------------"<<std::endl<<"P2 = "<<std::endl, F, N2, N2, C11, ldc);

                // U3 = P1 + P2 in C11
            DFElt U3Min, U3Max;
            if (Protected::NeedPreAddReduction (U3Min, U3Max, H1.Outmin, H1.Outmax, H2.Outmin, H2.Outmax, WH)){
                freduce(F,N2,N2,C11,ldc);
                freduce(F,N2,N2,T,ldt);
            }
            faddin (DF, uplo, N2, T, ldt, C11, ldc);
                //WriteMatrix (std::cerr<<"---------------"<<std::endl<<"U3 = "<<std::endl, F, N2, N2, C11, ldc);

            fflas_delete (T);
            WH.Outmin = min3 (U3Min, U4Min, U5Min);
            WH.Outmax = max3 (U3Max, U4Max, U5Max);
        }
         return C;
    }
//============================
            // version with accumulation and 3 temps
            // S1 = A11 - A21 in C12

            // A22' = A22 Y in T1

            // S2 = A21 + A22' in T2

            // P4 = S1 x S2^T in  T3

            // P1 = A11 x A11^T in C12

            // C11 = P1 + beta.C11 in C11
        
            // U3 = P2 + P1 = A12 x A12^T + P1 in C11

            // S3 = A11 - S2 in T2
        
            // U1 = P1 - P5 = -S3 x S3^T + P1 in C12

            // U2 = U1 - P4 in C12

            // U5 = U2 - P4^T + C22 in C22 (add)

            // S4 = S3  + A12 x Y in T2
            // A12xY on the fly

            // - P3 = - A22 Y x S4^T + C21 in C21

            // U4 = U2 - P3  in C21

////////////////////////////
            // version without accumulation but requires 1 temp


        
            // S1 = A11 - A21 in C12

            // A22' = A22 Y in C11

            // S2 = A21 + A22' in C21

            // P4 = S1 x S2^T in  C22

            // S3 = A11 - S2 in C21

            // -P5 = -S3 x S3^T  in T

            // S4 = S3 x Y + A12 in C12

            // -P3 = -A22' x S4^T in C21

            // P1 = A11 x A11^T in C11

            // U1 = P1 - P5 in T

            // U2 = U1 - P4 in T

            // U4 = U2 - P3 in C21

            // U5 = U2 - P4^T in C22
        
            // P2 = A12 x A12^T in C12

            // U3 = P1 + P2 in C11

//        ======================
            // algo du papier sans acc mais pas en place
            // S1 = A11 - A21

            // S2 = A21 + A22 x Y^T

            // S3 = A11 - S2

            // S4 = S3 x Y + A12

            // P1 = A11 x A11^T in C11

            // P2 = A12 x A12^T

            // P3 = A22 Y x S4^T

            // P4 = S1 x S2^T

            // P5 = S3 x S3^T

            // U1 = P1 - P5

            // U2 = U1 - P4

            // U3 = P1 + P2 in C11

            // U4 = U2 - P3 in C21

            // U5 = U2 - P4^T in C22
        
}


#endif // __FFLASFFPACK_fflas_fsyrk_strassen_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
