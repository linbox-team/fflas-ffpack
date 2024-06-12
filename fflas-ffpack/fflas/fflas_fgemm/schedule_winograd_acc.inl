/* Copyright (C) 2014 the LinBox group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

/** @file fflas/fflas_fgemm/winograd_acc.inl
 * @ingroup MMalgos
 * @brief Winograd implementation
 * @bib ISSAC09 Scheduling
 */

#ifndef __FFLASFFPACK_fgemm_winograd_acc_INL
#define __FFLASFFPACK_fgemm_winograd_acc_INL

namespace FFLAS { namespace BLAS3 {


    // 3 temps and 23 ops
    // TODO: Add check for modular reductions before final additions
    template < class Field,class FieldTrait >
    inline void WinogradAcc_3_23 (const Field& F,
                                  const FFLAS_TRANSPOSE ta,
                                  const FFLAS_TRANSPOSE tb,
                                  const size_t mr, const size_t nr, const size_t kr,
                                  const typename Field::Element alpha,
                                  typename Field::ConstElement_ptr A,const size_t lda,
                                  typename Field::ConstElement_ptr B,const size_t ldb,
                                  const typename Field::Element  beta,
                                  typename Field::Element_ptr C, const size_t ldc,
                                  MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > & WH
                                 )
    {
        MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > H = WH ;
        H.recLevel = H.recLevel - 1 ;

        FFLASFFPACK_check(!F.isZero(beta));

        typename Field::Element mbeta  ;
        F.neg(mbeta,beta);

        size_t lb, cb, la, ca;
        size_t x3rd = std::max(mr,kr);
        typename Field::ConstElement_ptr A11=A, A12, A21, A22;
        typename Field::ConstElement_ptr B11=B, B12, B21, B22;
        typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;



        size_t ldX3;
        // Three temporary submatrices are required

        if (ta == FflasTrans) {
            A21 = A + mr;
            A12 = A + kr*lda;
            A22 = A12 + mr;
            la = kr;
            ca = mr;
        }
        else { // ta == FflasNoTrans
            A12 = A + kr;
            A21 = A + mr*lda;
            A22 = A21 + kr;
            la = mr;
            ca = kr;
        }
        if (tb == FflasTrans) {
            B21 = B + kr;
            B12 = B + nr*ldb;
            B22 = B12 + kr;
            lb = nr;
            cb = kr;
            ldX3 = x3rd;
        }
        else { // ta == FflasNoTrans
            B12 = B + nr;
            B21 = B + kr*ldb;
            B22 = B21 + nr;
            lb = kr;
            ldX3 = cb = nr;
        }

        // P2 = alpha . A12 * B21 + beta . C11  in C11
        fgemm (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, beta, C11, ldc, H);

        typename Field::Element_ptr X3 = fflas_new (F, x3rd, nr);

        // T3 = B22 - B12 in X3
        fsub(F,lb,cb,B22,ldb,B12,ldb,X3,ldX3);

        typename Field::Element_ptr X2 = fflas_new (F, mr, kr);

        // S3 = A11 - A21 in X2
        fsub(F,la,ca,A11,lda,A21,lda,X2,ca);

        // C22 = C22 - C12 if beta != 0
        fsubin(F,mr,nr,C12,ldc,C22,ldc);

        // C21 = C21 - C22
        fsubin(F,mr,nr,C22,ldc,C21,ldc);

        // P7 = alpha . S3 * T3 + beta . C22 in C22
        fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, X3, ldX3, beta, C22, ldc, H);

        // T1 = B12 - B11 in X3
        fsub(F,lb,cb,B12,ldb,B11,ldb,X3,ldX3);

        // S1 = A21 + A22 in X2
        fadd(F,la,ca,A21,lda,A22,lda,X2,ca);

        // P5 = alpha . S1*T1 + beta . C12 in C12
        fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, X3, ldX3, beta, C12, ldc, H);

        // T2 = B22 - T1 in X3
        fsub(F,lb,cb,B22,ldb,X3,ldX3,X3,ldX3);

        // S2 = S1 - A11 in X2
        fsubin(F,la,ca,A11,lda,X2,ca);

        typename Field::Element_ptr X1 = fflas_new (F, mr, nr);

        // P6 = alpha . S2 * T2 in X1
        fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, X3, ldX3, F.zero, X1, nr, H);

        // T4 = T2 - B21 in X3
        fsubin(F,lb,cb,B21,ldb,X3,ldX3);

        // S4 = A12 -S2 in X2
        fsub(F,la,ca,A12,lda,X2,ca,X2,ca);

        // P4 = alpha . A22 * T4 - beta . C21 in C21
        fgemm (F, ta, tb, mr, nr, kr, alpha, A22, lda, X3, ldX3, mbeta, C21, ldc, H);

        // P1 = alpha . A11 * B11 in X3
        fgemm (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X3, nr, H);

        //  U1 = P2 + P1 in C11
        faddin(F,mr,nr,X3,nr,C11,ldc);

        // U2 = P1 + P6 in tmpU2/X1  and
        faddin(F, mr, nr, X3, nr, X1, nr);

        // U3 = P7 + U2 in tmpU3/X3  and
        fadd(F, mr, nr, X1, nr, C22, ldc, X3, nr);

        // U7 = P5 + U3 in C22    and
        fadd(F, mr, nr, C12, ldc, X3, nr, C22, ldc);

        // U4 = P5 + U2 in C12    and
        faddin(F, mr, nr, X1, nr, C12, ldc);

        fflas_delete (X1);

        // U6 = U3 - P4 in C21    and
        fsub(F, mr, nr, X3, nr, C21, ldc, C21, ldc);

        fflas_delete (X3);

        // P3 = alpha . S4*B22 in X1
        fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, B22, ldb, F.one, C12, ldc, H);

        fflas_delete (X2);

    } // WinogradAccOld

    // 3 temps and 21 ops
    template < class Field, class FieldTrait>
    inline void WinogradAcc_3_21 (const Field& F,
                                  const FFLAS_TRANSPOSE ta,
                                  const FFLAS_TRANSPOSE tb,
                                  const size_t mr, const size_t nr, const size_t kr,
                                  const typename Field::Element alpha,
                                  typename Field::ConstElement_ptr A,const size_t lda,
                                  typename Field::ConstElement_ptr B,const size_t ldb,
                                  const typename Field::Element  beta,
                                  typename Field::Element_ptr C, const size_t ldc,
                                  MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > & WH
                                 )
    {
        typedef MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > MMH_t;
        typedef typename  MMH_t::DelayedField::Element_ptr DFEptr;
        typedef typename  MMH_t::DelayedField::ConstElement_ptr DFCEptr;
        typedef typename  MMH_t::DelayedField::Element DFElt;

        const typename MMH_t::DelayedField & DF = WH.delayedField;

        FFLASFFPACK_check(!DF.isZero(beta));

        size_t lb, cb, la, ca;
        size_t x3rd = std::max(mr,kr);
        typename Field::ConstElement_ptr A11=A, A12, A21, A22;
        typename Field::ConstElement_ptr B11=B, B12, B21, B22;
        typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;

        typename Field::Element mbeta;
        F.neg(mbeta,beta);
        DFElt betadf;
        DF.init(betadf);
        if (F.isMOne(beta)) {
            DF.assign(betadf, DF.mOne);
        } else {
            DF.assign(betadf, beta);
        }

        size_t ldX3;

        if (ta == FflasTrans) {
            A21 = A + mr;
            A12 = A + kr*lda;
            A22 = A12 + mr;
            la = kr;
            ca = mr;
        } else { // ta == FflasNoTrans
            A12 = A + kr;
            A21 = A + mr*lda;
            A22 = A21 + kr;
            la = mr;
            ca = kr;
        }
        if (tb == FflasTrans) {
            B21 = B + kr;
            B12 = B + nr*ldb;
            B22 = B12 + kr;
            lb = nr;
            cb = kr;
            ldX3 = x3rd;
        } else { // ta == FflasNoTrans
            B12 = B + nr;
            B21 = B + kr*ldb;
            B22 = B21 + nr;
            lb = kr;
            ldX3 = cb = nr;
        }

        // Three temporary submatrices are required
        typename Field::Element_ptr X3 = fflas_new (F, x3rd, nr);

        // T1 = B12 - B11 in X3
        fsub(DF,lb,cb,(DFCEptr)B12,ldb,(DFCEptr)B11,ldb,(DFEptr)X3,ldX3);

        typename Field::Element_ptr X2 = fflas_new(F,mr,kr);

        // S1 = A21 + A22 in X2
        fadd(DF,la,ca,(DFCEptr)A21,lda,(DFCEptr)A22,lda,(DFEptr)X2,ca);

        typename Field::Element_ptr X1 = fflas_new(F,mr,nr);
        // P5 = alpha . S1*T1  in X1
        MMH_t H5(F, WH.recLevel-1,
                 2*WH.Amin, 2*WH.Amax,
                 -(WH.Bmax-WH.Bmin),
                 WH.Bmax-WH.Bmin,
                 0, 0);
        fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, X3, ldX3, F.zero, X1, nr, H5);

        DFElt C22Min, C22Max;
        DFElt C12Min, C12Max;
        // This test will be optimized out
        if (Protected::NeedDoublePreAddReduction (C12Min, C12Max, H5.Outmin, H5.Outmax, WH.Cmin, WH.Cmax, betadf, WH)){
            freduce(F,mr,nr,X1,nr);
            H5.initOut();
        }
        C22Min = C12Min; C22Max = C12Max;

        // C22 = P5 + beta C22 in C22
        fadd(DF,mr,nr,(DFCEptr)X1,nr,betadf,(DFCEptr)C22,ldc,(DFEptr)C22,ldc);

        // C12 = P5 + beta C12 in C12
        fadd(DF,mr,nr,(DFCEptr)X1,nr,betadf,(DFCEptr)C12,ldc,(DFEptr)C12,ldc);

        // P1 = alpha . A11 * B11 in X1
        MMH_t H1(F, WH.recLevel-1,
                 WH.Amin, WH.Amax,
                 WH.Bmin, WH.Bmax,
                 0, 0);
        fgemm (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X1, nr, H1);

        // P2 = alpha . A12 * B21 + beta . C11  in C11
        MMH_t H2(F, WH.recLevel-1,
                 WH.Amin, WH.Amax,
                 WH.Bmin, WH.Bmax,
                 WH.Cmin, WH.Cmax);
        fgemm (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, beta, C11, ldc, H2);

        //  U1 = P2 + P1 in C11
        DFElt U1Min, U1Max;
        if (Protected::NeedPreAddReduction (U1Min,U1Max, H1.Outmin, H1.Outmax, H2.Outmin,H2.Outmax, WH) ){
            freduce(F,mr,nr,X1,nr);
            freduce(F,mr,nr,C11,ldc);
        }
        faddin(DF,mr,nr,(DFCEptr)X1,nr,(DFEptr)C11,ldc);

        // T2 = B22 - T1 in X3
        fsub(DF,lb,cb,(DFCEptr)B22,ldb,(DFCEptr)X3,ldX3,(DFEptr)X3,ldX3);

        // S2 = S1 - A11 in X2
        fsubin(DF,la,ca,(DFCEptr)A11,lda,(DFEptr)X2,ca);

        // U2 = P6 + P1 = alpha . S2 * T2 + P1 in X1
        MMH_t H6(F, WH.recLevel-1,
                 2*WH.Amin-WH.Amax, 2*WH.Amax-WH.Amin,
                 2*WH.Bmin-WH.Bmax, 2*WH.Bmax-WH.Bmin,
                 H1.Outmin, H1.Outmax);
        fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, X3, ldX3, F.one, X1, nr, H6);

        // U4 = U2 + C12 in C12
        DFElt U4Min, U4Max;
        if (Protected::NeedPreAddReduction (U4Min, U4Max, H6.Outmin, H6.Outmax, C12Min, C12Max, WH)){
            freduce(F,mr,nr,C12,ldc);
            freduce(F,mr,nr,X1,nr);
        }
        faddin(DF,mr,nr,(DFCEptr)X1,nr,(DFEptr)C12,ldc);

        // T4 = T2 - B21 in X3
        fsubin(DF,lb,cb,(DFCEptr)B21,ldb,(DFEptr)X3,ldX3);

        // S4 = A12 -S2 in X2
        fsub(DF,la,ca,(DFCEptr)A12,lda,(DFCEptr)X2,ca,(DFEptr)X2,ca);

        // P4 = alpha . A22 * T4 - beta . C21 in C21
        MMH_t H4(F, WH.recLevel-1,
                 WH.Amin, WH.Amax,
                 2*WH.Bmin-2*WH.Bmax, 2*WH.Bmax-2*WH.Bmin,
                 WH.Cmin, WH.Cmax);
        fgemm (F, ta, tb, mr, nr, kr, alpha, A22, lda, X3, ldX3, mbeta, C21, ldc, H4);

        // U5 = P3 + U4 = alpha . S4*B22 + U4 in C12
        MMH_t H3(F, WH.recLevel-1,
                 2*WH.Amin-2*WH.Amax, 2*WH.Amax-2*WH.Amin,
                 WH.Bmin, WH.Bmax,
                 U4Min, U4Max);
        fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, B22, ldb, F.one, C12, ldc, H3);

        // T3 = B22 - B12 in X3
        fsub(DF,lb,cb,(DFCEptr)B22,ldb,(DFCEptr)B12,ldb,(DFEptr)X3,ldX3);

        // S3 = A11 - A21 in X2
        fsub(DF,la,ca,(DFCEptr)A11,lda,(DFCEptr)A21,lda,(DFEptr)X2,ca);

        // U3 = P7 + U2  = alpha . S3 * T3 + U2 in X1
        MMH_t H7(F, WH.recLevel-1,
                 WH.Amin-WH.Amax, WH.Amax-WH.Amin,
                 WH.Bmin-WH.Bmax, WH.Bmax-WH.Bmin,
                 H6.Outmin, H6.Outmax);
        fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, X3, ldX3, F.one, X1, nr, H7);

        fflas_delete (X2);
        fflas_delete (X3);

        // U7 =  U3 + C22 in C22
        DFElt U7Min, U7Max;
        if (Protected::NeedPreAddReduction (U7Min, U7Max, H7.Outmin, H7.Outmax, C22Min, C22Max, WH)){
            freduce(F,mr,nr,X1,nr);
            freduce(F,mr,nr,C22,ldc);
        }
        faddin(DF,mr,nr,(DFCEptr)X1,nr,(DFEptr)C22,ldc);

        // U6 = U3 - P4 in C21
        DFElt U6Min, U6Max;
        if (Protected::NeedPreSubReduction(U6Min, U6Max, H7.Outmin, H7.Outmax, H4.Outmin, H4.Outmax, WH)){
            freduce(F,mr,nr,X1,nr);
            freduce(F,mr,nr,C21,ldc);
        }
        fsub(DF,mr,nr,(DFCEptr)X1,nr,(DFCEptr)C21,ldc,(DFEptr)C21,ldc);

        fflas_delete (X1);

        // Updating WH with Outmin, Outmax of the result
        WH.Outmin = min4 (U1Min, H3.Outmin, U6Min, U7Min);
        WH.Outmax = max4 (U1Max, H3.Outmax, U6Max, U7Max);
    } // WinogradAcc

    // 3 temps and 21 ops
    template < class Field, class FieldTrait>
    inline void WinogradAcc_IP (const Field& F,
                                const FFLAS_TRANSPOSE ta,
                                const FFLAS_TRANSPOSE tb,
                                const size_t mr, const size_t nr, const size_t kr,
                                const typename Field::Element alpha,
                                typename Field::ConstElement_ptr A,const size_t lda,
                                typename Field::ConstElement_ptr B,const size_t ldb,
                                const typename Field::Element  beta,
                                typename Field::Element_ptr C, const size_t ldc,
                                MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > & WH
            )
    {
        typedef MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > MMH_t;
        typedef typename  MMH_t::DelayedField::Element_ptr DFEptr;
        typedef typename  MMH_t::DelayedField::ConstElement_ptr DFCEptr;
        typedef typename  MMH_t::DelayedField::Element DFElt;

        const typename MMH_t::DelayedField & DF = WH.delayedField;

        FFLASFFPACK_check(!DF.isZero(beta));

        size_t lb, cb, la, ca;
        typename Field::ConstElement_ptr A11=A, A12, A21, A22;
        typename Field::ConstElement_ptr B11=B, B12, B21, B22;
        typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;

        typename Field::Element mbeta,malpha;
        F.neg(mbeta,beta);
        F.neg(malpha,alpha);
        DFElt betadf,mbetadf;
        DFElt alphadf,malphadf;
        DF.init(betadf);
        DF.init(mbetadf);
        if (F.isMOne(beta)) {
            DF.assign(betadf, DF.mOne);
            DF.assign(mbetadf, DF.one);
        } else if (F.isOne(beta)){
            DF.assign(betadf,DF.one);    
            DF.assign(mbetadf,DF.mOne);    
        }else{
            DF.assign (betadf, beta);
            DF.neg (mbetadf, beta);
        }
        if (F.isMOne(alpha)) {
            DF.assign(alphadf, DF.mOne);
            DF.assign(malphadf, DF.one);
        } else if (F.isOne(alpha)){
            DF.assign(alphadf,DF.one);    
            DF.assign(malphadf,DF.mOne);    
        }else{
            DF.assign (alphadf, beta);
            DF.neg (malphadf, beta);
        }

        if (ta == FflasTrans) {
            A21 = A + mr;
            A12 = A + kr*lda;
            A22 = A12 + mr;
            la = kr;
            ca = mr;
        } else { // ta == FflasNoTrans
            A12 = A + kr;
            A21 = A + mr*lda;
            A22 = A21 + kr;
            la = mr;
            ca = kr;
        }
        if (tb == FflasTrans) {
            B21 = B + kr;
            B12 = B + nr*ldb;
            B22 = B12 + kr;
            lb = nr;
            cb = kr;
        } else { // ta == FflasNoTrans
            B12 = B + nr;
            B21 = B + kr*ldb;
            B22 = B21 + nr;
            lb = kr;
            cb = nr;
        }
        // std::cerr<<std::endl;
        // WriteMatrix(std::cerr<<"A = "<<std::endl, F, 2*mr, 2*kr, A, lda);
        // WriteMatrix(std::cerr<<"B = "<<std::endl, F, 2*kr, 2*nr, B, ldb);
        // WriteMatrix(std::cerr<<"C = "<<std::endl, F, 2*mr, 2*nr, C, ldc);

        // WriteMatrix(std::cerr<<"A = ", F, 2*mr, 2*kr, A, lda, FFLAS::FflasSageMath);
        // WriteMatrix(std::cerr<<"B = ", F, 2*kr, 2*nr, B, ldb, FFLAS::FflasSageMath);
        // WriteMatrix(std::cerr<<"C = ", F, 2*mr, 2*nr, C, ldc, FFLAS::FflasSageMath);
        
        // S1 = A21 - A11 in A21
        fsubin (DF, la, ca, (DFCEptr) A11, lda, (DFEptr) A21, lda);
        // WriteMatrix(std::cerr<<"S1 = A21-A11"<<std::endl, F, mr, kr, A21, lda);
        // T1 = B12 - B22 in B12
        fsubin (DF, lb, cb, (DFCEptr) B22, ldb, (DFEptr) B12, ldb);
        // WriteMatrix(std::cerr<<"T1 = B12-B22"<<std::endl, F, kr, nr, B12, ldb);
        // C21 = C21 - C22 in C21
        fsubin (DF, mr, nr, (DFCEptr) C22, ldc, (DFEptr) C21, ldc);
        // WriteMatrix(std::cerr<<"C21 = C21-C22"<<std::endl, F, mr, nr, C21, ldc);
        // P1 = alpha . S1*T1 + beta C22  in C22
        MMH_t H1(F, WH.recLevel-1,
                 -(WH.Amax-WH.Amin), WH.Amax-WH.Amin,
                 -(WH.Bmax-WH.Bmin), WH.Bmax-WH.Bmin,
                 -(WH.Cmax-WH.Cmin), WH.Cmax-WH.Cmin );
        fgemm (F, ta, tb, mr, nr, kr, alpha, A21, lda, B12, ldb, beta, C22, ldc, H1);
        // WriteMatrix(std::cerr<<"P1 = S1*T1 + C22"<<std::endl, F, mr, nr, C22, ldc);

        // S2 = S1 + A22 in A21
        faddin (DF, la, ca, (DFCEptr) A22, lda, (DFEptr) A21, lda);
        // WriteMatrix(std::cerr<<"S1="<<std::endl, F, mr, kr, A21, lda);
        // WriteMatrix(std::cerr<<"A22="<<std::endl, F, mr, kr, A22, lda);
        // WriteMatrix(std::cerr<<"S2 = S1+A22"<<std::endl, F, mr, kr, A21, lda);
        // T2 = T1 - B11 in B12
        fsubin (DF, lb, cb, (DFCEptr) B11, ldb, (DFEptr) B12, ldb);
        // WriteMatrix(std::cerr<<"T2 = T1-B11"<<std::endl, F, kr, nr, B12, ldb);

        //TODO
        // DFElt C22Min, C22Max;
        // DFElt C12Min, C12Max;
        // // This test will be optimized out
        //DFElt U1Min, U1Max;
        //if (Protected::NeedPreAddReduction (U1Min,U1Max, H1.Outmin, H1.Outmax, H2.Outmin,H2.Outmax, WH) ){
        //     freduce(F,mr,nr,X1,nr);
        //     freduce(F,mr,nr,C11,ldc);
        // }

        // C12 = C22 - beta C12  in C12
        fadd (DF,mr,nr,(DFCEptr)C22,ldc,mbetadf,(DFCEptr)C12,ldc,(DFEptr)C12,ldc);
        // WriteMatrix(std::cerr<<"C12 = C22 - beta C12"<<std::endl, F, mr, nr, C12, ldc);
      // P2 = - alpha . S2 * T2 +  C22  in C22
        MMH_t H2(F, WH.recLevel-1,
                 2*WH.Amin-WH.Amax, 2*WH.Amax-WH.Amin,
                 WH.Bmin-2*WH.Bmax, WH.Bmax-2*WH.Bmin,
                 H1.Outmin, H1.Outmax);
        fgemm (F, ta, tb, mr, nr, kr, malpha, A21, lda, B12, ldb, F.one, C22, ldc, H2);
        // WriteMatrix(std::cerr<<"P2 = -S2*T2 + C22"<<std::endl, F, mr, nr, C22, ldc);

        // C11 = C22 - beta .  C11  in C11
        fadd (DF, mr, nr, (DFCEptr) C22, ldc, mbetadf, (DFCEptr) C11, ldc, (DFEptr)C11,ldc);
        // WriteMatrix(std::cerr<<"C11 = C22- beta C11"<<std::endl, F, mr, nr, C11, ldc);
      // P3 = alpha . A11 * B11 +  C22  in C22
        MMH_t H3(F, WH.recLevel-1,
                 WH.Amin, WH.Amax,
                 WH.Bmin, WH.Bmax,
                 H2.Outmin, H2.Outmax);
        fgemm (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.one, C22, ldc, H3);
        // WriteMatrix(std::cerr<<"P3 = A11 *B11 + C22"<<std::endl, F, mr, nr, C22, ldc);
        // C11 =  C11 - C22 in C11
        fsubin (DF, mr, nr, (DFCEptr) C22, ldc, (DFEptr) C11, ldc);
        // WriteMatrix(std::cerr<<"C11 = C11 - C22"<<std::endl, F, mr, nr, C11, ldc);
       // C21 = C22 + beta .  C21  in C11
        fadd (DF, mr, nr, (DFCEptr) C22, ldc, betadf, (DFCEptr) C21, ldc, (DFEptr)C21,ldc);
        // WriteMatrix(std::cerr<<"C21 = C22 + beta C21"<<std::endl, F, mr, nr, C21, ldc);
        // T3 = T2 + B21 in B12
        faddin (DF, lb, cb, (DFCEptr) B21, ldb, (DFEptr) B12, ldb);
        // WriteMatrix(std::cerr<<"T3 = T2 + B21"<<std::endl, F, kr, nr, B12, ldb);
    // P4 = alpha . A22 * T3 +  C21  in C21
        MMH_t H4(F, WH.recLevel-1,
                 WH.Amin, WH.Amax,
                 2*WH.Bmin-2*2*WH.Bmax, 2*WH.Bmax-2*WH.Bmin,
                 H3.Outmin, H3.Outmax); // TODO
        fgemm (F, ta, tb, mr, nr, kr, alpha, A22, lda, B12, ldb, F.one, C21, ldc, H4);
        // WriteMatrix(std::cerr<<"P4 = A22 * T3 + C21"<<std::endl, F, mr, nr, C21, ldc);
        // T4 = T3 + B22 - B21 in B12
        faddin (DF, lb, cb, (DFCEptr) B22, ldb, (DFEptr) B12, ldb);
        fsubin (DF, lb, cb, (DFCEptr) B21, ldb, (DFEptr) B12, ldb);
        // S3 =  S2 - A12 in A21
        fsubin (DF, la, ca, (DFCEptr) A12, lda, (DFEptr) A21, lda);
        // P5 = malpha . S3 * B22 -  C12  in C12
        MMH_t H5(F, WH.recLevel-1,
                 WH.Amin, WH.Amax,
                 2*WH.Bmin-2*2*WH.Bmax, 2*WH.Bmax-2*WH.Bmin,
                 WH.Cmin, WH.Cmax); // TODO
        fgemm (F, ta, tb, mr, nr, kr, malpha, A21, lda, B22, ldb, F.mOne, C12, ldc, H5);
        // S4 =  S3 + A12 + A11 in A21
        faddin (DF, la, ca, (DFCEptr) A12, lda, (DFEptr) A21, lda);
        faddin (DF, la, ca, (DFCEptr) A11, lda, (DFEptr) A21, lda);
        // P6 = malpha . S4 * T4 +  C22  in C12
        MMH_t H6(F, WH.recLevel-1,
                 2*WH.Amin, 2*WH.Amax,
                 WH.Bmin-WH.Bmax, WH.Bmax-WH.Bmin,
                 H3.Outmin, H3.Outmax); // TODO
        fgemm (F, ta, tb, mr, nr, kr, alpha, A21, lda, B12, ldb, F.one, C22, ldc, H6);
        // C12 =  C12 + C22 in C12
        faddin (DF, mr, nr, (DFCEptr) C22, ldc, (DFEptr) C12, ldc);
        // B12 = T4 + B11 
        faddin (DF, lb, cb, (DFCEptr) B11, ldb, (DFEptr) B12, ldb);
        // A21 =  S4 - A22 in A21
        fsubin (DF, la, ca, (DFCEptr) A22, lda, (DFEptr) A21, lda);
        // P6 = alpha . A12 * B21 - C11  in C12
        MMH_t H7(F, WH.recLevel-1,
                 WH.Amin, WH.Amax,
                 WH.Bmin, WH.Bmax,
                 H3.Outmin, H3.Outmax); 
        fgemm (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.mOne, C11, ldc, H7);
        
        // // U7 =  U3 + C22 in C22
        // DFElt U7Min, U7Max;
        // if (Protected::NeedPreAddReduction (U7Min, U7Max, H7.Outmin, H7.Outmax, C22Min, C22Max, WH)){
        //     freduce(F,mr,nr,X1,nr);
        //     freduce(F,mr,nr,C22,ldc);
        // }
        // faddin(DF,mr,nr,(DFCEptr)X1,nr,(DFEptr)C22,ldc);

        // // U6 = U3 - P4 in C21
        // DFElt U6Min, U6Max;
        // if (Protected::NeedPreSubReduction(U6Min, U6Max, H7.Outmin, H7.Outmax, H4.Outmin, H4.Outmax, WH)){
        //     freduce(F,mr,nr,X1,nr);
        //     freduce(F,mr,nr,C21,ldc);
        // }
        // fsub(DF,mr,nr,(DFCEptr)X1,nr,(DFCEptr)C21,ldc,(DFEptr)C21,ldc);


        // Updating WH with Outmin, Outmax of the result
        //TODO
        // WH.Outmin = min4 (U1Min, H3.Outmin, U6Min, U7Min);
        // WH.Outmax = max4 (U1Max, H3.Outmax, U6Max, U7Max);
    } // WinogradAcc

    // 2 temps and 24 ops
    // TODO: Add check for modular reductions before final additions
    template < class Field, class FieldTrait >
    inline void WinogradAcc_2_24 (const Field& F,
                                  const FFLAS_TRANSPOSE ta,
                                  const FFLAS_TRANSPOSE tb,
                                  const size_t mr, const size_t nr, const size_t kr,
                                  const typename Field::Element alpha,
                                  const typename Field::Element_ptr A,const size_t lda,
                                  const typename Field::Element_ptr B,const size_t ldb,
                                  const typename Field::Element  beta,
                                  typename Field::Element_ptr C, const size_t ldc,
                                  MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > & WH
                                 )
    {
        MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > H = WH ;
        H.recLevel = H.recLevel - 1 ;

        FFLASFFPACK_check(!F.isZero(beta));

        typename Field::Element malpha ;
        F.neg(malpha,alpha);

        // A, B and c submatrices
        const typename Field::Element_ptr A11=A, A12, A21, A22;
        const typename Field::Element_ptr B11=B, B12, B21, B22;
        typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;



        size_t la, ca, lb, cb; // lines and columns in A,B sub matrices

        // Three temporary submatrices are required

        if (ta == FflasTrans) {
            A21 = A + mr;
            A12 = A + kr*lda;
            A22 = A12 + mr;
            la = kr ;
            ca = mr ;
        }
        else { // ta == FflasNoTrans
            A12 = A + kr;
            A21 = A + mr*lda;
            A22 = A21 + kr;
            la = mr ;
            ca = kr ;

        }
        if (tb == FflasTrans) {
            B21 = B + kr;
            B12 = B + nr*ldb;
            B22 = B12 + kr;
            lb = nr ;
            cb = kr ;

        }
        else { // ta == FflasNoTrans
            B12 = B + nr;
            B21 = B + kr*ldb;
            B22 = B21 + nr;
            lb = kr ;
            cb = nr ;
        }

        // Z1 = C22 - C12         in C22
        fsubin(F,mr,nr,C12,ldc,C22,ldc);
        // Z3 = C12-C21           in C12
        fsubin(F,mr,nr,C21,ldc,C12,ldc);
        // S1 = A21 + A22         in X
        typename Field::Element_ptr X = fflas_new(F,mr,std::max(nr,kr));
        fadd(F,la,ca,A21,lda,A22,lda,X,ca);
        // T1 = B12 - B11         in Y
        typename Field::Element_ptr Y = fflas_new(F,nr,kr);
        fsub(F,lb,cb,B12,ldb,B11,ldb,Y,cb);
        // P5 = a S1 T1 + b Z3    in C12
        fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, Y, cb, beta, C12, ldc, H);
        // S2 = S1 - A11          in X
        fsubin(F,la,ca,A11,lda,X,ca);
        // T2 = B22 - T1          in Y
        fsub(F,lb,cb,B22,ldb,Y,cb,Y,cb);
        // P6 = a S2 T2 + b C21   in C21
        fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, Y, cb, beta, C21, ldc, H);
        // S4 = A12 - S2          in X
        fsub(F,la,ca,A12,lda,X,ca,X,ca);
        // W1 = P5 + beta Z1      in C22
        fadd(F,mr,nr,C12,ldc,beta,C22,ldc,C22,ldc);
        // P3 = a S4 B22 + P5     in C12
        fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, B22, ldb, F.one, C12, ldc, H);
        // P1 = a A11 B11         in X
        fgemm (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X, nr, H);
        // U2 = P6 + P1           in C21
        faddin(F,mr,nr,X,nr,C21,ldc);
        // P2 = a A12 B21 + b C11 in C11
        fgemm (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, beta, C11, ldc, H);
        // U1 = P1 + P2           in C11
        faddin(F,mr,nr,X,nr,C11,ldc);
        // U5 = U2 + P3           in C12
        faddin(F,mr,nr,C21,ldc,C12,ldc);
        // S3 =  A11 - A21        in X ;
        fsub(F,la,ca,A11,lda,A21,lda,X,ca);
        // T3 = B22 - B12         in Y
        fsub(F,lb,cb,B22,ldb,B12,ldb,Y,cb);
        // U3 = a S3 T3 + U2      in C21
        fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, Y, cb, F.one, C21, ldc, H);
        fflas_delete (X);
        // U7 = U3 + W1           in C22
        faddin(F,mr,nr,C21,ldc,C22,ldc);
        // T1_ = B12 - B11        in Y
        fsub(F,lb,cb,B12,ldb,B11,ldb,Y,cb);
        // T2_ = B22 - T1_        in Y
        fsub(F,lb,cb,B22,ldb,Y,cb,Y,cb);
        // T4 = T2_ - B21         in Y
        fsub(F,lb,cb,Y,cb,B21,ldb,Y,cb);
        // U6 = -a A22 T4 + U3    in C21;
        fgemm (F, ta, tb, mr, nr, kr, malpha, A22, lda, Y, cb, F.one, C21, ldc, H);
        fflas_delete (Y);


    } // WinogradAccOld

    // 2 temps and 27 ops
    // TODO: Add check for modular reductions before final additions
    template < class Field, class FieldTrait >
    inline void WinogradAcc_2_27 (const Field& F,
                                  const FFLAS_TRANSPOSE ta,
                                  const FFLAS_TRANSPOSE tb,
                                  const size_t mr, const size_t nr, const size_t kr,
                                  const typename Field::Element alpha,
                                  const typename Field::Element_ptr A,const size_t lda,
                                  const typename Field::Element_ptr B,const size_t ldb,
                                  const typename Field::Element  beta,
                                  typename Field::Element_ptr C, const size_t ldc,
                                  MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > & WH)
    {
        MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > H = WH ;
        H.recLevel = H.recLevel - 1 ;

        FFLASFFPACK_check(!F.isZero(beta));

        typename Field::Element malpha ;
        F.neg(malpha,alpha);

        // A, B and c submatrices
        const typename Field::Element_ptr A11=A, A12, A21, A22;
        const typename Field::Element_ptr B11=B, B12, B21, B22;
        typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;



        size_t la, ca, lb, cb; // lines and columns in A,B sub matrices

        // Three temporary submatrices are required

        if (ta == FflasTrans) {
            A21 = A + mr;
            A12 = A + kr*lda;
            A22 = A12 + mr;
            la = kr ;
            ca = mr ;
        }
        else { // ta == FflasNoTrans
            A12 = A + kr;
            A21 = A + mr*lda;
            A22 = A21 + kr;
            la = mr ;
            ca = kr ;

        }
        if (tb == FflasTrans) {
            B21 = B + kr;
            B12 = B + nr*ldb;
            B22 = B12 + kr;
            lb = nr ;
            cb = kr ;

        }
        else { // ta == FflasNoTrans
            B12 = B + nr;
            B21 = B + kr*ldb;
            B22 = B21 + nr;
            lb = kr ;
            cb = nr ;
        }

        // Z1 = C22 - C12         in C22
        fsubin(F,mr,nr,C12,ldc,C22,ldc);
        // Z3 = C12-C21           in C12
        fsubin(F,mr,nr,C21,ldc,C12,ldc);
        // S1 = A21 + A22         in X
        typename Field::Element_ptr X = fflas_new(F,mr,std::max(nr,kr));
        fadd(F,la,ca,A21,lda,A22,lda,X,ca);
        // T1 = B12 - B11         in Y
        typename Field::Element_ptr Y = fflas_new(F,nr,std::max(kr,mr));
        fsub(F,lb,cb,B12,ldb,B11,ldb,Y,cb);
        // P5 = a S1 T1 + b Z3    in C12
        fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, Y, cb, beta, C12, ldc, H);
        // S2 = S1 - A11          in X
        fsubin(F,la,ca,A11,lda,X,ca);
        // T2 = B22 - T1          in Y
        fsub(F,lb,cb,B22,ldb,Y,cb,Y,cb);
        // P6 = a S2 T2 + b C21   in C21
        fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, Y, cb, beta, C21, ldc, H);
        // S4 = A12 - S2          in X
        fsub(F,la,ca,A12,lda,X,ca,X,ca);
        // W1 = P5 + beta Z1      in C22
        fadd(F,mr,nr,C12,ldc,beta,C22,ldc,C22,ldc);
        // P3 = a S4 B22 + P5     in C12
        fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, B22, ldb, F.zero, Y, nr, H);
        fadd(F,mr,nr,Y,nr,C12,ldc,C12,ldc);
        // P1 = a A11 B11         in X
        fgemm (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X, nr, H);
        // U2 = P6 + P1           in C21
        faddin(F,mr,nr,X,nr,C21,ldc);
        // P2 = a A12 B21 + b C11 in C11
        fgemm (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.zero, Y, nr, H);
        fadd(F,mr,nr,Y,nr,beta,C11,ldc,C11,ldc);
        // U1 = P1 + P2           in C11
        faddin(F,mr,nr,X,nr,C11,ldc);
        // U5 = U2 + P3           in C12
        faddin(F,mr,nr,C21,ldc,C12,ldc);
        // S3 =  A11 - A21        in X ;
        fsub(F,la,ca,A11,lda,A21,lda,X,ca);
        // T3 = B22 - B12         in Y
        fsub(F,lb,cb,B22,ldb,B12,ldb,Y,cb);
        // U3 = a S3 T3 + U2      in C21
        fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, Y, cb, F.one, C21, ldc, H);
        // U7 = U3 + W1           in C22
        faddin(F,mr,nr,C21,ldc,C22,ldc);
        // T1_ = B12 - B11        in Y
        fsub(F,lb,cb,B12,ldb,B11,ldb,Y,cb);
        // T2_ = B22 - T1_        in Y
        fsub(F,lb,cb,B22,ldb,Y,cb,Y,cb);
        // T4 = T2_ - B21         in Y
        fsub(F,lb,cb,Y,cb,B21,ldb,Y,cb);
        // U6 = -a A22 T4 + U3    in C21;
        fgemm (F, ta, tb, mr, nr, kr, alpha, A22, lda, Y, cb, F.zero, X, nr, H);
        fflas_delete (Y);
        fsub(F,mr,nr,C21,ldc,X,nr,C21,ldc);
        fflas_delete (X);


    } // WinogradAcc3


} // BLAS3

} // FFLAS

#endif // __FFLASFFPACK_fgemm_winograd_acc_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
