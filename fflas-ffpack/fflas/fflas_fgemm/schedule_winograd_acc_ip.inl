/* Copyright (C) 2014 the LinBox group
 *
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

/** @file fflas/fflas_fgemm/winograd_acc2.inl
 * @ingroup MMalgos
 * @brief Winograd implementation
 * @bib ISSAC09 Scheduling
 */

#ifndef __FFLASFFPACK_fgemm_winograd_acc_ip_INL
#define __FFLASFFPACK_fgemm_winograd_acc_ip_INL

namespace FFLAS { namespace BLAS3 {

    template < class Field, class FieldTrait >
    inline void WinogradAcc_LR (const Field& F,
                                const FFLAS_TRANSPOSE ta,
                                const FFLAS_TRANSPOSE tb,
                                const size_t mr, const size_t nr, const size_t kr,
                                const typename Field::Element alpha,
                                typename Field::Element_ptr A,const size_t lda,
                                typename Field::Element_ptr B,const size_t ldb,
                                const typename Field::Element  beta,
                                typename Field::Element_ptr C, const size_t ldc,
                                const MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > & WH
                               )
    {
        MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > H = WH ;
        H.recLevel = H.recLevel - 1 ;

        FFLASFFPACK_check(!F.isZero(beta));

        typename Field::Element malpha ;
        F.neg(malpha,alpha);

        // A, B and c submatrices
        typename Field::Element_ptr A11=A, A12, A21, A22;
        typename Field::Element_ptr B11=B, B12, B21, B22;
        typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;


        typename Field::Element mbeta ;
        F.neg(mbeta,beta);

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
        // S1 = A21 + A22         in X
        typename Field::Element_ptr X = fflas_new (F, std::max(std::max(mr*nr,kr*nr),mr*kr), 1);
        fadd(F,la,ca,A21,lda,A22,lda,X,ca);
        // T1 = B12 - B11         in Y
        typename Field::Element_ptr Y = fflas_new (F, std::max(mr,kr), nr);
        fsub(F,lb,cb,B12,ldb,B11,ldb,Y,cb);
        // Z2 = C21 - Z1          in C21
        fsubin(F,mr,nr,C22,ldc,C21,ldc);
        // T3 = B22 - B12         in B12 ;
        fsub(F,lb,cb,B22,ldb,B12,ldb,B12,ldb);
        // S3 =  A11 - A21        in A21
        fsub(F,la,ca,A11,lda,A21,lda,A21,lda);
        // P7 = a S3 T3 + b Z1    in C22
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A21, lda, B12, ldb, beta, C22, ldc, H);
        // S2 = S1 - A11          in A21
        fsub(F,la,ca,X,ca,A11,lda,A21,lda);
        // T2 = B22 - T1          in B12
        fsub(F,lb,cb,B22,ldb,Y,cb,B12,ldb);
        // P5 = a S1 T1 + b C12   in C12
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, X, ca, Y, cb, beta, C12, ldc, H);
        // T4 = T2 - B21          in X
        fsub(F,lb,cb,B12,ldb,B21,ldb,X,cb);
        // W1 = a A22 T4          in Y;
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A22, lda, X, cb, F.zero, Y, nr, H);
        // P4 = W1 - b Z2         in C21
        fadd(F,mr,nr,Y,nr,mbeta,C21,ldc,C21,ldc);
        // S4 = A12 - S2          in A22
        fsub(F,la,ca,A12,lda,A21,lda,A22,lda);
        // P6 = a S2 T2           in X
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A21, lda, B12, ldb, F.zero, X, nr, H);
        // W2 = a A12 B21         in Y
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.zero, Y, nr, H);
        // P2 = W2 + beta C11     in C11
        fadd(F,mr,nr,Y,nr,beta,C11,ldc,C11,ldc);
        // P1 = a A11 B11         in Y
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, Y, nr, H);
        // U1 = P1 + P2           in C11
        faddin(F,mr,nr,Y,nr,C11,ldc);
        // U2 = P6 + P1           in X
        faddin(F,mr,nr,Y,nr,X,nr);
        fflas_delete (Y);
        // U3 = U2 + P7           in C22
        faddin(F,mr,nr,X,nr,C22,ldc);
        // U4 = U2 + P5           in X
        faddin(F,mr,nr,C12,ldc,X,nr);
        // U6 = U3 - P4           in C21
        fsub(F,mr,nr,C22,ldc,C21,ldc,C21,ldc);
        // U7 = U3 + P5           in C22
        faddin(F,mr,nr,C12,ldc,C22,ldc);
        // P3 = a S4 B22          in C12
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A22, lda, B22, ldb, F.zero, C12, ldc, H);
        // U5 = U4 + P3           in C12
        faddin(F,mr,nr,X,nr,C12,ldc);

        fflas_delete (X);


    } // WinogradAccOld

    template < class Field, class FieldTrait >
    inline void WinogradAcc_R_S (const Field& F,
                                 const FFLAS_TRANSPOSE ta,
                                 const FFLAS_TRANSPOSE tb,
                                 const size_t mr, const size_t nr, const size_t kr,
                                 const typename Field::Element alpha,
                                 const typename Field::Element_ptr A,const size_t lda,
                                 typename Field::Element_ptr B,const size_t ldb,
                                 const typename Field::Element  beta,
                                 typename Field::Element_ptr C, const size_t ldc,
                                 const MMHelper<Field, MMHelperAlgo::Winograd,FieldTrait > & WH
                                )
    {
        MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > H = WH ;
        H.recLevel = H.recLevel - 1 ;

        FFLASFFPACK_check(!F.isZero(beta));

        typename Field::Element malpha ;
        F.neg(malpha,alpha);

        // A, B and c submatrices
        const typename Field::Element_ptr A11=A, A12, A21, A22;
        typename Field::Element_ptr B11=B, B12, B21, B22;
        typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;


        typename Field::Element mbeta ;
        F.neg(mbeta,beta);

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

        FFLASFFPACK_check(mr == nr && kr == nr);

        // Z1 = C22 - C12         in C22
        fsubin(F,mr,nr,C12,ldc,C22,ldc);
        // T1 = B12 - B11         in X
        // typename Field::Element_ptr X = fflas_new (F, std::max(mr,kr)*nr];
        typename Field::Element_ptr X = fflas_new (F, mr, nr);
        fsub(F,lb,cb,B12,ldb,B11,ldb,X,cb);
        // Z2 = C21 - Z1          in C21
        fsubin(F,mr,nr,C22,ldc,C21,ldc);
        // T3 = B22 - B12         in B12 ;
        fsub(F,lb,cb,B22,ldb,B12,ldb,B12,ldb);
        // S3 =  A11 - A21        in Y
        typename Field::Element_ptr Y = fflas_new (F, mr, kr);
        fsub(F,la,ca,A11,lda,A21,lda,Y,ca);
        // P7 = a S3 T3 + b Z1    in C22
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, Y, ca, B12, ldb, beta, C22, ldc, H);
        // S1 = A21 + A22         in Y
        fadd(F,la,ca,A21,lda,A22,lda,Y,ca);
        // T2 = B22 - T1          in B12
        fsub(F,lb,cb,B22,ldb,X,cb,B12,ldb);
        // P5 = a S1 T1 + b C12   in C12
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, Y, ca, X, cb, beta, C12, ldc, H);
        // T4 = T2 - B21          in X
        fsub(F,lb,cb,B12,ldb,B21,ldb,X,cb);
        // P4 = a A22 T4 - b Z2   in C21
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A22, lda, X, cb, mbeta, C21, ldc, H);
        // W1 = a A12 B21          in X;
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.zero, X, nr, H);
        // P2 = W1 + beta C11     in C11
        fadd(F,mr,nr,X,nr,beta,C11,ldc,C11,ldc);
        // S2 = S1 - A11          in Y
        fsubin(F,la,ca,A11,lda,Y,ca);
        // P6 = a S2 T2           in B21
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, Y, ca, B12, ldb, F.zero, B21, ldb, H);
        // S4 = A12 - S2          in Y
        fsub(F,la,ca,A12,lda,Y,ca,Y,ca);
        // P1 = a A11 B11         in X
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X, nr, H);
        // U2 = P6 + P1           in B21
        faddin(F,mr,nr,X,nr,B21,ldb);
        // U3 = U2 + P7           in C22
        faddin(F,mr,nr,B21,ldb,C22,ldc);
        // U4 = U2 + P5           in B21
        faddin(F,mr,nr,C12,ldc,B21,ldb);
        // U6 = U3 - P4           in C21
        fsub(F,mr,nr,C22,ldc,C21,ldc,C21,ldc);
        // U1 = P1 + P2           in C11
        faddin(F,mr,nr,X,nr,C11,ldc);
        fflas_delete (X);
        // U7 = U3 + P5           in C22
        faddin(F,mr,nr,C12,ldc,C22,ldc);
        // P3 = a S4 B22          in C12
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, Y, ca, B22, ldb, F.zero, C12, ldc, H);
        fflas_delete (Y);
        // U5 = U4 + P3           in C12
        faddin(F,mr,nr,B21,ldb,C12,ldc);




    } // WinogradAccOld


    template < class Field ,class FieldTrait>
    inline void WinogradAcc_L_S (const Field& F,
                                 const FFLAS_TRANSPOSE ta,
                                 const FFLAS_TRANSPOSE tb,
                                 const size_t mr, const size_t nr, const size_t kr,
                                 const typename Field::Element alpha,
                                 typename Field::Element_ptr A,const size_t lda,
                                 const typename Field::Element_ptr B,const size_t ldb,
                                 const typename Field::Element  beta,
                                 typename Field::Element_ptr C, const size_t ldc,
                                 const MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > & WH
                                )
    {
        MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > H = WH ;
        H.recLevel = H.recLevel - 1 ;

        FFLASFFPACK_check(!F.isZero(beta));

        typename Field::Element malpha ;
        F.neg(malpha,alpha);

        // A, B and c submatrices
        typename Field::Element_ptr A11=A, A12, A21, A22;
        const typename Field::Element_ptr B11=B, B12, B21, B22;
        typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;


        typename Field::Element mbeta ;
        F.neg(mbeta,beta);

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

        FFLASFFPACK_check(mr == nr && kr == nr);

        // Z1 = C22 - C12         in C22
        fsubin(F,mr,nr,C12,ldc,C22,ldc);
        // Z2 = C21 - Z1          in C21
        fsubin(F,mr,nr,C22,ldc,C21,ldc);
        // S3 =  A11 - A21        in X
        typename Field::Element_ptr X = fflas_new (F, mr, nr);
        fsub(F,la,ca,A11,lda,A21,lda,X,ca);
        // S1 = A21 + A22         in A21
        faddin(F,la,ca,A22,lda,A21,lda);
        // T3 = B22 - B12         in Y ;
        typename Field::Element_ptr Y = fflas_new (F, mr, kr);
        fsub(F,lb,cb,B22,ldb,B12,ldb,Y,cb);
        // P7 = a S3 T3 + b Z1    in C22
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, X, ca, Y, cb, beta, C22, ldc, H);
        // T1 = B12 - B11         in X
        fsub(F,lb,cb,B12,ldb,B11,ldb,X,cb);
        // T2 = B22 - T1          in Y
        fsub(F,lb,cb,B22,ldb,X,cb,Y,cb);
        // P5 = a S1 T1 + b C12   in C12
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A21, lda, X, cb, beta, C12, ldc, H);
        // S2 = S1 - A11          in A21
        fsubin(F,la,ca,A11,lda,A21,lda);
        // P1 = a A11 B11         in X
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X, nr, H);
        // S4 = A12 - S2          in A11
        fsub(F,la,ca,A12,lda,A21,lda,A11,lda);
        // P2 = a A12 B21 + b C11 in C11;
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, beta, C11, ldc, H);
        // U1 = P1 + P2           in C11
        faddin(F,mr,nr,X,nr,C11,ldc);
        // P6 = a S2 T2           in A12
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A21, lda, Y, cb, F.zero, A12, lda, H);
        // T4 = T2 - B21          in Y
        fsubin(F,lb,cb,B21,ldb,Y,cb);
        // W2 = a A22 T4          in A21
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A22, lda, Y, cb, F.zero, A21, lda, H);
        // P4 = W2 - beta Z2      in C21
        fadd(F,mr,nr,A21,lda,mbeta,C21,ldc,C21,ldc);
        // U2 = P6 + P1           in X
        faddin(F,mr,nr,A12,lda,X,nr);
        // U3 = U2 + P7           in C22
        faddin(F,mr,nr,X,nr,C22,ldc);
        // U6 = U3 - P4           in C21
        fsub(F,mr,nr,C22,ldc,C21,ldc,C21,ldc);
        // U7 = U3 + P5           in C22
        faddin(F,mr,nr,C12,ldc,C22,ldc);
        // U4 = U2 + P5           in C12
        faddin(F,mr,nr,X,nr,C12,ldc);
        fflas_delete (X);
        // W3 = a S4 B22          in Y
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A11, lda, B22, ldb, F.zero, Y, nr, H);
        // U5 = U4 + W3           in C12
        faddin(F,mr,nr,Y,nr,C12,ldc);
        fflas_delete (Y);


    } // WinogradAccOld

} // BLAS3
} // FFLAS

#endif // __FFLASFFPACK_fgemm_winograd_acc_ip_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
