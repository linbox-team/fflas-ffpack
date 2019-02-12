/*
 * Copyright (C) 2014 the LinBox group
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

/** @file fflas/fflas_fgemm/winograd_ip.inl
 * @ingroup MMalgos
 * @brief Winograd implementation
 * @bib ISSAC09 Scheduling
 */

#ifndef __FFLASFFPACK_fgemm_winograd_ip_INL
#define __FFLASFFPACK_fgemm_winograd_ip_INL

namespace FFLAS { namespace BLAS3 {

    template < class Field, class FieldTrait >
    inline void Winograd_LR_S (const Field& F,
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
        MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait> H = WH ;
        H.recLevel = H.recLevel - 1 ;

        FFLASFFPACK_check(F.isZero(beta));

        // FFLASFFPACK_check(mr == nr && mr == kr);
        FFLASFFPACK_check(kr == nr);

        size_t lb, cb, la, ca;
        typename Field::Element_ptr A11=A, A12, A21, A22;
        typename Field::Element_ptr B11=B, B12, B21, B22;
        typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;


        if (ta == FflasTrans) {
            A21  = A + mr;
            A12  = A + kr*lda;
            A22  = A12 + mr;
            la   = kr;
            ca   = mr;
        }
        else {
            A12  = A + kr;
            A21  = A + mr*lda;
            A22  = A21 + kr;
            la   = mr;
            ca   = kr;
        }
        if (tb == FflasTrans) {
            B21  = B + kr;
            B12  = B + nr*ldb;
            B22  = B12 + kr;
            lb   = nr;
            cb   = kr;
        }
        else {
            B12  = B + nr;
            B21  = B + kr*ldb;
            B22  = B21 + nr;
            lb   = kr;
            cb   = nr;
        }


        // S3 = A11 - A21         in C11
        fsub(F,la,ca,A11,lda,A21,lda,C11,ldc);
        // S1 =  A21 + A22        in A21
        faddin(F,la,ca,A22,lda,A21,lda);
        // T1 = B12 - B11         in C22
        fsub(F,lb,cb,B12,ldb,B11,ldb,C22,ldc);
        // T3 = B22 - B12         in B12
        fsub(F,lb,cb,B22,ldb,B12,ldb,B12,ldb);
        // P7 = S3 T3             in C21
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, C11, ldc, B12, ldb, F.zero, C21, ldc, H);
        // S2 = S1 - A11          in C12
        fsub(F,la,ca,A21,lda,A11,lda,C12,ldc);
        // P1 = A11 B11           in C11
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, C11, ldc, H);
        // T2 = B22 - T1          in B11
        fsub(F,lb,cb,B22,ldb,C22,ldc,B11,ldb);
        // P5 = S1 T1             in A11
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A21, lda, C22, ldc, F.zero, A11, lda, H);
        // T4 = T2 - B21          in C22
        fsub(F,lb,cb,B11,ldb,B21,ldb,C22,ldc);
        // P4 = A22 T4            in A21
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A22, lda, C22, ldc, F.zero, A21, lda, H);
        // S4 = A12 - S2          in A22
        fsub(F,la,ca,A12,lda,C12,ldc,A22,lda);
        // P6 = S2 T2             in C22
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, C12, ldc, B11, ldb, F.zero, C22, ldc, H);
        // U2 = P1 + P6           in C22
        faddin(F,mr,nr,C11,ldc,C22,ldc);
        // P2 = A12 B21           in C12
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.zero, C12, ldc, H);
        // U1 = P1 + P2           in C11
        faddin(F,mr,nr,C12,ldc,C11,ldc);
        // U4 = U2 + P5           in C12
        fadd(F,mr,nr,C22,ldc,A11,lda,C12,ldc);
        // U3 = U2 + P7           in C22
        faddin(F,mr,nr,C21,ldc,C22,ldc);
        // U6 = U3 - P4           in C21
        fsub(F,mr,nr,C22,ldc,A21,lda,C21,ldc);
        // U7 = U3 + P5           in C22
        faddin(F,mr,nr,A11,lda,C22,ldc);
        // P3 = S4 B22            in A12
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A22, lda, B22, ldb, F.zero, A12, lda, H);
        // U5 = U4 + P3           in C12
        faddin(F,mr,nr,A12,lda,C12,ldc);



    } // WinogradIP


    template < class Field, class FieldTrait >
    inline void Winograd_L_S(const Field& F,
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

        FFLASFFPACK_check(F.isZero(beta));

        FFLASFFPACK_check(kr == nr && kr <= mr);

        size_t lb, cb, la, ca;
        typename Field::Element_ptr A11=A, A12, A21, A22;
        const typename Field::Element_ptr B11=B, B12, B21, B22;
        typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;


        if (ta == FflasTrans) {
            A21  = A + mr;
            A12  = A + kr*lda;
            A22  = A12 + mr;
            la   = kr;
            ca   = mr;
        }
        else {
            A12  = A + kr;
            A21  = A + mr*lda;
            A22  = A21 + kr;
            la   = mr;
            ca   = kr;
        }
        if (tb == FflasTrans) {
            B21  = B + kr;
            B12  = B + nr*ldb;
            B22  = B12 + kr;
            lb   = nr;
            cb   = kr;
        }
        else {
            B12  = B + nr;
            B21  = B + kr*ldb;
            B22  = B21 + nr;
            lb   = kr;
            cb   = nr;
        }


        // S3 = A11 - A21         in C22
        fsub(F,la,ca,A11,lda,A21,lda,C22,ldc);
        // S1 =  A21 + A22        in A21
        fadd(F,la,ca,A22,lda,A21,lda,A21,lda);
        // S2 = S1 - A11          in C12
        fsub(F,la,ca,A21,lda,A11,lda,C12,ldc);
        // T1 = B12 - B11         in C21
        fsub(F,lb,cb,B12,ldb,B11,ldb,C21,ldc);
        // P1 = A11 B11           in C11
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, C11, ldc, H);
        // T3 = B22 - B12         in A11
        fsub(F,lb,cb,B22,ldb,B12,ldb,A11,lda);
        // P7 = S3 T3             in X
        typename Field::Element_ptr X = fflas_new (F, mr, nr);
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, C22, ldc, A11, lda, F.zero, X, nr, H);
        // T2 = B22 - T1          in A11
        fsub(F,lb,cb,B22,ldb,C21,ldc,A11,lda);
        // P5 = S1 T1             in C22
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A21, lda, C21, ldc, F.zero, C22, ldc, H);
        // S4 = A12 - S2          in C21
        fsub(F,la,ca,A12,lda,C12,ldc,C21,ldc);
        // P3 = S4 B22            in A21
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, C21, ldc, B22, ldb, F.zero, A21, lda, H);
        // P6 = S2 T2             in C21
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, C12, ldc, A11, lda, F.zero, C21, ldc, H);
        // T4 = T2 - B21          in A11
        fsubin(F,lb,cb,B21,ldb,A11,lda);
        // U2 = P1 + P6           in C21
        faddin(F,mr,nr,C11,ldc,C21,ldc);
        // U4 = U2 + P5           in C12
        fadd(F,mr,nr,C22,ldc,C21,ldc,C12,ldc);
        // U3 = U2 + P7           in C21
        faddin(F,mr,nr,X,nr,C21,ldc);
        // U7 = U3 + P5           in C22
        faddin(F,mr,nr,C21,ldc,C22,ldc);
        // U5 = U4 + P3           in C12
        faddin(F,la,ca,A21,lda,C12,ldc);
        // P2 = A12 B21           in X
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.zero, X, nr, H);
        // U1 = P1 + P2           in C11
        faddin(F,mr,nr,X,nr,C11,ldc);
        fflas_delete (X);
        // P4 = A22 T4            in A21
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A22, lda, A11, lda, F.zero, A21, lda, H);
        // U6 = U3 - P4           in C21
        fsubin(F,mr,nr,A21,lda,C21,ldc);


    } // WinogradIP

    template < class Field, class FieldTrait >
    inline void Winograd_R_S(const Field& F,
                             const FFLAS_TRANSPOSE ta,
                             const FFLAS_TRANSPOSE tb,
                             const size_t mr, const size_t nr, const size_t kr,
                             const typename Field::Element alpha,
                             const typename Field::Element_ptr A,const size_t lda,
                             typename Field::Element_ptr B,const size_t ldb,
                             const typename Field::Element  beta,
                             typename Field::Element_ptr C, const size_t ldc,
                             const MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > & WH
                            )
    {
        MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > H = WH ;
        H.recLevel = H.recLevel - 1 ;

        FFLASFFPACK_check(F.isZero(beta));

        FFLASFFPACK_check(kr == nr && kr <= mr);

        size_t lb, cb, la, ca;
        const typename Field::Element_ptr A11=A, A12, A21, A22;
        typename Field::Element_ptr B11=B, B12, B21, B22;
        typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;


        if (ta == FflasTrans) {
            A21  = A + mr;
            A12  = A + kr*lda;
            A22  = A12 + mr;
            la   = kr;
            ca   = mr;
        }
        else {
            A12  = A + kr;
            A21  = A + mr*lda;
            A22  = A21 + kr;
            la   = mr;
            ca   = kr;
        }
        if (tb == FflasTrans) {
            B21  = B + kr;
            B12  = B + nr*ldb;
            B22  = B12 + kr;
            lb   = nr;
            cb   = kr;
        }
        else {
            B12  = B + nr;
            B21  = B + kr*ldb;
            B22  = B21 + nr;
            lb   = kr;
            cb   = nr;
        }


        // S3 = A11 - A21         in C22
        fsub(F,la,ca,A11,lda,A21,lda,C22,ldc);
        // S1 =  A21 + A22        in C21
        fadd(F,la,ca,A22,lda,A21,lda,C21,ldc);
        // T1 = B12 - B11         in C12
        fsub(F,lb,cb,B12,ldb,B11,ldb,C12,ldc);
        // P1 = A11 B11           in C11
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, C11, ldc, H);
        // S2 = S1 - A11          in B11
        fsub(F,la,ca,C21,ldc,A11,lda,B11,ldb);
        // T3 = B22 - B12         in B12
        fsub(F,lb,cb,B22,ldb,B12,ldb,B12,ldb);
        // P7 = S3 T3             in X
        typename Field::Element_ptr X = fflas_new (F, mr, nr);
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, C22, ldc, B12, ldb, F.zero, X, nr, H);
        // T2 = B22 - T1          in B12
        fsub(F,lb,cb,B22,ldb,C12,ldc,B12,ldb);
        // P5 = S1 T1             in C22
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, C21, ldc, C12, ldc, F.zero, C22, ldc, H);
        // T4 = T2 - B21          in C12
        fsub(F,lb,cb,B12,ldb,B21,ldb,C12,ldc);
        // P6 = S2 T2             in C21
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, B11, ldb, B12, ldb, F.zero, C21, ldc, H);
        // P4 = A22 T4            in B12
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A22, lda, C12, ldc, F.zero, B12, ldb, H);
        // S4 = A12 - S2          in B11
        fsub(F,la,ca,A12,lda,B11,ldb,B11,ldb);
        // U2 = P1 + P6           in C21
        faddin(F,mr,nr,C11,ldc,C21,ldc);
        // U4 = U2 + P5           in C12
        fadd(F,mr,nr,C22,ldc,C21,ldc,C12,ldc);
        // U3 = U2 + P7           in C21
        faddin(F,mr,nr,X,nr,C21,ldc);
        fflas_delete (X);
        // U7 = U3 + P5           in C22
        faddin(F,mr,nr,C21,ldc,C22,ldc);
        // U6 = U3 - P4           in C21
        fsubin(F,mr,nr,B12,ldb,C21,ldc);
        // P3 = S4 B22            in B12
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, B11, ldb, B22, ldb, F.zero, B12, ldb, H);
        // U5 = U4 + P3           in C12
        faddin(F,la,ca,B12,ldb,C12,ldc);
        // P2 = A12 B21           in B12
        fgemm2 (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.zero, B12, ldb, H);
        // U1 = P1 + P2           in C11
        faddin(F,mr,nr,B12,ldb,C11,ldc);


    } // WinogradIP

} // BLAS3


} // FFLAS

#endif // __FFLASFFPACK_fgemm_winograd_ip_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
