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

/** @file fflas/fflas_fgemm/schedule_bini.inl
 * @ingroup MMalgos
 * @brief Bini implementation
 */

#ifndef __FFLASFFPACK_fgemm_bini_INL
#define __FFLASFFPACK_fgemm_bini_INL

namespace FFLAS { namespace BLAS3 {

    template < class Field >
    inline void Bini (const Field& F,
                      const FFLAS_TRANSPOSE ta,
                      const FFLAS_TRANSPOSE tb,
                      const size_t mr, const size_t nr, const size_t kr,
                      const typename Field::Element alpha,
                      const typename Field::Element_ptr A,const size_t lda,
                      const typename Field::Element_ptr B,const size_t ldb,
                      const typename Field::Element  beta,
                      typename Field::Element_ptr C, const size_t ldc,
                      const size_t kmax, const size_t w, const FFLAS_BASE base,
                      const size_t rec_level)
    {

        FFLASFFPACK_check(F.isZero(beta));
        FFLASFFPACK_check(rec_level>0);

        size_t imaxb, jmaxb, imaxa, jmaxa, ldx2;
        // size_t x3rd = std::max(mr,kr);
        const typename Field::Element_ptr d11,d12,d21,d22;
        typename Field::Element_ptr d11c,d12c,d21c,d22c,dx1,dx2;
        const typename Field::Element_ptr A11=A, A12, A21, A22;
        const typename Field::Element_ptr B11=B, B12, B21, B22;
        typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C+nr+mr*ldc;


        size_t x1rd = std::max(nr,kr);
        size_t ldx1;
        if (ta == FflasTrans) {
            A21 = A + mr;
            A12 = A + kr*lda;
            A22 = A12 + mr;
            imaxa = kr;
            jmaxa = mr;
            ldx1 = mr;
        }
        else {
            A12 = A + kr;
            A21 = A + mr*lda;
            A22 = A21 + kr;
            imaxa = mr;
            jmaxa = kr;
            ldx1  = x1rd;
        }
        if (tb == FflasTrans) {
            B21 = B + kr;
            B12 = B + nr*ldb;
            B22 = B12 + kr;
            imaxb = nr;
            jmaxb = kr;
            ldx2 = kr;
        }
        else {
            B12 = B + nr;
            B21 = B + kr*ldb;
            B22 = B21 + nr;
            imaxb = kr;
            ldx2 = jmaxb = nr;
        }



    } // Bini

} // BLAS3


} // FFLAS

#endif // __FFLASFFPACK_fgemm_bini_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
