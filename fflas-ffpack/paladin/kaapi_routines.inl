/* fflas/fflas_pftrsm.inl
 * Copyright (C) 2013 Ziad Sultan
 *
 * Written by Ziad Sultan  < Ziad.Sultan@imag.fr >
 * Time-stamp: <17 Jun 14 14:32:29 Jean-Guillaume.Dumas@imag.fr>
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


#ifndef __FFLASFFPACK_KAAPI_ROUTINES_INL
#define __FFLASFFPACK_KAAPI_ROUTINES_INL





#ifdef __FFLASFFPACK_USE_KAAPI
namespace FFLAS {

    template<class Field, class Helper>
    struct Taskfgemm15 : public ka::Task<15>::Signature<
                         Field,
                         FFLAS_TRANSPOSE,
                         FFLAS_TRANSPOSE,
                         size_t ,
                         size_t ,
                         size_t ,
                         typename Field::Element,
                         ka::R<typename Field::Element>,
                         size_t ,
                         ka::R<typename Field::Element>,
                         size_t ,
                         typename Field::Element,
                         ka::RW<typename Field::Element>,
                         size_t ,
                         Helper
                         //size_t // winograd
                         >{};
    /*
       template<class Field>
       struct Taskfgemm14 : public ka::Task<14>::Signature<
       Field,
       FFLAS_TRANSPOSE,
       FFLAS_TRANSPOSE,
       size_t ,
       size_t ,
       size_t ,
       typename Field::Element,
       ka::R<typename Field::Element>,
       size_t ,
       ka::R<typename Field::Element>,
       size_t ,
       typename Field::Element,
       ka::RW<typename Field::Element>,
       size_t
       >{};


*/
    template<class Field>
    struct Taskftrsm12: public ka::Task<12>::Signature<
                        Field , /* Field F */
                        FFLAS::FFLAS_SIDE ,
                        FFLAS::FFLAS_UPLO ,
                        FFLAS::FFLAS_TRANSPOSE ,
                        FFLAS::FFLAS_DIAG ,
                        size_t ,   /* size : M */
                        size_t ,   /* size : N */
                        typename Field::Element ,
                        ka::R<typename Field::Element >, /* Matrix A */
                        size_t , /* lda */
                        ka::RW<typename Field::Element >, /* Matrix B */
                        size_t  /* ldb */
                        >{};



    template<class Field, class Helper>
    void spawnerfgemm(const Field& F,
                      const FFLAS::FFLAS_TRANSPOSE ta,
                      const FFLAS::FFLAS_TRANSPOSE tb,
                      size_t BlockRowDim,
                      size_t BlockColDim,
                      size_t k,
                      const typename Field::Element alpha,
                      ka::pointer_r<typename Field::Element> A,
                      const size_t lda,
                      ka::pointer_r<typename Field::Element> B,
                      const size_t ldb,
                      const typename Field::Element beta,
                      ka::pointer_rw<typename Field::Element> C, const size_t ldc,
                      Helper WH){
        ka::Spawn<Taskfgemm15<Field, Helper> >()( F, ta, tb, BlockRowDim, BlockColDim, k, alpha, A.ptr(), lda, B.ptr() , ldb,
                                                  beta, C.ptr(), ldc, WH);
    }

}

template<class Field, class Helper>
struct TaskBodyCPU<FFLAS::Taskfgemm15<Field, Helper> >{
    void operator()(const Field& F,
                    const FFLAS::FFLAS_TRANSPOSE ta,
                    const FFLAS::FFLAS_TRANSPOSE tb,
                    size_t BlockRowDim,
                    size_t BlockColDim,
                    size_t k,
                    const typename Field::Element alpha,
                    ka::pointer_r<typename Field::Element> A,
                    const size_t lda,
                    ka::pointer_r<typename Field::Element> B,
                    const size_t ldb,
                    const typename Field::Element beta,
                    ka::pointer_rw<typename Field::Element> C, const size_t ldc,
                    Helper WH
                    //	Helper & WH
                    //	size_t w
                   )
    {
        FFLAS::MMHelper<Field, FFLAS::MMHelperAlgo::Winograd, typename FFLAS::FieldTraits<Field>::value> W(WH);
        /*
           FFLAS::MMHelper<Field, FFLAS::MMHelperAlgo::Winograd> WH;
           WH(F,w);*/
        FFLAS::fgemm( F, ta, tb, BlockRowDim, BlockColDim, k, alpha, A.ptr(), lda, B.ptr() , ldb,
                      beta, C.ptr(), ldc, W);
    }
};


/*
   template<class Field>
   struct TaskBodyCPU<FFLAS::Taskfgemm14<Field> >{
   void operator()(const Field& F,
   const FFLAS::FFLAS_TRANSPOSE ta,
   const FFLAS::FFLAS_TRANSPOSE tb,
   size_t BlockRowDim,
   size_t BlockColDim,
   size_t k,
   const typename Field::Element alpha,
   ka::pointer_r<typename Field::Element> A,
   const size_t lda,
   ka::pointer_r<typename Field::Element> B,
   const size_t ldb,
   const typename Field::Element beta,
   ka::pointer_rw<typename Field::Element> C, const size_t ldc)
   {
   FFLAS::fgemm( F, ta, tb, BlockRowDim, BlockColDim, k, alpha, A.ptr(), lda, B.ptr() , ldb,
   beta, C.ptr(), ldc);
   }
   };
   */
template<class Field>
struct TaskBodyCPU<FFLAS::Taskftrsm12<Field> > {
    void operator()(const Field & F, const FFLAS::FFLAS_SIDE Side,
                    const FFLAS::FFLAS_UPLO Uplo,
                    const FFLAS::FFLAS_TRANSPOSE TransA,
                    const FFLAS::FFLAS_DIAG Diag,
                    const size_t M, const size_t N,
                    const typename Field::Element alpha,
                    ka::pointer_r<typename Field::Element > A, const size_t lda,
                    ka::pointer_rw<typename Field::Element > B, const size_t ldb )
    {

        FFLAS::ftrsm(F, Side, Uplo, TransA, Diag, M, N, alpha, A.ptr(), lda, B.ptr(), ldb);
    }
};


#endif


#endif //  __FFLASFFPACK_KAAPI_ROUTINES_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
