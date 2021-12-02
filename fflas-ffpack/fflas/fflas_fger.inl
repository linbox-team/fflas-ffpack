/* fflas/fflas_fger.inl
 * Copyright (C) 2005 Clement Pernet
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

#ifndef __FFLASFFPACK_fger_INL
#define __FFLASFFPACK_fger_INL

namespace FFLAS {

    template<class Field>
    inline void
    fger (const Field& F, const size_t M, const size_t N,
          const typename Field::Element alpha,
          typename Field::ConstElement_ptr x, const size_t incx,
          typename Field::ConstElement_ptr y, const size_t incy,
          typename Field::Element_ptr A, const size_t lda)
    {
        MMHelper<Field, MMHelperAlgo::Classic> H(F,0);
        fger (F, M, N, alpha, x, incx, y, incy, A, lda, H);
        freduce (F, M, N, A, lda);
    }
} //FFLAS

namespace FFLAS { namespace Protected {

    template<class FloatElement, class Field>
    inline void
    fger_convert (const Field& F, const size_t M, const size_t N,
                  const typename Field::Element alpha,
                  typename Field::ConstElement_ptr x, const size_t incx,
                  typename Field::ConstElement_ptr y, const size_t incy,
                  typename Field::Element_ptr A, const size_t lda)
    {
        Givaro::ModularBalanced<FloatElement> G((FloatElement) F.characteristic());
        FloatElement alphaf;
        F.convert (alphaf, alpha);
        G.reduce(alphaf);
        FloatElement* Af = fflas_new (G,M,N);
        FloatElement* Xf = fflas_new (G,M);
        FloatElement* Yf = fflas_new (G,N);

        fconvert(F, M, N, Af, N, A, lda);
        freduce(G, M, N, Af, N);
        fconvert(F, M, Xf, 1, x, incx);
        freduce(G, M, Xf, 1);
        fconvert(F, N, Yf, 1, y, incy);
        freduce(G, N, Yf, 1);
        fger (G, M, N, alphaf, Xf, 1, Yf, 1, Af, N);

        finit (F, M, N, Af, N, A, lda);

        fflas_delete (Af);
        fflas_delete (Xf);
        fflas_delete (Yf);
    }
}// Protected
}// FFLAS

namespace FFLAS{

    template<class Field>
    inline void
    fger (const Field& F, const size_t M, const size_t N,
          const typename Field::Element alpha,
          typename Field::ConstElement_ptr x, const size_t incx,
          typename Field::ConstElement_ptr y, const size_t incy,
          typename Field::Element_ptr A, const size_t lda,
          MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> > & H)
    {
        if (F.isZero(alpha)) { return ; }

        if (F.cardinality() < DOUBLE_TO_FLOAT_CROSSOVER && F.cardinality() > 2){
            return Protected::fger_convert<float,Field>(F,M,N,alpha,x, incx, y,incy, A, lda);
        } else if  (16*F.cardinality() < Givaro::ModularBalanced<double>::maxCardinality()){
            return Protected::fger_convert<double,Field>(F,M,N,alpha,x, incx, y,incy, A, lda);
        } else {
            FFPACK::failure()(__func__,__LINE__,"Invalid ConvertTo Mode for this field");
        }
    }


    template<class Field,class AnyTag>
    inline void
    fger (const Field& F, const size_t M, const size_t N,
          const typename Field::Element alpha,
          typename Field::ConstElement_ptr x, const size_t incx,
          typename Field::ConstElement_ptr y, const size_t incy,
          typename Field::Element_ptr A, const size_t lda,
          //	      MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DefaultTag> & H)
         MMHelper<Field, MMHelperAlgo::Classic, AnyTag> & H)
         {
             if (F.isZero(alpha)) { return ; }

             typename Field::Element tmp;
             typename Field::ConstElement_ptr xi=x, yj=y;
             typename Field::Element_ptr Ai=A;


             if ( M < N ){
                 if ( F.isOne( alpha ) )
                     for ( ; Ai < A+M*lda; Ai+=lda, xi+=incx ){
                         yj = y;
                         for (size_t j = 0; j < N; ++j, yj+=incy )
                             F.axpyin( *(Ai+j), *xi, *yj );
                     }
                 else if ( F.isMOne( alpha ) )
                     for ( ; Ai < A+M*lda; Ai+=lda, xi+=incx ){
                         F.neg( tmp, *xi );
                         yj = y;
                         for (size_t j = 0; j < N; ++j, yj+=incy )
                             F.axpyin( *(Ai+j), tmp, *yj );
                     }
                 else
                     for ( ; Ai < A+M*lda; Ai+=lda, xi+=incx ){
                         F.mul( tmp, alpha, *xi );
                         yj = y;
                         for (size_t j = 0; j < N; ++j, yj+=incy )
                             F.axpyin( *(Ai+j), tmp, *yj );
                     }
             } else {
                 if ( F.isOne( alpha ) ){
                     for ( ; Ai < A+N; ++Ai, yj+=incy ){
                         xi = x;
                         for (size_t i = 0; i < M; ++i, xi+=incx )
                             F.axpyin( *(Ai+i*lda), *xi, *yj );
                     }
                 }
                 else if ( F.isMOne( alpha ) )
                     for ( ; Ai < A+N; ++Ai, yj+=incy ){
                         F.neg( tmp, *yj );
                         xi = x;
                         for (size_t i = 0; i < M; ++i, xi+=incx )
                             F.axpyin( *(Ai+i*lda), *xi, tmp );
                     }
                 else
                     for ( ; Ai < A+N; ++Ai, yj+=incy ){
                         F.mul( tmp, alpha, *yj );
                         xi = x;
                         for (size_t i = 0; i < M; ++i, xi+=incx )
                             F.axpyin( *(Ai+i*lda), *xi, tmp );
                     }
             }

         }

    inline void
    fger( const Givaro::DoubleDomain& F, const size_t M, const size_t N,
          const Givaro::DoubleDomain::Element alpha,
          const Givaro::DoubleDomain::ConstElement_ptr x, const size_t incx,
          const Givaro::DoubleDomain::ConstElement_ptr y, const size_t incy,
          Givaro::DoubleDomain::Element_ptr A, const size_t lda,
          MMHelper<Givaro::DoubleDomain, MMHelperAlgo::Classic, ModeCategories::DefaultTag> & H)
    {
        if (F.isZero(alpha)) { return ; }
        FFLASFFPACK_check(lda);
#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_dger( CblasRowMajor, (int)M, (int)N, alpha, x, (int)incx, y, (int)incy, A, (int)lda );
    }

    template<class Field>
    inline void
    fger(const Field& F, const size_t M, const size_t N,
         const typename Field::Element alpha,
         const typename Field::ConstElement_ptr x, const size_t incx,
         const typename Field::ConstElement_ptr y, const size_t incy,
         typename Field::Element_ptr A, const size_t lda,
         MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DefaultBoundedTag> & H)
    {
        H.setOutBounds (1, alpha, 1.0);
        MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DefaultTag> Hd(F,0);
        fger (F, M, N, alpha, x, incx, y, incy, A, lda, Hd);
    }

    inline void
    fger( const Givaro::FloatDomain& F, const size_t M, const size_t N,
          const Givaro::FloatDomain::Element alpha,
          const Givaro::FloatDomain::ConstElement_ptr x, const size_t incx,
          const Givaro::FloatDomain::ConstElement_ptr y, const size_t incy,
          Givaro::FloatDomain::Element_ptr A, const size_t lda,
          MMHelper<Givaro::FloatDomain, MMHelperAlgo::Classic, ModeCategories::DefaultTag> & H)
    {
        if (F.isZero(alpha)) { return ; }

        FFLASFFPACK_check(lda);
#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        cblas_sger( CblasRowMajor, (int)M, (int)N, alpha, x, (int)incx, y, (int)incy, A, (int)lda );
    }




    template<class Field>
    inline void
    fger (const Field& F, const size_t M, const size_t N,
          const typename Field::Element alpha,
          typename Field::ConstElement_ptr x, const size_t incx,
          typename Field::ConstElement_ptr y, const size_t incy,
          typename Field::Element_ptr A, const size_t lda,
          MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::LazyTag> & H)
    {
        if (F.isZero(alpha)) { return ; }

        typedef MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::LazyTag> HelperType;
        typedef typename HelperType::DelayedField delayedField;
        typedef typename HelperType::DelayedField::Element DFElt;
        typedef typename HelperType::DelayedField::ConstElement_ptr DFCElt_ptr;
        typedef typename HelperType::DelayedField::Element_ptr DFElt_ptr;
        typedef typename Field::Element_ptr					Element_ptr;
        typedef MMHelper<delayedField, MMHelperAlgo::Classic, ModeCategories::DefaultBoundedTag> DelayedHelperType;

        DelayedHelperType Hfp(H);

        if (Hfp.MaxDelayedDim(1.0) < 1){
            if (Hfp.Amin < H.FieldMin || Hfp.Amax>H.FieldMax){
                Hfp.initA();
                freduce_constoverride (F, M, x, incx);
            }
            if (Hfp.Bmin < H.FieldMin || Hfp.Bmax>H.FieldMax){
                Hfp.initB();
                freduce_constoverride (F, N, y, incy);
            }
            if (Hfp.Cmin < H.FieldMin || Hfp.Cmax>H.FieldMax){
                Hfp.initC();
                freduce (F, M, N, A, lda);
            }
        }
        Hfp.Outmin = Hfp.FieldMin;
        Hfp.Outmax = Hfp.FieldMax;
        if (F.isOne(alpha) || F.isMOne(alpha)){
            DFElt alphadf;
            if (F.isMOne( alpha)) alphadf = -F.one;
            else alphadf = F.one;

            fger (H.delayedField, M, N, alphadf, (DFCElt_ptr)x, incx, (DFCElt_ptr)y, incy, (DFElt_ptr)A, lda, Hfp);

            H.Outmin = Hfp.Outmin;
            H.Outmax = Hfp.Outmax;
        } else {
            Element_ptr sY  = FFLAS::fflas_new (F, N);
            fscal(F, N, alpha, y, incy, sY, 1);

            fger (H.delayedField, M, N, 1.0,  (DFCElt_ptr)x, incx, (DFCElt_ptr) sY, 1,  (DFElt_ptr)A, lda, Hfp);

            FFLAS::fflas_delete(sY);

            H.setOutBounds (1, alpha, 1.0);

        }

    }

    template<class Field>
    inline void
    fger (const Field& F, const size_t M, const size_t N,
          const typename Field::Element alpha,
          typename Field::ConstElement_ptr x, const size_t incx,
          typename Field::ConstElement_ptr y, const size_t incy,
          typename Field::Element_ptr A, const size_t lda,
          MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DelayedTag> & H)
    {
        if (F.isZero(alpha)) { return ; }

        if (Protected::AreEqual<Field, Givaro::Modular<int64_t> >::value ||
            Protected::AreEqual<Field, Givaro::ModularBalanced<int64_t> >::value){
            if (F.cardinality() < Givaro::ModularBalanced<double>::maxCardinality()&& F.cardinality() > 2)
                return Protected::fger_convert<double,Field>(F,M,N,alpha,x,incx,y,incy, A,lda);
            else{
                // Stay over int64_t
                MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::LazyTag, ParSeqHelper::Sequential> HG(H);
                HG.recLevel = 0;
                typename Field::Element_ptr sY  = FFLAS::fflas_new (F, N);
                fscal(F, N, alpha, y, incy, sY, 1);
                fgemm(F,FflasNoTrans,FflasNoTrans,M,N,1,F.one,x,incx,sY,1,F.one,A,lda,HG);
                FFLAS::fflas_delete(sY);
                freduce(F,M,N,A,lda);
                H.initOut();
                return;
            }
        }
        typedef MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DelayedTag> ModularHelperType;
        typedef typename ModularHelperType::DelayedField	delayedField;
        typedef typename delayedField::Element	DFElt;
        typedef typename delayedField::ConstElement_ptr	DFCElt_ptr;
        typedef typename delayedField::Element_ptr	DFElt_ptr;
        typedef typename Field::Element_ptr					Element_ptr;
        typedef MMHelper<delayedField, MMHelperAlgo::Classic, ModeCategories::DefaultBoundedTag> DelayedHelperType;

        DelayedHelperType Hfp(H);

        if (F.isOne(alpha) || F.isMOne(alpha)){
            DFElt alphadf;
            if (F.isMOne( alpha)) alphadf = -F.one;
            else alphadf = F.one;

            fger (H.delayedField, M, N, alphadf, (DFCElt_ptr)x, incx, (DFCElt_ptr)y, incy, (DFElt_ptr)A, lda, Hfp);

        } else {
            Element_ptr sY  = FFLAS::fflas_new (F, N);
            fscal(F, N, alpha, y, incy, sY, 1);

            fger (H.delayedField, M, N, H.delayedField.one, (DFCElt_ptr)x, incx, (DFCElt_ptr)sY, (size_t)1, (DFElt_ptr)A, lda, Hfp);

            FFLAS::fflas_delete(sY);

        }

        H.initOut();
    }

} // FFLAS
//#include "fflas-ffpack/fflas/fflas_fger_mp.inl" moved to fflas.h
#endif // __FFLASFFPACK_fger_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
