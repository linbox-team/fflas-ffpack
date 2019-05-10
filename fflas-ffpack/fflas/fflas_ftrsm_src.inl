/* Copyright (C) 2005 C. Pernet
 * Written by C. Pernet
 *  ========LICENCE========
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
 */

#define Mjoin(pre, nam) my_join(pre, nam)
#define my_join(pre, nam) pre ## nam

#ifdef __FFLAS__TRANSPOSE
#define __FFLAS__Acolinc lda
#define __FFLAS__Arowinc 1
#ifdef __FFLAS__LOW
#define __FFLAS__UPPER
#else
#define __FFLAS__LOWER
#endif
#else
#ifdef __FFLAS__LOW
#define __FFLAS__LOWER
#else
#define __FFLAS__UPPER
#endif
#define __FFLAS__Acolinc 1
#define __FFLAS__Arowinc lda
#endif

#ifdef __FFLAS__LEFT
#define __FFLAS__SIDE Left
#define __FFLAS__Na M
#define __FFLAS__Nb N
#ifdef __FFLAS__TRANSPOSE
#define __FFLAS__Acopcolinc __FFLAS__Na
#define __FFLAS__Acoprowinc 1
#else // __FFLAS__NOTRANSPOSE
#define __FFLAS__Acopcolinc 1
#define __FFLAS__Acoprowinc __FFLAS__Na
#endif
#define __FFLAS__Mb nsplit
#define __FFLAS__Nb2 N
#define __FFLAS__Mb2 M-nsplit
#define __FFLAS__Mbrest nrestsplit
#define __FFLAS__Nbrest N
#define __FFLAS__Mupdate M-(i+1)*nsplit
#define __FFLAS__Nupdate N
#define __FFLAS__Anorminc __FFLAS__Acolinc
#define __FFLAS__Acopnorminc __FFLAS__Acopcolinc
#define __FFLAS__Bnorminc 1
#define __FFLAS__Bnormnext ldb
#define __FFLAS__Bdim N
#ifdef __FFLAS__LOWER
#define __FFLAS__Atriang A + i * nsplit * (lda + 1)
#define __FFLAS__Aupdate A + i * nsplit * (lda + 1) + nsplit*__FFLAS__Arowinc
#define __FFLAS__Arest A + (__FFLAS__Na - nrestsplit) * (lda + 1)
#define __FFLAS__Anormnext __FFLAS__Arowinc
#define __FFLAS__Acopnormnext __FFLAS__Acoprowinc
#define __FFLAS__Bupdate B + (i+1)*nsplit*ldb
#define __FFLAS__Brec B + i * nsplit * ldb
#define __FFLAS__Brest B + (M - nrestsplit) * ldb
#define __FFLAS__A1 A
#define __FFLAS__A2 A + nsplit * __FFLAS__Arowinc
#define __FFLAS__A3 A + nsplit * (lda + 1)
#define __FFLAS__B1 B
#define __FFLAS__B2 B + nsplit * ldb
#define __FFLAS__Normdim i
#else // __FFLAS__UPPER
#define __FFLAS__Atriang A + (__FFLAS__Na - (i + 1) * nsplit) * (lda + 1)
#define __FFLAS__Aupdate A + (__FFLAS__Na  - (i + 1) * nsplit) * __FFLAS__Acolinc
#define __FFLAS__Arest A
#define __FFLAS__Anormnext lda + 1
#define __FFLAS__Acopnormnext __FFLAS__Na + 1
#define __FFLAS__Bupdate B
#define __FFLAS__Brec B + (M - (i + 1) * nsplit) * ldb
#define __FFLAS__Brest B
#define __FFLAS__A1 A + (__FFLAS__Na - nsplit) * (lda + 1)
#define __FFLAS__A2 A + (__FFLAS__Na - nsplit) * __FFLAS__Acolinc
#define __FFLAS__A3 A
#define __FFLAS__B1 B + (M - nsplit)*ldb
#define __FFLAS__B2 B
#define __FFLAS__Normdim __FFLAS__Na-i-1
#endif
#else  // __FFLAS__RIGHT
#define __FFLAS__SIDE Right
#define __FFLAS__Na N
#define __FFLAS__Nb nsplit
#ifdef __FFLAS__TRANSPOSE
#define __FFLAS__Acopcolinc __FFLAS__Na
#define __FFLAS__Acoprowinc 1
#else // __FFLAS__NOTRANSPOSE
#define __FFLAS__Acopcolinc 1
#define __FFLAS__Acoprowinc __FFLAS__Na
#endif
#define __FFLAS__Mb M
#define __FFLAS__Mb2 M
#define __FFLAS__Nb2 N-nsplit
#define __FFLAS__Mbrest M
#define __FFLAS__Nbrest nrestsplit
#define __FFLAS__Mupdate M
#define __FFLAS__Nupdate N - (i + 1) * nsplit
#define __FFLAS__Anorminc __FFLAS__Arowinc
#define __FFLAS__Acopnorminc __FFLAS__Acoprowinc
#define __FFLAS__Bnorminc ldb
#define __FFLAS__Bnormnext 1
#define __FFLAS__Bdim M
#ifdef __FFLAS__UPPER
#define __FFLAS__Atriang A + i * nsplit * (lda + 1)
#define __FFLAS__Aupdate A + i * nsplit * (lda + 1) + nsplit * __FFLAS__Acolinc
#define __FFLAS__Arest A + (__FFLAS__Na - nrestsplit) * (lda + 1)
#define __FFLAS__Anormnext __FFLAS__Acolinc
#define __FFLAS__Acopnormnext __FFLAS__Acopcolinc
#define __FFLAS__Bupdate B + (i + 1) * nsplit
#define __FFLAS__Brec B + i * nsplit
#define __FFLAS__Brest B + (N - nrestsplit)
#define __FFLAS__A1 A
#define __FFLAS__A2 A + nsplit * __FFLAS__Acolinc
#define __FFLAS__A3 A + nsplit * (lda + 1)
#define __FFLAS__B1 B
#define __FFLAS__B2 B + nsplit
#define __FFLAS__Normdim i
#else // __FFLAS__LOWER
#define __FFLAS__Atriang A + (__FFLAS__Na - (i + 1) * nsplit) * (lda + 1)
#define __FFLAS__Aupdate A + (__FFLAS__Na - (i + 1) * nsplit) * __FFLAS__Arowinc
#define __FFLAS__Arest A
#define __FFLAS__Anormnext lda + 1
#define __FFLAS__Acopnormnext __FFLAS__Na + 1
#define __FFLAS__Bupdate B
#define __FFLAS__Brec B + (N - (i + 1) * nsplit)
#define __FFLAS__Brest B
#define __FFLAS__A1 A + (__FFLAS__Na - nsplit) * (lda + 1)
#define __FFLAS__A2 A + (__FFLAS__Na - nsplit) * __FFLAS__Arowinc
#define __FFLAS__A3 A
#define __FFLAS__B1 B + (N - nsplit)
#define __FFLAS__B2 B
#define __FFLAS__Normdim __FFLAS__Na - i -1
#endif
#endif

#ifdef __FFLAS__UP
#define __FFLAS__UPLO Upper
#else
#define __FFLAS__UPLO Lower
#endif

#ifdef __FFLAS__UNIT
#define __FFLAS__DIAG Unit
#else
#define __FFLAS__DIAG NonUnit
#endif

#ifdef __FFLAS__TRANSPOSE
#define __FFLAS__TRANS Trans
#else
#define __FFLAS__TRANS NoTrans
#endif

#ifdef __FFLAS__DOUBLE
#define __FFLAS__ELEMENT double
#define __FFLAS__DOMAIN Givaro::DoubleDomain
#define __FFLAS__BLAS_PREFIX d
#endif

#ifdef __FFLAS__FLOAT
#define __FFLAS__ELEMENT float
#define __FFLAS__DOMAIN Givaro::FloatDomain
#define __FFLAS__BLAS_PREFIX s
#endif

#ifdef __FFLAS__GENERIC
#define __FFLAS__ELEMENT Element
#endif

#ifdef __FFLAS_MULTIPRECISION
#define __FFLAS__ELEMENT FFPACK::rns_double_elt
#define __FFLAS__DOMAIN FFPACK::RNSInteger<FFPACK::rns_double>
#define __FFLAS__BLAS_PREFIX imp
#endif

#ifndef __FFLAS__GENERIC
template <>
class Mjoin(ftrsm, Mjoin(__FFLAS__SIDE, Mjoin(__FFLAS__UPLO, Mjoin(__FFLAS__TRANS, __FFLAS__DIAG))))<__FFLAS__ELEMENT>{
public:

    // TRSM with delayed updates: assumes input in Zp and ensures output in Zp.
    // The multiple MatMul updates (recursive sequence) are done over Z
    template<class Field, class ParSeqTrait>
    void delayed (const Field& F, const size_t M, const size_t N,
#ifdef __FFLAS__TRSM_READONLY
                  typename Field::ConstElement_ptr
#else //__FFLAS__TRSM_READONLY
                  typename Field::Element_ptr
#endif //__FFLAS__TRSM_READONLY
                  A, const size_t lda,
                  typename Field::Element_ptr B, const size_t ldb,
                  const size_t nblas, size_t nbblocsblas,
                  TRSMHelper<StructureHelper::Recursive, ParSeqTrait> & H)

    {
        //static __FFLAS__DOMAIN D(F); // is this safe ??
        __FFLAS__DOMAIN D(F); // is this safe ??

        if ( __FFLAS__Na <= nblas ){
            freduce (F, M, N, B, ldb);

#define __FFLAS__Atrsm A
#define __FFLAS__Atrsm_lda lda
#ifndef __FFLAS__UNIT
#ifdef __FFLAS__TRSM_READONLY
            //! @warning this is C99 (-Wno-vla)
            //typename Field::Element Acop[__FFLAS__Na*__FFLAS__Na];
            typename Field::Element_ptr Acop = FFLAS::fflas_new(F,__FFLAS__Na,__FFLAS__Na);
            typename Field::Element_ptr Acopi = Acop;
#undef __FFLAS__Atrsm
#undef __FFLAS__Atrsm_lda
#define __FFLAS__Atrsm Acop
#define __FFLAS__Atrsm_lda __FFLAS__Na
#endif //__FFLAS__TRSM_READONLY
            typename Field::Element inv;
#ifdef __FFLAS__TRSM_READONLY
            typename Field::ConstElement_ptr
#else //__FFLAS__TRSM_READONLY
            typename Field::Element_ptr
#endif //__FFLAS__TRSM_READONLY
            Ai = A;
            typename Field::Element_ptr Bi = B;
#ifdef __FFLAS__LEFT
#ifdef __FFLAS__UP
            Ai += __FFLAS__Acolinc;
#ifdef __FFLAS__TRSM_READONLY
            Acopi += __FFLAS__Acopcolinc;
#endif //__FFLAS__TRSM_READONLY
#endif //__FFLAS__UP
#endif //__FFLAS__LEFT
#ifdef __FFLAS__RIGHT
#ifdef __FFLAS__LOW
            Ai += __FFLAS__Arowinc;
#ifdef __FFLAS__TRSM_READONLY
            Acopi += __FFLAS__Acoprowinc;
#endif //__FFLAS__TRSM_READONLY
#endif //__FFLAS__LOW
#endif //__FFLAS__RIGHT
            for (size_t i = 0; i < __FFLAS__Na; ++i){
#ifdef _FF_DEBUG
                if ( F.isZero(*(A+i*(lda+1))) ) throw PreconditionFailed(__func__,__FILE__,__LINE__,"Triangular matrix not invertible");
#endif //_FF_DEBUG
                F.inv (inv, *(A + i * (lda+1)));
#ifndef __FFLAS_MULTIPRECISION
#ifdef __FFLAS__TRSM_READONLY
                fscal (F, __FFLAS__Normdim, inv, Ai, __FFLAS__Anorminc, Acopi, __FFLAS__Acopnorminc);
                Acopi += __FFLAS__Acopnormnext;
#else //__FFLAS__TRSM_READONLY
                fscalin (F, __FFLAS__Normdim, inv, Ai, __FFLAS__Anorminc);
#endif //__FFLAS__TRSM_READONLY
#endif //__FFLAS_MULTIPRECISION
                FFLAS::fscalin (F, __FFLAS__Bdim, inv, Bi, __FFLAS__Bnorminc);
                Ai += __FFLAS__Anormnext;
                Bi += __FFLAS__Bnormnext;

            }
#endif // __FFLAS__UNIT
#ifndef __FFLAS_MULTIPRECISION
#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
            openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
            Mjoin(cblas_,Mjoin(__FFLAS__BLAS_PREFIX,trsm))
            (CblasRowMajor,
             Mjoin (Cblas, __FFLAS__SIDE),
             Mjoin (Cblas, __FFLAS__UPLO),
             Mjoin (Cblas, __FFLAS__TRANS),
             CblasUnit,
             (int)M, (int)N, D.one, __FFLAS__Atrsm, (int)__FFLAS__Atrsm_lda, B, (int)ldb );
            freduce (F, M, N, B, ldb);
#endif //__FFLAS_MULTIPRECISION

#ifndef __FFLAS__UNIT
            Ai = A;
#ifdef __FFLAS__LEFT
#ifdef __FFLAS__UP
            Ai += __FFLAS__Acolinc;
#endif //__FFLAS__UP
#endif //__FFLAS__LEFT
#ifdef __FFLAS__RIGHT
#ifdef __FFLAS__LOW
            Ai += __FFLAS__Arowinc;
#endif //__FFLAS__LOW
#endif //__FFLAS__RIGHT

#ifndef __FFLAS__TRSM_READONLY
#ifndef __FFLAS_MULTIPRECISION
            for (size_t i = 0; i < __FFLAS__Na; ++i){
                fscalin( F, __FFLAS__Normdim, *(A + i * (lda+1)) , Ai, __FFLAS__Anorminc);
                Ai += __FFLAS__Anormnext;
            }
#endif //__FFLAS_MULTIPRECISION
#endif //__FFLAS__TRSM_READONLY

#ifdef __FFLAS__TRSM_READONLY
            FFLAS::fflas_delete(Acop);
#endif //__FFLAS__TRSM_READONLY
#endif // __FFLAS__UNIT
        } else { // __FFLAS__Na <= nblas
            size_t nbblocsup = (nbblocsblas + 1) / 2;
            size_t nsplit = nbblocsup * nblas;

            this->delayed (F, __FFLAS__Mb, __FFLAS__Nb,
                           __FFLAS__A1, lda, __FFLAS__B1, ldb, nblas, nbblocsup, H);



#ifdef __FFLAS__RIGHT
            fgemm (D, FflasNoTrans, Mjoin (Fflas, __FFLAS__TRANS),
                   __FFLAS__Mb2, __FFLAS__Nb2, nsplit, D.mOne,
                   __FFLAS__B1, ldb, __FFLAS__A2, lda,
                   F.one, __FFLAS__B2, ldb, H.parseq);
#else
            fgemm (D, Mjoin (Fflas, __FFLAS__TRANS), FflasNoTrans,
                   __FFLAS__Mb2, __FFLAS__Nb2, nsplit, D.mOne,
                   __FFLAS__A2, lda, __FFLAS__B1, ldb,
                   F.one, __FFLAS__B2, ldb, H.parseq);
#endif //__FFLAS__RIGHT

            this->delayed (F, __FFLAS__Mb2, __FFLAS__Nb2,
                           __FFLAS__A3, lda, __FFLAS__B2, ldb, nblas, nbblocsblas - nbblocsup, H);
        }
    }
    template <class Field, class ParSeqTrait>
    void operator () (const Field& F, const size_t M, const size_t N,
#ifdef __FFLAS__TRSM_READONLY
                      typename Field::ConstElement_ptr
#else
                      typename Field::Element_ptr
#endif //__FFLAS__TRSM_READONLY
                      A, const size_t lda,
                      typename Field::Element_ptr B, const size_t ldb,
                      TRSMHelper<StructureHelper::Recursive, ParSeqTrait> & H)
    {
#if defined(__FFLAS_MULTIPRECISION) && defined(BENCH_PERF_FTRSM_MP)
        FFLAS::Timer chrono;chrono.start();
#endif


        if (!M || !N ) return;

        //static __FFLAS__DOMAIN D(F);
        __FFLAS__DOMAIN D(F);
        size_t nblas = TRSMBound<Field> (F);
        size_t ndel = DotProdBoundClassic (F, F.one);
        ndel = (ndel / nblas)*nblas;
        size_t nsplit = ndel;
        size_t nbblocsplit = (__FFLAS__Na-1) / nsplit;
        size_t nrestsplit = ((__FFLAS__Na-1) % nsplit) +1;
        for ( size_t  i = 0; i < nbblocsplit; ++i) {
            this->delayed (F, __FFLAS__Mb, __FFLAS__Nb,
                           __FFLAS__Atriang, lda, __FFLAS__Brec, ldb, nblas, nsplit / nblas, H);

#ifdef __FFLAS__RIGHT
            fgemm (F, FflasNoTrans, Mjoin (Fflas, __FFLAS__TRANS),
                   __FFLAS__Mupdate, __FFLAS__Nupdate, nsplit, F.mOne,
                   __FFLAS__Brec, ldb, __FFLAS__Aupdate, lda,
                   F.one, __FFLAS__Bupdate, ldb, H.parseq);
#else
            fgemm (F, Mjoin (Fflas, __FFLAS__TRANS),  FflasNoTrans,
                   __FFLAS__Mupdate, __FFLAS__Nupdate, nsplit, F.mOne,
                   __FFLAS__Aupdate, lda, __FFLAS__Brec, ldb,
                   F.one, __FFLAS__Bupdate, ldb, H.parseq);
#endif //__FFLAS__RIGHT
        }
        if (nrestsplit)
            this->delayed (F, __FFLAS__Mbrest, __FFLAS__Nbrest,
                           __FFLAS__Arest, lda, __FFLAS__Brest, ldb, nblas, nrestsplit / nblas, H);

#if defined(__FFLAS_MULTIPRECISION) && defined(BENCH_PERF_FTRSM_MP)
        chrono.stop();
        F.t_trsm+=chrono.usertime();
#endif

    }


}; //class ftrsm....

#else // __FFLAS__GENERIC

template <class Element>
class Mjoin(ftrsm, Mjoin(__FFLAS__SIDE, Mjoin(__FFLAS__UPLO, Mjoin(__FFLAS__TRANS, __FFLAS__DIAG)))) {
public:

    template<class Field, class ParSeqTrait>
    void operator()	(const Field& F, const size_t M, const size_t N,
#ifdef __FFLAS__TRSM_READONLY
                     typename Field::ConstElement_ptr
#else
                     typename Field::Element_ptr
#endif
                     A, const size_t lda,
                     typename Field::Element_ptr B, const size_t ldb,
                     TRSMHelper<StructureHelper::Recursive, ParSeqTrait> & H)
    {
        if (__FFLAS__Na == 1){

#ifndef __FFLAS__UNIT
            typename Field::Element inv;
            F.init(inv);
#ifdef _FF_DEBUG
            if ( F.isZero(*A) ) throw PreconditionFailed(__func__,__FILE__,__LINE__,"Triangular matrix not invertible");
#endif //_FF_DEBUG
            F.inv(inv, *A);
            FFLAS::fscalin(F, __FFLAS__Bdim, inv, B, __FFLAS__Bnorminc);

#endif //__FFLAS__UNIT
        } else { // __FFLAS__Na > 1
            size_t nsplit = __FFLAS__Na >> 1;
            this->operator() (F, __FFLAS__Mb, __FFLAS__Nb, __FFLAS__A1, lda, __FFLAS__B1, ldb, H);
#ifdef __FFLAS__RIGHT
            fgemm (F, FflasNoTrans , Mjoin (Fflas, __FFLAS__TRANS),
                   __FFLAS__Mb2, __FFLAS__Nb2, nsplit, F.mOne,
                   __FFLAS__B1, ldb, __FFLAS__A2, lda,
                   F.one, __FFLAS__B2, ldb, H.parseq);
#else //__FFLAS__RIGHT
            fgemm (F, Mjoin (Fflas, __FFLAS__TRANS), FFLAS::FflasNoTrans,
                   __FFLAS__Mb2, __FFLAS__Nb2, nsplit, F.mOne,
                   __FFLAS__A2, lda, __FFLAS__B1, ldb,
                   F.one, __FFLAS__B2, ldb, H.parseq);
#endif //__FFLAS__RIGHT
            this->operator() (F, __FFLAS__Mb2, __FFLAS__Nb2, __FFLAS__A3, lda, __FFLAS__B2, ldb, H);
        }
    }
};

#endif // __FFLAS__GENERIC

#ifdef __FFLAS__LOWER
#undef __FFLAS__LOWER
#else
#undef __FFLAS__UPPER
#endif
#undef __FFLAS__UPLO
#undef __FFLAS__DIAG
#undef __FFLAS__SIDE
#undef __FFLAS__TRANS
#undef __FFLAS__Na
#undef __FFLAS__Mb
#undef __FFLAS__Nb
#undef __FFLAS__Mbrest
#undef __FFLAS__Nbrest
#undef __FFLAS__Mupdate
#undef __FFLAS__Nupdate
#undef __FFLAS__Atriang
#undef __FFLAS__Aupdate
#undef __FFLAS__Arest
#undef __FFLAS__Bupdate
#undef __FFLAS__Brec
#undef __FFLAS__Brest
#undef __FFLAS__Bnorminc
#undef __FFLAS__Bnormnext
#undef __FFLAS__Anormnext
#undef __FFLAS__Acopnormnext
#undef __FFLAS__Anorminc
#undef __FFLAS__Acopnorminc
#undef __FFLAS__ELEMENT
#undef __FFLAS__BLAS_PREFIX
#undef __FFLAS__DOMAIN
#undef __FFLAS__A1
#undef __FFLAS__A2
#undef __FFLAS__A3
#undef __FFLAS__B1
#undef __FFLAS__B2
#undef __FFLAS__Nb2
#undef __FFLAS__Mb2
#undef __FFLAS__Bdim
#undef __FFLAS__Normdim
#undef __FFLAS__Acolinc
#undef __FFLAS__Arowinc
#undef __FFLAS__Acopcolinc
#undef __FFLAS__Acoprowinc
#undef __FFLAS__Atrsm_lda
#undef __FFLAS__Atrsm
#undef Mjoin
#undef my_join
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
