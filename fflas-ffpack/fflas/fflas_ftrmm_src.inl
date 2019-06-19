/* Copyright (C) 2005 C. Pernet
 * Written by C. Pernet
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
 *
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
#define __FFLAS__Mb nsplit
#define __FFLAS__Nb2 N
#define __FFLAS__Mb2 M-nsplit
#define __FFLAS__Mbrest nrestsplit
#define __FFLAS__Nbrest N
#define __FFLAS__Mupdate nrestsplit + i * nsplit
#define __FFLAS__Nupdate N
#define __FFLAS__Bdim N
#define __FFLAS__Bnorminc 1
#define __FFLAS__Cnorminc 1
#ifdef __FFLAS__LOWER
#define __FFLAS__Atriang  A + (nbblocsplit - (i + 1)) * nsplit * (lda + 1)
#define __FFLAS__Aupdate __FFLAS__Atriang + nsplit * __FFLAS__Arowinc
#define __FFLAS__Arest A + nbblocsplit * nsplit * (lda+1)
#define __FFLAS__Brec B + (nbblocsplit - (i+1)) * nsplit * ldb
#define __FFLAS__Crec C + (nbblocsplit - (i+1)) * nsplit * ldc
#define __FFLAS__Bupdate B + (nbblocsplit - i) * nsplit * ldb
#define __FFLAS__Cupdate C + (nbblocsplit - i) * nsplit * ldc
#define __FFLAS__Brest B + nbblocsplit * nsplit * ldb
#define __FFLAS__Crest C + nbblocsplit * nsplit * ldc
#define __FFLAS__A1 A + (nsplit) * (lda + 1)
#define __FFLAS__A2 A + (nsplit) * __FFLAS__Arowinc
#define __FFLAS__A3 A
#define __FFLAS__B1 B + (nsplit) * ldb
#define __FFLAS__C1 C + (nsplit) * ldc
#define __FFLAS__B2 B
#define __FFLAS__C2 C
#else
#define __FFLAS__Atriang A + (nrestsplit + i * nsplit) * (lda + 1)
#define __FFLAS__Aupdate A + (nrestsplit + i * nsplit) * __FFLAS__Acolinc
#define __FFLAS__Arest A
#define __FFLAS__Brec B + (nrestsplit + i * nsplit) * ldb
#define __FFLAS__Crec C + (nrestsplit + i * nsplit) * ldc
#define __FFLAS__Bupdate B
#define __FFLAS__Cupdate C
#define __FFLAS__Brest B
#define __FFLAS__Crest C
#define __FFLAS__A1 A
#define __FFLAS__A2 A + (M-nsplit) * __FFLAS__Acolinc
#define __FFLAS__A3 A + (M-nsplit) * (lda + 1)
#define __FFLAS__B1 B
#define __FFLAS__C1 C
#define __FFLAS__B2 B + (M-nsplit) * ldb
#define __FFLAS__C2 C + (M-nsplit) * ldc
#endif
#else
#define __FFLAS__SIDE Right
#define __FFLAS__Na N
#define __FFLAS__Nb nsplit
#define __FFLAS__Mb M
#define __FFLAS__Mb2 M
#define __FFLAS__Nb2 N-nsplit
#define __FFLAS__Mbrest M
#define __FFLAS__Nbrest nrestsplit
#define __FFLAS__Mupdate M
#define __FFLAS__Nupdate nrestsplit + i * nsplit
#define __FFLAS__Bdim M
#define __FFLAS__Bnorminc ldb
#define __FFLAS__Cnorminc ldc
#ifdef __FFLAS__UPPER
#define __FFLAS__Atriang A + (nbblocsplit - (i + 1)) * nsplit * (lda + 1)
#define __FFLAS__Aupdate __FFLAS__Atriang + nsplit * __FFLAS__Acolinc
#define __FFLAS__Arest A + nbblocsplit * nsplit * (lda+1)
#define __FFLAS__Brec B + (nbblocsplit - (i+1)) * nsplit
#define __FFLAS__Crec C + (nbblocsplit - (i+1)) * nsplit
#define __FFLAS__Bupdate B + (nbblocsplit - i) * nsplit
#define __FFLAS__Cupdate C + (nbblocsplit - i) * nsplit
#define __FFLAS__Brest B + nbblocsplit * nsplit
#define __FFLAS__Crest C + nbblocsplit * nsplit
#define __FFLAS__A1 A + (nsplit) * (lda + 1)
#define __FFLAS__A2 A + (nsplit) * __FFLAS__Acolinc
#define __FFLAS__A3 A
#define __FFLAS__B1 B + nsplit
#define __FFLAS__C1 C + nsplit
#define __FFLAS__B2 B
#define __FFLAS__C2 C
#else
#define __FFLAS__Atriang A + (nrestsplit + i * nsplit) * (lda + 1)
#define __FFLAS__Aupdate A + (nrestsplit + i * nsplit) * __FFLAS__Arowinc
#define __FFLAS__Arest A
#define __FFLAS__Brec B + (nrestsplit + i * nsplit)
#define __FFLAS__Crec C + (nrestsplit + i * nsplit)
#define __FFLAS__Bupdate B
#define __FFLAS__Cupdate C
#define __FFLAS__Brest B
#define __FFLAS__Crest C
#define __FFLAS__A1 A
#define __FFLAS__A2 A + (N-nsplit) * __FFLAS__Arowinc
#define __FFLAS__A3 A + (N-nsplit) * (lda + 1)
#define __FFLAS__B1 B
#define __FFLAS__C1 C
#define __FFLAS__B2 B + N-nsplit
#define __FFLAS__C2 C + N-nsplit
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



#ifndef __FFLAS__GENERIC
template <>
class Mjoin(ftrmm, Mjoin(__FFLAS__SIDE, Mjoin(__FFLAS__UPLO, Mjoin(__FFLAS__TRANS, __FFLAS__DIAG))))<__FFLAS__ELEMENT>{
public:

    template <class Field>
    void delayed (const Field& F, const size_t M, const size_t N,
                  typename Field::ConstElement_ptr A, const size_t lda,
                  typename Field::Element_ptr B, const size_t ldb)
    {
#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif
        Mjoin(cblas_,Mjoin(__FFLAS__BLAS_PREFIX,trmm))
        (CblasRowMajor,
         Mjoin (Cblas, __FFLAS__SIDE),
         Mjoin (Cblas, __FFLAS__UPLO),
         Mjoin (Cblas, __FFLAS__TRANS),
         Mjoin (Cblas, __FFLAS__DIAG),
         (int)M, (int)N, 1.0, A, (int)lda, B, (int)ldb );
        freduce(F, M, N, B, ldb);
    }

    template <class Field>
    void delayed (const Field& F, const size_t M, const size_t N,
                  typename Field::ConstElement_ptr A, const size_t lda,
                  typename Field::ConstElement_ptr B, const size_t ldb,
                  const typename Field::Element beta,
                  typename Field::Element_ptr C, const size_t ldc)
    {
#if defined(__FFLASFFPACK_OPENBLAS_NUM_THREADS) and not defined (__FFLASFFPACK_OPENBLAS_NT_ALREADY_SET)
        openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

        // copy before call to BLAS
        typename Field::Element_ptr tmp=NULL;
        if (!F.isZero(beta)){
            tmp=FFLAS::fflas_new(F,M,N);
            fscal (F, M, N, beta, C, ldc, tmp, N);
        }
        fassign (F, M, N, B, ldb, C, ldc);

        Mjoin(cblas_,Mjoin(__FFLAS__BLAS_PREFIX,trmm))
        (CblasRowMajor,
         Mjoin (Cblas, __FFLAS__SIDE),
         Mjoin (Cblas, __FFLAS__UPLO),
         Mjoin (Cblas, __FFLAS__TRANS),
         Mjoin (Cblas, __FFLAS__DIAG),
         (int)M, (int)N, 1.0, A, (int)lda, C, (int)ldc);
        freduce(F, M, N, C, ldc);
        // TODO: reduce number of modulo by delayed some adds / scal, etc
        if (!F.isZero(beta)){
            faddin (F, M, N, tmp,N, C, ldc);
            fflas_delete(tmp);
        }
    }

    template <class Field>
    void operator () (const Field& F, const size_t M, const size_t N,
                      typename Field::ConstElement_ptr A, const size_t lda,
                      typename Field::Element_ptr B, const size_t ldb)
    {

        if (!M || !N ) return;

        size_t nsplit = DotProdBoundClassic (F,  F.one);
        size_t nbblocsplit = (__FFLAS__Na-1) / nsplit;
        size_t nrestsplit = ((__FFLAS__Na-1) % nsplit) +1;
        FFLASFFPACK_check(__FFLAS__Na == nsplit*nbblocsplit+nrestsplit);

        if (nrestsplit)
            this->delayed (F, __FFLAS__Mbrest, __FFLAS__Nbrest,
                           __FFLAS__Arest, lda, __FFLAS__Brest, ldb);

        for ( size_t  i = 0; i < nbblocsplit; ++i) {

#ifdef __FFLAS__RIGHT
            fgemm (F, FflasNoTrans, Mjoin (Fflas, __FFLAS__TRANS),
                   __FFLAS__Mupdate, __FFLAS__Nupdate, nsplit, F.one,
                   __FFLAS__Brec, ldb, __FFLAS__Aupdate, lda, F.one, __FFLAS__Bupdate, ldb);
#else
            fgemm (F, Mjoin (Fflas, __FFLAS__TRANS),  FflasNoTrans,
                   __FFLAS__Mupdate, __FFLAS__Nupdate, nsplit, F.one,
                   __FFLAS__Aupdate, lda, __FFLAS__Brec, ldb, F.one, __FFLAS__Bupdate, ldb);
#endif

            this->delayed (F, __FFLAS__Mb, __FFLAS__Nb,
                           __FFLAS__Atriang, lda, __FFLAS__Brec, ldb);


        }
    }

    template <class Field>
    void operator () (const Field& F, const size_t M, const size_t N,
                      typename Field::ConstElement_ptr A, const size_t lda,
                      typename Field::ConstElement_ptr B, const size_t ldb,
                      const typename Field::Element beta,
                      typename Field::Element_ptr C, const size_t ldc)
    {

        if (!M || !N ) return;

        size_t nsplit = DotProdBoundClassic (F,  F.one);
        size_t nbblocsplit = (__FFLAS__Na-1) / nsplit;
        size_t nrestsplit = ((__FFLAS__Na-1) % nsplit) +1;
        FFLASFFPACK_check(__FFLAS__Na == nsplit*nbblocsplit+nrestsplit);

        if (nrestsplit)
            this->delayed (F, __FFLAS__Mbrest, __FFLAS__Nbrest,
                           __FFLAS__Arest, lda, __FFLAS__Brest, ldb, beta, __FFLAS__Crest, ldc);

        for ( size_t  i = 0; i < nbblocsplit; ++i) {

#ifdef __FFLAS__RIGHT
            fgemm (F, FflasNoTrans, Mjoin (Fflas, __FFLAS__TRANS),
                   __FFLAS__Mupdate, __FFLAS__Nupdate, nsplit, F.one,
                   __FFLAS__Brec, ldb, __FFLAS__Aupdate, lda, F.one, __FFLAS__Cupdate, ldc);
#else
            fgemm (F, Mjoin (Fflas, __FFLAS__TRANS),  FflasNoTrans,
                   __FFLAS__Mupdate, __FFLAS__Nupdate, nsplit, F.one,
                   __FFLAS__Aupdate, lda, __FFLAS__Brec, ldb, F.one, __FFLAS__Cupdate, ldc);
#endif

            this->delayed (F, __FFLAS__Mb, __FFLAS__Nb,
                           __FFLAS__Atriang, lda, __FFLAS__Brec, ldb, beta, __FFLAS__Crec, ldc);


        }
    }

}; //class ftrmm....

#else // __FFLAS__GENERIC

template <class Element>
class Mjoin(ftrmm, Mjoin(__FFLAS__SIDE, Mjoin(__FFLAS__UPLO, Mjoin(__FFLAS__TRANS, __FFLAS__DIAG)))) {
public:

    template<class Field>
    void operator() (const Field& F, const size_t M, const size_t N,
                     typename Field::ConstElement_ptr A, const size_t lda,
                     typename Field::Element_ptr B, const size_t ldb)
    {

        if (__FFLAS__Na == 1)
#ifdef __FFLAS__NONUNIT
            fscalin(F, __FFLAS__Bdim, *A, B, __FFLAS__Bnorminc);
#else
        ;
#endif

        else { // __FFLAS__Na > 1
            size_t nsplit = __FFLAS__Na >> 1;
            this->operator() (F, __FFLAS__Mb2, __FFLAS__Nb2, __FFLAS__A1, lda, __FFLAS__B1, ldb);

#ifdef __FFLAS__RIGHT
            fgemm (F, FflasNoTrans , Mjoin (Fflas, __FFLAS__TRANS),
                   __FFLAS__Mb2, __FFLAS__Nb2, nsplit, F.one,
                   __FFLAS__B2, ldb, __FFLAS__A2, lda, F.one, __FFLAS__B1, ldb);
#else
            fgemm (F, Mjoin (Fflas, __FFLAS__TRANS), FflasNoTrans,
                   __FFLAS__Mb2, __FFLAS__Nb2, nsplit, F.one,
                   __FFLAS__A2, lda, __FFLAS__B2, ldb, F.one, __FFLAS__B1, ldb);
#endif
            this->operator() (F, __FFLAS__Mb, __FFLAS__Nb, __FFLAS__A3, lda, __FFLAS__B2, ldb);
        }
    }
    template<class Field>
    void operator() (const Field& F, const size_t M, const size_t N,
                     typename Field::ConstElement_ptr A, const size_t lda,
                     typename Field::ConstElement_ptr B, const size_t ldb,
                     const typename Field::Element beta,
                     typename Field::Element_ptr C, const size_t ldc)
    {

        if (__FFLAS__Na == 1){
            if (!F.isZero(beta))
                fscalin (F, __FFLAS__Bdim, beta, C, __FFLAS__Cnorminc);
#ifdef __FFLAS__NONUNIT
            //TODO: replace by faxpby when available
            faxpy (F, __FFLAS__Bdim, *A, B, __FFLAS__Bnorminc, C, __FFLAS__Cnorminc);
#else
            faddin (F, __FFLAS__Bdim, B, __FFLAS__Bnorminc, C, __FFLAS__Cnorminc);
#endif
        }
        else { // __FFLAS__Na > 1
            size_t nsplit = __FFLAS__Na >> 1;
            this->operator() (F, __FFLAS__Mb2, __FFLAS__Nb2, __FFLAS__A1, lda, __FFLAS__B1, ldb, beta, __FFLAS__C1, ldc);

#ifdef __FFLAS__RIGHT
            fgemm (F, FflasNoTrans , Mjoin (Fflas, __FFLAS__TRANS),
                   __FFLAS__Mb2, __FFLAS__Nb2, nsplit, F.one,
                   __FFLAS__B2, ldb, __FFLAS__A2, lda, F.one, __FFLAS__C1, ldc);
#else
            fgemm (F, Mjoin (Fflas, __FFLAS__TRANS), FflasNoTrans,
                   __FFLAS__Mb2, __FFLAS__Nb2, nsplit, F.one,
                   __FFLAS__A2, lda, __FFLAS__B2, ldb, F.one, __FFLAS__C1, ldc);
#endif
            this->operator() (F, __FFLAS__Mb, __FFLAS__Nb, __FFLAS__A3, lda, __FFLAS__B2, ldb, beta, __FFLAS__C2, ldc);
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
#undef __FFLAS__Cupdate
#undef __FFLAS__Brec
#undef __FFLAS__Crec
#undef __FFLAS__Brest
#undef __FFLAS__Crest
#undef __FFLAS__ELEMENT
#undef __FFLAS__BLAS_PREFIX
#undef __FFLAS__DOMAIN
#undef __FFLAS__A1
#undef __FFLAS__A2
#undef __FFLAS__A3
#undef __FFLAS__B1
#undef __FFLAS__C1
#undef __FFLAS__B2
#undef __FFLAS__C2
#undef __FFLAS__Nb2
#undef __FFLAS__Mb2
#undef __FFLAS__Bdim
#undef __FFLAS__Acolinc
#undef __FFLAS__Arowinc
#undef __FFLAS__Bnorminc
#undef __FFLAS__Cnorminc
#undef Mjoin
#undef my_join
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
