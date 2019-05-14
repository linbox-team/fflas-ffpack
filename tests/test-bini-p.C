/*
 * Copyright (C) 2015 the FFLAS-FFPACK group
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *
 * This file is Free Software and part of FFLAS-FFPACK.
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


#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "assert.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/flimits.h"

#include <givaro/udl.h>

// using namespace FFPACK;
#define NEWWINO
// #define NOTRANDOM
//
#define DIVIDE_INTO(x,y) (((x) + (y) - 1)/(y))

const int algos = 6 ;
const int algos_k = 2 ;

using Givaro::Modular;
using Givaro::ModularBalanced;
using Givaro::Timer;
using FFLAS::FieldTraits;
typedef std::vector<Timer> time_v ;
typedef std::vector<int> int_v ;

const int selec[] = {
    0
    ,1
    ,2
    ,3
    ,4
    ,5
};

const int selec_k[] = {
    0
    ,1
};

const char * descr[] = {
    "322 low mem"
    , "322 first 1"
    , "322 4 tmp  "
    , "223 low mem"
    , "232 first 1"
    , "232 all tmp"
    , "comp left  "
    , "comp right "
    // , "322 sqrt   "
};

const char * descr_k[] = {
    "comp left  "
    , "comp right "
};

namespace FFLAS { /*  compression */

    template<class Elem, int Num>
    struct Packer ;

    template<>
    struct Packer<double,2> {
        uint64_t bits = (limits<double>::digits()/2) ;
        double   base = (double) (1_ui64 << bits) ;
        uint64_t mask = (1_ui64 << bits) - 1_ui64 ;

        template<class T>
        void accu(double * p, T * w) {
            *p *= base ;
            *p += (double)*w ;
        }
    } ;


    /* ****** */
    /*  pack  */
    /* ****** */

    /*  pack nb words (a,b,c) -> [a|b|c] */
    template<class wide_T, class pack_T, int Nb>
    void pack_word( pack_T * packed,
                    const wide_T * words, int32_t stride,
                    Packer<pack_T,Nb> & packer) ;


    template<class wide_T>
    void pack_word/*<wide_T,double,2>*/( double * packed,
                                         const wide_T * words, int32_t stride,
                                         Packer<double,2> & packer)
    {
        // std::cout << "pack " << *words << '+' << *(words+stride) << " * " << (uint64_t) packer.base << " = ";
        // words += stride ;
        *packed = (double) *words ;
        words += stride ;
        packer.accu(packed,words);
        // std::cout << (uint64_t) *packed << std::endl;
    }

    /*  pack nb words (a,b) -> [a|b|0]  filling with zeros */
    template<class wide_T, class pack_T, int Nb>
    void pack_word_part( pack_T * packed, int32_t nb,
                         const wide_T * words, int32_t stride,
                         Packer<pack_T,Nb> & packer) ;

    template<class wide_T>
    void pack_word_part/* <wide_T,double,2> */( double * packed, int32_t nb,
                                                const wide_T * words, int32_t stride,
                                                Packer<double,2> & packer)
    {
        assert(nb == 1);
        *packed = (double) *words ;
        // words += stride ;
        // packer.accu(packed,words);
        *packed *= packer.base ;
    }

    /* ****** */
    /* unpack */
    /* ****** */

    template<class wide_T, class pack_T, int Nb>
    void unpack_word( wide_T * words, int32_t stride,
                      const pack_T * packed,
                      Packer<pack_T,Nb> & packer);

    template<class wide_T>
    void unpack_word/* <wide_T,double,2> */( wide_T * words, int32_t stride,
                                             const double * packed,
                                             Packer<double ,2> & packer)
    {
        uint64_t pck = (uint64_t) *packed ;
        words += stride ;
        *words = (double) (pck & packer.mask) ;
        words -= stride ;
        pck >>= packer.bits ;
        *words = (double) pck /*  & packer.mask  */ ;
    }


    template<class wide_T, class pack_T, int Nb>
    void unpack_word_part( wide_T * words, int32_t stride,
                           const pack_T * packed, int32_t nb,
                           Packer<pack_T,Nb> & packer);

    template<class wide_T>
    void unpack_word_part/* <wide_T,double,2> */( wide_T * words, int32_t stride,
                                                  const double * packed, int32_t nb,
                                                  Packer<double,2> & packer)
    {
        assert(nb == 1);
        words += stride ;
        *words = 0 ;
        words -= stride ;
        uint64_t pck = (uint64_t) *packed ;
        pck >>= packer.bits ;
        *words =  (double)pck /*  & packer.mask  */ ;
    }

    /* ****** */
    /*  pack  */
    /* ****** */

    template<class wide_T, class pack_T, int Nb, bool row_packed>
    void pack_matrix( pack_T * packed, int32_t row_p, int32_t col_p, int32_t ldm_p,
                      const wide_T * elemts, int32_t row_e, int32_t col_e, int32_t ldm_e,
                      Packer<pack_T,Nb> & packer)
    {
        if (row_packed == true) {
            for (int32_t i = 0 ; i < row_e ;  i++ ) {
                const wide_T * e_p = elemts + i * ldm_e ;
                pack_T * p_p = packed + i * ldm_p ;
                int32_t j = 0 ;
                for ( ; j < col_e/Nb*Nb ;  j+=Nb, e_p+=Nb, p_p++) {
                    pack_word<wide_T>(p_p,e_p,1,packer);

                }
                if (j < col_e)
                    pack_word_part<wide_T>(p_p,col_e-j,e_p,1,packer);
            }
        }
        else { /*  col_packed */
            int32_t i = 0 ;
            int32_t ii = 0 ;
            for ( ; i < row_e/Nb*Nb ;  i += Nb , ii++) {
                const wide_T * e_p = elemts + i * ldm_e ;
                pack_T * p_p = packed + ii * ldm_p ;
                for (int32_t j = 0 ; j < col_e ;  j++, e_p++, p_p++) {
                    pack_word<wide_T>(p_p,e_p,ldm_e,packer);

                }
            }
            if (i < row_e)
                pack_word_part<wide_T>(packed+i*ldm_p,row_e-i,elemts+ii*ldm_e,ldm_e,packer);

        }
    }

    /* ****** */
    /* unpack */
    /* ****** */

    template<class wide_T, class pack_T, int Nb, bool row_packed>
    void unpack_matrix( wide_T * elemts, int32_t row_e, int32_t col_e, int32_t ldm_e,
                        const pack_T * packed, int32_t row_p, int32_t col_p, int32_t ldm_p,
                        Packer<pack_T,Nb> & packer)
    {
        if (row_packed == true) {
            for (int32_t i = 0 ; i < row_e ;  i++ ) {
                wide_T * e_p = elemts + i * ldm_e ;
                const pack_T * p_p = packed + i * ldm_p ;
                int32_t j = 0 ;
                for ( ; j < col_e/Nb*Nb ;  j+=Nb, e_p+=Nb, p_p++) {
                    unpack_word<wide_T>(e_p,1,p_p,packer);

                }
                if (j < col_e)
                    unpack_word_part<wide_T>(e_p,1,p_p,col_e-j,packer);
            }
        }
        else { /*  col_packed */
            int32_t i = 0 ;
            int32_t ii = 0 ;
            for ( ; i < row_e/Nb*Nb ;  i += Nb , ii++) {
                wide_T * e_p = elemts + i * ldm_e ;
                const pack_T * p_p = packed + ii * ldm_p ;
                for (int32_t j = 0 ; j < col_e ;  j++, e_p++, p_p++) {
                    unpack_word<wide_T>(e_p,ldm_e,p_p,packer);

                }
            }
            if (i < row_e)
                unpack_word_part<wide_T>(elemts+i*ldm_e,ldm_e,packed+ii*ldm_p,row_e-i,packer);

        }
    }

    /*  compress A */
    template<class Field, bool left_compress >
    void
    fgemm_compressed(const Field & F,
                     int m, int n, int k,
                     const typename Field::Element * A, int lda,
                     const typename Field::Element * B, int ldb,
                     typename Field::Element * C, int ldc
                    )
    {
        Givaro::ZRing<double>   NoField;
        double * A_k, * B_k, * C_k ;

        typedef typename Field::Element elem_t ;
        Packer<elem_t,2> packer ;

        int m_k = m , n_k = n , lda_k = lda, ldb_k = ldb, ldc_k = ldc ;
        if (left_compress) {
            m_k = DIVIDE_INTO(m,2)*2 ;
            lda_k = m_k ;
            ldc_k = n ;

            A_k =  FFLAS::fflas_new<double>(m_k*k) ;
            //!@bug don't zero all, just the "border"
            FFLAS::fzero(NoField,m_k,k,A_k,k);

            B_k = const_cast<typename Field::Element *>(B) ;

            pack_matrix<elem_t,elem_t,2,false>(A_k,m_k,k,lda_k,
                                               A,m,k,lda,
                                               packer);
        }
        else {
            n_k = DIVIDE_INTO(n,2)*2 ;
            ldb_k = n_k ;
            ldc_k = n_k ;

            A_k = const_cast<typename Field::Element *>(A) ;
            B_k = FFLAS::fflas_new<double>(k*n_k)  ;
            //!@bug don't zero all, just the "border"
            FFLAS::fzero(NoField,k,n_k,B_k,n_k);

            pack_matrix<elem_t,elem_t,2,true>(B_k,k,n_k,ldb_k,
                                              B,k,n,ldb,
                                              packer);
        }

        C_k = FFLAS::fflas_new<double>(m_k*n_k) ;
        //!@bug don't zero all, just the "border"
        FFLAS::fzero(NoField,m_k,n_k,C_k,n_k);

        pack_matrix<elem_t,elem_t,2,!left_compress>(C_k,m_k,n_k,ldc_k,
                                                    C,m,n,ldc,
                                                    packer);

#if 0
        double * C_e = FFLAS::fflas_new<double>(m*ldc);
        unpack_matrix<elem_t,elem_t,2,!left_compress>(C_e,m,n,ldc,
                                                      C_k,m_k,n_k,ldc_k,
                                                      packer);

        int	faux = 0 ;
        for (int i = 0 ; i < m ; ++i) {
            for (int j = 0 ; j < n ; ++j) {
                if (! (C[i*ldc+j] == C_e[i*ldc+j]) ) {
                    ++faux ;
                }
            }
        }
        if (faux) {
            std::cout << "bad pack/unpack ; bad/all = " << faux << '/' << m*n << " ~~  " << (double)faux/(double)(m*n) << std::endl;
        }

        if (faux && (n<20)) {
            std::cout << "IN " << std::endl;
            for (int i = 0 ; i < m ; ++i) {
                for (int j = 0 ; j < n ; ++j)
                    std::cout << C[i*ldc+j] << ' ';
                std::cout << std::endl;
            }
            std::cout << "OUT" << std::endl;
            for (int i = 0 ; i < m ; ++i) {
                for (int j = 0 ; j < n ; ++j)
                    std::cout << C_e[i*ldc+j] << ' ';
                std::cout << std::endl;
            }


        }

        if (faux)
            exit(-1);
#endif




        Givaro::DoubleDomain G ;

        fgemm(G,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
              m_k,n_k,k, 1, A_k,lda_k, B_k,ldb_k, 0, C_k, ldc_k);

        // cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans,
        // m_k,n_k,k, 1, A_k,lda_k, B_k,ldb_k, 0, C_k, ldc_k);


        unpack_matrix<elem_t,elem_t,2,!left_compress>(C,m,n,ldc,
                                                      C_k,m_k,n_k,ldc_k,
                                                      packer);

        if (left_compress)
            FFLAS::fflas_delete(A_k);
        else
            FFLAS::fflas_delete(B_k);
        FFLAS::fflas_delete(C_k);
    }

}

namespace FFLAS { /*  tools */


    template<class Field>
    void finit_fuzzy(Field & F, size_t m, size_t n, double * C, size_t ldc)
    {


        if (n == ldc)
            // FFLAS::vectorised::modp<true,true>(C,C,m*n,p,invp,0,p-1);
            FFLAS::vectorised::modp<Field,true>(F,C,m*n,C);
        else
            for (size_t i = 0 ; i < m ; ++i)
                // FFLAS::vectorised::modp<true,true>(C+i*ldc,C+i*ldc,n,p,invp,0,p-1);
                FFLAS::vectorised::modp<Field,true>(F,C+i*ldc,n,C+i*ldc);
    }


    // C = a*A + B
    void add(const size_t m, const size_t n,
             double a,
             const double *A, const size_t lda,
             const double *B, const size_t ldb,
             double *C, const size_t ldc)
    {
        const double *Ai = A,*Bi = B;
        double *Ci       = C;
        for (;Ai < A+m*lda ; Ai+=lda,Bi+=ldb,Ci+=ldc)
            for (size_t j = 0 ; j < n ; ++j)
                Ci[j] = a * Ai[j] + Bi[j];
    }

    // C = C-(A+B)
    void subadd(const size_t m, const size_t n,
                const double *A, const size_t lda,
                const double *B, const size_t ldb,
                double *C, const size_t ldc)
    {
        const double *Ai = A,*Bi = B;
        double *Ci       = C;
        for (;Ai < A+m*lda ; Ai+=lda,Bi+=ldb,Ci+=ldc)
            for (size_t j = 0 ; j < n ; ++j) {
                Ci[j] = Ci[j] - Ai[j] - Bi[j] ;
            }

    }

    // C = -(A+B)
    void negadd(const size_t m, const size_t n,
                const double *A, const size_t lda,
                const double *B, const size_t ldb,
                double *C, const size_t ldc)
    {
        const double *Ai = A,*Bi = B;
        double *Ci       = C;
        for (;Ai < A+m*lda ; Ai+=lda,Bi+=ldb,Ci+=ldc)
            for (size_t j = 0 ; j < n ; ++j) {
                Ci[j] =  - Ai[j] - Bi[j] ;
            }

    }


    // C = C+A-B
    void addsub(const size_t m, const size_t n,
                const double *A, const size_t lda,
                const double *B, const size_t ldb,
                double *C, const size_t ldc)
    {
        const double *Ai = A,*Bi = B;
        double *Ci       = C;
        for (;Ai < A+m*lda ; Ai+=lda,Bi+=ldb,Ci+=ldc)
            for (size_t j = 0 ; j < n ; ++j) {
                Ci[j] = Ci[j] + Ai[j] - Bi[j] ;
            }

    }


    // C = (C+B)/e
    template<class Field>
    void addscalinf(const Field & F, const size_t m, const size_t n,
                    const double *B, const size_t ldb,
                    double e,
                    double *C, const size_t ldc)
    {
        const double * Bi = B;
        double * Ci = C;
        for (;Bi < B+m*ldb ; Ci+=ldc, Bi += ldb)
            for (size_t j = 0 ; j < n ; ++j)
                Ci[j]= (Ci[j]+Bi[j])*e ;
        // F.init( Ci[j], (Ci[j]+Bi[j])/e );

    }

    // C = (C-B)/e
    template<class Field>
    void subscalinf(const Field & F, const size_t m, const size_t n,
                    const double *B, const size_t ldb,
                    double e,
                    double *C, const size_t ldc)
    {
        const double * Bi = B;
        double * Ci = C;
        for (;Bi < B+m*ldb ; Ci+=ldc, Bi += ldb)
            for (size_t j = 0 ; j < n ; ++j)
                Ci[j]= (Ci[j]-Bi[j])*e ;
        // F.init( Ci[j], (Ci[j]-Bi[j])/e );

    }

    // C = (D-B)/e
    template<class Field>
    void subscal(const Field & F, const size_t m, const size_t n,
                 const double *D, const size_t ldd,
                 const double *B, const size_t ldb,
                 double e,
                 double *C, const size_t ldc)
    {
        const double * Bi = B;
        const double * Di = D;
        double * Ci = C;
        for (;Bi < B+m*ldb ; Ci+=ldc, Bi += ldb, Di += ldd)
            for (size_t j = 0 ; j < n ; ++j)
                Ci[j] = (Di[j]-Bi[j])*e ;

    }

    // C = (D+B)/e
    template<class Field>
    void addscal(const Field & F, const size_t m, const size_t n,
                 const double *D, const size_t ldd,
                 const double *B, const size_t ldb,
                 double e,
                 double *C, const size_t ldc)
    {
        const double * Bi = B;
        const double * Di = D;
        double * Ci = C;
        for (;Bi < B+m*ldb ; Ci+=ldc, Bi += ldb, Di += ldd)
            for (size_t j = 0 ; j < n ; ++j)
                Ci[j] = (Di[j]+Bi[j])*e ;

    }

    // C = C + (D-B)/e
    template<class Field>
    void subscalacc(const Field & F, const size_t m, const size_t n,
                    const double *D, const size_t ldd,
                    const double *B, const size_t ldb,
                    double e,
                    double *C, const size_t ldc)
    {
        const double * Bi = B;
        const double * Di = D;
        double * Ci = C;
        for (;Bi < B+m*ldb ; Ci+=ldc, Bi += ldb, Di += ldd)
            for (size_t j = 0 ; j < n ; ++j)
                Ci[j] += (Di[j]-Bi[j])*e ;

    }

#ifndef TRE
    // #ifndef NDEBUG
    // #define TRE 1
    // #else
#define TRE (size_t)(__FFLASFFPACK_WINOTHRESHOLD)
    // #define TRE (size_t)(__FFLASFFPACK_WINOTHRESHOLD*0.9)
    // #endif
#endif
    template<class Field>
    double * gemm_fflas(const Field & F,
                        const size_t m, const size_t n, const size_t k,
                        const double *A, size_t lda,
                        const double *B, size_t ldb,
                        double * C, size_t ldc,
                        int rec =  0)
    {
        Givaro::DoubleDomain R;
        FFLAS::fgemm(R,
                     FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
                     m,n,k,
                     1,
                     A,lda, B,ldb,
                     0,
                     C, ldc);

        // cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans,
        // m,n,k,1,A,lda,B,ldb,0,C,ldc);

        return C;
    }
} // FFLAS

namespace FFLAS { namespace Protected { namespace Rec {

    // Field must be Givaro::Modular<double>
    template<class Field>
    double *
    gemm_bini_322_0(const Field & F
                    , const size_t m
                    , const size_t n
                    , const size_t k
                    , const double *A , const size_t lda
                    , const double *B , const size_t ldb
                    , double *C , const size_t ldc
                    , int rec
                    , const double & epsilon
                   )
    {
        Givaro::ZRing<double>   NoField;
        // const double p = (double)F.characteristic();
        size_t M = (n>m)?std::min(k,m):std::min(k,n);
        // std::cout << rec << ',' <<  M  << std::endl;
        // Field G(p*p);

        if ( M < TRE  || rec <= 0) {
            return gemm_fflas(F, m,n,k, A,lda, B,ldb, C, ldc);
        }

        assert(k/2*2==k); // k divisible par 2
        assert(n/2*2==n); // n divisible par 2
        assert(m/3*3==m); // m divisible par 3


        size_t n2 = n/2;
        size_t k2 = k/2;
        size_t m3 = m/3;

        // std::cout << "€ = " << epsilon << std::endl;

        // sub matrices in A
        const double * A11 = A;
        const double * A12 = A   +k2;
        const double * A21 = A   +lda*m3;
        const double * A22 = A21 +k2;
        const double * A31 = A21 +lda*m3;
        const double * A32 = A31 +k2;

        // sub matrices in C
        double * C11 = C;
        double * C12 = C   +n2;
        double * C21 = C   +ldc*m3;
        double * C22 = C21 +n2;
        double * C31 = C21 +ldc*m3;
        double * C32 = C31 +n2;

        // sub matrices in B
        const double * B11 = B;
        const double * B12 = B   +n2;
        const double * B21 = B   +ldb*k2;
        const double * B22 = B21 +n2;

        FFLAS::fzero(NoField,m,n,C,ldc);

        /*
         * Algo :
         * S1  := A11  +A22;
         * S4  := e*A12+A22;
         * S5  := A11  +e*A12;
         * S6  := A21  +A32;
         * S9  := A21  +e*A31;
         * S10 := e*A31+A32;
         *
         * T1  := e*B11 +B22;
         * T2  := B21   +B22;
         * T4  := -e*B11+B21;
         * T5  := e*B12 +B22;
         * T6  := B11   +e*B22;
         * T7  := B11   +B12;
         * T9  := B12   -e*B22;
         * T10 := B11   +e*B21;
         *
         * P1 := S1 *T1;
         * P2 := A22*T2;
         * P3 := A11*B22;
         * P4 := S4 *T4;
         * P5 := S5 *T5;
         * P6 := S6 *T6;
         * P7 := A21*T7;
         * P8 := A32*B11;
         * P9 := S9 *T9;
         * P10:= S10*T10;
         *
         * C11 := (P1-P2-P3+P4)/e;
         * C12 := (P3-P5)/(-e) ;
         * C21 := P4+P6-P10 ;
         * C22 := P1-P5+P9;
         * C31 := (-P8+P10)/e;
         * C32 := (P6-P7-P8+P9)/e;
         *
         */

        double * S1 = FFLAS::fflas_new<double>(m3*k2) ;
        // double * C11t = FFLAS::fflas_new<double>(n2*m3) ;
        // S1  := A11  +A22;
        FFLAS::fadd(NoField,m3,k2,A11,lda,A22,lda,S1,k2);
        // T1  := e*B11 +B22;
        double * T1 = FFLAS::fflas_new<double>(n2*k2) ; // ou aire
        add(k2,n2,epsilon,B11,ldb,B22,ldb,T1,n2);
        // P1 := S1 *T1; (dans C22)
        gemm_bini_322_0(F,m3,n2,k2,S1,k2,T1,n2,C22,ldc,rec-1,epsilon);
        // S4  := e*A12+A22;
        double * eA12 = FFLAS::fflas_new<double >(m3*k2);
        FFLAS::fscal(NoField,m3,k2,epsilon,A12,lda,eA12,k2) ;
        FFLAS::fadd(NoField,m3,k2,eA12,k2,A22,lda,S1,k2);
        // T4  := -e*B11+B21;
        add(k2,n2,-epsilon,B11,ldb,B21,ldb,T1,n2);
        // P4 := S4 *T4; (dans C21)
        gemm_bini_322_0(F,m3,n2,k2,S1,k2,T1,n2,C21,ldc,rec-1,epsilon);
        // C11 = P1+P4
        FFLAS::fadd(NoField,m3,n2,C21,ldc,C22,ldc,C11,ldc);
        // T2  := B21  +B22;
        FFLAS::fadd(NoField,k2,n2,B21,ldb,B22,ldb,T1,n2);
        // P2 := A22*T2;
        double * P1 = FFLAS::fflas_new<double>(n2*m3) ; // ou aire
        gemm_bini_322_0(F,m3,n2,k2,A22,lda,T1,n2,P1,n2,rec-1,epsilon);
        // P3 := A11*B22; (dans C12)
        gemm_bini_322_0(F,m3,n2,k2,A11,lda,B22,ldb,C12,ldc,rec-1,epsilon);
        // C11 -= (P2+P3)
        subadd(m3,n2,P1,n2,C12,ldc,C11,ldc);
        // S5  := A11  +e*A12;
        FFLAS::fadd(NoField,m3,k2,eA12,k2,A11,lda,S1,k2);
        // T5  := e*B12 +B22;
        add(k2,n2,epsilon,B12,ldb,B22,ldb,T1,n2);
        // P5 := S5 *T5;
        double * P2 = FFLAS::fflas_new<double>(n2*m3) ; // ou aire
        gemm_bini_322_0(F,m3,n2,k2,S1,k2,T1,n2,P2,n2,rec-1,epsilon);
        // C12 -= P5
        subscalinf(NoField,m3,n2,P2,n2,-(double)1/epsilon,C12,ldc);
        // S6  := A21  +A32;
        FFLAS::fadd(NoField,m3,k2,A21,lda,A32,lda,S1,k2);
        // T6  := B11   +e*B22;
        add(k2,n2,epsilon,B22,ldb,B11,ldb,T1,n2);
        // P6 := S6 *T6; dans C32
        gemm_bini_322_0(F,m3,n2,k2,S1,k2,T1,n2,C32,ldc,rec-1,epsilon);
        // C21+= P6
        FFLAS::faddin(NoField,m3,n2,C32,ldc,C21,ldc);
        // T7  := B11  +B12;
        FFLAS::fadd(NoField,k2,n2,B11,ldb,B12,ldb,T1,n2);
        // P7 := A21*T7; !signe
        gemm_bini_322_0(F,m3,n2,k2,A21,lda,T1,n2,P1,n2,rec-1,epsilon);
        // P8 := A32*B11; dans C31 !signe
        gemm_bini_322_0(F,m3,n2,k2,A32,lda,B11,ldb,C31,ldc,rec-1,epsilon);
        // C32 -= P8+P7
        subadd(m3,n2,P1,n2,C31,ldc,C32,ldc);
        // S9  := A21  +e*A31;
        double * eA31 = eA12 ;
        FFLAS::fscal(NoField,m3,k2,epsilon,A31,lda,eA31,k2);
        FFLAS::fadd(NoField,m3,k2,eA31,k2,A21,lda,S1,k2);
        // T9  := B12   -e*B22;
        add(k2,n2,-epsilon,B22,ldb,B12,ldb,T1,n2);
        // P9 := S9 *T9;
        gemm_bini_322_0(F,m3,n2,k2,S1,k2,T1,n2,P1,n2,rec-1,epsilon);
        // C32= (C32+P9)/p
        addscalinf(NoField,m3,n2,P1,n2,(double)1/epsilon,C32,ldc);
        // C22+= P9-P5
        addsub(m3,n2,P1,n2,P2,n2,C22,ldc);
        FFLAS::fflas_delete( P2);
        // S10 := e*A31+A32;
        FFLAS::fadd(NoField,m3,k2,eA31,k2,A32,lda,S1,k2);
        FFLAS::fflas_delete( eA12 );
        // T10 := B11   +e*B21;
        add(k2,n2,epsilon,B21,ldb,B11,ldb,T1,n2);
        // P10:= S10*T10;
        gemm_bini_322_0(F,m3,n2,k2,S1,k2,T1,n2,P1,n2,rec-1,epsilon);
        FFLAS::fflas_delete( S1);
        FFLAS::fflas_delete( T1);
        // C21-= P10
        FFLAS::fsubin(NoField,m3,n2,P1,n2,C21,ldc);
        // C31= (C31-P10)/(-epsilon)
        subscalinf(NoField,m3,n2,P1,n2,-(double)1/epsilon,C31,ldc);
        FFLAS::fflas_delete( P1);
        // C11 := (P1+P-P3+P4)/e;
        FFLAS::fscalin(NoField,m3,n2,(double)1/epsilon,C11,ldc);

        return C;

    }

    // Field must be Givaro::Modular<double>
    template<class Field>
    double *
    gemm_bini_322_mem(const Field & F
                      , const size_t m
                      , const size_t n
                      , const size_t k
                      , const double *A , const size_t lda
                      , const double *B , const size_t ldb
                      , double *C , const size_t ldc
                      , int rec
                      , const double & epsilon
                     )
    {
        Givaro::ZRing<double>   NoField;
        // const double p = (double)F.characteristic();
        size_t M = (n>m)?std::min(k,m):std::min(k,n);
        // std::cout << rec << ',' <<  M  << std::endl;
        // Field G(p*p);

        if ( M < TRE  || rec <= 0) {
            // std::cout << "ffw" << std::endl;
            return gemm_fflas(F, m,n,k, A,lda, B,ldb, C, ldc);
            // return gemm_fflas(NoField, m,n,k, A,lda, B,ldb, C, ldc);
        }

        assert(k/2*2==k); // k divisible par 2
        assert(n/2*2==n); // n divisible par 2
        assert(m/3*3==m); // m divisible par 3

        // std::cout << "tested" << std::endl;

        size_t n2 = n/2;
        size_t k2 = k/2;
        size_t m3 = m/3;

        // std::cout << "€ = " << epsilon << std::endl;

        // sub matrices in A
        const double * A11 = A;
        const double * A12 = A   +k2;
        const double * A21 = A   +lda*m3;
        const double * A22 = A21 +k2;
        const double * A31 = A21 +lda*m3;
        const double * A32 = A31 +k2;

        // sub matrices in C
        double * C11 = C;
        double * C12 = C   +n2;
        double * C21 = C   +ldc*m3;
        double * C22 = C21 +n2;
        double * C31 = C21 +ldc*m3;
        double * C32 = C31 +n2;

        // sub matrices in B
        const double * B11 = B;
        const double * B12 = B   +n2;
        const double * B21 = B   +ldb*k2;
        const double * B22 = B21 +n2;

        FFLAS::fzero(F,m,n,C,ldc);

        /*
         * Algo :
         * S1  := A11  +A22;
         * S4  := e*A12+A22;
         * S5  := A11  +e*A12;
         * S6  := A21  +A32;
         * S9  := A21  +e*A31;
         * S3 := e*A31+A32;
         *
         * T1  := e*B11 +B22;
         * T2  := B21   +B22;
         * T4  := -e*B11+B21;
         * T5  := e*B12 +B22;
         * T6  := B11   +e*B22;
         * T7  := B11   +B12;
         * T9  := B12   -e*B22;
         * T3 := B11   +e*B21;
         *
         * P1 := S1 *T1;
         * P2 := A22*T2;
         * P10 := A11*B22;
         * P4 := S4 *T4;
         * P5 := S5 *T5;
         * P6 := S6 *T6;
         * P7 := A21*T7;
         * P8 := A32*B11;
         * P9 := S9 *T9;
         * P3:= S3*T3;
         *
         * C11 := (P1-P2-P10+P4)/e;
         * C12 := (P10-P5)/(-e) ;
         * C21 := P4+P6-P3 ;
         * C22 := P1-P5+P9;
         * C31 := (-P8+P3)/e;
         * C32 := (P6-P7-P8+P9)/e;
         *
         */


        // P10
        gemm_bini_322_mem(F,m3,n2,k2,A11,lda,B22,ldb,C11,ldc,rec-1,epsilon);
        // S5
        double * X = FFLAS::fflas_new<double>(m3*k2);
        add(m3,k2,epsilon,A12,lda,A11,lda,X,k2);
        // T5
        // double * Y = FFLAS::fflas_new<double>(std::max(k2,m3)*n2);
        double * Y = FFLAS::fflas_new<double>(k2*n2);
        add(k2,n2,epsilon,B12,ldb,B22,ldb,Y,n2);
        // P5
        gemm_bini_322_mem(F,m3,n2,k2,X,k2,Y,n2,C22,ldc,rec-1,epsilon);
        // C12
        subscal(NoField,m3,n2,C22,ldc,C11,ldc,(double)1/epsilon,C12,ldc);
        // T2
        FFLAS::fadd(NoField,k2,n2,B21,ldb,B22,ldb,Y,n2);
        // P2
        gemm_bini_322_mem(F,m3,n2,k2,A22,lda,Y,n2,C31,ldc,rec-1,epsilon);
        // C11
        FFLAS::faddin(NoField,m3,n2,C31,ldc,C11,ldc);
        // S1
        FFLAS::fadd(NoField,m3,k2,A11,lda,A22,lda,X,k2);
        // T1
        add(k2,n2,epsilon,B11,ldb,B22,ldb,Y,n2);
        // P1
        gemm_bini_322_mem(F,m3,n2,k2,X,k2,Y,n2,C21,ldc,rec-1,epsilon);
        // C22
        FFLAS::fsub(NoField,m3,n2,C21,ldc,C22,ldc,C22,ldc);
        // C11
        FFLAS::fsub(NoField,m3,n2,C21,ldc,C11,ldc,C11,ldc);
        // S4
        add(m3,k2,epsilon,A12,lda,A22,lda,X,k2);
        // T4
        add(k2,n2,-epsilon,B11,ldb,B21,ldb,Y,n2);
        // P4
        gemm_bini_322_mem(F,m3,n2,k2,X,k2,Y,n2,C21,ldc,rec-1,epsilon);
        // C11
        addscalinf(NoField,m3,n2,C21,ldc,(double)1/epsilon,C11,ldc);
        // S9
        add(m3,k2,epsilon,A31,lda,A21,lda,X,k2);
        // T9
        add(k2,n2,-epsilon,B22,ldb,B12,ldb,Y,n2);
        // P9
        gemm_bini_322_mem(F,m3,n2,k2,X,k2,Y,n2,C32,ldc,rec-1,epsilon);
        //  C22
        FFLAS::faddin(NoField,m3,n2,C32,ldc,C22,ldc);
        // S6
        FFLAS::fadd(NoField,m3,k2,A21,lda,A32,lda,X,k2);
        // T6
        add(k2,n2,epsilon,B22,ldb,B11,ldb,Y,n2);
        // P6
        gemm_bini_322_mem(F,m3,n2,k2,X,k2,Y,n2,C31,ldc,rec-1,epsilon);
        // C21
        FFLAS::faddin(NoField,m3,n2,C31,ldc,C21,ldc);
        // C32
        FFLAS::faddin(NoField,m3,n2,C31,ldc,C32,ldc);
        // T7
        FFLAS::fadd(NoField,k2,n2,B11,ldb,B12,ldb,Y,n2);
        // P7
        gemm_bini_322_mem(F,m3,n2,k2,A21,lda,Y,n2,C31,ldc,rec-1,epsilon);
        // if (epsilon > 1 && rec == 2) { FFLAS::finit(G,m3,n2,C31,ldc);}
        // C32
        FFLAS::fsubin(NoField,m3,n2,C31,ldc,C32,ldc);
        // S3
        add(m3,k2,epsilon,A31,lda,A32,lda,X,k2);
        // T3
        add(k2,n2,epsilon,B21,ldb,B11,ldb,Y,n2);
        // P3
        gemm_bini_322_mem(F,m3,n2,k2,X,k2,Y,n2,C31,ldc,rec-1,epsilon);
        FFLAS::fflas_delete( X);
        FFLAS::fflas_delete( Y );
        // C21
        FFLAS::fsubin(NoField,m3,n2,C31,ldc,C21,ldc);
        // P8
        Y = FFLAS::fflas_new<double>(m3*n2);
        gemm_bini_322_mem(F,m3,n2,k2,A32,lda,B11,ldb,Y,n2,rec-1,epsilon);
        // C31
        subscalinf(NoField,m3,n2,Y,n2,(double)1/epsilon,C31,ldc);
        // FFLAS::fsubin(NoField,m3,n2,Y,n2,C31,ldc);
        // C32
        subscalinf(NoField,m3,n2,Y,n2,(double)1/epsilon,C32,ldc);
        // FFLAS::fsubin(NoField,m3,n2,Y,n2,C32,ldc);
        // FFLAS::fscalin(NoField,m3,n,(double)1/epsilon,C31,ldc);
        FFLAS::fflas_delete( Y );


        return C;

    }

    // Field must be Givaro::Modular<double>
    template<class Field>
    double *
    gemm_bini_223_mem(const Field & F
                      , const size_t m
                      , const size_t n
                      , const size_t k
                      , const double *A , const size_t lda
                      , const double *B , const size_t ldb
                      , double *C , const size_t ldc
                      , int rec
                      , const double & epsilon
                     )
    {
        Givaro::ZRing<double>   NoField;
        // const double p = (double)F.characteristic();
        size_t M = (n>m)?std::min(k,m):std::min(k,n);
        // std::cout << rec << ',' <<  M  << std::endl;
        // Field G(p*p);

        if ( M < TRE  || rec <= 0) {
            // std::cout << "ffw" << std::endl;
            return gemm_fflas(F, m,n,k, A,lda, B,ldb, C, ldc);
            // return gemm_fflas(NoField, m,n,k, A,lda, B,ldb, C, ldc);
        }

        assert(k/2*2==k); // k divisible par 2
        assert(n/3*3==n); // n divisible par 2
        assert(m/2*2==m); // m divisible par 3

        // std::cout << "tested" << std::endl;

        size_t m2 = m/2;
        size_t k2 = k/2;
        size_t n3 = n/3;

        // std::cout << "€ = " << epsilon << std::endl;

        // sub matrices in A
        const double * A11 = A;
        const double * A12 = A   +k2;
        const double * A21 = A   +lda*m2;
        const double * A22 = A21 +k2;

        // sub matrices in C
        double * C11 = C;
        double * C12 = C   +n3;
        double * C13 = C   +2*n3;
        double * C21 = C   +ldc*m2;
        double * C22 = C21 +n3;
        double * C23 = C21 +2*n3;



        // sub matrices in B
        const double * B11 = B;
        const double * B12 = B   +n3;
        const double * B13 = B   +2*n3;
        const double * B21 = B   +ldb*k2;
        const double * B22 = B21 +n3;
        const double * B23 = B21 +2*n3;


        FFLAS::fzero(F,m,n,C,ldc);

        /*
         * Algo :
         * S1  := B11  +B22;
         * S4  := e*B21+B22;
         * S5  := B11  +e*B21;
         * S6  := B12  +B23;
         * S9  := B12  +e*B13;
         * S3 := e*B13+B23;
         *
         * T1  := e*A11 +A22;
         * T2  := A12   +A22;
         * T4  := -e*A11+A12;
         * T5  := e*A21 +A22;
         * T6  := A11   +e*A22;
         * T7  := A11   +A21;
         * T9  := A21   -e*A22;
         * T3  := A11   +e*A12;
         *
         * P1 := S1 *T1;
         * P2 := T2 * B22;
         * P10 := A22 * B11;
         * P4 := S4 *T4;
         * P5 := S5 *T5;
         * P6 := S6 *T6;
         * P7 := T7*B12;
         * P8 := A11*B23;
         * P9 := S9 *T9;
         * P3 := S3*T3;
         *
         * C11 := (P1-P2-P10+P4)/e;
         * C21 := (P10-P5)/(-e) ;
         * C12 := P4+P6-P3 ;
         * C22 := P1-P5+P9;
         * C13 := (-P8+P3)/e;
         * C23 := (P6-P7-P8+P9)/e;
         *
         */


        // P10
        gemm_bini_223_mem(F,m2,n3,k2,A22,lda,B11,ldb,C11,ldc,rec-1,epsilon);
        // S5
        double * Y = FFLAS::fflas_new<double>(k2*n3);
        add(k2,n3,epsilon,B21,ldb,B11,ldb,Y,n3);
        // T5
        double * X = FFLAS::fflas_new<double>(m2*k2);
        add(m2,k2,epsilon,A21,lda,A22,lda,X,k2);
        // P5
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C22,ldc,rec-1,epsilon);
        // C12
        subscal(NoField,m2,n3,C22,ldc,C11,ldc,(double)1/epsilon,C21,ldc);
        // T2
        FFLAS::fadd(NoField,m2,k2,A12,lda,A22,lda,X,k2);
        // P2
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,B22,ldb,C13,ldc,rec-1,epsilon);
        // C11
        FFLAS::faddin(NoField,m2,n3,C13,ldc,C11,ldc);
        // S1
        FFLAS::fadd(NoField,k2,n3,B11,ldb,B22,ldb,Y,n3);
        // T1
        add(m2,k2,epsilon,A11,lda,A22,lda,X,k2);
        // P1
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C12,ldc,rec-1,epsilon);
        // C22
        FFLAS::fsub(NoField,m2,n3,C12,ldc,C22,ldc,C22,ldc);
        // C11
        FFLAS::fsub(NoField,m2,n3,C12,ldc,C11,ldc,C11,ldc);
        // S4
        add(k2,n3,epsilon,B21,ldb,B22,ldb,Y,n3);
        // T4
        add(m2,k2,-epsilon,A11,lda,A12,lda,X,k2);
        // P4
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C12,ldc,rec-1,epsilon);
        // C11
        addscalinf(NoField,m2,n3,C12,ldc,(double)1/epsilon,C11,ldc);
        // S9
        add(k2,n3,epsilon,B13,ldb,B12,ldb,Y,n3);
        // T9
        add(m2,k2,-epsilon,A22,lda,A21,lda,X,k2);
        // P9
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C23,ldc,rec-1,epsilon);
        //  C22
        FFLAS::faddin(NoField,m2,n3,C23,ldc,C22,ldc);
        // S6
        FFLAS::fadd(NoField,k2,n3,B12,ldb,B23,ldb,Y,n3);
        // T6
        add(m2,k2,epsilon,A22,lda,A11,lda,X,k2);
        // P6
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C13,ldc,rec-1,epsilon);
        // C21
        FFLAS::faddin(NoField,m2,n3,C13,ldc,C12,ldc);
        // C32
        FFLAS::faddin(NoField,m2,n3,C13,ldc,C23,ldc);
        // T7
        FFLAS::fadd(NoField,m2,k2,A11,lda,A21,lda,X,k2);
        // P7
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,B12,ldb,C13,ldc,rec-1,epsilon);
        // if (epsilon > 1 && rec == 2) { FFLAS::finit(G,m2,n3,C31,ldc);}
        // C32
        FFLAS::fsubin(NoField,m2,n3,C13,ldc,C23,ldc);
        // S3
        add(k2,n3,epsilon,B13,ldb,B23,ldb,Y,n3);
        // T3
        add(m2,k2,epsilon,A12,lda,A11,lda,X,k2);
        // P3
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C13,ldc,rec-1,epsilon);
        FFLAS::fflas_delete( Y );
        FFLAS::fflas_delete( X );
        // C21
        FFLAS::fsubin(NoField,m2,n3,C13,ldc,C12,ldc);
        // P8
        Y = FFLAS::fflas_new<double>(m2*n3);
        gemm_bini_223_mem(F,m2,n3,k2,A11,lda,B23,ldb,Y,n3,rec-1,epsilon);
        // C31
        subscalinf(NoField,m2,n3,Y,n3,(double)1/epsilon,C13,ldc);
        // C32
        subscalinf(NoField,m2,n3,Y,n3,(double)1/epsilon,C23,ldc);
        FFLAS::fflas_delete( Y );


        return C;

    }

    // Field must be Givaro::Modular<double>
    template<class Field>
    double *
    gemm_bini_322_2(const Field & F
                    , const size_t m
                    , const size_t n
                    , const size_t k
                    , const double *A , const size_t lda
                    , const double *B , const size_t ldb
                    , double *C , const size_t ldc
                    , int rec
                    , const double & epsilon
                   )
    {
        Givaro::ZRing<double>   NoField;
        // const double p = (double)F.characteristic();
        size_t M = (n>m)?std::min(k,m):std::min(k,n);
        // std::cout << rec << ',' <<  M  << std::endl;
        // Field G(p*p);

        if ( M < TRE  || rec <= 0) {
            // std::cout << "ffw" << std::endl;
            return gemm_fflas(F, m,n,k, A,lda, B,ldb, C, ldc);
            // return gemm_fflas(NoField, m,n,k, A,lda, B,ldb, C, ldc);
        }

        assert(k/2*2==k); // k divisible par 2
        assert(n/2*2==n); // n divisible par 2
        assert(m/3*3==m); // m divisible par 3

        // std::cout << "tested" << std::endl;

        size_t n2 = n/2;
        size_t k2 = k/2;
        size_t m3 = m/3;

        // std::cout << "€ = " << epsilon << std::endl;

        // sub matrices in A
        const double * A11 = A;
        const double * A12 = A   +k2;
        const double * A21 = A   +lda*m3;
        const double * A22 = A21 +k2;
        const double * A31 = A21 +lda*m3;
        const double * A32 = A31 +k2;

        // sub matrices in C
        double * C11 = C;
        double * C12 = C   +n2;
        double * C21 = C   +ldc*m3;
        double * C22 = C21 +n2;
        double * C31 = C21 +ldc*m3;
        double * C32 = C31 +n2;

        // sub matrices in B
        const double * B11 = B;
        const double * B12 = B   +n2;
        const double * B21 = B   +ldb*k2;
        const double * B22 = B21 +n2;

        FFLAS::fzero(F,m,n,C,ldc);

        /*
         * Algo :
         * S1 := A11  +A22;
         * S4 := e*A12+A22;
         * S5 := A11  +e*A12;
         * S3 := e*A31+A32;
         * S6 := A21  +A32;
         * S9 := A21  +e*A31;
         *
         * T1 := e*B11 +B22;
         * T2 := B21   +B22;
         * T3 := B11   +e*B21;
         * T4 := -e*B11+B21;
         * T5 := e*B12 +B22;
         * T6 := B11   +e*B22;
         * T7 := B11   +B12;
         * T9 := B12   -e*B22;
         *
         * P1 := S1 *T1;
         * P2 := A22*T2;
         * P10 := A11*B22;
         * P4 := S4 *T4;
         * P5 := S5 *T5;
         * P6 := S6 *T6;
         * P7 := A21*T7;
         * P8 := A32*B11;
         * P9 := S9 *T9;
         * P3:= S3*T3;
         *
         * C11 := (P1-P2-P10+P4)/e;
         * C12 := (P10-P5)/(-e) ;
         * C21 := P4+P6-P3 ;
         * C22 := P1-P5+P9;
         * C31 := (-P8+P3)/e;
         * C32 := (P6-P7-P8+P9)/e;
         *
         */

        double * U = FFLAS::fflas_new<double>(m3*n2);
        double * V = FFLAS::fflas_new<double>(m3*n2);
        double * X = FFLAS::fflas_new<double>(m3*std::max(k2,n2));
        double * Y = FFLAS::fflas_new<double>(std::max(k2,m3)*n2);

        // S4
        add(m3,k2,epsilon,A12,lda,A22,lda,X,k2);
        // T4
        add(k2,n2,-epsilon,B11,ldb,B21,ldb,Y,n2);
        // P4
        gemm_bini_322_2(F,m3,n2,k2,X,k2,Y,n2,U,n2,rec-1,epsilon);
        // S9
        add(m3,k2,epsilon,A31,lda,A21,lda,X,k2);
        // T9
        add(k2,n2,-epsilon,B22,ldb,B12,ldb,Y,n2);
        // P9
        gemm_bini_322_2(F,m3,n2,k2,X,k2,Y,n2,V,n2,rec-1,epsilon);
        // S5
        add(m3,k2,epsilon,A12,lda,A11,lda,X,k2);
        // T5
        add(k2,n2,epsilon,B12,ldb,B22,ldb,Y,n2);
        // P5
        gemm_bini_322_2(F,m3,n2,k2,X,k2,Y,n2,C12,ldc,rec-1,epsilon);
        // S3
        add(m3,k2,epsilon,A31,lda,A32,lda,X,k2);
        // T3
        add(k2,n2,epsilon,B21,ldb,B11,ldb,Y,n2);
        // P3
        gemm_bini_322_2(F,m3,n2,k2,X,k2,Y,n2,C31,ldc,rec-1,epsilon);
        // C22 = P9-P5
        FFLAS::fsub(NoField,m3,n2,V,n2,C12,ldc,C22,ldc);
        // C21 = P4-P3
        FFLAS::fsub(NoField,m3,n2,U,n2,C31,ldc,C21,ldc);
        // T2
        FFLAS::fadd(NoField,k2,n2,B21,ldb,B22,ldb,Y,n2);
        // P2
        gemm_bini_322_2(F,m3,n2,k2,A22,lda,Y,n2,X,n2,rec-1,epsilon);
        // XXX approximate
        // C11 = (P4 - P2) / e
        subscal(NoField,m3,n2,U,n2,X,n2,1./epsilon,C11,ldc);
        // T7
        FFLAS::fadd(NoField,k2,n2,B11,ldb,B12,ldb,Y,n2);
        // P7
        gemm_bini_322_2(F,m3,n2,k2,A21,lda,Y,n2,X,n2,rec-1,epsilon);
        // XXX approximate
        // C32 = (P9-P7) / e
        subscal(NoField,m3,n2,V,n2,X,n2,1./epsilon,C32,ldc);
        // S1
        FFLAS::fadd(NoField,m3,k2,A11,lda,A22,lda,X,k2);
        // T1
        add(k2,n2,epsilon,B11,ldb,B22,ldb,Y,n2);
        // P1
        gemm_bini_322_2(F,m3,n2,k2,X,k2,Y,n2,U,n2,rec-1,epsilon);
        // C22 += P1
        FFLAS::faddin(NoField,m3,n2,U,n2,C22,ldc);
        // P10
        gemm_bini_322_2(F,m3,n2,k2,A11,lda,B22,ldb,V,n2,rec-1,epsilon);
        // C12 = (P5-P10)/e
        subscalinf(NoField,m3,n2,V,n2,1./epsilon,C12,ldc);
        // XXX approximate
        // C11 = C11 + (P1-P10)/e
        subscalacc(NoField,m3,n2,U,n2,V,n2,1./epsilon,C11,ldc);
        // S6
        FFLAS::fadd(NoField,m3,k2,A21,lda,A32,lda,X,k2);
        // T6
        add(k2,n2,epsilon,B22,ldb,B11,ldb,Y,n2);
        // P6
        gemm_bini_322_2(F,m3,n2,k2,X,k2,Y,n2,U,n2,rec-1,epsilon);
        // C21 += P6
        FFLAS::faddin(NoField,m3,n2,U,n2,C21,ldc);
        // P8
        gemm_bini_322_2(F,m3,n2,k2,A32,lda,B11,ldb,V,n2,rec-1,epsilon);
        // C31 = (P3-P8)/2
        subscalinf(NoField,m3,n2,V,n2,1./epsilon,C31,ldc);
        // XXX approximate
        // C32 = C32 + (P6-P8)/e
        subscalacc(NoField,m3,n2,U,n2,V,n2,1./epsilon,C32,ldc);


        FFLAS::fflas_delete( X);
        FFLAS::fflas_delete( Y );
        FFLAS::fflas_delete( U);
        FFLAS::fflas_delete( V);


        return C;

    }


    // Field must be Givaro::Modular<double>
    template<class Field>
    double *
    gemm_bini_232_2(const Field & F
                    , const size_t m
                    , const size_t n
                    , const size_t k
                    , const double *A , const size_t lda
                    , const double *B , const size_t ldb
                    , double *C , const size_t ldc
                    , int rec
                    , const double & epsilon
                   )
    {
        Givaro::ZRing<double>   NoField;
        // const double p = (double)F.characteristic();
        size_t M = (n>m)?std::min(k,m):std::min(k,n);
        // Field G(p*p);

        if ( M < TRE  || rec <= 0) {
            // std::cout << "ffw" << std::endl;
            return gemm_fflas(F, m,n,k, A,lda, B,ldb, C, ldc);
            // return gemm_fflas(NoField, m,n,k, A,lda, B,ldb, C, ldc);
        }

        assert(k/3*3==k); // k divisible par 3
        assert(n/2*2==n); // n divisible par 2
        assert(m/2*2==m); // m divisible par 2

        // std::cout << "tested" << std::endl;

        size_t n2 = n/2;
        size_t k3 = k/3;
        size_t m2 = m/2;

        // std::cout << "€ = " << epsilon << std::endl;

        // sub matrices in B
        const double * B11 = B;
        const double * B12 = B   +n2;
        const double * B21 = B   +ldb*k3;
        const double * B22 = B21 +n2;
        const double * B31 = B21 +ldb*k3;
        const double * B32 = B31 +n2;

        // sub matrices in C
        double * C11 = C;
        double * C12 = C   +n2;
        double * C21 = C   +ldc*m2;
        double * C22 = C21 +n2;

        // sub matrices in A

        const double * A11 = A;
        const double * A12 = A   +k3;
        const double * A13 = A   +2*k3;
        const double * A21 = A   +lda*m2;
        const double * A22 = A21 +k3;
        const double * A23 = A21 +2*k3;


        FFLAS::fzero(F,m,n,C,ldc);

        /*
         * Algo :
         *
         * S1  := A11  +A22*e;
         * S3  := -(A11+A21);
         * S4  := A11+A12*e;
         * S5  := A21 - A22*e;
         * S6  := A12*e  + A23;
         * S8 := -(A13+A23):
         * S9  := A22*e  + A23;
         * S10 := -A12*e+A13;
         *
         * T1  := B11 +B22;
         * T4  := e*B12+B22;
         * T5  := B11 +e*B12;
         * T6  := B21   +B32;
         * T9  := B21   + e*B31;
         * T10 := e*B31   +B32;
         *
         * P1 := Bini232(S1,T1 ,e);
         * P2 := Bini232(A11,B22 ,e);
         * P3 := Bini232(S3,B11,e);
         * P4 := Bini232(S4,T4 ,e);
         * P5 := Bini232(S5,T5 ,e);
         * P6 := Bini232(S6,T6 ,e);
         * P7 := Bini232(A23,B21 ,e);
         * P8 := Bini232(S8,B32,e);
         * P9 := Bini232(S9,T9 ,e);
         * P10:= Bini232(S10,T10,e);
         *
         *
         * C11 := evalm(P1-P4+(P6-P7+P8+P10)/e);
         * C12 := evalm((-P2+P4)/e+P10) ;
         * C21 := evalm(P5+(-P7+P9)/e) ;
         * C22 := evalm((P1-P2+P3+P5)/e+P6-P9);
         *
         */

        double * U = FFLAS::fflas_new<double>(m2*n2);
        double * V = FFLAS::fflas_new<double>(m2*n2);
        double * X = FFLAS::fflas_new<double>(m2*k3);
        double * Y = FFLAS::fflas_new<double>(k3*n2);

        // S1
        add(m2,k3,epsilon,A22,lda,A11,lda,X,k3);
        // T1
        FFLAS::fadd(NoField,k3,n2,B11,ldb,B22,ldb,Y,n2);
        // P1 (in U)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
        // S3
        negadd(m2,k3,A11,lda,A21,lda,X,k3);
        // P3 (in V)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,B11,ldb,V,n2,rec-1,epsilon);
        // C22 = (P1+P3)/e
        // FFLAS::fadd(NoField,m2,n2,U,n2,V,n2,C22,ldc); // XXX acc
        addscal(NoField,m2,n2,U,n2,V,n2,(double)1/epsilon,C22,ldc);
        // S6
        add(m2,k3,epsilon,A12,lda,A23,lda,X,k3);
        // T6
        FFLAS::fadd(NoField,k3,n2,B21,ldb,B32,ldb,Y,n2);
        // P6 (in V)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,V,n2,rec-1,epsilon);
        // C22 += P6
        FFLAS::faddin(NoField,m2,n2,V,n2,C22,ldc);
        // S8
        negadd(m2,k3,A13,lda,A23,lda,X,k3);
        // P8 (in C11)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,B32,ldb,C11,ldc,rec-1,epsilon);
        // C11 = (P8+P6)/e
        addscalinf(NoField,m2,n2,V,n2,(double)1/epsilon,C11,ldc);
        // C11 += P1
        FFLAS::faddin(NoField,m2,n2,U,n2,C11,ldc);
        // S4
        add(m2,k3,epsilon,A12,lda,A11,lda,X,k3);
        // T4
        add(k3,n2,epsilon,B12,ldb,B22,ldb,Y,n2);
        // P4 (in U)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
        // C11 -= P4
        FFLAS::fsubin(NoField,m2,n2,U,n2,C11,ldc);
        // P2 (in C12)
        gemm_bini_232_2(F,m2,n2,k3,A11,lda,B22,ldb,C12,ldc,rec-1,epsilon);
        // S5
        add(m2,k3,-epsilon,A22,lda,A21,lda,X,k3);
        // T5
        add(k3,n2,epsilon,B12,ldb,B11,ldb,Y,n2);
        // P5 (in V)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,V,n2,rec-1,epsilon);
        // C22 += (P5-P2)/e
        subscalacc(NoField,m2,n2,V,n2,C12,ldc,(double)1/epsilon,C22,ldc);
        // C12 = (P4-P2)/e
        subscalinf(NoField,m2,n2,U,n2,-(double)1/epsilon,C12,ldc);
        // S9
        add(m2,k3,epsilon,A22,lda,A23,lda,X,k3);
        // T9
        add(k3,n2,epsilon,B31,ldb,B21,ldb,Y,n2);
        // P9 (in U)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
        // C22 -= P9
        FFLAS::fsubin(NoField,m2,n2,U,n2,C22,ldc);
        // P7 (in C21)
        gemm_bini_232_2(F,m2,n2,k3,A23,lda,B21,ldb,C21,ldc,rec-1,epsilon);
        // C11 = C11  - P7/e
        add(m2,n2,-(double)1/epsilon,C21,ldc,C11,ldc,C11,ldc);
        // C21 =  (P9-P7)/e
        subscalinf(NoField,m2,n2,U,n2,-(double)1/epsilon,C21,ldc);
        // C21 += P5
        FFLAS::faddin(NoField,m2,n2,V,n2,C21,ldc);
        // S10
        add(m2,k3,-epsilon,A12,lda,A13,lda,X,k3);
        // T10
        add(k3,n2,epsilon,B31,ldb,B32,ldb,Y,n2);
        // P10 (in U)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
        // C12 += P10
        FFLAS::faddin(NoField,m2,n2,U,n2,C12,ldc);
        // C11 += P10/e
        add(m2,n2,(double)1/epsilon,U,n2,C11,ldc,C11,ldc);


        FFLAS::fflas_delete( X );
        FFLAS::fflas_delete( Y );
        FFLAS::fflas_delete( U );
        FFLAS::fflas_delete( V );


        return C;

    }

    template<class Field>
    double *
    gemm_bini_232_3_acc(const Field & F
                        , const size_t m
                        , const size_t n
                        , const size_t k
                        , const double *A , const size_t lda
                        , const double *B , const size_t ldb
                        , double *C , const size_t ldc
                        , int rec
                        , const double & epsilon
                       )
    {
        if (rec != 0)
            exit(-1);
        Givaro::DoubleDomain R;
        FFLAS::fgemm(R,
                     FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
                     m,n,k,
                     1,
                     A,lda, B,ldb,
                     1,
                     C, ldc);


    }

    template<class Field>
    double *
    gemm_bini_232_3(const Field & F
                    , const size_t m
                    , const size_t n
                    , const size_t k
                    , const double *A , const size_t lda
                    , const double *B , const size_t ldb
                    , double *C , const size_t ldc
                    , int rec
                    , const double & epsilon
                   )
    {
        Givaro::ZRing<double>   NoField;
        // const double p = (double)F.characteristic();
        size_t M = (n>m)?std::min(k,m):std::min(k,n);
        // Field G(p*p);

        if ( M < TRE  || rec <= 0) {
            // std::cout << "ffw" << std::endl;
            return gemm_fflas(F, m,n,k, A,lda, B,ldb, C, ldc);
            // return gemm_fflas(NoField, m,n,k, A,lda, B,ldb, C, ldc);
        }

        assert(k/3*3==k); // k divisible par 3
        assert(n/2*2==n); // n divisible par 2
        assert(m/2*2==m); // m divisible par 2

        // std::cout << "tested" << std::endl;

        size_t n2 = n/2;
        size_t k3 = k/3;
        size_t m2 = m/2;

        // std::cout << "€ = " << epsilon << std::endl;

        // sub matrices in B
        const double * B11 = B;
        const double * B12 = B   +n2;
        const double * B21 = B   +ldb*k3;
        const double * B22 = B21 +n2;
        const double * B31 = B21 +ldb*k3;
        const double * B32 = B31 +n2;

        // sub matrices in C
        double * C11 = C;
        double * C12 = C   +n2;
        double * C21 = C   +ldc*m2;
        double * C22 = C21 +n2;

        // sub matrices in A

        const double * A11 = A;
        const double * A12 = A   +k3;
        const double * A13 = A   +2*k3;
        const double * A21 = A   +lda*m2;
        const double * A22 = A21 +k3;
        const double * A23 = A21 +2*k3;


        FFLAS::fzero(F,m,n,C,ldc);

        /*
         * Algo :
         *
         * S1  := A11  +A22*e;
         * S3  := -(A11+A21);
         * S4  := A11+A12*e;
         * S5  := A21 - A22*e;
         * S6  := A12*e  + A23;
         * S8 := -(A13+A23):
         * S9  := A22*e  + A23;
         * S10 := -A12*e+A13;
         *
         * T1  := B11 +B22;
         * T4  := e*B12+B22;
         * T5  := B11 +e*B12;
         * T6  := B21   +B32;
         * T9  := B21   + e*B31;
         * T10 := e*B31   +B32;
         *
         * P1 := Bini232(S1,T1 ,e);
         * P2 := Bini232(A11,B22 ,e);
         * P3 := Bini232(S3,B11,e);
         * P4 := Bini232(S4,T4 ,e);
         * P5 := Bini232(S5,T5 ,e);
         * P6 := Bini232(S6,T6 ,e);
         * P7 := Bini232(A23,B21 ,e);
         * P8 := Bini232(S8,B32,e);
         * P9 := Bini232(S9,T9 ,e);
         * P10:= Bini232(S10,T10,e);
         *
         *
         * C11 := evalm(P1-P4+(P6-P7+P8+P10)/e);
         * C12 := evalm((-P2+P4)/e+P10) ;
         * C21 := evalm(P5+(-P7+P9)/e) ;
         * C22 := evalm((P1-P2+P3+P5)/e+P6-P9);
         *
         */

        // could be just one band for the scalings



        double * U = FFLAS::fflas_new<double>(m2*n2);
        double * V = FFLAS::fflas_new<double>(std::max(k3,m2)*n2);
        double * X = FFLAS::fflas_new<double>(m2*k3);
        double * Y = FFLAS::fflas_new<double>(k3*n2);

        // S1
        double * eA22 = FFLAS::fflas_new<double>(std::max(m2,n2)*k3);
        FFLAS::fscal(NoField,m2,k3,epsilon,A22,lda,eA22,k3);
        FFLAS::fadd(NoField,m2,k3,eA22,k3,A11,lda,X,k3);
        // T1
        FFLAS::fadd(NoField,k3,n2,B11,ldb,B22,ldb,Y,n2);
        // P1 (in U)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
        // S3
        negadd(m2,k3,A11,lda,A21,lda,X,k3);
        // P3 (in V)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,B11,ldb,V,n2,rec-1,epsilon);
        // C22 = (P1+P3)/e
        addscal(NoField,m2,n2,U,n2,V,n2,(double)1/epsilon,C22,ldc);
        // S6
        double * eA12 = FFLAS::fflas_new<double>(m2*k3);
        FFLAS::fscal(NoField,m2,k3,epsilon,A12,lda,eA12,k3);
        FFLAS::fadd(NoField,m2,k3,eA12,k3,A23,lda,X,k3);
        // T6
        FFLAS::fadd(NoField,k3,n2,B21,ldb,B32,ldb,Y,n2);
        // P6 (in V)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,V,n2,rec-1,epsilon);
        // C22 += P6
        FFLAS::faddin(NoField,m2,n2,V,n2,C22,ldc);
        // S8
        negadd(m2,k3,A13,lda,A23,lda,X,k3);
        // P8 (in C11)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,B32,ldb,C11,ldc,rec-1,epsilon);
        // C11 = (P8+P6)/e
        addscalinf(NoField,m2,n2,V,n2,(double)1/epsilon,C11,ldc);
        // C11 += P1
        FFLAS::faddin(NoField,m2,n2,U,n2,C11,ldc);
        // S4
        FFLAS::fadd(NoField,m2,k3,eA12,k3,A11,lda,X,k3);
        // T4
        double * eB12 = V ; // FFLAS::fflas_new<double>(n2*k3);
        FFLAS::fscal(NoField,k3,n2,epsilon,B12,ldb,eB12,n2);
        FFLAS::fadd(NoField,k3,n2,eB12,n2,B22,ldb,Y,n2);
        // P4 (in U)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
        // C11 -= P4
        FFLAS::fsubin(NoField,m2,n2,U,n2,C11,ldc);
        // P2 (in C12)
        gemm_bini_232_2(F,m2,n2,k3,A11,lda,B22,ldb,C12,ldc,rec-1,epsilon);
        // S5
        FFLAS::fsub(NoField,m2,k3,A21,lda,eA22,k3,X,k3);
        // T5
        FFLAS::fadd(NoField,k3,n2,eB12,n2,B11,ldb,Y,n2);
        // FFLAS::fflas_delete( eB12);
        // P5 (in V)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,V,n2,rec-1,epsilon);
        // C22 += (P5-P2)/e
        subscalacc(NoField,m2,n2,V,n2,C12,ldc,(double)1/epsilon,C22,ldc);
        // C12 = (P4-P2)/e
        subscalinf(NoField,m2,n2,U,n2,-(double)1/epsilon,C12,ldc);
        // S9
        FFLAS::fadd(NoField,m2,k3,eA22,k3,A23,lda,X,k3);
        double * eB31 = eA22 ;
        FFLAS::fscal(NoField,k3,n2,epsilon,B31,ldb,eB31,n2);
        // T9
        FFLAS::fadd(NoField,k3,n2,eB31,n2,B21,ldb,Y,n2);
        // P9 (in U)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
        // C22 -= P9
        FFLAS::fsubin(NoField,m2,n2,U,n2,C22,ldc);
        // P7 (in C21)
        gemm_bini_232_2(F,m2,n2,k3,A23,lda,B21,ldb,C21,ldc,rec-1,epsilon);
        // C11 = C11  - P7/e
        add(m2,n2,-(double)1/epsilon,C21,ldc,C11,ldc,C11,ldc);
        // C21 =  (P9-P7)/e
        subscalinf(NoField,m2,n2,U,n2,-(double)1/epsilon,C21,ldc);
        // C21 += P5
        FFLAS::faddin(NoField,m2,n2,V,n2,C21,ldc);
        // S10
        FFLAS::fsub(NoField,m2,k3,A13,lda,eA12,k3,X,k3);
        FFLAS::fflas_delete( eA12);
        // T10
        FFLAS::fadd(NoField,k3,n2,eB31,n2,B32,ldb,Y,n2);
        FFLAS::fflas_delete( eA22);
        // P10 (in U)
        gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
        // C12 += P10
        FFLAS::faddin(NoField,m2,n2,U,n2,C12,ldc);
        // C11 += P10/e
        add(m2,n2,(double)1/epsilon,U,n2,C11,ldc,C11,ldc);


        FFLAS::fflas_delete( X );
        FFLAS::fflas_delete( Y );
        FFLAS::fflas_delete( U );
        FFLAS::fflas_delete( V );


        return C;

    }

#if 0
    template<class Field>
    double *
    gemm_bini_322_sqrt(const Field & F
                       , const size_t m
                       , const size_t n
                       , const size_t k
                       , const double *A , const size_t lda
                       , const double *B , const size_t ldb
                       , double *C , const size_t ldc
                       , int rec
                       , const double & epsilon
                      )
    {
        Givaro::ZRing<double>   NoField;
        // const double p = (double)F.characteristic();
        size_t M = (n>m)?std::min(k,m):std::min(k,n);
        // std::cout << rec << ',' <<  M  << std::endl;
        // Field G(p*p);

        if ( M < TRE  || rec <= 0) {
            // std::cout << "ffw" << std::endl;
            return gemm_fflas(F, m,n,k, A,lda, B,ldb, C, ldc);
            // return gemm_fflas(NoField, m,n,k, A,lda, B,ldb, C, ldc);
        }

        assert(k/2*2==k); // k divisible par 2
        assert(n/3*3==n); // n divisible par 2
        assert(m/2*2==m); // m divisible par 3

        // std::cout << "tested" << std::endl;

        size_t m2 = m/2;
        size_t k2 = k/2;
        size_t n3 = n/3;

        // std::cout << "€ = " << epsilon << std::endl;

        // sub matrices in A
        const double * A11 = A;
        const double * A12 = A   +k2;
        const double * A21 = A   +lda*m2;
        const double * A22 = A21 +k2;

        // sub matrices in C
        double * C11 = C;
        double * C12 = C   +n3;
        double * C13 = C   +2*n3;
        double * C21 = C   +ldc*m2;
        double * C22 = C21 +n3;
        double * C23 = C21 +2*n3;



        // sub matrices in B
        const double * B11 = B;
        const double * B12 = B   +n3;
        const double * B13 = B   +2*n3;
        const double * B21 = B   +ldb*k2;
        const double * B22 = B21 +n3;
        const double * B23 = B21 +2*n3;


        FFLAS::fzero(F,m,n,C,ldc);

        /*
         * Algo :
         * S1  := B11  +B22;
         * S4  := e*B21+B22;
         * S5  := B11  +e*B21;
         * S6  := B12  +B23;
         * S9  := B12  +e*B13;
         * S3 := e*B13+B23;
         *
         * T1  := e*A11 +A22;
         * T2  := A12   +A22;
         * T4  := -e*A11+A12;
         * T5  := e*A21 +A22;
         * T6  := A11   +e*A22;
         * T7  := A11   +A21;
         * T9  := A21   -e*A22;
         * T3  := A11   +e*A12;
         *
         * P1 := S1 *T1;
         * P2 := T2 * B22;
         * P10 := A22 * B11;
         * P4 := S4 *T4;
         * P5 := S5 *T5;
         * P6 := S6 *T6;
         * P7 := T7*B12;
         * P8 := A11*B23;
         * P9 := S9 *T9;
         * P3 := S3*T3;
         *
         * C11 := (P1-P2-P10+P4)/e;
         * C21 := (P10-P5)/(-e) ;
         * C12 := P4+P6-P3 ;
         * C22 := P1-P5+P9;
         * C13 := (-P8+P3)/e;
         * C23 := (P6-P7-P8+P9)/e;
         *
         */


        // P10
        gemm_bini_223_mem(F,m2,n3,k2,A22,lda,B11,ldb,C11,ldc,rec-1,epsilon);
        // S5
        double * Y = FFLAS::fflas_new<double>(k2*n3);
        add(k2,n3,epsilon,B21,ldb,B11,ldb,Y,n3);
        // T5
        double * X = FFLAS::fflas_new<double>(m2*k2);
        add(m2,k2,epsilon,A21,lda,A22,lda,X,k2);
        // P5
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C22,ldc,rec-1,epsilon);
        // C12
        subscal(NoField,m2,n3,C22,ldc,C11,ldc,(double)1/epsilon,C21,ldc);
        // T2
        FFLAS::fadd(NoField,m2,k2,A12,lda,A22,lda,X,k2);
        // P2
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,B22,ldb,C13,ldc,rec-1,epsilon);
        // C11
        FFLAS::faddin(NoField,m2,n3,C13,ldc,C11,ldc);
        // S1
        FFLAS::fadd(NoField,k2,n3,B11,ldb,B22,ldb,Y,n3);
        // T1
        add(m2,k2,epsilon,A11,lda,A22,lda,X,k2);
        // P1
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C12,ldc,rec-1,epsilon);
        // C22
        FFLAS::fsub(NoField,m2,n3,C12,ldc,C22,ldc,C22,ldc);
        // C11
        FFLAS::fsub(NoField,m2,n3,C12,ldc,C11,ldc,C11,ldc);
        // S4
        add(k2,n3,epsilon,B21,ldb,B22,ldb,Y,n3);
        // T4
        add(m2,k2,-epsilon,A11,lda,A12,lda,X,k2);
        // P4
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C12,ldc,rec-1,epsilon);
        // C11
        addscalinf(NoField,m2,n3,C12,ldc,(double)1/epsilon,C11,ldc);
        // S9
        add(k2,n3,epsilon,B13,ldb,B12,ldb,Y,n3);
        // T9
        add(m2,k2,-epsilon,A22,lda,A21,lda,X,k2);
        // P9
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C23,ldc,rec-1,epsilon);
        //  C22
        FFLAS::faddin(NoField,m2,n3,C23,ldc,C22,ldc);
        // S6
        FFLAS::fadd(NoField,k2,n3,B12,ldb,B23,ldb,Y,n3);
        // T6
        add(m2,k2,epsilon,A22,lda,A11,lda,X,k2);
        // P6
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C13,ldc,rec-1,epsilon);
        // C21
        FFLAS::faddin(NoField,m2,n3,C13,ldc,C12,ldc);
        // C32
        FFLAS::faddin(NoField,m2,n3,C13,ldc,C23,ldc);
        // T7
        FFLAS::fadd(NoField,m2,k2,A11,lda,A21,lda,X,k2);
        // P7
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,B12,ldb,C13,ldc,rec-1,epsilon);
        // if (epsilon > 1 && rec == 2) { FFLAS::finit(G,m2,n3,C31,ldc);}
        // C32
        FFLAS::fsubin(NoField,m2,n3,C13,ldc,C23,ldc);
        // S3
        add(k2,n3,epsilon,B13,ldb,B23,ldb,Y,n3);
        // T3
        add(m2,k2,epsilon,A12,lda,A11,lda,X,k2);
        // P3
        gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C13,ldc,rec-1,epsilon);
        FFLAS::fflas_delete( Y );
        FFLAS::fflas_delete( X );
        // C21
        FFLAS::fsubin(NoField,m2,n3,C13,ldc,C12,ldc);
        // P8
        Y = FFLAS::fflas_new<double>(m2*n3);
        gemm_bini_223_mem(F,m2,n3,k2,A11,lda,B23,ldb,Y,n3,rec-1,epsilon);
        // C31
        subscalinf(NoField,m2,n3,Y,n3,(double)1/epsilon,C13,ldc);
        // C32
        subscalinf(NoField,m2,n3,Y,n3,(double)1/epsilon,C23,ldc);
        FFLAS::fflas_delete( Y );


        return C;

    }
#endif


} // Rec
} // Protected
} // FFLAS

namespace FFLAS { namespace Protected {

    template<class Field>
    typename Field::Element *
    gemm_bini_p(const Field &F
                , const size_t m
                , const size_t n
                , const size_t k
                , const typename Field::Element *A
                , const size_t lda
                , const typename Field::Element *B
                , const size_t ldb
                , typename Field::Element *C
                , const size_t ldc
                , int rec
                , size_t algo
               )
    {

        assert(k/6*6==k); // k divisible par 6
        assert(n/6*6==n); // n divisible par 6
        assert(m/6*6==m); // m divisible par 6

        // e-formule
        double epsilon = (double) F.characteristic() ;
        switch(algo) {
        case 0 :
            Rec::gemm_bini_322_mem(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
            FFLAS::finit_fuzzy(F,m,n,C,ldc);
            // FFLAS::finit(F,m,n,C,ldc);
            break;
        case 1 :
            Rec::gemm_bini_322_0(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
            FFLAS::finit_fuzzy(F,m,n,C,ldc);
            // FFLAS::finit(F,m,n,C,ldc);
            break;
        case 2 :
            Rec::gemm_bini_322_2(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
            FFLAS::finit_fuzzy(F,m,n,C,ldc);
            break;
        case 3 :
            Rec::gemm_bini_223_mem(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
            FFLAS::finit_fuzzy(F,m,n,C,ldc);
            // FFLAS::finit(F,m,n,C,ldc);
            break;
        case 4 :
            Rec::gemm_bini_232_2(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
            FFLAS::finit_fuzzy(F,m,n,C,ldc);
            break;
        case 5 :
            Rec::gemm_bini_232_3(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
            FFLAS::finit_fuzzy(F,m,n,C,ldc);
            break;
#if 0
        case 8 : {
                     double epsilon2 = sqrt((double)epsilon);
                     std::cout << epsilon2 << std::endl;
                     Rec::gemm_bini_322_sqrt(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon2);
                     // FFLAS::finit_fuzzy(F,m,n,C,ldc);
                     for(size_t i = 0  ; i < m ; ++i) {
                         for(size_t j = 0  ; j < n ; ++j)
                             C[i*ldc+j] = rint(fmod(C[i*ldc+j],epsilon2));
                     }
                     break;
                 }
#endif
        default :
                 std::cout << " not an algo :" << algo << std::endl;;
                 exit(-1);
        }



        return C;

    }

    template<class Field>
    typename Field::Element *
    gemm_bini_e(const Field &F
                , const size_t m
                , const size_t n
                , const size_t k
                , const typename Field::Element *A
                , const size_t lda
                , const typename Field::Element *B
                , const size_t ldb
                , typename Field::Element *C
                , const size_t ldc
                , int rec
                , size_t algo
               )
    {

        assert(k/2*2==k); // k divisible par 2
        assert(n/2*2==n); // n divisible par 2
        assert(m/3*3==m); // m divisible par 3

        // e-formule
        double epsilon = 1./(1<<27);
        switch(algo) {
        case 0 :
            Rec::gemm_bini_322_mem(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
            break;
        case 1 :
            Rec::gemm_bini_322_0(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
            break;
        case 2 :
            Rec::gemm_bini_322_2(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
            break;
        case 3 :
            Rec::gemm_bini_223_mem(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
            break;
        case 4 :
            Rec::gemm_bini_232_2(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
            break;
        case 5 :
            Rec::gemm_bini_232_3(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
            break;
        default :
            std::cout << " not an algo :" << algo << std::endl;;
            exit(-1);
        }


        // vire les e.
        // FFLAS::finit_fuzzy(F,m,n,C,ldc);
        FFLAS::finit_fuzzy(F,m,n,C,ldc);

        return C;

    }

    template<class Field>
    typename Field::Element *
    gemm_compress(const Field &F
                  , const size_t m
                  , const size_t n
                  , const size_t k
                  , const typename Field::Element *A
                  , const size_t lda
                  , const typename Field::Element *B
                  , const size_t ldb
                  , typename Field::Element *C
                  , const size_t ldc
                  , int rec
                  , size_t algo
                 )
    {

        assert(k/6*6==k); // k divisible par 6
        assert(n/6*6==n); // n divisible par 6
        assert(m/6*6==m); // m divisible par 6

        switch(algo) {
        case 0 :
            fgemm_compressed<Field,true>(F,(int)m,(int)n,(int)k,A,(int)lda,B,(int)ldb,C,(int)ldc);
            FFLAS::freduce(F,m,n,C,ldc);
            break;
        case 1 :
            fgemm_compressed<Field,false>(F,(int)m,(int)n,(int)k,A,(int)lda,B,(int)ldb,C,(int)ldc);
            FFLAS::freduce(F,m,n,C,ldc);
            break;
        default :
            std::cout << " not an algo :" << algo << std::endl;;
            exit(-1);
        }



        return C;

    }

} // Protected
} // FFLAS

template<class Field>
void check_equal(const Field & F,int m,int n,
                 typename Field::Element * D,int ldd,
                 typename Field::Element * E,int lde,
                 const char * nomalgo, const char * madescr, int & ok_p)
{
    int faux = 0 ;
    for (int i = 0 ; i < m ; ++i) {
        for (int j = 0 ; j < n ; ++j) {
            if (!F.areEqual(D[i*ldd+j],E[i*lde+j])) {
                ++faux ;
            }
        }
    }
    if (faux) {
        std::cout << nomalgo << " " << madescr << " : bad/all = " << faux << '/' << m*n << " ~~  " << (double)faux/(double)(m*n) << std::endl;
    }
    else ok_p ++ ;


#if 1
    if (faux && (n<20)) {
        std::cout << "OK" << std::endl;
        for (int i = 0 ; i < m ; ++i) {
            for (int j = 0 ; j < n ; ++j)
                std::cout << D[i*ldd+j] << ' ';
            std::cout << std::endl;
        }
        std::cout << "KO" << std::endl;
        for (int i = 0 ; i < m ; ++i) {
            for (int j = 0 ; j < n ; ++j)
                std::cout << E[i*lde+j] << ' ';
            std::cout << std::endl;
        }


        std::cout << "Diff" << std::endl;
        for (int i = 0 ; i < m ; ++i) {
            for (int j = 0 ; j < n ; ++j)
                std::cout << D[i*ldd+j]-E[i*lde+j] << ' ';
            std::cout << std::endl;
        }
    }
#endif
}


template<class Field>
void test_algos(const Field &F, int m, int n, int k
                , const typename Field::Element * A, int lda
                , const typename Field::Element * B, int ldb
                , int r
                , time_v & tim_k, time_v & tim_e , time_v & tim_p
                , int_v & ok_k, int_v & ok_e, int_v & ok_p
                , FFLAS::Timer & tim_wd, int & nb_wd
                , bool with_e
                , bool with_k
               )
{
    FFLAS::Timer tmp ;
    typedef typename Field::Element Element;

    Element * D = FFLAS::fflas_new<Element>(m*n);
    Element * C = FFLAS::fflas_new<Element>(m*n);

    tmp.clear();tmp.start();
    fgemm(F,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
          m,n,k, 1, A,k, B,n, 0, D, n);
    tmp.stop(); tim_wd += tmp ; nb_wd ++;

    /*  bini_p */
    if (with_e) {
        for (int algo = 0 ; algo < algos ; ++algo) {
            tmp.clear();tmp.start();
            FFLAS::Protected::gemm_bini_e(F,m,n,k,A,k,B,n,C,n,r,selec[algo]);
            tmp.stop(); tim_e[algo] += tmp ;

            /*  checking  */
            check_equal(F,m,n,D,n,C,n,"bini_e",descr[algo],ok_e[algo]) ;
        }
    }

    /*  compress */
    if (with_k && std::is_same<typename FieldTraits<Field>::category,FFLAS::FieldCategories::ModularTag>::value && (! FieldTraits<Field>::balanced)) {
        for (int algo = 0 ; algo < algos_k ; ++algo) {
            tmp.clear();tmp.start();
            FFLAS::Protected::gemm_compress(F,m,n,k,A,k,B,n,C,n,r,selec_k[algo]);
            tmp.stop(); tim_k[algo] += tmp ;

            /*  checking  */
            check_equal(F,m,n,D,n,C,n,"compress",descr_k[algo],ok_k[algo]) ;


        }
    }

    /*  bini_p */
    for (int algo = 0 ; algo < algos ; ++algo) {
        tmp.clear();tmp.start();
        FFLAS::Protected::gemm_bini_p(F,m,n,k,A,k,B,n,C,n,r,selec[algo]);
        tmp.stop(); tim_p[algo] += tmp ;

        /*  checking  */
        check_equal(F,m,n,D,n,C,n,"bini_p",descr[algo],ok_p[algo]) ;


    }

    FFLAS::fflas_delete(C);
    FFLAS::fflas_delete(D);
}

template<class T>
struct changeField {
    typedef T other ;
};

template<>
struct changeField<Modular<double> > {
    typedef Givaro::Modular<float> other;
};

template<>
struct changeField<ModularBalanced<double> > {
    typedef ModularBalanced<float> other;
};

double descrip(int algo, int_v & ok_e, time_v & tim_e, int iters, const char ** madescr, const char * nom)
{
    int min_e = -1 ;
    double bini_e = -1 ;
    for (int i = 0 ; i < algo ; ++i){
        if (ok_e[i] == (int)iters) {
            double bini1 = tim_e[i].usertime()/(double)ok_e[i] ;
            if (bini_e <  0) {
                bini_e = bini1;
                min_e = (int) i ;
            }
            else if (bini1 < bini_e) {
                min_e  = (int)i ;
                bini_e = bini1 ;
            }
        }
    }
    for (int i = 0 ; i < algo ; ++i){
        if (ok_e[i] == (int)iters) {
            double bini1 = tim_e[i].usertime()/(double)ok_e[i] ;
            std::cout << nom << " ( " << madescr[i] << " ) : " ;
            if ((int)i == min_e) std::cout << " * " ;
            else std::cout << "   ";
            std::cout << bini1  << 's'<<  std::endl;
        }
    }

    return bini_e ;
}


template<class Field>
void test(int m, int k, int n, int p, int r, bool with_e, bool with_k, 	int iters = 4, uint64_t seed=0)
{

    typedef typename Field::Element Element;

    Element * A = FFLAS::fflas_new<Element>(m*k);
    Element * B = FFLAS::fflas_new<Element>(n*k);


    Field F(p);
    typename Field::RandIter G(F,seed);
    F.write(std::cout<< " * Field " ) << std::endl;

    typedef typename changeField<Field>::other Field_f  ;
    typedef typename Field_f::Element Element_f ;
    Field_f F_f(p);
    Element_f * A_f = FFLAS::fflas_new<Element_f>(m*k);
    Element_f * B_f = FFLAS::fflas_new<Element_f>(n*k);
    Element_f * C_f = FFLAS::fflas_new<Element_f>(m*n);

#if defined(NOTRANDOM)
    int i0 ;
    int j0 ;
    Element p2 ; F.init(p2,(int)F.mOne/2);
    std::cout << p2 << std::endl;
#warning "not random"
    for (int i = 0 ; i < m ; ++i)
        for (int j = 0 ; j < k ; ++j) {
            i0 = i/(m/3);
            j0 = j/(k/2);
            if      (i0 == 0 and j0 == 0) A[i*k+j] = F.mOne ;
            else if (i0 == 0 and j0 == 1) A[i*k+j] = F.zero ;
            else if (i0 == 1 and j0 == 0) A[i*k+j] = F.mOne ;
            else if (i0 == 1 and j0 == 1) A[i*k+j] = F.mOne ;
            else if (i0 == 2 and j0 == 0) A[i*k+j] = F.mOne ;
            else if (i0 == 2 and j0 == 1) A[i*k+j] = F.mOne ;
            else A[i*k+j] = F.mOne ;
        }
    for (int i = 0 ; i < k ; ++i)
        for (int j = 0 ; j < n ; ++j) {
            i0 = i/(k/2);
            j0 = j/(n/2);
            if      (i0 == 0 and j0 == 0) B[i*n+j] = F.mOne ;
            else if (i0 == 0 and j0 == 1) B[i*n+j] = F.mOne ;
            else if (i0 == 1 and j0 == 0) B[i*n+j] = F.mOne ;
            else if (i0 == 1 and j0 == 1) B[i*n+j] = F.zero ;
            else B[i*n+j] = F.mOne ;

        }
#endif

    time_v tim_e(algos), tim_p(algos), tim_k(algos_k);
    FFLAS::Timer tim_wd; tim_wd.clear();
    FFLAS::Timer tim_wf; tim_wf.clear();
    FFLAS::Timer tmp;
    for (int i = 0 ; i < algos ; ++i) {
        tim_e[i].clear();
        tim_p[i].clear();
    }
    for (int i = 0 ; i < algos_k ; ++i) {
        tim_k[i].clear();
    }

    int_v ok_p(algos,0),  ok_e(algos,0), ok_k(algos_k,0);
    int nb_wd = 0 , nb_wf = 0 ;

    for (int b = 0 ; b < iters ; ++b) {
        std::cout << "iter " << b+1 << " of " << iters << std::endl;
#if not defined(NOTRANDOM)
        FFPACK::RandomMatrix(F, m, k, A, k, G);
        FFPACK::RandomMatrix(F, k, n, B, n, G);
#endif
        FFLAS::finit(F_f,m,k,A,k,A_f,k);
        FFLAS::finit(F_f,k,n,B,n,B_f,n);

        tmp.clear();tmp.start();
        fgemm(F_f,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
              m,n,k, 1, A_f,k, B_f,n, 0, C_f, n);
        tmp.stop(); tim_wf += tmp ; nb_wf ++ ;

        test_algos(F,m,n,k,A,k,B,n,r,
                   tim_k,tim_e,tim_p,
                   ok_k,ok_e,ok_p,
                   tim_wd,nb_wd,
                   with_e,with_k);
    }
    std::cout  << std::endl << "results" << std::endl;

    double bini_e = descrip(algos,ok_e,tim_e,iters,descr,"Bini_e");
    double bini_p = descrip(algos,ok_p,tim_p,iters,descr,"Bini_p");
    double bini_k = descrip(algos_k,ok_k,tim_k,iters,descr_k,"Bini_k");


    double t_wd = tim_wd.usertime()/(double)(nb_wd);
    double t_wf = tim_wf.usertime()/(double)(nb_wf);

    std::cout << "Wino d : " << t_wd << 's'<<  std::endl;
    std::cout << "Wino f : " << t_wf << 's'<<  std::endl;
    double wino =  std::min(t_wd,t_wf) ;
    if (bini_e>=0)
        std::cout << "Gain e: " << ((bini_e-wino)/wino)*100 << '%' << std::endl;
    if (bini_p>=0)
        std::cout << "Gain p: " << ((bini_p-wino)/wino)*100 << '%' << std::endl;
    if (bini_k>=0)
        std::cout << "Gain k: " << ((bini_k-wino)/wino)*100 << '%' << std::endl;




    FFLAS::fflas_delete( A );
    FFLAS::fflas_delete( B);

    FFLAS::fflas_delete( A_f );
    FFLAS::fflas_delete( B_f);
    FFLAS::fflas_delete( C_f);
}



int main(int ac, char **av) {
    static int m = 36 ;
    static int n = 12 ;
    static int k = 18 ;
    static int p = 101;
    bool  eps = false ;
    bool  kom = false ;
    int r = 1 ;
    uint64_t seed = getSeed();
    int iters = 4;

    static Argument as[] = {
        { 'p', "-p P", "Set the field characteristic.",  TYPE_INT , &p },
        { 'n', "-n N", "Set the number of cols in C.",   TYPE_INT , &n },
        { 'm', "-m N", "Set the number of rows in C.",   TYPE_INT , &m },
        { 'k', "-k N", "Set the number of rows in B.",   TYPE_INT , &k },
        { 'r', "-k N", "Set the recursive number Bini.", TYPE_INT , &r },
        { 'i', "-i N", "Set the numebr of iterations.", TYPE_INT , &iters },
        { 's', "-s N", "Set the seed                 .", TYPE_UINT64 , &seed },
        { 'e', "-e " , "epsilon                 .", TYPE_NONE , &eps },
        { 'c', "-c " , "compress                .", TYPE_NONE , &kom},
        END_OF_ARGUMENTS
    };
    FFLAS::parseArguments(ac,av,as);

    srand(seed);
    srand48(seed);

    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "size: " << m << ',' << k << ',' << n << std::endl;
    std::cout << "mod : " << p << std::endl;
    std::cout << "rec : " << r << std::endl;
    std::cout << "seed: " << seed << std::endl;
    std::cout << "thre: " << TRE << std::endl;
    std::cout << "=====================================================" << std::endl;
    test<Modular<double> > (m,k,n,p,r,eps,kom,iters,seed);
    std::cout << "=====================================================" << std::endl;
    test<ModularBalanced<double> > (m,k,n,p,r,eps,kom,iters,seed);
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
