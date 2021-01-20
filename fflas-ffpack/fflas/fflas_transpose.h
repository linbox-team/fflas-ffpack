/*
 * Copyright (C) 2020 the FFLAS-FFPACK group
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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
/** @file fflas/fflas_transpose.h
 * @brief transpose the storage of the matrix (switch between row and col major mode)
 */


#ifndef __FFLASFFPACK_transpose_H
#define __FFLASFFPACK_transpose_H
#include "fflas-ffpack/utils/debug.h"
#include "fflas-ffpack/fflas/fflas.h"

#ifndef FFLAS_TRANSPOSE_BLOCKSIZE 
#define FFLAS_TRANSPOSE_BLOCKSIZE 32 // MUST BE A POWER OF TWO
#endif

#include "fflas-ffpack/utils/fflas_memory.h"
#include "fflas-ffpack/fflas/fflas_simd.h"

namespace FFLAS {

    /**************************************************************************/
    template <typename Field, typename Simd,
                typename std::enable_if<Simd::template is_same_element<Field>::value>::type* = nullptr>
    struct BlockTransposeSIMD {
    private:
        using Element          = typename Field::Element;
        using Element_ptr      = typename Field::Element_ptr;
        using ConstElement_ptr = typename Field::ConstElement_ptr;
        using vect_t           = typename Simd::vect_t;

        template <bool B, class T = void>
        using enable_if_t = typename std::enable_if<B, T>::type;

        template <size_t s1, size_t s2, class T = void>
        using IsSimdSize = enable_if_t<s1 == s2 && Simd::vect_size == s1, T>;
    public:
        /* */
        static inline constexpr size_t size() { return Simd::vect_size; }

        /* */
        static inline const std::string info() {
            return "transpose blocks with " + Simd::type_string()
                    + " [vect_size=" + std::to_string(Simd::vect_size) + "]";
        }

        /* transpose for vect_size == 1 */
        template <size_t _s = Simd::vect_size, IsSimdSize<_s, 1>* = nullptr>
        void transpose (const Field& F, ConstElement_ptr A, size_t lda,
                                            Element_ptr B, size_t ldb) const {
            F.assign (*B, *A);
        }

#define LD(i) R##i=Simd::loadu(A+lda*i)
#define ST(i) Simd::storeu(B+ldb*i,R##i)
#define PCK(i,j) Simd::unpacklohi(R##i,R##j,R##i,R##j);

        /* transpose for vect_size == 2 */
        template <size_t _s = Simd::vect_size, IsSimdSize<_s, 2>* = nullptr>
        void transpose (const Field& F, ConstElement_ptr A, size_t lda,
                                            Element_ptr B, size_t ldb) const {
            vect_t R0,R1;
            LD(0);LD(1);
            PCK(0,1);
            ST(0);ST(1);
        }

        /* transpose for vect_size == 4 */
        template <size_t _s = Simd::vect_size, IsSimdSize<_s, 4>* = nullptr>
        void transpose (const Field& F, ConstElement_ptr A, size_t lda,
                                            Element_ptr B, size_t ldb) const {
            vect_t R0,R1,R2,R3;
            LD(0);LD(1);LD(2);LD(3);
            PCK(0,2); PCK(1,3);
            PCK(0,1); PCK(2,3);
            ST(0);ST(1);ST(2);ST(3);
        }

        /* transpose for vect_size == 8 */
        template <size_t _s = Simd::vect_size, IsSimdSize<_s, 8>* = nullptr>
        void transpose (const Field& F, ConstElement_ptr A, size_t lda,
                                            Element_ptr B, size_t ldb) const {
            vect_t R0,R1,R2,R3,R4,R5,R6,R7;
            LD(0);LD(1);LD(2);LD(3);LD(4);LD(5);LD(6);LD(7);
            PCK(0,4); PCK(1,5); PCK(2,6); PCK(3,7);
            PCK(0,2); PCK(1,3); PCK(4,6); PCK(5,7);
            PCK(0,1); PCK(2,3); PCK(4,5); PCK(6,7);
            ST(0);ST(1);ST(2);ST(3);ST(4);ST(5);ST(6);ST(7);
        }

        /* transpose for vect_size == 16 */
        template <size_t _s = Simd::vect_size, IsSimdSize<_s, 16>* = nullptr>
        void transpose (const Field& F, ConstElement_ptr A, size_t lda,
                                            Element_ptr B, size_t ldb) const {
            vect_t R0,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15;
            LD(0);LD(1);LD(2);LD(3);LD(4);LD(5);LD(6);LD(7);LD(8);LD(9);LD(10);LD(11);LD(12);LD(13);LD(14);LD(15);
            PCK(0,8); PCK(1,9); PCK(2,10); PCK(3,11); PCK(4,12); PCK(5,13);  PCK(6,14);  PCK(7,15);
            PCK(0,4); PCK(1,5); PCK(2,6);  PCK(3,7);  PCK(8,12); PCK(9,13);  PCK(10,14); PCK(11,15);
            PCK(0,2); PCK(1,3); PCK(4,6);  PCK(5,7);  PCK(8,10); PCK(9,11);  PCK(12,14); PCK(13,15);
            PCK(0,1); PCK(2,3); PCK(4,5);  PCK(6,7);  PCK(8,9);  PCK(10,11); PCK(12,13); PCK(14,15);
            ST(0);ST(1);ST(2);ST(3);ST(4);ST(5);ST(6);ST(7);ST(8);ST(9);ST(10);ST(11);ST(12);ST(13);ST(14);ST(15);
        }

        /* transpose for vect_size == 32 */
        template <size_t _s = Simd::vect_size, IsSimdSize<_s, 32>* = nullptr>
        void transpose (const Field& F, ConstElement_ptr A, size_t lda,
                                            Element_ptr B, size_t ldb) const {
            vect_t R0,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15,R16,R17,R18,R19,R20,R21,R22,R23,R24,R25,R26,R27,R28,R29,R30,R31;
            LD(0);LD(1);LD(2);LD(3);LD(4);LD(5);LD(6);LD(7);LD(8);LD(9);LD(10);LD(11);LD(12);LD(13);LD(14);LD(15);
            LD(16);LD(17);LD(18);LD(19);LD(20);LD(21);LD(22);LD(23);LD(24);LD(25);LD(26);LD(27);LD(28);LD(29);LD(30);LD(31);
            PCK(0,16); PCK(1,17); PCK(2,18); PCK(3,19); PCK(4,20); PCK(5,21);  PCK(6,22);  PCK(7,23);  PCK(8,24);  PCK(9,25);  PCK(10,26); PCK(11,27); PCK(12,28); PCK(13,29); PCK(14,30); PCK(15,31);
            PCK(0,8);  PCK(1,9);  PCK(2,10); PCK(3,11); PCK(4,12); PCK(5,13);  PCK(6,14);  PCK(7,15);  PCK(16,24); PCK(17,25); PCK(18,26); PCK(19,27); PCK(20,28); PCK(21,29); PCK(22,30); PCK(23,31);
            PCK(0,4);  PCK(1,5);  PCK(2,6);  PCK(3,7);  PCK(8,12); PCK(9,13);  PCK(10,14); PCK(11,15); PCK(16,20); PCK(17,21); PCK(18,22); PCK(19,23); PCK(24,28); PCK(25,29); PCK(26,30); PCK(27,31);
            PCK(0,2);  PCK(1,3);  PCK(4,6);  PCK(5,7);  PCK(8,10); PCK(9,11);  PCK(12,14); PCK(13,15); PCK(16,18); PCK(17,19); PCK(20,22); PCK(21,23); PCK(24,26); PCK(25,27); PCK(28,30); PCK(29,31);
            PCK(0,1);  PCK(2,3);  PCK(4,5);  PCK(6,7);  PCK(8,9);  PCK(10,11); PCK(12,13); PCK(14,15); PCK(16,17); PCK(18,19); PCK(20,21); PCK(22,23); PCK(24,25); PCK(26,27); PCK(28,29); PCK(30,31);
            ST(0);ST(1);ST(2);ST(3);ST(4);ST(5);ST(6);ST(7);ST(8);ST(9);ST(10);ST(11);ST(12);ST(13);ST(14);ST(15);
            ST(16);ST(17);ST(18);ST(19);ST(20);ST(21);ST(22);ST(23);ST(24);ST(25);ST(26);ST(27);ST(28);ST(29);ST(30);ST(31);
        }

#undef LD
#undef ST
#undef PCK
    };


    /**************************************************************************/
    namespace _ftranspose_impl {
        /* It is the caller responsability to call those function with the
         * correct template parameters and function parameters.
         * In particular, bs (the block size) must always be a multiple of
         * Simd::vect_size.
         */

        /* m is the number of rows
         * n is the number of cols
         * lda must be >= n
         * ldb must be >= m
         * memory area of A and B must not overlap
         */
        template <size_t bs, typename Field, typename BTSimd>
        void not_inplace (const Field& F, const BTSimd& BTS, const size_t m,
                            const size_t n, typename Field::ConstElement_ptr A,
                            const size_t lda, typename Field::Element_ptr B,
                            const size_t ldb) {
            const size_t vs = BTS.size();

            for (size_t ib = 0; ib < m; ib+=bs) {
                const size_t ibend = std::min (m, ib+bs);
                for (size_t jb = 0; jb < n; jb+=bs) {
                    const size_t jbend = std::min (n, jb+bs);
                    size_t iv;
                    for (iv = ib; iv+vs <= ibend; iv+=vs) {
                        size_t jv;
                        for (jv = jb; jv+vs <= jbend; jv+=vs) {
                            BTS.transpose(F, A+iv*lda+jv, lda, B+jv*ldb+iv, ldb);
                        }
                        /* remaining cols that cannot be handled with Simd */
                        for (size_t j = jv; j < jbend; j++) {
                            for (size_t i = iv; i < iv+vs; i++) {
                                F.assign (*(B+j*ldb+i), *(A+i*lda+j));
                            }
                        }
                    }
                    /* remaining rows that cannot be handled with Simd */
                    for (size_t i = iv; i < ibend; i++) {
                        for (size_t j = jb; j < jbend; j++) {
                            F.assign (*(B+j*ldb+i), *(A+i*lda+j));
                        }
                    }
                }
            }
        }

        /* m is the number of rows and cols (A is a square matrice)
         * lda must be >= m
         */
        template <size_t bs, typename Field, typename BTSimd>
        void square_inplace_v1 (const Field& F, const BTSimd& BTS,
                                const size_t m, typename Field::Element_ptr A,
                                const size_t lda) {
            const size_t n = m;
            const size_t vs = BTS.size();

            typename Field::Element TMP[bs*bs];
            finit (F, bs, bs, TMP, bs);

            for (size_t ib = 0; ib < m; ib+=bs) {
                /* First, we tranpose diagonal blocks [ib..ib+bs, ib..ib+bs] */
                const size_t l = std::min (bs, m-ib);
                const size_t ibend = std::min (m, ib+bs);
                size_t iv;

                fassign(F, l, l, A+ib*lda+ib, lda, TMP, bs);

                for (iv = ib; iv+vs <= ibend; iv+=vs) {
                    const size_t ivend = std::min (m, iv+vs);
                    size_t jv;

                    BTS.transpose(F, A+iv*lda+iv, lda, A+iv*lda+iv, lda);

                    for (jv = iv+vs; jv+vs <= ibend; jv+=vs){
                        BTS.transpose(F, A+jv*lda+iv, lda, A+iv*lda+jv, lda);
                        BTS.transpose(F, TMP+(iv-ib)*bs+(jv-ib), bs, A+jv*lda+iv, lda);
                    }
                    /* remaining cols that cannot be handled with Simd */
                    for (size_t i = iv; i < ivend; i++) {
                        for (size_t j = jv; j < ibend; j++) {
                            F.assign (*(A+i*lda+j), *(A+j*lda+i));
                            F.assign (*(A+j*lda+i), *(TMP+(i-ib)*bs+(j-ib)));
                        }
                    }
                }
                /* remaining rows that cannot be handled with Simd */
                for (size_t i = iv; i < ibend; i++) {
                    for (size_t j = i+1; j < ibend; j++){
                        F.assign (*(A+i*lda+j), *(A+j*lda+i));
                        F.assign (*(A+j*lda+i), *(TMP+(i-ib)*bs+(j-ib)));
                    }
                }
                
                /* Then, we transpose and swap diagonal blocks
                 * [i..i+bs, j..j+bs] and [j..j+bs,i..i+bs]
                 */
                for (size_t jb = ib+bs; jb < n; jb+=bs) {
                    const size_t jbend = std::min (n, jb+bs);

                    /* copy only the first block */
                    fassign(F, bs, std::min (bs, n-jb), A+ib*lda+jb, lda, TMP, bs);
                    for (size_t iv = ib; iv+vs <= ibend; iv+=vs) {
                        const size_t ivend = std::min (m, iv+vs);
                        size_t jv;
                        for (jv = jb; jv+vs <= jbend; jv+=vs) {
                            BTS.transpose(F, A+jv*lda+iv, lda, A+iv*lda+jv, lda);
                            BTS.transpose(F, TMP+(iv-ib)*bs+(jv-jb), bs, A+jv*lda+iv, lda);
                        }
                        /* remaining cols that cannot be handled with Simd */
                        for (size_t i = iv; i < ivend; i++) {
                            for (size_t j = jv; j < jbend; j++) {
                                F.assign (*(A+i*lda+j), *(A+j*lda+i));
                                F.assign (*(A+j*lda+i), *(TMP+(i-ib)*bs+(j-jb)));
                            }
                        }
                    }
                }
            }
        }

        /* m is the number of rows and cols (A is a square matrice)
         * lda must be >= m
         *
         * This variant does not separate diagonal block ffrom the other ones.
         * remark: diagonal blocks are transposed and written twice
         */
        template <size_t bs, typename Field, typename BTSimd>
        void square_inplace_v2 (const Field& F, const BTSimd& BTS,
                                const size_t m, typename Field::Element_ptr A,
                                const size_t lda) {
            const size_t n = m;
            const size_t vs = BTS.size();

            typename Field::Element TMP1[bs*bs], TMP2[bs*bs];
            finit (F, bs, bs, TMP1, bs);
            finit (F, bs, bs, TMP2, bs);

            for (size_t ib = 0; ib < m; ib+=bs) {
                const size_t l = std::min (bs, m-ib);
                const size_t ibend = std::min (m, ib+bs);

                for (size_t jb = ib; jb < n; jb+=bs) {
                    const size_t k = std::min (bs, n-jb);
                    const size_t jbend = std::min (n, jb+bs);
                    size_t iv;

	                fassign (F, l, k, A+ib*lda+jb, lda, TMP1, bs);
	                fassign (F, k, l, A+jb*lda+ib, lda, TMP2, bs);       
	                for (iv = ib; iv+vs <= ibend; iv+=vs) {
                        size_t jv;
	                    for (jv = jb; jv+vs <= jbend; jv+=vs) {
                            /* Since each block are copied in TMP1 and TMP2, we
                             * put back their transposed at the right position
                             * in the result
                             */
                            BTS.transpose(F, TMP2+(jv-jb)*bs+(iv-ib), bs, A+iv*lda+jv, lda);
	                        BTS.transpose(F, TMP1+(iv-ib)*bs+(jv-jb), bs, A+jv*lda+iv, lda);
	                    }
                        /* remaining cols that cannot be handled with Simd */
                        for (size_t j = jv; j < jbend; j++) {
                            for (size_t i = iv; i < iv+vs; i++) {
                                F.assign(*(A+i*lda+j),*(TMP2+(j-jb)*bs+(i-ib)));
                                F.assign(*(A+j*lda+i),*(TMP1+(i-ib)*bs+(j-jb)));
                            }
                        }
                    }
                    /* remaining rows that cannot be handled with Simd */
                    for (size_t i = iv; i < ibend; i++) {
                        for (size_t j = i+1; j < jbend; j++) {
                            F.assign(*(A+i*lda+j),*(TMP2+(j-jb)*bs+(i-ib)));
                            F.assign(*(A+j*lda+i),*(TMP1+(i-ib)*bs+(j-jb)));
                        }
                    }
                }
            }
        }

        /* m is the number of rows
         * n is the number of cols
         * lda is assumed to be == n
         * ldb is assumed to be == m
         *
         * Copy A and use the not_inplace function
         */
        template <size_t bs, typename Field, typename BTSimd>
        void nonsquare_inplace_v1 (const Field& F, const BTSimd& BTS,
                                    const size_t m, const size_t n,
                                    typename Field::Element_ptr A) {
            const size_t lda = n;
            const size_t ldb = m;
            typename Field::Element_ptr TMP = fflas_new (F, m, n);
            fassign (F, m, n, A, lda, TMP, lda);
            not_inplace<bs> (F, BTS, m, n, TMP, lda, A, ldb);
            fflas_delete (TMP);
        }

        /* m is the number of rows
         * n is the number of cols
         * lda is assumed to be == n
         * ldb is assumed to be == m
         *
         * Use algorithm from
         *  Bryan Catanzaro, Alexander Keller, Michael Garland, “A Decomposition
         *  for In-place Matrix Transposition”. Principles and Practices of
         *  Parallel Programming 2014, pages 193-206, Orlando, Florida.
         */
        template <size_t bs, typename Field, typename BTSimd>
        void nonsquare_inplace_v2 (const Field& F, const BTSimd& BTS,
                                    const size_t m, const size_t n,
                                    typename Field::Element_ptr A) {
            const size_t lda = n;

            size_t u, v, d;
            size_t c = Givaro::gcdext (d, u, v, m, n);
            size_t a = m / c;
            size_t b = n / c;
            typename Field::Element_ptr tmp = fflas_new (F, std::max (m, n));

            /* Column rotate */
            if (c > 1)
            {
                for (size_t j = 0; j < n; j++)
                {
                    for (size_t i = 0; i < m; i++)
                        tmp[i] = A[((i+j/b) % m)*lda + j];
                    for (size_t i = 0; i < m; i++)
                        A[i*lda + j] = tmp[i];
                }
            }
            /* Row shuffle */
            for (size_t i = 0; i < m; i++)
            {
                for (size_t j = 0; j < n; j++)
                    tmp[(((i+j/b) % m) + j*m) % n] = A[i*lda + j];
                for (size_t j = 0; j < n; j++)
                    A[i*lda + j] = tmp[j];
            }
            /* Column shuffle */
            for (size_t j = 0; j < n; j++)
            {
                for (size_t i = 0; i < m; i++)
                    tmp[i] = A[((j+i*n - i/a) % m)*lda + j];
                for (size_t i = 0; i < m; i++)
                    A[i*lda + j] = tmp[i];
            }
            fflas_delete (tmp);
        }
    }

    /**************************************************************************/
    /*
     * Perfom transposition on the matrix A and store the result in matrix B.
     *      B[j,i] = A[i,j]     for 0 <= i < m, 0 <= j < n
     *
     * Requirements:
     *  - m (=the number of rows of A) must be less than or equal to lda
     *  - n (=the number of cols of A) must be less than or equal to ldb
     *  - if A and B are distinct pointers, the memory areas of A and B must not
     *  overlap.
     *  - if A and A are the same pointer and m == n, lda and ldb must be equal.
     *  - if A and A are the same pointer and m != n, lda must be equal to n and
     *  ldb must be equal to m.
     */
    template <typename Field, typename Simd = Simd<typename Field::Element>,
                size_t bs=FFLAS_TRANSPOSE_BLOCKSIZE,
                typename std::enable_if<Simd::template is_same_element<Field>::value>::type* = nullptr,
                typename std::enable_if<bs >= 1 && bs % Simd::vect_size == 0>::type* = nullptr>
    inline typename Field::Element_ptr
    ftranspose (const Field& F, const size_t m, const size_t n,
                    typename Field::ConstElement_ptr A, const size_t lda,
                    typename Field::Element_ptr B, const size_t ldb)
    {
        /* stride of matrix A must be at least the number of columns */
        FFLASFFPACK_check (n <= lda);
        /* stride of matrix B must be at least the number of rows */
        FFLASFFPACK_check (m <= ldb);

        BlockTransposeSIMD<Field, Simd> BTS;
        //std::cout << "ftranspose: " << BTS.info() << std::endl;

        if (A != B) { /* not in place */
            /* Memory area must not overlap (not easy to check with
             * standard-compliant code, see https://stackoverflow.com/questions/51699331/how-to-check-that-two-arbitrary-memory-ranges-are-not-overlapped-in-c-c)
             */

            _ftranspose_impl::not_inplace<bs> (F, BTS, m, n, A, lda, B, ldb);
        }
        else if (m == n) { /* in place with square matrices */
            /* both strides must be identical */
            FFLASFFPACK_check (lda == ldb);

            _ftranspose_impl::square_inplace_v1<bs> (F, BTS, m, B, ldb);
            //_ftranspose_impl::square_inplace_v2<bs> (F, BTS, m, B, ldb);
        }
        else { /* in place with non square matrices */
            /* no strides are allowed in this case */
            FFLASFFPACK_check (n == lda);
            FFLASFFPACK_check (m == ldb);

            _ftranspose_impl::nonsquare_inplace_v1<bs> (F, BTS, m, n, B);
            //_ftranspose_impl::nonsquare_inplace_v2<bs> (F, BTS, m, n, B);
        }
        return B;
    }
} // end of FFLAS namespace


#endif // __FFLASFFPACK_transpoose_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
