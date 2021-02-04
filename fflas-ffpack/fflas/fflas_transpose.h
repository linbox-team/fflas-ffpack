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

        /* transpose for vect_size == 2 */
        template <size_t _s = Simd::vect_size, IsSimdSize<_s, 2>* = nullptr>
        void transpose (const Field& F, ConstElement_ptr A, size_t lda,
                                            Element_ptr B, size_t ldb) const {
            vect_t R0,R1;
            LD(0);LD(1);
            Simd::transpose (R0, R1);
            ST(0);ST(1);
        }

        /* transpose for vect_size == 4 */
        template <size_t _s = Simd::vect_size, IsSimdSize<_s, 4>* = nullptr>
        void transpose (const Field& F, ConstElement_ptr A, size_t lda,
                                            Element_ptr B, size_t ldb) const {
            vect_t R0,R1,R2,R3;
            LD(0);LD(1);LD(2);LD(3);
            Simd::transpose (R0, R1, R2, R3);
            ST(0);ST(1);ST(2);ST(3);
        }

        /* transpose for vect_size == 8 */
        template <size_t _s = Simd::vect_size, IsSimdSize<_s, 8>* = nullptr>
        void transpose (const Field& F, ConstElement_ptr A, size_t lda,
                                            Element_ptr B, size_t ldb) const {
            vect_t R0,R1,R2,R3,R4,R5,R6,R7;
            LD(0);LD(1);LD(2);LD(3);LD(4);LD(5);LD(6);LD(7);
            Simd::transpose (R0, R1, R2, R3, R4, R5, R6, R7);
            ST(0);ST(1);ST(2);ST(3);ST(4);ST(5);ST(6);ST(7);
        }

        /* transpose for vect_size == 16 */
        template <size_t _s = Simd::vect_size, IsSimdSize<_s, 16>* = nullptr>
        void transpose (const Field& F, ConstElement_ptr A, size_t lda,
                                            Element_ptr B, size_t ldb) const {
            vect_t R0,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15;
            LD(0);LD(1);LD(2);LD(3);LD(4);LD(5);LD(6);LD(7);LD(8);LD(9);LD(10);LD(11);LD(12);LD(13);LD(14);LD(15);
            Simd::transpose (R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11,
                             R12, R13, R14, R15);
            ST(0);ST(1);ST(2);ST(3);ST(4);ST(5);ST(6);ST(7);ST(8);ST(9);ST(10);ST(11);ST(12);ST(13);ST(14);ST(15);
        }

#undef LD
#undef ST
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
                    for (size_t iv = ib; iv+vs <= ibend; iv+=vs) {
                        for (size_t jv = jb; jv+vs <= jbend; jv+=vs) {
                            BTS.transpose(F, A+iv*lda+jv, lda, B+jv*ldb+iv, ldb);
                        }
                    }
                }
            }

            size_t last_i = (m / vs) * vs;
            /* remaining cols that cannot be handled with Simd */
            for (size_t j = (n / vs) * vs; j < n; j++) {
                fassign (F, last_i, A+j, lda, B+j*ldb, 1);
            }
            /* remaining rows that cannot be handled with Simd */
            for (size_t i = last_i; i < m; i++) {
                fassign (F, n, A+i*lda, 1, B+i, ldb);
            }
        }

        /* m is the number of rows and cols (A is a square matrice)
         * lda must be >= m
         *
         * This variant does not separate diagonal block ffrom the other ones.
         * remark: diagonal blocks are transposed and written twice
         */
        template <size_t bs, typename Field, typename BTSimd>
        void square_inplace (const Field& F, const BTSimd& BTS, const size_t m,
                             typename Field::Element_ptr A, const size_t lda) {
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
                            fassign (F, vs, TMP2+(j-jb)*bs+(iv-ib), 1, A+iv*lda+j, lda);
                            fassign (F, vs, TMP1+(iv-ib)*bs+(j-jb), bs, A+j*lda+iv, 1);
                        }
                    }
                    /* remaining rows that cannot be handled with Simd */
                    for (size_t i = iv; i < ibend; i++) {
                        size_t j = iv;
                        fassign (F, jbend-j, TMP2+(j-jb)*bs+(i-ib), bs, A+i*lda+j, 1);
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
            /* Memory area must not overlap: not easy to check with
             * standard-compliant code, see https://stackoverflow.com/questions/51699331/how-to-check-that-two-arbitrary-memory-ranges-are-not-overlapped-in-c-c
             */

            _ftranspose_impl::not_inplace<bs> (F, BTS, m, n, A, lda, B, ldb);
        }
        else if (m == n) { /* in place with square matrices */
            /* both strides must be identical */
            FFLASFFPACK_check (lda == ldb);

            _ftranspose_impl::square_inplace<bs> (F, BTS, m, B, ldb);
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
