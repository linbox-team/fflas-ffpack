/*
 * Copyright (C) the FFLAS-FFPACK group
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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
 *.
 */
#include <iomanip>
#include <iostream>
#include <random>

#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/fflas_memory.h"

#include <givaro/modular.h>
#include <recint/rint.h>

#include "fflas-ffpack/fflas/fflas_transpose.h"

using namespace std;
using namespace FFLAS;
using namespace FFPACK;

using Givaro::Modular;
using Givaro::ModularBalanced;

/******************************************************************************/
template <typename Elt>
class Test {
    public:
        using Field = Modular<Elt>;
        using Elt_ptr = typename Field::Element_ptr;
        using Residu = typename Field::Residu_t;
        template <bool B, class T = void>
        using enable_if_t = typename std::enable_if<B, T>::type;
        template <typename Simd>
        using is_same_element = typename Simd::template is_same_element<Field>;
        template <typename E>
        using enable_if_no_simd_t = enable_if_t<Simd<E>::vect_size == 1>;
        template <typename E>
        using enable_if_simd128_t = enable_if_t<sizeof(E)*Simd<E>::vect_size == 16>;
        template <typename E>
        using enable_if_simd256_t = enable_if_t<sizeof(E)*Simd<E>::vect_size == 32>;
        template <typename E>
        using enable_if_simd512_t = enable_if_t<sizeof(E)*Simd<E>::vect_size == 64>;

        /* ctor */
        Test () : F(cardinality()) {
        }

        /* */
        template <typename _E = Elt,
                  enable_if_t<!is_same<_E, Givaro::Integer>::value>* = nullptr>
        static Residu
        cardinality () {
            return Field::maxCardinality();
        }

        template <typename _E = Elt,
                  enable_if_t<is_same<_E, Givaro::Integer>::value>* = nullptr>
        static Residu
        cardinality () {
            /* Test Givaro::Integer with a large (=more than a word) modulus */
            return Givaro::Integer ("0x100000000000000000000000000000000");
        }

        /* main test function */
        template <typename Simd = NoSimd<Elt>,
                  enable_if_t<is_same_element<Simd>::value>* = nullptr>
        bool test_ftranspose (size_t m, size_t n, Elt_ptr A, size_t lda,
                                                    Elt_ptr B, size_t ldb) {
            Elt_ptr M, Mt;
            if (A == B) {
                M = fflas_new (F, ldb, lda);
                finit (F, ldb, lda, M, lda);
                fassign (F, m, n, A, lda, M, lda);
                Mt = M;
            }
            else {
                M = A;
                Mt = B;
            }
            ftranspose<Field, Simd> (F, m, n, M, lda, Mt, ldb);

            bool ok = true;
            for (size_t i = 0; i < m; i++)
                for (size_t j = 0; j < n; j++)
                    ok = ok && ((*(Mt+j*ldb+i)) == (*(A+i*lda+j)));

            if (!ok) {
                cerr << "Error, with " << m << "x" << n << " matrix"
                     << (lda == n ? "" : " with stride " + to_string(lda))
                     << (A == B ? " (inplace variant)" : "") << " over "
                     << F.type_string() << " using " << Simd::type_string()
                     << endl;
                for (size_t i = 0; i < m; i++)
                    for (size_t j = 0; j < n; j++)
                        if ((*(Mt+j*ldb+i)) != (*(A+i*lda+j))) {
                            cerr << "   at index (" << j << ", " << i
                                 << ") expected " << *(A+i*lda+j) << " got "
                                 << *(Mt+j*ldb+i) << endl;
                        }
            }

            if (A == B) {
                fflas_delete (M);
            }
            return ok;
        }

        /* perform the tests for all sizes and types of matrices */
        template <typename Simd = NoSimd<Elt>,
                  enable_if_t<is_same_element<Simd>::value>* = nullptr>
        bool doTests () {

            /* square matrices */
            bool ok_sq_not_inplace = true, ok_sq_inplace = true;
            size_t nrows[] = { 3*FFLAS_TRANSPOSE_BLOCKSIZE,
                               3*FFLAS_TRANSPOSE_BLOCKSIZE+Simd::vect_size,
                               3*FFLAS_TRANSPOSE_BLOCKSIZE+Simd::vect_size+3 };
            for (auto m: nrows) {
                Elt_ptr M = fflas_new (F, m, m);
                Elt_ptr Mt = fflas_new (F, m, m);

                finit (F, m, m, M, m);
                for (size_t i = 0; i < m; i++)
                    for (size_t j = 0; j < m; j++)
                        F.assign (M[i*m+j], (i << 8) + j);
                finit (F, m, m, Mt, m);

                /* not inplace full matrix */
                bool b1 = test_ftranspose<Simd> (m, m, M, m, Mt, m);

                /* not inplace submatrix */
                size_t s = m - FFLAS_TRANSPOSE_BLOCKSIZE;
                bool b2 = test_ftranspose<Simd> (s, s, M, m, Mt, m);

                ok_sq_not_inplace = ok_sq_not_inplace && b1 && b2;

                /* inplace full matrix */
                bool b3 = test_ftranspose<Simd> (m, m, M, m, M, m);

                /* inplace submatrix */
                s = m - Simd::vect_size;
                bool b4 = test_ftranspose<Simd> (s, s, M, m, M, m);

                ok_sq_inplace = ok_sq_inplace && b3 && b4;

                fflas_delete (M);
                fflas_delete (Mt);
            }

            /* non square matrices */
            bool ok_nsq_not_inplace = true, ok_nsq_inplace = true;
            size_t ncols[] = { 2*FFLAS_TRANSPOSE_BLOCKSIZE,
                               4*FFLAS_TRANSPOSE_BLOCKSIZE+Simd::vect_size,
                               3*FFLAS_TRANSPOSE_BLOCKSIZE+2*Simd::vect_size+1};
            for (size_t i = 0; i < 3; i++) {
                size_t m = nrows[i];
                size_t n = ncols[i];

                Elt_ptr M = fflas_new (F, m, n);
                Elt_ptr Mt = fflas_new (F, n, m);

                finit (F, m, n, M, n);
                for (size_t i = 0; i < m; i++)
                    for (size_t j = 0; j < n; j++)
                        F.assign (M[i*n+j], (i << 8) + j);
                finit (F, n, m, Mt, m);

                /* not inplace full matrix */
                bool b1 = test_ftranspose<Simd> (m, n, M, n, Mt, m);

                /* not inplace submatrix */
                size_t s = m - Simd::vect_size;
                size_t t = n - 2*Simd::vect_size+1;
                bool b2 = test_ftranspose<Simd> (s, t, M, n, Mt, m);

                ok_nsq_not_inplace = ok_nsq_not_inplace && b1 && b2;

                /* inplace full matrix */
                bool b3 = test_ftranspose<Simd> (m, n, M, n, M, m);

                ok_nsq_inplace = ok_nsq_not_inplace && b3;

                fflas_delete (M);
                fflas_delete (Mt);
            }

            /* print results */
            bool ok = ok_sq_not_inplace && ok_sq_inplace
                        && ok_nsq_not_inplace && ok_nsq_inplace;
            std::cout << F.type_string()
                      << string (36-F.type_string().size(), '.') << " "
                      << Simd::type_string()
                      << string (36-Simd::type_string().size(), '.') << " "
                      << (ok ? "PASSED" : "FAILED")
                      << endl;

            return ok;
        }

        /* run tests: call doTests for all available Simd structs */
        template <typename _E = Elt,
                  enable_if_t<is_same<_E, Elt>::value>* = nullptr,
                  enable_if_no_simd_t<_E>* = nullptr>
        bool run () {
            doTests ();
            return true;
        }

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
        template <typename _E = Elt,
                  enable_if_t<is_same<_E, Elt>::value>* = nullptr,
                  enable_if_t<Simd<_E>::vect_size != 1>* = nullptr,
                  enable_if_simd128_t<_E>* = nullptr>
        bool run () {
            doTests ();
            doTests<Simd128<Elt>> ();
            return true;
        }
#endif

#ifdef __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS
        template <typename _E = Elt,
                  enable_if_t<is_same<_E, Elt>::value>* = nullptr,
                  enable_if_t<Simd<_E>::vect_size != 1>* = nullptr,
                  enable_if_simd256_t<_E>* = nullptr>
        bool run () {
            doTests ();
            doTests<Simd128<Elt>> ();
            doTests<Simd256<Elt>> ();
            return true;
        }
#endif

#ifdef __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS
        template <typename _E = Elt,
                  enable_if_t<is_same<_E, Elt>::value>* = nullptr,
                  enable_if_t<Simd<_E>::vect_size != 1>* = nullptr,
                  enable_if_simd512_t<_E>* = nullptr>
        bool run () {
            doTests ();
            doTests<Simd128<Elt>> ();
            doTests<Simd256<Elt>> ();
            doTests<Simd512<Elt>> ();
            return true;
        }
#endif

    protected:
        Field F;
};

/******************************************************************************/
int main(int argc, char** argv)
{
    std::cout << std::setprecision(17);
    std::cerr << std::setprecision(17);

    bool ok = true;

    ok = ok && Test<float>().run();
    ok = ok && Test<double>().run();
    ok = ok && Test<uint64_t>().run();
    ok = ok && Test<int64_t>().run();
    ok = ok && Test<uint32_t>().run();
    ok = ok && Test<int32_t>().run();
    ok = ok && Test<uint16_t>().run();
    ok = ok && Test<int16_t>().run();
    ok = ok && Test<Givaro::Integer>().run();
    ok = ok && Test<RecInt::rint<6>>().run();
    ok = ok && Test<RecInt::ruint<6>>().run();
    ok = ok && Test<RecInt::rint<7>>().run();
    ok = ok && Test<RecInt::ruint<7>>().run();
    ok = ok && Test<RecInt::rint<8>>().run();
    ok = ok && Test<RecInt::ruint<8>>().run();

    return !ok;
}
