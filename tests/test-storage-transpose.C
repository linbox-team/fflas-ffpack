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
    Test (size_t mm, size_t nn) : F(cardinality()), _mm(mm), _nn(nn) {
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
        bool ok = true;

            /* square matrices */
        size_t nrows[] = { 3*FFLAS_TRANSPOSE_BLOCKSIZE,
                           3*FFLAS_TRANSPOSE_BLOCKSIZE+Simd::vect_size,
                           3*FFLAS_TRANSPOSE_BLOCKSIZE+Simd::vect_size+3,
                           _mm};
        for (auto m: nrows) {
            Elt_ptr M = fflas_new (F, m, m);
            Elt_ptr Mt = fflas_new (F, m, m);

            finit (F, m, m, M, m);
            for (size_t i = 0; i < m; i++)
                for (size_t j = 0; j < m; j++)
                    F.assign (M[i*m+j], (uint64_t)(i << 8) + j);
            finit (F, m, m, Mt, m);

                /* not inplace full matrix */
            ok &= test_ftranspose<Simd> (m, m, M, m, Mt, m);

                /* not inplace submatrix */
            size_t s = m - FFLAS_TRANSPOSE_BLOCKSIZE;
            ok &= test_ftranspose<Simd> (s, s, M, m, Mt, m);

                /* inplace full matrix */
            ok &= test_ftranspose<Simd> (m, m, M, m, M, m);

                /* inplace submatrix */
            s = m - Simd::vect_size;
            ok &= test_ftranspose<Simd> (s, s, M, m, M, m);

            fflas_delete (M);
            fflas_delete (Mt);
        }

            /* non square matrices */
        size_t ncols[] = { 2*FFLAS_TRANSPOSE_BLOCKSIZE,
                           4*FFLAS_TRANSPOSE_BLOCKSIZE+Simd::vect_size,
                           3*FFLAS_TRANSPOSE_BLOCKSIZE+2*Simd::vect_size+1,
                           _nn};
        for (size_t i = 0; i < 3; i++) {
            size_t m = nrows[i];
            size_t n = ncols[i];

            Elt_ptr M = fflas_new (F, m, n);
            Elt_ptr Mt = fflas_new (F, n, m);

            finit (F, m, n, M, n);
            for (size_t i = 0; i < m; i++)
                for (size_t j = 0; j < n; j++)
                    F.assign (M[i*n+j], (uint64_t)(i << 8) + j);
            finit (F, n, m, Mt, m);

                /* not inplace full matrix */
            ok &= test_ftranspose<Simd> (m, n, M, n, Mt, m);

                /* not inplace submatrix */
            size_t s = m - Simd::vect_size;
            size_t t = n - 2*Simd::vect_size+1;
            ok &= test_ftranspose<Simd> (s, t, M, n, Mt, m);

                /* inplace full matrix */
            ok &= test_ftranspose<Simd> (m, n, M, n, M, m);

            fflas_delete (M);
            fflas_delete (Mt);
        }

            /* print results */
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
        return doTests();
    }

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
    template <typename _E = Elt,
              enable_if_t<is_same<_E, Elt>::value>* = nullptr,
              enable_if_t<Simd<_E>::vect_size != 1>* = nullptr,
              enable_if_simd128_t<_E>* = nullptr>
    bool run () {
        return doTests() & doTests<Simd128<Elt>>();
    }
#endif

#ifdef __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS
    template <typename _E = Elt,
              enable_if_t<is_same<_E, Elt>::value>* = nullptr,
              enable_if_t<Simd<_E>::vect_size != 1>* = nullptr,
              enable_if_simd256_t<_E>* = nullptr>
    bool run () {
        return doTests() & doTests<Simd128<Elt>>()
            & doTests<Simd256<Elt>>();
    }
#endif

#ifdef __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS
    template <typename _E = Elt,
              enable_if_t<is_same<_E, Elt>::value>* = nullptr,
              enable_if_t<Simd<_E>::vect_size != 1>* = nullptr,
              enable_if_simd512_t<_E>* = nullptr>
    bool run () {
        return doTests() & doTests<Simd128<Elt>>()
            & doTests<Simd256<Elt>>() & doTests<Simd512<Elt>>();
    }

        /* Workaround for Simd512<(u)int32/16_t> which does not exist */
    template <typename _E = Elt,
              enable_if_t<is_same<_E, uint32_t>::value
                          || is_same<_E, int32_t>::value
                          || is_same<_E, uint16_t>::value
                          || is_same<_E, int16_t>::value
                          >* = nullptr>
    bool run () {
        return doTests() & doTests<Simd128<Elt>>()
            & doTests<Simd256<Elt>>();
    }
#endif

protected:
    Field F;
    size_t _mm,_nn;
};

/******************************************************************************/
int main(int argc, char** argv)
{
    std::cout << std::setprecision(17);
    std::cerr << std::setprecision(17);

    size_t m = 1000;
    size_t n = 2000;
    
    Argument as[] = {
        { 'm', "-m M", "Set the dimension m",         TYPE_INT , &m },
        { 'n', "-n N", "Set the dimension n",         TYPE_INT , &n },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);
    bool ok = true;

    ok &= Test<float>(m,n).run();
    ok &= Test<double>(m,n).run();
    ok &= Test<uint64_t>(m,n).run();
    ok &= Test<int64_t>(m,n).run();
    ok &= Test<uint32_t>(m,n).run();
    ok &= Test<int32_t>(m,n).run();
    ok &= Test<uint16_t>(m,n).run();
    ok &= Test<int16_t>(m,n).run();
    ok &= Test<Givaro::Integer>(m,n).run();
    ok &= Test<RecInt::rint<6>>(m,n).run();
    ok &= Test<RecInt::ruint<6>>(m,n).run();
    ok &= Test<RecInt::rint<7>>(m,n).run();
    ok &= Test<RecInt::ruint<7>>(m,n).run();
    ok &= Test<RecInt::rint<8>>(m,n).run();
    ok &= Test<RecInt::ruint<8>>(m,n).run();

    return !ok;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
