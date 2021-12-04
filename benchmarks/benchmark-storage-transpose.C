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
class Bench {
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
        Bench (size_t m, size_t n, size_t iters, bool inplace)
            : F(cardinality()), m(m), n(n), iters(iters), inplace(inplace) {
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
            /* Bench Givaro::Integer with a large (=more than a word) modulus */
            return Givaro::Integer ("0x100000000000000000000000000000000");
        }

        /* main test function */
        template <typename Simd = NoSimd<Elt>,
                  enable_if_t<is_same_element<Simd>::value>* = nullptr>
        void doBenchs () {
            Elt_ptr A = nullptr, B = nullptr;
            const size_t lda = n;
            const size_t ldb = m;

            if (inplace) {
                A = B = fflas_new (F, m, lda);
            }
            else {
                A = fflas_new (F, m, lda);
                B = fflas_new (F, n, ldb);
                finit (F, n, m, B, ldb);
            }

            finit (F, m, n, A, lda);
            for (size_t i = 0; i < m; i++)
                for (size_t j = 0; j < n; j++)
		  F.assign (A[i*lda+j], (uint64_t) (i << 16) + j);

            Givaro::Timer chrono;
            chrono.clear(); chrono.start();
            for (size_t i = 0; i < iters; i++) {
                ftranspose<Field, Simd> (F, m, n, A, lda, B, ldb);
            }
            chrono.stop();
            double time = chrono.usertime();

            fflas_delete (A);
            if (!inplace) {
                fflas_delete (B);
            }

            /* print results */
            /* mul by 2: 1 read + 1 write per element of the matrix */
            double throughput = (2*iters*m*n*sizeof(Elt)) / (time * (1 << 30));
            std::cout << F.type_string()
                      << string (27-F.type_string().size(), '.')
                      << Simd::type_string()
                      << string (26-Simd::type_string().size(), '.')
                      << (inplace ? "inplace " : "....... ")
                      << setw (7) << fixed << setprecision(2)
                      << 1000*time / ((double) iters) << "ms "
                      << setw (5) << fixed << setprecision(2)
                      << throughput << "GB/s" << endl;
        }

        /* run tests: call doBenchs for best Simd structs or all available Simd
         * structs, depending on parameter.
         */
        template <typename _E = Elt,
                  enable_if_t<is_same<_E, Elt>::value>* = nullptr,
                  enable_if_no_simd_t<_E>* = nullptr>
        void run (bool allsimd) {
            doBenchs ();
        }

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
        template <typename _E = Elt,
                  enable_if_t<is_same<_E, Elt>::value>* = nullptr,
                  enable_if_t<Simd<_E>::vect_size != 1>* = nullptr,
                  enable_if_simd128_t<_E>* = nullptr>
        void run (bool allsimd) {
            if (allsimd) {
                doBenchs ();
            }
            doBenchs<Simd128<Elt>> ();
        }
#endif

#ifdef __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS
        template <typename _E = Elt,
                  enable_if_t<is_same<_E, Elt>::value>* = nullptr,
                  enable_if_t<Simd<_E>::vect_size != 1>* = nullptr,
                  enable_if_simd256_t<_E>* = nullptr>
        void run (bool allsimd) {
            if (allsimd) {
                doBenchs ();
                doBenchs<Simd128<Elt>> ();
            }
            doBenchs<Simd256<Elt>> ();
        }
#endif

#ifdef __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS
        template <typename _E = Elt,
                  enable_if_t<is_same<_E, Elt>::value>* = nullptr,
                  enable_if_t<Simd<_E>::vect_size != 1>* = nullptr,
                  enable_if_simd512_t<_E>* = nullptr>
        void run (bool allsimd) {
            if (allsimd) {
                doBenchs ();
                doBenchs<Simd128<Elt>> ();
                doBenchs<Simd256<Elt>> ();
            }
            doBenchs<Simd512<Elt>> ();
        }

        /* Workaround for Simd512<(u)int32/16_t> which does not exist */
        template <typename _E = Elt,
                  enable_if_t<is_same<_E, uint32_t>::value
                              || is_same<_E, int32_t>::value
                              || is_same<_E, uint16_t>::value
                              || is_same<_E, int16_t>::value
                              >* = nullptr>
        void run (bool allsimd) {
            if (allsimd) {
                doBenchs ();
                doBenchs<Simd128<Elt>> ();
            }
            doBenchs<Simd256<Elt>> ();
        }
#endif

    protected:
        Field F;
        const size_t m;
        const size_t n;
        const size_t iters;
        const bool inplace;
};

/******************************************************************************/
int main(int argc, char** argv)
{
    size_t iters = 20;
    size_t m = 6144;
    size_t n = 4096;
    bool allsimd = false;
    Argument as[] = {
        { 'm', "-m M", "Set the row dimension of the matrix.",
                                                                TYPE_INT, &m },
        { 'n', "-n N", "Set the column dimension of the matrix.",
                                                                TYPE_INT, &n },
        { 'i', "-i R", "Set number of repetitions.", TYPE_INT , &iters },
        { 'a', "-a Y/N", "benchmarks all Simd structs (default: no).",
                                                        TYPE_BOOL, &allsimd },

        END_OF_ARGUMENTS
    };

    /* parse command-line */
    FFLAS::parseArguments (argc, argv, as);
    cout << "# To rerun this bench: " << argv[0] << " -m " << m << " -n " << n
         << " -i " << iters << (allsimd ? " -a" : "") <<endl;

    bool inplace = true;
    do {
        inplace = not inplace;
        cout << endl << "Benchmarking transpose on matrices of size " << m
             << "x" << n << (inplace ? " [inplace variant]" : "") << ":"
             << endl;

        Bench<double>(m,n,iters,inplace).run(allsimd);
        Bench<uint64_t>(m,n,iters,inplace).run(allsimd);
        Bench<float>(m,n,iters,inplace).run(allsimd);
        Bench<uint32_t>(m,n,iters,inplace).run(allsimd);
        Bench<uint16_t>(m,n,iters,inplace).run(allsimd);

        /* Bench<Givaro::Integer>(m,n,iters,inplace).run(allsimd); */
        /* Bench<RecInt::ruint<6>>(m,n,iters,inplace).run(allsimd); */
        /* Bench<RecInt::ruint<7>>(m,n,iters,inplace).run(allsimd); */
        /* Bench<RecInt::ruint<8>>(m,n,iters,inplace).run(allsimd); */
    }
    while (not inplace) ;

    return 0;
}
