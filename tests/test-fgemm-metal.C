/*
 * Copyright (C) 2026 FFLAS-FFPACK group
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

//--------------------------------------------------------------------------
//     Test for the Metal backend of the single precision fgemm (fflas_metal.h)
//--------------------------------------------------------------------------

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <givaro/modular-balanced.h>
#include <givaro/zring.h>

#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/args-parser.h"

static uint64_t rng;
static float entry (int bound) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return (float) ((int) ((rng >> 20) % (2 * bound + 1)) - bound); }

// C = alpha A B + beta C by FFLAS::metal::sgemm (the GPU and the BLAS in strips) against one
// cblas_sgemm, on integer-valued matrices small enough for both to be exact, with leading
// dimensions larger than the matrices (operands inside larger arrays, as in a recursion).
static bool test_sgemm (size_t m, size_t k, size_t n, float alpha, float beta, bool& taken)
{
    const size_t lda = k + 37, ldb = n + 5, ldc = n + 64;
    std::vector<float> A (m * lda), B (k * ldb), C (m * ldc), R;
    for (auto& x : A) x = entry (7);
    for (auto& x : B) x = entry (7);
    for (auto& x : C) x = entry (7);
    R = C;
    cblas_sgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, (int)m, (int)n, (int)k, alpha, A.data (), (int)lda,
                 B.data (), (int)ldb, beta, R.data (), (int)ldc);
    taken = FFLAS::metal::sgemm (m, n, k, alpha, A.data (), lda, B.data (), ldb, beta, C.data (), ldc);
    if (!taken) return true;        // no GPU here, or the product is below the threshold: nothing to compare
    for (size_t i = 0; i < m * ldc; i++)
        if (C[i] != R[i]) {
            std::cerr << "FAIL: metal::sgemm " << m << "x" << k << "x" << n << " alpha " << alpha << " beta " << beta
                      << " differs at row " << i / ldc << ", column " << i % ldc << std::endl;
            return false;
        }
    return true;
}

// fgemm over a single precision field with the backend in the way, against the same field's product
// applied to a vector (Freivalds): C v = A (B v)
template <class Field>
static bool test_fgemm (const Field& F, size_t m, size_t k, size_t n, int w, int bound)
{
    typedef typename Field::Element E;
    E *A = FFLAS::fflas_new (F, m, k), *B = FFLAS::fflas_new (F, k, n), *C = FFLAS::fflas_new (F, m, n);
    E *v = FFLAS::fflas_new (F, n), *Bv = FFLAS::fflas_new (F, k), *ABv = FFLAS::fflas_new (F, m), *Cv = FFLAS::fflas_new (F, m);
    for (size_t i = 0; i < m * k; i++) F.init (A[i], entry (bound));
    for (size_t i = 0; i < k * n; i++) F.init (B[i], entry (bound));
    for (size_t i = 0; i < n; i++) F.init (v[i], entry (bound));
    FFLAS::MMHelper<Field, FFLAS::MMHelperAlgo::Winograd> H (F, w);
    FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, k, F.one, A, k, B, n, F.zero, C, n, H);
    FFLAS::fgemv (F, FFLAS::FflasNoTrans, k, n, F.one, B, n, v, 1, F.zero, Bv, 1);
    FFLAS::fgemv (F, FFLAS::FflasNoTrans, m, k, F.one, A, k, Bv, 1, F.zero, ABv, 1);
    FFLAS::fgemv (F, FFLAS::FflasNoTrans, m, n, F.one, C, n, v, 1, F.zero, Cv, 1);
    bool ok = FFLAS::fequal (F, m, ABv, 1, Cv, 1);
    if (!ok) F.write (std::cerr << "FAIL: fgemm " << m << "x" << k << "x" << n << ", " << w << " Winograd levels, over ") << std::endl;
    FFLAS::fflas_delete (A, B, C);
    FFLAS::fflas_delete (v, Bv, ABv, Cv);
    return ok;
}

int main (int argc, char** argv)
{
    uint64_t seed = 1234567;
    Argument as[] = {
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };
    FFLAS::parseArguments (argc, argv, as);
    rng = seed | 1;
    // small products, so that the test is quick: the backend reads its threshold at its first call
    setenv ("FFLAS_METAL_MIN", "256", 0);

    bool ok = true, taken = false, any = false;
    const size_t shapes[][3] = { {600, 600, 600}, {433, 700, 512}, {1000, 300, 777}, {300, 1100, 300}, {1299, 257, 259} };
    for (auto& s : shapes) {
        ok = ok && test_sgemm (s[0], s[1], s[2], 1.0f, 0.0f, taken);
        any = any || taken;
        ok = ok && test_sgemm (s[0], s[1], s[2], -2.0f, 3.0f, taken);
    }
    ok = ok && test_sgemm (100, 600, 600, 1.0f, 0.0f, taken);   // below the threshold: must be declined
    if (taken) { std::cerr << "FAIL: a product below the threshold was taken" << std::endl; ok = false; }

    Givaro::ZRing<float> Z;
    Givaro::ModularBalanced<float> F (101);
    Givaro::Modular<float> G (1009);
    for (int w = 0; w <= 1; w++) {
        ok = ok && test_fgemm (Z, 640, 600, 700, w, 3);    // entries small enough for the check to be exact in float
        ok = ok && test_fgemm (F, 1024, 1024, 1024, w, 50);
        ok = ok && test_fgemm (G, 800, 1200, 600, w, 500);
    }

    unsigned long products, gpu, blas;
    fflas_metal_stats (&products, &gpu, &blas);
    std::cout << (any ? "Metal backend: " : "no GPU for the Metal backend here (nothing compared with it): ")
              << products << " products shared, " << gpu << " strips on the GPU, " << blas << " on the BLAS" << std::endl;
    if (any && gpu == 0) { std::cerr << "FAIL: the GPU computed no strip" << std::endl; ok = false; }
    return ok ? 0 : -1;
}
