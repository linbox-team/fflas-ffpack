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
//          Test for the Arm SME kernels of fadd and fsub (fflas_sme.h)
//--------------------------------------------------------------------------

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <givaro/zring.h>

#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/args-parser.h"

// fadd and fsub over ZRing<T> on vectors of every alignment and of lengths around the SVE
// vector, out of place and in place, against a plain loop.  The entries are small integers, so
// the comparison is exact.
template <class T>
static bool test (uint64_t seed, size_t& kernels)
{
    Givaro::ZRing<T> F;
    const size_t maxn = 1200, pad = 64;
    T* a = FFLAS::fflas_new<T> (maxn + pad, Alignment::CACHE_LINE);
    T* b = FFLAS::fflas_new<T> (maxn + pad, Alignment::CACHE_LINE);
    T* c = FFLAS::fflas_new<T> (maxn + pad, Alignment::CACHE_LINE);
    std::vector<T> ref (maxn);
    uint64_t s = seed | 1;
    auto next = [&s] () { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return (T) ((int) ((s >> 20) % 2001) - 1000); };
    const size_t lengths[] = { 0, 1, 15, 16, 17, 63, 64, 65, 128, 496, 500, 512, 1000, 1024, 1200 };
    bool ok = true;
    for (size_t n : lengths)
        for (size_t oc = 0; oc < 17 && ok; oc += (oc < 2 ? 1 : 5))
            for (size_t oa = 0; oa < 17 && ok; oa += 8)
                for (int op = 0; op < 4 && ok; op++) {
                    T *A = a + oa, *B = b + (oa + oc) % 16, *C = c + oc;
                    for (size_t i = 0; i < n; i++) { A[i] = next (); B[i] = next (); C[i] = next (); }
                    const bool add = op % 2 == 0, inplace = op >= 2;
                    const T* X = inplace ? C : A;     // C <- C +/- B when in place
                    for (size_t i = 0; i < n; i++) ref[i] = add ? X[i] + B[i] : X[i] - B[i];
                    if (FFLAS::sme::whole_vectors (C, X, (const T*) B, n)) kernels++;
                    if (add) FFLAS::fadd (F, n, X, 1, B, 1, C, 1);
                    else FFLAS::fsub (F, n, X, 1, B, 1, C, 1);
                    for (size_t i = 0; i < n; i++)
                        if (C[i] != ref[i]) {
                            std::cerr << "FAIL: " << (add ? "fadd" : "fsub") << (inplace ? " in place" : "")
                                      << " n=" << n << " offsets " << oc << ", " << oa << " at " << i << std::endl;
                            ok = false;
                            break;
                        }
                }
    FFLAS::fflas_delete (a);
    FFLAS::fflas_delete (b);
    FFLAS::fflas_delete (c);
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
    size_t kf = 0, kd = 0;
    bool ok = test<float> (seed, kf) && test<double> (seed, kd);
    std::cout << "SME vector length " << FFLAS::sme::vector_bytes () << " bytes"
              << (FFLAS::sme::vector_bytes () ? "" : " (no SME on this processor, or FFLAS_SME=0: plain loops tested)")
              << "; cases that went to the kernels: " << kf << " (float), " << kd << " (double)" << std::endl;
    return ok ? 0 : -1;
}
