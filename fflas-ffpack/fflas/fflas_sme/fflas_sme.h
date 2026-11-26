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

/** @file fflas/fflas_sme/fflas_sme.h
 * @brief Additions and subtractions of float and double vectors in streaming SVE, the
 * vector mode of Arm's Scalable Matrix Extension (SME): 512-bit vectors on Apple's M4.
 *
 * The SIMD layer of the library covers SSE, AVX and AVX-512; on an Arm machine its
 * elementwise loops are plain C++ that the compiler vectorises for NEON at best.  With
 * __FFLASFFPACK_HAVE_SME (configure --enable-sme) fadd and fsub over an unparametric
 * float or double domain (the additions of the Winograd recursion over a delayed field)
 * send a vector to the kernels below when it starts on a vector boundary and holds whole
 * vectors, the case in which they are faster; every other vector takes the plain loop,
 * as before.  The kernels are header-only and need no compiler flag: each carries the
 * target attribute and enters and leaves streaming mode itself.  Whether the processor
 * has SME is checked once at run time; FFLAS_SME=0 in the environment switches the
 * kernels off.
 *
 * Measured on an Apple M4 Pro (clang 21): inside Winograd's recursion, on rows of 432 to
 * 864 floats that are aligned whole vectors, 1.4 to 1.6 times faster than the plain
 * loop; on rows that are not, slower (three to eight times with a masked last vector,
 * 1.2 to 1.4 times with the head and tail in plain C), hence the rule above.  On arrays
 * too large for the cache both are bound by memory and the plain loop is about 10%
 * faster.  The modular addition and the modular reduction were also tried and are half
 * the speed of the plain loops in streaming SVE, so they are not here.  The matrix
 * products are the BLAS's: Apple's Accelerate runs them on the SME units with or without
 * this option.
 */

#ifndef __FFLASFFPACK_fflas_sme_H
#define __FFLASFFPACK_fflas_sme_H

#include <arm_sme.h>
#include <cstddef>
#include <cstdlib>
#if defined(__APPLE__)
#include <sys/sysctl.h>
#elif defined(__linux__)
#include <sys/auxv.h>
#endif

#ifndef __FFLASFFPACK_SME_MIN
#define __FFLASFFPACK_SME_MIN 64   /* vectors shorter than this take the plain loop */
#endif

namespace FFLAS { namespace sme {

#define __FFLASFFPACK_SME_TARGET __attribute__((target("sme")))

    /* the streaming vector length in bytes (64 on the M4); 0 if the processor has no SME or the
     * environment says FFLAS_SME=0 */
    __FFLASFFPACK_SME_TARGET inline size_t query_vector_bytes () { return svcntsb(); }

    inline size_t vector_bytes ()
    {
        static const size_t bytes = [] () -> size_t {
            const char* e = std::getenv ("FFLAS_SME");
            if (e && std::atoi (e) == 0) return 0;
            bool have = false;
#if defined(__APPLE__)
            int v = 0;
            size_t len = sizeof (v);
            have = sysctlbyname ("hw.optional.arm.FEAT_SME", &v, &len, NULL, 0) == 0 && v != 0;
#elif defined(__linux__) && defined(HWCAP2_SME)
            have = (getauxval (AT_HWCAP2) & HWCAP2_SME) != 0;
#endif
            return have ? query_vector_bytes () : 0;
        } ();
        return bytes;
    }

    /* C <- A + B and C <- A - B on n elements, n a multiple of the vector length, all three
     * arrays on a vector boundary (C may be A or B) */
#define __FFLASFFPACK_SME_KERNEL(NAME, T, CNT, BITS, OP)                                              \
    __arm_locally_streaming __FFLASFFPACK_SME_TARGET __attribute__((noinline))                         \
    inline void NAME (T* C, const T* A, const T* B, size_t n)                                          \
    {                                                                                                  \
        const size_t vl = CNT ();                                                                      \
        const svbool_t g = svptrue_b##BITS ();                                                         \
        for (size_t i = 0; i < n; i += vl)                                                             \
            svst1_f##BITS (g, C + i, OP##_f##BITS##_x (g, svld1_f##BITS (g, A + i), svld1_f##BITS (g, B + i))); \
    }

    __FFLASFFPACK_SME_KERNEL (add_kernel, float, svcntw, 32, svadd)
    __FFLASFFPACK_SME_KERNEL (sub_kernel, float, svcntw, 32, svsub)
    __FFLASFFPACK_SME_KERNEL (add_kernel, double, svcntd, 64, svadd)
    __FFLASFFPACK_SME_KERNEL (sub_kernel, double, svcntd, 64, svsub)

#undef __FFLASFFPACK_SME_KERNEL
#undef __FFLASFFPACK_SME_TARGET

    /* the processor has SME, the three arrays start on a vector boundary and the n elements are
     * whole vectors */
    template<class T>
    inline bool whole_vectors (const T* C, const T* A, const T* B, size_t n)
    {
        const size_t v = vector_bytes ();
        return v && n >= __FFLASFFPACK_SME_MIN && (n * sizeof (T)) % v == 0
            && (size_t) C % v == 0 && (size_t) A % v == 0 && (size_t) B % v == 0;
    }

    /* C <- A +/- B by the kernels if the arrays qualify; false if the caller must do it.  Only
     * float and double have kernels. */
    template<class T>
    inline bool addsub (bool, T*, const T*, const T*, size_t) { return false; }

    inline bool addsub (bool add, float* C, const float* A, const float* B, size_t n)
    {
        if (!whole_vectors (C, A, B, n)) return false;
        if (add) add_kernel (C, A, B, n); else sub_kernel (C, A, B, n);
        return true;
    }

    inline bool addsub (bool add, double* C, const double* A, const double* B, size_t n)
    {
        if (!whole_vectors (C, A, B, n)) return false;
        if (add) add_kernel (C, A, B, n); else sub_kernel (C, A, B, n);
        return true;
    }

} // sme
} // FFLAS

#endif // __FFLASFFPACK_fflas_sme_H
