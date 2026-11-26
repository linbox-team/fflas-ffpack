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

/** @file fflas/fflas_metal/fflas_metal.h
 * @brief The GPU of an Apple-silicon Mac as a second engine for the single precision
 * matrix product (Metal Performance Shaders), shared with the BLAS.
 *
 * Enabled by configure --enable-metal (__FFLASFFPACK_HAVE_METAL), which also builds the
 * one compiled file, fflas_metal.m, into libfflas_metal.a; a program then links
 * -lfflas_metal -framework Foundation -framework Metal -framework MetalPerformanceShaders
 * (fflas-ffpack-config --libs and the pkg-config file say so).
 *
 * The classical fgemm over Givaro::FloatDomain, where every product over a single
 * precision field ends (ZRing<float>, Modular<float>, ModularBalanced<float>, and the
 * fields that are converted to them), offers its product to FFLAS::metal::sgemm before
 * calling cblas_sgemm.  A product whose dimensions are all large is cut into strips of
 * rows of C; the GPU takes strips from the top and the BLAS from the bottom until they
 * meet.  On unified memory the GPU reads and writes the operands where they are: the
 * pages that hold them are wrapped as Metal buffers for the call, nothing is copied.
 * The GPU's product is, like Accelerate's sgemm, a chain of fused multiply-adds, so on
 * the integer-valued matrices of a delayed-reduction product, which no rounding
 * touches, the result is the same.  Metal has no double precision: products over double
 * are not changed.
 *
 * Measured on an Apple M4 Pro over Accelerate, square products: 1.2x faster than the
 * BLAS alone at n = 2000 to 2304, 1.7 to 2x at 3456 and 2.2 to 2.6x from 4608 to 13824.
 */

#ifndef __FFLASFFPACK_fflas_metal_H
#define __FFLASFFPACK_fflas_metal_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** computes rows first_row .. first_row + rows - 1 of the product with the BLAS */
typedef void (*fflas_metal_strip_fn) (void* ctx, size_t first_row, size_t rows);

/** C = alpha A B + beta C in single precision, row-major, A m x k, B k x n, C m x n,
 * shared between the GPU and the caller's BLAS: blas_strip (ctx, i, r) must compute rows
 * i .. i+r-1 of the result, and is called on the calling thread while the GPU works.
 * @return 0 if C was computed; nonzero if the product was not taken (a dimension below
 * the threshold, no GPU, a second call while one is running, or pages that Metal cannot
 * wrap), in which case C is untouched and the caller computes it.
 *
 * Environment: FFLAS_METAL=0 never takes a product; FFLAS_METAL_MIN=d sets the smallest
 * dimension taken (default 2000); FFLAS_METAL_STRIPS=s the number of strips (default:
 * strips of at least 432 rows, 16 at most); FFLAS_METAL_ONLY=1 gives every strip to the
 * GPU.
 */
int fflas_metal_sgemm (size_t m, size_t n, size_t k, float alpha,
                       const float* A, size_t lda, const float* B, size_t ldb,
                       float beta, float* C, size_t ldc,
                       fflas_metal_strip_fn blas_strip, void* ctx);

/** the number of products taken, and of strips computed by the GPU and by the BLAS, since
 * the start of the program (any pointer may be NULL) */
void fflas_metal_stats (unsigned long* products, unsigned long* gpu_strips, unsigned long* blas_strips);

#ifdef __cplusplus
}

#include "fflas-ffpack/config-blas.h"

namespace FFLAS { namespace metal {

    /** C = alpha A B + beta C shared between the GPU and cblas_sgemm; false if the product
     * was not taken and the caller must call the BLAS itself. */
    inline bool sgemm (const size_t m, const size_t n, const size_t k, const float alpha,
                       const float* A, const size_t lda, const float* B, const size_t ldb,
                       const float beta, float* C, const size_t ldc)
    {
        struct Product {
            size_t n, k, lda, ldb, ldc;
            float alpha, beta;
            const float *A, *B;
            float *C;
        } P = { n, k, lda, ldb, ldc, alpha, beta, A, B, C };
        return fflas_metal_sgemm (m, n, k, alpha, A, lda, B, ldb, beta, C, ldc,
                                  [] (void* ctx, size_t i, size_t rows) {
                                      const Product* p = static_cast<const Product*> (ctx);
                                      cblas_sgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
                                                   (int)rows, (int)p->n, (int)p->k, p->alpha,
                                                   p->A + i * p->lda, (int)p->lda, p->B, (int)p->ldb,
                                                   p->beta, p->C + i * p->ldc, (int)p->ldc);
                                  }, &P) == 0;
    }

} // metal
} // FFLAS
#endif

#endif // __FFLASFFPACK_fflas_metal_H
