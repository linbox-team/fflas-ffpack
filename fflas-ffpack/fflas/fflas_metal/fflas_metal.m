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

/* fflas_metal.m: see fflas_metal.h.  Objective-C, built with -fobjc-arc (configure
 * --enable-metal does it); link with the Foundation, Metal and MetalPerformanceShaders
 * frameworks.  It calls no BLAS itself: the caller's strip function does. */
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <MetalPerformanceShaders/MetalPerformanceShaders.h>
#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>

#include "fflas_metal.h"

#define NKERN 8
static id<MTLDevice> g_dev;
static id<MTLCommandQueue> g_cq;
static MPSMatrixMultiplication *g_kern[NKERN];
static size_t g_km[NKERN], g_kk[NKERN], g_kn[NKERN];
static float g_ka[NKERN], g_kb[NKERN];
static int g_nkern, g_knext;
static int g_init;              /* 0 not tried, 1 ready, -1 no GPU or switched off */
static size_t g_min = 2000, g_strips = 0;
static int g_only;
static unsigned long g_products, g_gpu_strips, g_blas_strips;
static pthread_mutex_t g_lock = PTHREAD_MUTEX_INITIALIZER;   /* one product at a time */

static void metal_init(void) {
    const char *e = getenv("FFLAS_METAL");
    g_init = -1;
    if (e && atoi(e) == 0) return;
    if ((e = getenv("FFLAS_METAL_MIN")) && atol(e) > 0) g_min = (size_t)atol(e);
    if ((e = getenv("FFLAS_METAL_STRIPS")) && atol(e) > 0) g_strips = (size_t)atol(e);
    if ((e = getenv("FFLAS_METAL_ONLY"))) g_only = atoi(e) != 0;
    g_dev = MTLCreateSystemDefaultDevice();
    if (g_dev) g_cq = [g_dev newCommandQueue];
    if (g_dev && g_cq) g_init = 1;
}

/* the kernel for strips of rows x n with inner dimension k (cached: a product uses at most two) */
static MPSMatrixMultiplication *kernel(size_t rows, size_t k, size_t n, float alpha, float beta) {
    for (int i = 0; i < g_nkern; i++)
        if (g_km[i] == rows && g_kk[i] == k && g_kn[i] == n && g_ka[i] == alpha && g_kb[i] == beta) return g_kern[i];
    MPSMatrixMultiplication *mm = [[MPSMatrixMultiplication alloc] initWithDevice:g_dev transposeLeft:NO transposeRight:NO
                                                                       resultRows:rows resultColumns:n interiorColumns:k
                                                                            alpha:alpha beta:beta];
    int slot = g_nkern < NKERN ? g_nkern++ : (g_knext++ % NKERN);
    g_kern[slot] = mm;
    g_km[slot] = rows; g_kk[slot] = k; g_kn[slot] = n; g_ka[slot] = alpha; g_kb[slot] = beta;
    return mm;
}

/* the pages that hold [p, p + bytes) as a shared Metal buffer, no copy; *off: p's offset in it */
static id<MTLBuffer> wrap(const void *p, size_t bytes, size_t *off) {
    const uintptr_t page = (uintptr_t)getpagesize();
    uintptr_t lo = (uintptr_t)p & ~(page - 1), hi = ((uintptr_t)p + bytes + page - 1) & ~(page - 1);
    *off = (size_t)((uintptr_t)p - lo);
    return [g_dev newBufferWithBytesNoCopy:(void *)lo length:(NSUInteger)(hi - lo) options:MTLResourceStorageModeShared deallocator:nil];
}

/* what the two engines share during one product (a plain local: the GPU's block holds a pointer) */
typedef struct {
    pthread_mutex_t mu;
    size_t lo, hi;          /* strips lo .. hi-1 are not taken yet */
    size_t *failed, nfailed; /* strips the GPU could not compute */
    unsigned long ngpu;
} share_t;

int fflas_metal_sgemm(size_t m, size_t n, size_t k, float alpha, const float *A, size_t lda, const float *B, size_t ldb,
                      float beta, float *C, size_t ldc, fflas_metal_strip_fn blas_strip, void *ctx) {
    if (g_init < 0 || !blas_strip) return 1;
    if (pthread_mutex_trylock(&g_lock) != 0) return 1;
    int rc = 1;
    @autoreleasepool {
        if (!g_init) metal_init();
        if (g_init == 1 && m >= g_min && n >= g_min && k >= g_min) {
            /* rows * rowBytes of each operand must lie inside its buffer */
            size_t offA, offB, offC;
            id<MTLBuffer> bA = wrap(A, m * lda * sizeof(float), &offA);
            id<MTLBuffer> bB = wrap(B, k * ldb * sizeof(float), &offB);
            id<MTLBuffer> bC = wrap(C, m * ldc * sizeof(float), &offC);
            if (bA && bB && bC) {
                /* few, tall strips: the BLAS is inefficient on thin ones and every strip of the GPU's is
                 * a command buffer; measured on an M4 Pro, strips of at least 432 rows, 16 at most */
                size_t strips = g_strips ? g_strips : m / 432;
                if (!g_strips && strips > 16) strips = 16;
                if (strips > m) strips = m;
                if (strips < 2) strips = g_strips == 1 ? 1 : 2;
                const size_t h = (m + strips - 1) / strips;
                strips = (m + h - 1) / h;
                MPSMatrixMultiplication *kfull = kernel(h, k, n, alpha, beta);
                const size_t hlast = m - (strips - 1) * h;
                MPSMatrixMultiplication *klast = hlast == h ? kfull : kernel(hlast, k, n, alpha, beta);
                MPSMatrixDescriptor *dB = [MPSMatrixDescriptor matrixDescriptorWithRows:k columns:n rowBytes:ldb * 4 dataType:MPSDataTypeFloat32];
                MPSMatrix *xB = [[MPSMatrix alloc] initWithBuffer:bB offset:offB descriptor:dB];
                share_t sh;
                pthread_mutex_init(&sh.mu, NULL);
                sh.lo = 0; sh.hi = strips; sh.nfailed = 0; sh.ngpu = 0;
                sh.failed = calloc(strips, sizeof(size_t));
                share_t *shp = &sh;
                dispatch_group_t grp = dispatch_group_create();
                dispatch_group_async(grp, dispatch_get_global_queue(QOS_CLASS_USER_INTERACTIVE, 0), ^{
                    for (;;) {
                        @autoreleasepool {
                            pthread_mutex_lock(&shp->mu);
                            size_t s = shp->lo < shp->hi ? shp->lo++ : (size_t)-1;
                            pthread_mutex_unlock(&shp->mu);
                            if (s == (size_t)-1) break;
                            size_t rows = s == strips - 1 ? hlast : h;
                            MPSMatrixDescriptor *dA = [MPSMatrixDescriptor matrixDescriptorWithRows:rows columns:k rowBytes:lda * 4 dataType:MPSDataTypeFloat32];
                            MPSMatrixDescriptor *dC = [MPSMatrixDescriptor matrixDescriptorWithRows:rows columns:n rowBytes:ldc * 4 dataType:MPSDataTypeFloat32];
                            MPSMatrix *xA = [[MPSMatrix alloc] initWithBuffer:bA offset:offA + s * h * lda * 4 descriptor:dA];
                            MPSMatrix *xC = [[MPSMatrix alloc] initWithBuffer:bC offset:offC + s * h * ldc * 4 descriptor:dC];
                            id<MTLCommandBuffer> cb = [g_cq commandBuffer];
                            [(rows == h ? kfull : klast) encodeToCommandBuffer:cb leftMatrix:xA rightMatrix:xB resultMatrix:xC];
                            [cb commit];
                            [cb waitUntilCompleted];
                            if (cb.status != MTLCommandBufferStatusCompleted) {     /* the BLAS redoes it below */
                                shp->failed[shp->nfailed++] = s;
                                break;
                            }
                            shp->ngpu++;
                        }
                    }
                });
                unsigned long nblas = 0;
                if (!g_only)
                    for (;;) {
                        pthread_mutex_lock(&sh.mu);
                        size_t s = sh.lo < sh.hi ? --sh.hi : (size_t)-1;
                        pthread_mutex_unlock(&sh.mu);
                        if (s == (size_t)-1) break;
                        size_t rows = s == strips - 1 ? hlast : h;
                        blas_strip(ctx, s * h, rows);
                        nblas++;
                    }
                dispatch_group_wait(grp, DISPATCH_TIME_FOREVER);
                /* a strip the GPU failed on, and those it left when it stopped */
                for (size_t i = 0; i < sh.nfailed; i++) {
                    size_t s = sh.failed[i], rows = s == strips - 1 ? hlast : h;
                    blas_strip(ctx, s * h, rows);
                    nblas++;
                }
                for (size_t s = sh.lo; s < sh.hi; s++) {
                    size_t rows = s == strips - 1 ? hlast : h;
                    blas_strip(ctx, s * h, rows);
                    nblas++;
                }
                free(sh.failed);
                pthread_mutex_destroy(&sh.mu);
                g_products++;
                g_gpu_strips += sh.ngpu;
                g_blas_strips += nblas;
                rc = 0;
            }
        }
    }
    pthread_mutex_unlock(&g_lock);
    return rc;
}

void fflas_metal_stats(unsigned long *products, unsigned long *gpu_strips, unsigned long *blas_strips) {
    if (products) *products = g_products;
    if (gpu_strips) *gpu_strips = g_gpu_strips;
    if (blas_strips) *blas_strips = g_blas_strips;
}
