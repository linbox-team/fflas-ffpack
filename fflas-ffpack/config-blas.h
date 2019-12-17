/* config-blas.h
 * Copyright (C) 2005  Pascal Giorgi
 *               2007  Clement Pernet
 * Written by Pascal Giorgi <pgiorgi@uwaterloo.ca>
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
 *
 */


#ifndef __FFLASFFPACK_config_blas_H
#define __FFLASFFPACK_config_blas_H

// #include "fflas-ffpack/utils/fflas_memory.h"
// #ifndef __FFLASFFPACK_CONFIGURATION
// #include "fflas-ffpack/fflas-ffpack-config.h"
// #endif

// #ifdef OPTIMISATION_MODE
// #include "fflas-ffpack/config.h"
// #endif

#ifdef HAVE_MKL
#define __FFLASFFPACK_HAVE_MKL
#endif

#ifdef __FFLASFFPACK_HAVE_MKL
#include <mkl.h>

#endif


#ifndef CBLAS_INT
#ifdef blasint /*  openblas */
#define CBLAS_INT blasint
#elif defined( MKL_INT )
#define CBLAS_INT MKL_INT
#else
#define CBLAS_INT int
#endif /* blasint */
#endif /*  CBLAS_INT */

#ifdef CUDA_BLAS

#define sgemv_ cublas_sgemv
#define sgemm_ cublas_sgemm
#define strsm_ cublas_strsm
#define strmm_ cublas_strmm

#endif // CUDA_BLAS


#ifndef __FFLASFFPACK_HAVE_MKL

#define CBLAS_ENUM_DEFINED_H
enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, AtlasConj=114};
enum CBLAS_UPLO  {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG  {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE  {CblasLeft=141, CblasRight=142};

// #define CBLAS_INDEX int


#ifndef __FFLASFFPACK_HAVE_CBLAS

// CBLAS are not available define our own wrapper

// define external link to BLAS function
extern "C" {

#define CBLAS_EXTERNALS
    static const char* EXT_BLAS_TRANSPOSE    (CBLAS_TRANSPOSE t) { if (t == CblasNoTrans) return "N"; else if (t == CblasTrans) return "T"; else return "";}
    static const char* EXT_BLAS_TRANSPOSE_tr (CBLAS_TRANSPOSE t) { if (t == CblasNoTrans) return "T"; else if (t == CblasTrans) return "N"; else return "";}

    static const char* EXT_BLAS_UPLO         (CBLAS_UPLO t)      { if (t == CblasUpper) return "U"; else return "L";}
    static const char* EXT_BLAS_UPLO_tr      (CBLAS_UPLO t)      { if (t == CblasUpper) return "L"; else return "U";}

    static const char* EXT_BLAS_DIAG         (CBLAS_DIAG t)      { if (t == CblasUnit)  return "U"; else return "N";}

    static const char* EXT_BLAS_SIDE         (CBLAS_SIDE t)      { if (t == CblasLeft)  return "L"; else return "R";}
    static const char* EXT_BLAS_SIDE_tr      (CBLAS_SIDE t)      { if (t == CblasLeft)  return "R"; else return "L";}


    // level 1 routines
    void   daxpy_   (const int*, const double*, const double*, const int*, double*, const int*);
    void   saxpy_   (const int*, const float*, const float*, const int*, float*, const int*);
    double ddot_    (const int*, const double*, const int*, const double*, const int*);
    float  sdot_    (const int*, const float*, const int*, const float*, const int*);
    double dasum_   (const int*, const double*, const int*);
    int    idamax_  (const int*, const double*, const int*);
    double dnrm2_   (const int*, const double*, const int*);

    // level 2 routines
    void dgemv_ (const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
    void sgemv_ (const char*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
    void dger_  (const int*, const int*, const double*, const double*, const int*, const double*, const int*, double*, const int*);
    void sger_  (const int*, const int*, const float*, const float*, const int*, const float*, const int*, float*, const int*);

    void dcopy_  (const int *, const double *, const int *, double *, const int *);
    void scopy_  (const int *, const float  *, const int *, float  *, const int *);

    void dscal_  (const int *, const double *, double *, const int *);
    void sscal_  (const int *, const float  *, float  *, const int *);



    // level 3 routines
    void dtrsm_ (const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);
    void strsm_ (const char*, const char*, const char*, const char*, const int*, const int*, const float*, const float*, const int*, float*, const int*);
    void dtrmm_ (const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);
    void strmm_ (const char*, const char*, const char*, const char*, const int*, const int*, const float*, const float*, const int*, float*, const int*);
    void sgemm_ (const char*, const char*, const int*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
    void dgemm_ (const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
}

// define C wrappers
extern "C" {


    // level 1 routines

    static inline void cblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY)
    {
        daxpy_ (&N,&alpha, X, &incX, Y, &incY);
    }

    static inline void cblas_saxpy(const int N, const float alpha, const float *X, const int incX, float *Y, const int incY)
    {
        saxpy_ (&N,&alpha, X, &incX, Y, &incY);
    }

    static inline double cblas_ddot(const int N, const double *X, const int incX, const double *Y, const int incY)
    {
        return ddot_ (&N, X, &incX, Y, &incY);
    }

    static inline float cblas_sdot(const int N, const float *X, const int incX, const float *Y, const int incY)
    {
        return sdot_ (&N, X, &incX, Y, &incY);
    }


    static inline double cblas_dasum(const int N, const double *X, const int incX){
        return dasum_ (&N, X, &incX);
    }

    static inline int cblas_idamax(const int N, const double *X, const int incX){
        return idamax_ (&N, X, &incX);
    }

    static inline double cblas_dnrm2(const int N, const double *X, const int incX){
        return dnrm2_(&N, X, &incX);
    }


    // level 2 routines

    static inline void cblas_dgemv(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const double alpha,
                                   const double *A, const int lda, const double *X, const int incX, const double beta, double *Y, const int incY)
    {
        if (Order == CblasRowMajor)
            dgemv_ ( EXT_BLAS_TRANSPOSE_tr(TransA), &N, &M, &alpha, A, &lda, X, &incX, &beta, Y, &incY);
        else
            dgemv_ ( EXT_BLAS_TRANSPOSE(TransA), &M, &N, &alpha, A, &lda, X, &incX, &beta, Y, &incY);
    }
    static inline void cblas_sgemv(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const float alpha,
                                   const float *A, const int lda, const float *X, const int incX, const float beta, float *Y, const int incY)
    {
        if (Order == CblasRowMajor)
            sgemv_ ( EXT_BLAS_TRANSPOSE_tr(TransA), &N, &M, &alpha, A, &lda, X, &incX, &beta, Y, &incY);
        else
            sgemv_ ( EXT_BLAS_TRANSPOSE(TransA), &M, &N, &alpha, A, &lda, X, &incX, &beta, Y, &incY);
    }

    static inline void cblas_dger(const enum CBLAS_ORDER Order, const int M, const int N, const double alpha, const double *X, const int incX,
                                  const double *Y, const int incY, double *A, const int lda)
    {
        if (Order == CblasRowMajor)
            dger_ (&N, &M, &alpha, Y, &incY, X, &incX, A, &lda);
        else
            dger_ (&M, &N, &alpha, X, &incX, Y, &incY, A, &lda);
    }

    static inline void cblas_sger(const enum CBLAS_ORDER Order, const int M, const int N, const float alpha, const float *X, const int incX,
                                  const float *Y, const int incY, float *A, const int lda)
    {
        if (Order == CblasRowMajor)
            sger_ (&N, &M, &alpha, Y, &incY, X, &incX, A, &lda);
        else
            sger_ (&M, &N, &alpha, X, &incX, Y, &incY, A, &lda);
    }

    static inline void cblas_dcopy(const int N, const double *X, const int incX, double *Y, const int incY)
    {
        dcopy_(&N,X,&incX,Y,&incY);
    }


    static inline void cblas_scopy(const int N, const float *X, const int incX, float *Y, const int incY)
    {
        scopy_(&N,X,&incX,Y,&incY);
    }

    static inline void cblas_dscal(const int N, const double alpha,  double *Y, const int incY)
    {
        dscal_(&N,&alpha,Y,&incY);
    }

    static inline void cblas_sscal(const int N, const float alpha,  float *Y, const int incY)
    {
        sscal_(&N,&alpha,Y,&incY);
    }


    // level 3 routines

    static inline void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                                   const enum CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double *A, const int lda,
                                   double *B, const int ldb)
    {
        if (Order == CblasRowMajor)
            dtrsm_ ( EXT_BLAS_SIDE_tr(Side), EXT_BLAS_UPLO_tr(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &N, &M, &alpha, A, &lda, B, &ldb);
        else
            dtrsm_ ( EXT_BLAS_SIDE(Side), EXT_BLAS_UPLO(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &M, &N, &alpha, A, &lda, B, &ldb);
    }
    static inline void cblas_strsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                                   const enum CBLAS_DIAG Diag, const int M, const int N, const float alpha, const float *A, const int lda,
                                   float *B, const int ldb)
    {
        if (Order == CblasRowMajor)
            strsm_ ( EXT_BLAS_SIDE_tr(Side), EXT_BLAS_UPLO_tr(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &N, &M, &alpha, A, &lda, B, &ldb);
        else
            strsm_ ( EXT_BLAS_SIDE(Side), EXT_BLAS_UPLO(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &M, &N, &alpha, A, &lda, B, &ldb);
    }

    static inline void cblas_dtrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                                   const enum CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double *A, const int lda,
                                   double *B, const int ldb)
    {
        if (Order == CblasRowMajor)
            dtrmm_ ( EXT_BLAS_SIDE_tr(Side), EXT_BLAS_UPLO_tr(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &N, &M, &alpha, A, &lda, B, &ldb);
        else
            dtrmm_ ( EXT_BLAS_SIDE(Side), EXT_BLAS_UPLO(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &M, &N, &alpha, A, &lda, B, &ldb);
    }
    static inline void cblas_strmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                                   const enum CBLAS_DIAG Diag, const int M, const int N, const float alpha, const float *A, const int lda,
                                   float *B, const int ldb)
    {
        if (Order == CblasRowMajor)
            strmm_ ( EXT_BLAS_SIDE_tr(Side), EXT_BLAS_UPLO_tr(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &N, &M, &alpha, A, &lda, B, &ldb);
        else
            strmm_ ( EXT_BLAS_SIDE(Side), EXT_BLAS_UPLO(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &M, &N, &alpha, A, &lda, B, &ldb);
    }

    static inline void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                                   const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb,
                                   const double beta, double *C, const int ldc)
    {
        if (Order == CblasRowMajor)
            dgemm_ ( EXT_BLAS_TRANSPOSE(TransB), EXT_BLAS_TRANSPOSE(TransA), &N, &M, &K, &alpha, B, &ldb, A, &lda, &beta, C, &ldc);
        else
            dgemm_ ( EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_TRANSPOSE(TransB), &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
    }
    static inline void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                                   const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb,
                                   const float beta, float *C, const int ldc)
    {
        if (Order == CblasRowMajor)
            sgemm_ ( EXT_BLAS_TRANSPOSE(TransB), EXT_BLAS_TRANSPOSE(TransA), &N, &M, &K, &alpha, B, &ldb, A, &lda, &beta, C, &ldc);
        else
            sgemm_ ( EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_TRANSPOSE(TransB), &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
    }

    static inline cblas_ssyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                              const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                              const float alpha, const float *A, const int lda,
                              const float beta, float *C, const int ldc){
       if (Order == CblasRowMajor)
           ssryk_ (EXT_BLAS_UPLO_tr(Uplo), EXT_BLAS_TRANSPOSE(Trans), N, K, alpha, A, lda, beta, C, ldc); // @TODO check this
       else
           ssryk_ (EXT_BLAS_UPLO (Uplo), EXT_BLAS_TRANSPOSE(Trans), N, K, alpha, A, lda, beta, C, ldc); 
    }
    void cblas_dsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const double *A, const int lda,
                     const double beta, double *C, const int ldc){
        if (Order == CblasRowMajor)
            dsryk_ (EXT_BLAS_UPLO_tr(Uplo), EXT_BLAS_TRANSPOSE(Trans), N, K, alpha, A, lda, beta, C, ldc); // @TODO check this
        else
            dsryk_ (EXT_BLAS_UPLO (Uplo), EXT_BLAS_TRANSPOSE(Trans), N, K, alpha, A, lda, beta, C, ldc); 
    }

}

#else // CBLAS PRESENT
extern "C" {

#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    void openblas_set_num_threads(int num_threads);
#endif

    int cblas_errprn(int ierr, int info, char *form, ...);

    // level 1 routines

    void   cblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY);
    void   cblas_saxpy(const int N, const float alpha, const float *X, const int incX, float *Y, const int incY);

    double cblas_ddot(const int N, const double *X, const int incX, const double *Y, const int incY);
    float  cblas_sdot(const int N, const float *X, const int incX, const float *Y, const int incY);

    double cblas_dasum(const int N, const double *X, const int incX);

    int    cblas_idamax(const int N, const double *X, const int incX);

    double cblas_dnrm2(const int N, const double *X, const int incX);


    // level 2 routines

    void cblas_dgemv(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const double alpha,
                     const double *A, const int lda, const double *X, const int incX, const double beta, double *Y, const int incY);

    void cblas_sgemv(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const float alpha,
                     const float *A, const int lda, const float *X, const int incX, const float beta, float *Y, const int incY);

    void cblas_dger(const enum CBLAS_ORDER Order, const int M, const int N, const double alpha, const double *X, const int incX,
                    const double *Y, const int incY, double *A, const int lda);

    void cblas_sger(const enum CBLAS_ORDER Order, const int M, const int N, const float alpha, const float *X, const int incX,
                    const float *Y, const int incY, float *A, const int lda);

    void cblas_dcopy(const int N, const double *X, const int incX,
                     double *Y, const int incY);

    void cblas_scopy(const int N, const float *X, const int incX,
                     float *Y, const int incY);

    void cblas_dscal(const int N, const double alpha,
                     double *Y, const int incY);

    void cblas_sscal(const int N, const float alpha,
                     float *Y, const int incY);


    // level 3 routines

    void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                     const enum CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double *A, const int lda,
                     double *B, const int ldb);

    void cblas_strsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                     const enum CBLAS_DIAG Diag, const int M, const int N, const float alpha, const float *A, const int lda,
                     float *B, const int ldb);

    void cblas_dtrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                     const enum CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double *A, const int lda,
                     double *B, const int ldb);

    void cblas_strmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                     const enum CBLAS_DIAG Diag, const int M, const int N, const float alpha, const float *A, const int lda,
                     float *B, const int ldb);

    void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                     const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb,
                     const double beta, double *C, const int ldc) ;
    void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                     const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb,
                     const float beta, float *C, const int ldc) ;
    void cblas_ssyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                     const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                     const float alpha, const float *A, const int lda,
                     const float beta, float *C, const int ldc);
    void cblas_dsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const double *A, const int lda,
                 const double beta, double *C, const int ldc);
}
#endif // CBLAS ?

#endif // __FFLASFFPACK_HAVE_MKL

#ifdef __FFLASFFPACK_HAVE_MKL
#define blas_enum
#else
#define blas_enum enum
#endif

#ifdef __FFLASFFPACK_HAVE_LAPACK

#ifndef __FFLASFFPACK_HAVE_CLAPACK

#ifndef CBLAS_EXTERNALS
#define CBLAS_EXTERNALS
// static const char* EXT_BLAS_TRANSPOSE    (CBLAS_TRANSPOSE t) { if (t == CblasNoTrans) return "N"; else if (t == CblasTrans) return "T"; else return "";}
// static const char* EXT_BLAS_TRANSPOSE_tr (CBLAS_TRANSPOSE t) { if (t == CblasNoTrans) return "T"; else if (t == CblasTrans) return "N"; else return "";}

static const char* EXT_BLAS_UPLO         (CBLAS_UPLO t)      { if (t == CblasUpper) return "U"; else return "L";}
static const char* EXT_BLAS_UPLO_tr      (CBLAS_UPLO t)      { if (t == CblasUpper) return "L"; else return "U";}

static const char* EXT_BLAS_DIAG         (CBLAS_DIAG t)      { if (t == CblasUnit)  return "U"; else return "N";}

// static const char* EXT_BLAS_SIDE         (CBLAS_SIDE t)      { if (t == CblasLeft)  return "L"; else return "R";}
// static const char* EXT_BLAS_SIDE_tr      (CBLAS_SIDE t)      { if (t == CblasLeft)  return "R"; else return "L";}
#endif // CBLAS_EXTERNALS


// define external link to LAPACK routines
extern "C" {
    //!@bug we should also allow lapacke from MLK
    void dgetrf_ (const CBLAS_INT *, const CBLAS_INT *, double *, const CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
    void dgetri_ (const CBLAS_INT *, double *, const CBLAS_INT *, const CBLAS_INT *, double *, const CBLAS_INT *, CBLAS_INT *);
    void dtrtri_ (const char *, const char *, const CBLAS_INT *, double *, const CBLAS_INT *, CBLAS_INT *);
    void dswap_ (const CBLAS_INT *, double *, const CBLAS_INT *, double *, const CBLAS_INT *);
}

extern "C" {
    //!@bug we should check if lapacke is present
    int LAPACKE_dsytrf( int matrix_layout, char uplo, int n, double* a, int lda, int* ipiv );
    int LAPACKE_dsytrf_rook( int matrix_layout, char uplo, int n, double* a, int lda, int* ipiv );
    int LAPACKE_dsytrf_rk( int matrix_layout, char uplo, int n, double* a, int lda, double* e, int* ipiv );
    int LAPACKE_dsytrf_aa( int matrix_layout, char uplo, int n, double* a, int lda, int* ipiv );
}


// define C wrappers
extern "C" {
    // LAPACK routines

    // return A=P.L.U (L unitary) with ColMajor
    // return A=L.U.P (U unitary) with RowMajor
    //! @bug Order is not used. we should use ATLAS/interfaces/lapack/C/src/clapack_dgetrf.c or similar
    inline CBLAS_INT clapack_dgetrf(const blas_enum CBLAS_ORDER, const CBLAS_INT M, const CBLAS_INT N,
                                    double *A, const CBLAS_INT lda, CBLAS_INT *ipiv)
    {
        CBLAS_INT info;
        dgetrf_ ( &M, &N, A, &lda, ipiv, &info);
        return info;
    }

    inline CBLAS_INT clapack_dgetri(const blas_enum CBLAS_ORDER, const CBLAS_INT N, double *A,
                                    const CBLAS_INT lda, const CBLAS_INT *ipiv)
    {
        CBLAS_INT info;
        double *work;

#ifndef __FFLASFFPACK_AUTOIMPLEMENT_DGETRI
        // the optimum size of work can be determCBLAS_INTed via the
        // Lapack function ilaenv.
        work= new double[N];
        dgetri_ (&N, A, &lda, ipiv, work, &N,  &info);
        delete[] work;
#else
        work= new double[N*N];
        dtrtri_("U","N", &N, A, &lda, &info);
        if (info > 0)
            return 0;

        for (CBLAS_INT i=0;i<N;++i){
            for(CBLAS_INT j=i;j<N;++j){
                work[i*N+j]=A[i*N+j];
                if (j>i) A[i*N+j]=0.0;
            }
            work[i*N+i]=1.;
        }

        double cst=1.;
        dtrsm_ ("R", "L", "N", "U", &N, &N, &cst, work, &N, A, &N);

        CBLAS_INT ip;
        const CBLAS_INT incr=1;
        for (CBLAS_INT i=0; i<N; ++i){
            ip = ipiv[i]-1;
            if (ip != i)
                dswap_ (&N, &A[i*lda],&incr , &A[ip*lda], &incr);
        }

        delete[] work;
#endif
        return info;
    }

    inline CBLAS_INT clapack_dtrtri(const blas_enum CBLAS_ORDER Order,const blas_enum CBLAS_UPLO Uplo,
                                    const blas_enum CBLAS_DIAG Diag,const CBLAS_INT N, double *A, const CBLAS_INT lda)
    {
        CBLAS_INT info;
        if (Order == CblasRowMajor)
            dtrtri_ (EXT_BLAS_UPLO_tr(Uplo), EXT_BLAS_DIAG(Diag), &N, A, &lda, &info);
        else
            dtrtri_ (EXT_BLAS_UPLO(Uplo), EXT_BLAS_DIAG(Diag), &N, A, &lda, &info);

        return info;
    }


}
#else // CLAPACK PRESENT

// define external link to CLAPACK routines
extern "C" {
    // LAPACK routines

    CBLAS_INT clapack_dgetrf(const blas_enum CBLAS_ORDER Order, const CBLAS_INT M, const CBLAS_INT N,
                             double *A, const CBLAS_INT lda, CBLAS_INT *ipiv);
    CBLAS_INT clapack_dgetri(const blas_enum CBLAS_ORDER Order, const CBLAS_INT N, double *A,
                             const CBLAS_INT lda, const CBLAS_INT *ipiv);
    CBLAS_INT clapack_dtrtri(const blas_enum CBLAS_ORDER Order,const blas_enum CBLAS_UPLO Uplo,
                             const blas_enum CBLAS_DIAG Diag,const CBLAS_INT N, double *A, const CBLAS_INT lda);

}
#endif // CLAPACK ?

#endif // LAPACK ?


#endif //__FFLASFFPACK_config_blas_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
