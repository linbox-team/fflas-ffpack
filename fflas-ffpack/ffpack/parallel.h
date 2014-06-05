/* fflas/parallel.h
 * Copyright (C) 2013 Jean Guillaume Dumas Clement Pernet Ziad Sultan
 *
 * Written by Jean Guillaume Dumas Clement Pernet Ziad Sultan
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

#ifndef __FFLASFFPACK_fflas_parallel_H
#define __FFLASFFPACK_fflas_parallel_H




#ifdef __FFLASFFPACK_USE_OPENMP
#include "omp.h"
#endif

#ifdef __FFLASFFPACK_USE_KAAPI
#include "kaapi++"
//#include "kaapi-routines.h"
#endif


//OMP
#define COMMA ,
#define READ(Args...) COMMA Args
#define WRITE(Args...) COMMA Args
#define READWRITE(Args...) COMMA Args
#define NOREAD(void) 
#define NOWRITE(void) 
#define NOREADWRITE(void) 


// KAAPI

//#define DSLREADARRAY(A,pointer, range) A,pointer,range 

#ifdef __FFLASFFPACK_USE_OPENMP
#define GLOBALSHARED(a, Args...) shared(Args)
#define PRAGMA_OMP_TASK_IMPL( F ) _Pragma( #F )
#endif


#ifdef __FFLASFFPACK_USE_OPENMP
    #define TASK(r, w, rw, f, Args...) \
    PRAGMA_OMP_TASK_IMPL( omp task GLOBALSHARED(x  r w rw) )	\
    f(Args)
#else 
   #ifdef __FFLASFFPACK_USE_KAAPI
      #define TASK(r, w, rw, f,Args...)  ka::Spawn<Task ## f <Field> >()(Args)
   #else
      #define TASK(r,w,rw,f,Args...) f(Args)
   #endif
#endif

#define STR_IMPL(a) #a
#define STR(a) STR_IMPL(a)

// Macro PLUQ
#ifdef __FFLASFFPACK_USE_OPENMP 
#define PPLUQ(Rank, F, Diag, M, N, A, lda, P, Q) do{			\
    Rank = pPLUQ(F, Diag, M, N, A, lda, P, Q);			\
}while(0)
#else 
   #ifdef __FFLASFFPACK_USE_KAAPI
   #define PPLUQ(Rank, F, Diag, M, N, A, lda, P, Q) do{			\
    Rank = PPLUQ(F, Diag, M, N, A, lda, P, Q);				\
   }while(0)
   #else
   #define PPLUQ(Rank, F, Diag, M, N, A, lda, P, Q) do{			\
    Rank = PPLUQ(F, Diag, M, N, A, lda, P, Q);			\
   }while(0)
   #endif
#endif 

// Macro pragma omp taskwait
#ifdef __FFLASFFPACK_USE_OPENMP 
   #define WAIT PRAGMA_OMP_TASK_IMPL( omp taskwait )
#else
   #ifdef __FFLASFFPACK_USE_KAAPI
   #define WAIT do{	\
   }while(0)
   #else
   #define WAIT do{	\
   }while(0)
   #endif
#endif 

// Macro barrier
#ifdef __FFLASFFPACK_USE_OPENMP 
   #define BARRIER do{\
  }while(0)
#else
#ifdef __FFLASFFPACK_USE_KAAPI
   #define BARRIER do{	\
     ka::Sync(); \
  }while(0)
#else
   #define BARRIER do{	\
  }while(0)
#endif
#endif 

// Get number of threads
#ifdef __FFLASFFPACK_USE_OPENMP
# define HPAC_NUM_THREADS omp_get_num_threads()
#else
# ifdef __FFLASFFPACK_USE_KAAPI
#  define HPAC_NUM_THREADS kaapi_getconcurrency_cpu()
# else
#  define HPAC_NUM_THREADS 1
# endif
#endif


#endif
