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
#endif



//OMP
#define COMMA ,
#define READ(Args...) COMMA Args
#define WRITE(Args...) COMMA Args
#define READWRITE(Args...) COMMA Args
#define NOREAD(void) 
#define NOWRITE(void) 
#define NOREADWRITE(void) 
#ifdef __FFLASFFPACK_USE_OPENMP
#define GLOBALSHARED(a, Args...) shared(Args)
#define PRAGMA_OMP_TASK_IMPL( F ) _Pragma( #F )
#endif


// KAAPI
#define xstr(s) str(s)
#define str(s) #s
#define foo 4
#define SPAWN(f,N) CONC(ka::Spawn<Task ## f, N)
#define CONC(f, N) f ## N

// Macro computes number of Arguments
//#define NUMARGS(...)  (sizeof((string[]){__VA_ARGS__})/sizeof(string))
#define NUMARGS(...)				\
  PP_NARG_(__VA_ARGS__,PP_RSEQ_N())
#define PP_NARG_(...)				\
  PP_ARG_N(__VA_ARGS__)
#define PP_ARG_N(						\
		 _1, _2, _3, _4, _5, _6, _7, _8, _9,_10,	\
		 _11,_12,_13,_14,_15,_16,_17,_18,_19,_20,	\
		 _21,_22,_23,_24,_25,_26,_27,_28,_29,_30,	\
		 _31,_32,_33,_34,_35,_36,_37,_38,_39,_40,	\
		 _41,_42,_43,_44,_45,_46,_47,_48,_49,_50,	\
		 _51,_52,_53,_54,_55,_56,_57,_58,_59,_60,	\
		 _61,_62,_63,N,...) N
#define PP_RSEQ_N()		   \
    63,62,61,60,		   \
    59,58,57,56,55,54,53,52,51,50, \
    49,48,47,46,45,44,43,42,41,40, \
    39,38,37,36,35,34,33,32,31,30, \
    29,28,27,26,25,24,23,22,21,20, \
    19,18,17,16,15,14,13,12,11,10, \
    9,8,7,6,5,4,3,2,1,0

// Macro Task
#ifdef __FFLASFFPACK_USE_OPENMP
#define TASK(r, w, rw, f, Args...)				\
  PRAGMA_OMP_TASK_IMPL( omp task GLOBALSHARED(x  r w rw) )	\
  f(Args)
#else 
#ifdef __FFLASFFPACK_USE_KAAPI
#define TASK(r, w, rw, f, Args...) SPAWN(f, NUMARGS(Args)) <Field> >()(Args)
#else
#define TASK(r,w,rw,f,Args...) f(Args)
#endif
#endif

#define STR_IMPL(a) #a
#define STR(a) STR_IMPL(a)

// Macro PLUQ
#ifdef __FFLASFFPACK_USE_OPENMP 
#define PPLUQ(Rank, F, Diag, M, N, A, lda, P, Q) do{		\
    Rank = pPLUQ(F, Diag, M, N, A, lda, P, Q);			\
  }while(0)
#else 
#ifdef __FFLASFFPACK_USE_KAAPI
#define PPLUQ(Rank, F, Diag, M, N, A, lda, P, Q) do{			\
    Rank = PPLUQ(F, Diag, M, N, A, lda, P, Q);				\
  }while(0)
#else
#define PPLUQ(Rank, F, Diag, M, N, A, lda, P, Q) do{		\
    Rank = PPLUQ(F, Diag, M, N, A, lda, P, Q);			\
  }while(0)
#endif
#endif 

// Macro pragma omp taskwait
#ifdef __FFLASFFPACK_USE_OPENMP 
#define WAIT PRAGMA_OMP_TASK_IMPL( omp taskwait )
#else
#ifdef __FFLASFFPACK_USE_KAAPI
#define WAIT do{				\
  }while(0)
#else
#define WAIT do{				\
 }while(0)
#endif
#endif 

// Macro barrier
#ifdef __FFLASFFPACK_USE_OPENMP 
#define BARRIER do{				\
  }while(0)
#else
#ifdef __FFLASFFPACK_USE_KAAPI
#define BARRIER do{				\
    ka::Sync();					\
  }while(0)
#else
#define BARRIER do{				\
  }while(0)
#endif
#endif 

// Parallel Region
#ifdef __FFLASFFPACK_USE_OPENMP
#define HPAC_PAR_REGION  PRAGMA_OMP_TASK_IMPL( omp parallel )  \
  PRAGMA_OMP_TASK_IMPL( omp single ) 
#else
#define HPAC_PAR_REGION 
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


// Parallel Main

#ifdef __FFLASFFPACK_USE_OPENMP
#define BEGIN_PARALLEL_MAIN(Args...) int main(Args)  {
#endif
#ifdef __FFLASFFPACK_USE_KAAPI
#define BEGIN_PARALLEL_MAIN(Args...)			\
  struct doit	{					\
  void operator()(int argc, char** argv) 
#endif


#ifdef __FFLASFFPACK_USE_OPENMP
#define END_PARALLEL_MAIN(void)  return 0; }
#endif
#ifdef __FFLASFFPACK_USE_KAAPI
#define END_PARALLEL_MAIN(void)						\
  }; int main(int argc, char** argv) {					\
    try { ka::Community com = ka::System::join_community( argc, argv ); \
      ka::SpawnMain<doit>()(argc, argv);				\
      com.leave();							\
      ka::System::terminate();}						\
    catch (const std::exception& E) { ka::logfile() << "Catch : " << E.what() << std::endl;} \
    catch (...) { ka::logfile() << "Catch unknown exception: " << std::endl;} \
    return 0;}
#endif

#endif //__FFLASFFPACK_fflas_parallel_H
