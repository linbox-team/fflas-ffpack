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


#define __FFLASFFPACK_SEQUENTIAL

#ifdef __FFLASFFPACK_USE_OPENMP
#undef __FFLASFFPACK_SEQUENTIAL
#include "omp.h"
#elif defined (__FFLASFFPACK_USE_KAAPI)
#undef __FFLASFFPACK_SEQUENTIAL
#include "kaapi++"
#endif

#ifdef __FFLASFFPACK_FORCE_SEQ

#undef __FFLASFFPACK_USE_OPENMP
#undef __FFLASFFPACK_USE_KAAPI
#undef __FFLASFFPACK_USE_TBB
#define __FFLASFFPACK_SEQUENTIAL

#endif

/*********************************************************/
/*********************** SEQUENTIAL***********************/
/*********************************************************/

#ifdef __FFLASFFPACK_SEQUENTIAL // MACRO for sequential execution

// TASK is a function call
#define TASK(r,w,rw,f,Args...) f(Args)

#define WAIT
#define BARRIER
#define PAR_REGION
#define NUM_THREADS 1
#define MAX_THREADS 1

#define BEGIN_PARALLEL_MAIN(Args...) int main(Args)  {
#define END_PARALLEL_MAIN(void)  return 0; }

#endif // Macro for sequential


/*********************************************************/
/************************* OPENMP ************************/
/*********************************************************/

#ifdef __FFLASFFPACK_USE_OPENMP //OpenMP macros

#define GLOBALSHARED(a, Args...) shared(Args)
#define PRAGMA_OMP_TASK_IMPL( F ) _Pragma( #F )

// Task definition with OpenMP
#define TASK(r, w, rw, f, Args...)				\
  PRAGMA_OMP_TASK_IMPL( omp task GLOBALSHARED(x  r w rw) )	\
  f(Args)
// macro omp taskwait (waits for all childs of current task)
#define WAIT PRAGMA_OMP_TASK_IMPL( omp taskwait )

#define BARRIER

// parallel region
#define PAR_REGION  PRAGMA_OMP_TASK_IMPL( omp parallel )  \
  PRAGMA_OMP_TASK_IMPL( omp single )
// get number of threads in the parallel region
# define NUM_THREADS omp_get_num_threads()
// get number of threads specified with the global variable OMP_NUM_THREADS
# define MAX_THREADS omp_get_max_threads()

#define BEGIN_PARALLEL_MAIN(Args...) int main(Args)  {
#define END_PARALLEL_MAIN(void)  return 0; }

#endif // OpenMP macros


/*********************************************************/
/************************* KAAPI *************************/
/*********************************************************/

#ifdef __FFLASFFPACK_USE_KAAPI // KAAPI

#define SPAWN(f,N) CONCATENATE_ARGS(ka::Spawn<Task ## f, N)
#define CONCATENATE_ARGS(f, N) f ## N 

// TASK definition
#define TASK(r, w, rw, f, Args...) SPAWN(f, NUMARGS(Args)) <Field> >()(Args)

// WAIT do nothing in kaapi
#define WAIT

// BARRIER means synchronization in kaapi (waits for the execution of all tasks)
#define BARRIER do{				\
    ka::Sync();					\
  }while(0)

#define PAR_REGION

// Number of threads
#  define NUM_THREADS kaapi_getconcurrency_cpu()

#  define MAX_THREADS kaapi_getconcurrency_cpu()

// Begin parallel main
#define BEGIN_PARALLEL_MAIN(Args...)			\
  struct doit	{					\
  void operator()(int argc, char** argv)
//end parallel main
#define END_PARALLEL_MAIN(void)						\
  }; int main(int argc, char** argv) {					\
    try { ka::Community com = ka::System::join_community( argc, argv ); \
      ka::SpawnMain<doit>()(argc, argv);				\
      com.leave();							\
      ka::System::terminate();}						\
    catch (const std::exception& E) { ka::logfile() << "Catch : " << E.what() << std::endl;} \
    catch (...) { ka::logfile() << "Catch unknown exception: " << std::endl;} \
    return 0;}



// Macro computes number of Arguments

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
#endif // KAAPI macros




/*********************************************************/
/******************* TBB(coming sooon) *******************/ 
/*********************************************************/
#ifdef __FFLASFFPACK_USE_TBB



#endif

// common macros
#define COMMA ,
#define READ(Args...) COMMA Args
#define WRITE(Args...) COMMA Args
#define READWRITE(Args...) COMMA Args
#define NOREAD(void)
#define NOWRITE(void)
#define NOREADWRITE(void)

#define PPLUQ(Rank, F, Diag, M, N, A, lda, P, Q) do{\
  Rank = pPLUQ(F, Diag, M, N, A, lda, P, Q);\
  }while(0)

#endif //__FFLASFFPACK_fflas_parallel_H
