/* fflas/parallel.h
 * Copyright (C) 2013 Jean Guillaume Dumas Clement Pernet Ziad Sultan
 *
 * Written by Jean Guillaume Dumas Clement Pernet Ziad Sultan
 *<ziad.sultan@imag.fr>
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


#include "fflas-ffpack/config.h"

#ifndef __FFLASFFPACK_USE_OPENMP
#define  __FFLASFFPACK_SEQUENTIAL
#else
#include "omp.h"
#endif

#ifdef __FFLASFFPACK_SEQUENTIAL
#undef __FFLASFFPACK_USE_OPENMP
#undef __FFLASFFPACK_USE_TBB
#undef __FFLASFFPACK_USE_KAAPI
#elif defined (__FFLASFFPACK_USE_KAAPI)
#undef __FFLASFFPACK_SEQUENTIAL
#undef __FFLASFFPACK_USE_TBB
#undef __FFLASFFPACK_USE_OPENMP
#include "kaapi++"
#include "fflas-ffpack/fflas/kaapi_routines.inl"
#elif defined __FFLASFFPACK_USE_TBB
#undef __FFLASFFPACK_SEQUENTIAL
#undef __FFLASFFPACK_USE_OPENMP
#undef __FFLASFFPACK_USE_KAAPI
#include </usr/include/tbb/tbb.h>
#include </usr/include/tbb/task.h>
#include </usr/include/tbb/parallel_for.h>
#include </usr/include/tbb/task_group.h>
/*
      extern "C"
      {
      tbb::task_group g;
      }
      */
#ifdef __FFLASFFPACK_HAVE_MKL
#ifndef _MKL_H_ // temporary
#error "MKL (mkl.h) not present, while you have MKL enabled"
#endif
#undef index_t
#define index_t MKL_INT
#endif // __FFLASFFPACK_HAVE_MKL


#endif

#ifdef __FFLASFFPACK_FORCE_SEQ

#undef __FFLASFFPACK_USE_OPENMP
#undef __FFLASFFPACK_USE_KAAPI
#undef __FFLASFFPACK_USE_TBB
#define __FFLASFFPACK_SEQUENTIAL

#endif

#ifndef index_t
#define index_t size_t
#endif



/*********************************************************/
/*********************** SEQUENTIAL***********************/
/*********************************************************/

#ifdef __FFLASFFPACK_SEQUENTIAL // MACRO for sequential execution

// TASK is a function call
#define TASK(M, I) {I;}

#define WAIT
#define CHECK_DEPENDENCIES
#define BARRIER
#define PAR_BLOCK


#define SYNCH_GROUP(Args...) {{Args};}

#define THREAD_INDEX 0
#define NUM_THREADS 1
#define SET_THREADS(num_threads) {}
#define MAX_THREADS 1

#define READ(Args...)
#define WRITE(Args...)
#define READWRITE(Args...)
#define CONSTREFERENCE(...)
#define VALUE(...)

#define BEGIN_PARALLEL_MAIN(Args...) int main(Args)  {
#define END_PARALLEL_MAIN(void)  return 0; }

// for 1D with iterator control and range access through iter (strategy 1D)
#define FORBLOCK1D(iter, m, Helper, Args...)                            \
{ FFLAS::ForStrategy1D<std::remove_const<decltype(m)>::type, typename decltype(Helper)::Cut, typename  decltype(Helper)::Param> iter(m, Helper); \
    for(iter.initialize(); !iter.isTerminated(); ++iter)            \
    {Args;} }

// for strategy 1D
#define FOR1D(i, m, Helper, Args...)                                    \
FORBLOCK1D(_internal_iterator, m, Helper,                           \
           for(auto i=_internal_iterator.begin(); i!=_internal_iterator.end(); ++i) \
           { Args; })

// PARFOR1D does normal execution of the loop
#define PARFORBLOCK1D(iter,  m, Helper, Args...)			       \
for(std::remove_const<decltype(m)>::type iter=0; iter<m; ++iter)     \
{ Args; }

// PARFOR1D does normal execution of the loop
#define PARFOR1D(iter,  m, Helper, Args...)				       \
for(std::remove_const<decltype(m)>::type iter=0; iter<m; ++iter)     \
{ Args; }


////////////////////   CUTTING LOOP MACROS 2D //////////////////////

// for strategy 2D with access to the range and control of iterator
#define FORBLOCK2D(iter, m, n, Helper, Args...)                         \
{ FFLAS::ForStrategy2D<std::remove_const<decltype(m)>::type, typename decltype(Helper)::Cut, typename  decltype(Helper)::Param> iter(m,n,Helper); \
    for(iter.initialize(); !iter.isTerminated(); ++iter)             \
    { Args; } }

// for strategy 2D
#define FOR2D(i, j, m, n, Helper, Args...)                              \
FORBLOCK2D(_internal_iterator, m, n, Helper,                        \
           for(auto i=_internal_iterator.ibegin(); i!=_internal_iterator.iend(); ++i) \
           for(auto j=_internal_iterator.jbegin(); j!=_internal_iterator.jend(); ++j) \
           { Args; })

// parallel for strategy 2D with access to the range and control of iterator
#define PARFORBLOCK2D(iter, m, n, Helper, Args...)      \
FORBLOCK2D(iter, m, n, Helper, Args)

// parallel for strategy 2D
#define PARFOR2D(i, j, m, n, Helper, Args...)   \
FOR2D(i, j, m, n, Helper, Args)


#endif // Macro for sequential


/*********************************************************/
/************************* OPENMP ************************/
/*********************************************************/

#ifdef __FFLASFFPACK_USE_OPENMP //OpenMP macros

#define PRAGMA_OMP_IMPL(Args...) _Pragma(#Args)

#define TASK(M, I)                             \
PRAGMA_OMP_IMPL(omp task M)           \
{I;}


#define SYNCH_GROUP(Args...)     {{Args};}\
WAIT;


// macro omp taskwait (waits for all childs of current task)
#define WAIT PRAGMA_OMP_IMPL(omp taskwait)
#define GLOBALSHARED(a, Args...) shared(Args)
#define CONSTREFERENCE(Args...) shared(Args)
#define VALUE(Args...) firstprivate(Args)
#define BARRIER PRAGMA_OMP_IMPL(omp barrier)

////////////////////   CUTTING LOOP MACROS 1D //////////////////////


// for with iterator control and range access through iter (strategy 1D)
// Warning: by default we assume that there is no dependency between each iteration, hence we pass an empty MODE() to the tasks.
// TODO: add an optional MODE argument to the parameter list of FORBLOCK1D
#define FORBLOCK1D(iter, m, Helper, Args...)                                       \
{ FFLAS::ForStrategy1D<std::remove_const<decltype(m)>::type, typename decltype(Helper)::Cut, typename  decltype(Helper)::Param > iter(m, Helper); \
    for(iter.initialize(); !iter.isTerminated(); ++iter){ {Args;}  } }



// for strategy 1D
// WARNING: the inner code Args should not contain any coma outside parenthesis (e.g. declaration lists, and template param list)
#define FOR1D_1(i, m, Helper, Args ...)                             \
FORBLOCK1D(_internal_iterator, m, Helper,                           \
           const auto internal_iter_begin(_internal_iterator.begin());      \
           const auto internal_iter_end(_internal_iterator.end());      \
           TASK(VALUE(internal_iter_begin, internal_iter_end), \
                 {for(auto i=internal_iter_begin; i!=internal_iter_end; ++i) \
                 { Args; } });)                                         \
                 WAIT;


// WARNING: @fixme, the passed mode should not contain another VALUE mode
#define FOR1D_2(i, m, Helper, mode, Args ...)                             \
FORBLOCK1D(_internal_iterator, m, Helper,                           \
           const auto internal_iter_begin(_internal_iterator.begin());      \
           const auto internal_iter_end(_internal_iterator.end());      \
           TASK(VALUE(internal_iter_begin, internal_iter_end) mode, \
                 {for(auto i=internal_iter_begin; i!=internal_iter_end; ++i) \
                 { Args; } });)                                         \
                 WAIT;


#define CONCATVAARGS(a, ...) CONCATVAARGS_(a, __VA_ARGS__)
#define CONCATVAARGS_(a, ...) a##__VA_ARGS__
#define FOR1D(i, m, Helper, ...) CONCATVAARGS(FOR1D_, NUMARGS(__VA_ARGS__))(i, m, Helper, __VA_ARGS__)

/*
#define PARFORBLOCK1D(iter, m, Helper, Args...)                         \
{ FFLAS::ForStrategy1D<std::remove_const<decltype(m)>::type > iter(m, Helper); \
PRAGMA_OMP_IMPL(omp parallel for num_threads(iter.numblocks()) schedule(runtime)) \
for(iter.initialize(); !iter.isTerminated(); ++iter)        \
{Args;} }
*/

//parallel for 1D with iterator control and range access cannot be done with openmp: syntax of openmp does not allow the use of the iterator syntax
// Thus, PARFORBLOCK1D and PARFOR1D have the same implementation with no cutting. If using OpenMP the user can specify the cutting in runtime using the environmental variable: (see OpenMP spec for more details)
// export OMP_SCHEDULE="DYNAMIC"
// or export OMP_SCHEDULE="GUIDED,4"
// or export OMP_SCHEDULE="STATIC"
// or export OMP_SCHEDULE="AUTO"

#define PARFORBLOCK1D(iter, m, Helper, Args...)                       \
{ FFLAS::ForStrategy1D<std::remove_const<decltype(m)>::type, typename decltype(Helper)::Cut, typename  decltype(Helper)::Param > OMPstrategyIterator(m, Helper); \
    PRAGMA_OMP_IMPL(omp parallel for num_threads(OMPstrategyIterator.numblocks()) schedule(runtime)) \
    for(std::remove_const<decltype(m)>::type iter=0; iter<m; ++iter) \
    { Args; } }


// parallel for 1D
#define PARFOR1D(iter, m, Helper, Args...)                \
{ auto h = Helper; \
    PARFORBLOCK1D(iter, m, h, { Args; } ) \
}


////////////////////   CUTTING LOOP MACROS 2D //////////////////////

// for strategy 2D with access to the range and control of iterator
#define FORBLOCK2D(iter, m, n, Helper, Args...)                         \
{ auto h=Helper; \
    FFLAS::ForStrategy2D<std::remove_const<decltype(m)>::type, typename decltype(h)::Cut, typename  decltype(h)::Param  > iter(m,n,h); \
    for(iter.initialize(); !iter.isTerminated(); ++iter)            \
    {Args;} }

// for strategy 2D
// WARNING: the inner code Args should not contain any coma outside parenthesis (e.g. declaration lists, and template param list)
#define FOR2D(i, j, m, n, Helper, Args...)                              \
FORBLOCK2D(_internal_iterator, m, n, Helper,                        \
           TASK(,                                             \
                for(auto i=_internal_iterator.ibegin(); i!=_internal_iterator.iend(); ++i) \
                for(auto j=_internal_iterator.jbegin(); j!=_internal_iterator.jend(); ++j) \
                { Args; });) \
                WAIT;

// parallel for strategy 2D with access to the range and control of iterator
// WARNING: This is not doable : OMP requires an iteration over an interval of ints.
/* #define PARFORBLOCK2D(iter, m, n, Helper, Args...)                   \
 *      { FFLAS::ForStrategy2D<std::remove_const<decltype(m)>::type, typename decltype(Helper)::Cut, typename  decltype(Helper)::Param  > iter(m,n,Helper); \
 *          PRAGMA_OMP_IMPL(omp parallel for num_threads(iter.rownumblocks()*iter.colnumblocks()) schedule(runtime)) \
 *              for(iter.initialize(); !iter.isTerminated(); ++iter)    \
 *              {Args;} }
 */

// parallel for strategy 2D
/* #define PARFOR2D(i, j, m, n, Helper, Args...)                        \
 *     PARFORBLOCK2D(_internal_iterator, m, n, Helper,                  \
 *                   for(auto i=_internal_iterator.ibegin(); i!=_internal_iterator.iend(); ++i) \
 *                       for(auto j=_internal_iterator.jbegin(); j!=_internal_iterator.jend(); ++j) \
 *                       { Args; })
 */

// parallel region
#define PAR_BLOCK  PRAGMA_OMP_IMPL(omp parallel)   \
PRAGMA_OMP_IMPL(omp single)

# define THREAD_INDEX omp_get_thread_num()
// get the number of threads in the parallel region
# define NUM_THREADS omp_get_num_threads()
// set the number of threads in the parallel region
# define SET_THREADS(numthreads) omp_set_num_threads(numthreads)
// get the number of threads specified with the global variable OMP_NUM_THREADS
# define MAX_THREADS omp_get_max_threads()

#define BEGIN_PARALLEL_MAIN(Args...) int main(Args)  {
#define END_PARALLEL_MAIN(void)  return 0; }


//////////////////////////////////////////////
/////////////// dataflow macros //////////////
#ifdef __FFLASFFPACK_USE_DATAFLOW // OMP dataflow synch DSL features

#define READ(Args...) depend(in: Args)
#define WRITE(Args...) depend(out: Args)
#define READWRITE(Args...) depend(inout: Args)
//computes dependencies (no wait here)
#define CHECK_DEPENDENCIES

#else // OPENMP3.1 (explicit synch mode)

#define CHECK_DEPENDENCIES PRAGMA_OMP_IMPL(omp taskwait)

#define READ(Args...)
#define WRITE(Args...)
#define READWRITE(Args...)

#endif // end DATAFLOW FLAG
///////////////////////////////////////////////
///////////////////////////////////////////////



#endif // OpenMP macros


/*********************************************************/
/*************************** TBB  ************************/
/*********************************************************/
#ifdef __FFLASFFPACK_USE_TBB

// workaround to overload macro CONSTREFERENCE

// CONSTREFERENCE macro
/* #define REF1(a) =,&a */
/* #define REF2(a,b) =,&a, &b */
/* #define REF3(a,b,c) =,&a,&b,&c */
/* #define REF4(a,b,c,d) =,&a,&b,&c,&d */
/* #define REF5(a,b,c,d,e) =,&a,&b,&c,&d,&e */
/* #define REF6(a,b,c,d,e,f) =,&a,&b,&c,&d,&e,&f */
/* #define REF7(a,b,c,d,e,f,g) =,&a,&b,&c,&d,&e,&f,&g */
/* #define REF8(a,b,c,d,e,f,g,h) =,&a,&b,&c,&d,&e,&f,&g,&h */
/* #define REF9(a,b,c,d,e,f,g,h,i) =,&a,&b,&c,&d,&e,&f,&g,&h,&i */
/* #define REF10(a,b,c,d,e,f,g,h,i,enough) =,&a,&b,&c,&d,&e,&f,&g,&h,&i,&enough */
/* #define GET_REF(_1,_2,_3,_4,_5,_6,_7,_8,_9,_10, NAME,...) NAME */
/* #define CONSTREFERENCE(...) GET_REF(__VA_ARGS__, REF10,REF9,REF8,REF7,REF6,REF5,REF4,REF3,REF2,REF1)(__VA_ARGS__) */


#define REF1(a) ,&a
#define REF2(a,b) ,&a, &b
#define REF3(a,b,c) ,&a,&b,&c
#define REF4(a,b,c,d) ,&a,&b,&c,&d
#define REF5(a,b,c,d,e) ,&a,&b,&c,&d,&e
#define REF6(a,b,c,d,e,f) ,&a,&b,&c,&d,&e,&f
#define REF7(a,b,c,d,e,f,g) ,&a,&b,&c,&d,&e,&f,&g
#define REF8(a,b,c,d,e,f,g,h) ,&a,&b,&c,&d,&e,&f,&g,&h
#define REF9(a,b,c,d,e,f,g,h,i) ,&a,&b,&c,&d,&e,&f,&g,&h,&i
#define REF10(a,b,c,d,e,f,g,h,i,enough) ,&a,&b,&c,&d,&e,&f,&g,&h,&i,&enough
#define GET_REF(_1,_2,_3,_4,_5,_6,_7,_8,_9,_10, NAME,...) NAME
#define CONSTREFERENCE(...) GET_REF(__VA_ARGS__, REF10,REF9,REF8,REF7,REF6,REF5,REF4,REF3,REF2,REF1)(__VA_ARGS__)


// workaround to overload macro VALUE
#define VAL1(a) ,a
#define VAL2(a,b) ,a, b
#define VAL3(a,b,c) ,a,b,c
#define VAL4(a,b,c,d) ,a,b,c,d
#define VAL5(a,b,c,d,e) ,a,b,c,d,e

#define GET_VAL(_1,_2,_3,_4,_5, NAME,...) NAME
#define VALUE(...) GET_VAL(__VA_ARGS__, VAL5,VAL4,VAL3,VAL2,VAL1)(__VA_ARGS__)

// need task_group to lunch a group of tasks in parallel
#define SYNCH_GROUP(Args...) \
{tbb::task_group g;  \
    {{Args};}             \
    g.wait();}


// TBB task
#define TASK(M, I)                              \
{                                           \
    g.run([=M](){I;});			\
}

//#define MODE(Args...) Args
#define WAIT g.wait()
#define CHECK_DEPENDENCIES g.wait()
#define BARRIER
#define PAR_BLOCK

#define THREAD_INDEX tbb::this_task_arena::current_thread_index()
#define NUM_THREADS tbb::task_scheduler_init::default_num_threads()
#define SET_THREADS(numthreads) tbb::task_scheduler_init::initialize(numthreads)
#define MAX_THREADS tbb::task_scheduler_init::default_num_threads()
#define READ(Args...)
#define WRITE(Args...)
#define READWRITE(Args...)

#define BEGIN_PARALLEL_MAIN(Args...) int main(Args)  {
#define END_PARALLEL_MAIN(void)  return 0; }
#define CAPTURE(Args...) [Args]

// for strategy 1D with access to the iterator
#define FORBLOCK1D(iter, m, Helper, Args...)                            \
{ FFLAS::ForStrategy1D<std::remove_const<decltype(m)>::type, typename decltype(Helper)::Cut, typename  decltype(Helper)::Param  > iter(m, Helper); \
    for(iter.initialize(); !iter.isTerminated(); ++iter)            \
    {Args;} }

// for strategy 1D
#define FOR1D(i, m, Helper, Args...)                                    \
FORBLOCK1D(_internal_iterator, m, Helper,                           \
           for(auto i=_internal_iterator.begin(); i!=_internal_iterator.end(); ++i) \
           { Args; } )


// tbb parallel for 1D
#define PARFORBLOCK1D(iter,  m, Helper, Args...)			\
{ FFLAS::ForStrategy1D<std::remove_const<decltype(m)>::type, typename decltype(Helper)::Cut, typename  decltype(Helper)::Param> iter(m, Helper); \
    tbb::parallel_for(							\
                                                tbb::blocked_range<std::remove_const<decltype(m)>::type >(0, m, iter.blocksize() ), \
                                                [=, &iter](const tbb::blocked_range<std::remove_const<decltype(m)>::type > &iter) { \
                                                {Args;} });					\
}

// tbb parallel for 1D
/*
#define PARFOR1D(i,  m, Helper, Args...)                                \
PARFORBLOCK1D(_internal_iterator, m, Helper,                        \
for(auto i=_internal_iterator.begin(); i!=_internal_iterator.end(); ++i) \
{ Args; } )
*/

#define PARFOR1D(i,  m, Helper, Args...)				\
{ FFLAS::ForStrategy1D<std::remove_const<decltype(m)>::type, typename decltype(Helper)::Cut, typename  decltype(Helper)::Param> TBBstrategyIterator(m, Helper);	\
    tbb::parallel_for(							\
                                                tbb::blocked_range<std::remove_const<decltype(m)>::type >(0, m, TBBstrategyIterator.blocksize() ), \
                                                [=](const tbb::blocked_range<std::remove_const<decltype(m)>::type > &TBBblockrangeIterator) { \
                                                for(auto i = TBBblockrangeIterator.begin();		\
                                                    i < TBBblockrangeIterator.end() ; ++i){	\
                                                    {Args;} }});					\
}


// for strategy 2D with access to the iterator
#define FORBLOCK2D(iter, m, n, Helper, Args...)                         \
{ FFLAS::ForStrategy2D<std::remove_const<decltype(m)>::type, typename decltype(Helper)::Cut, typename  decltype(Helper)::Param> iter(m,n,Helper); \
    for(iter.initialize(); !iter.isTerminated(); ++iter)		\
    {Args;} }

// for strategy 2D
#define FOR2D(i, j, m, n, Helper, Args...)                              \
FORBLOCK2D(_internal_iterator, m, n, Helper,				\
           for(auto i=_internal_iterator.ibegin(); i!=_internal_iterator.iend(); ++i) \
           for(auto j=_internal_iterator.jbegin(); j!=_internal_iterator.jend(); ++j) \
           { Args; })

// parallel for strategy 2D with access to the range and control of iterator
#define PARFORBLOCK2D(iter, m, n, Helper, Args...)                      \
{ FFLAS::ForStrategy2D<std::remove_const<decltype(m)>::type, typename decltype(Helper)::Cut, typename  decltype(Helper)::Param> iter(m,n,Helper); \
    tbb::parallel_for(							\
                                                tbb::blocked_range2d<std::remove_const<decltype(m)>::type >(0, m, iter.rowblocksize(), 0, n, iter.colblocksize() ), \
                                                [=, &i](const tbb::blocked_range2d<std::remove_const<decltype(m)>::type > &iter) { \
                                                {Args;} });					\
}

// parallel for strategy 2D
#define PARFOR2D(i, j, m, n, Helper, Args...)                           \
PARFORBLOCK2D(_internal_iterator, m, n, Helper,                     \
              for(auto i=_internal_iterator.ibegin(); i!=_internal_iterator.iend(); ++i) \
              for(auto j=_internal_iterator.jbegin(); j!=_internal_iterator.jend(); ++j) \
              { Args; })

#endif // end TBB macros

/*********************************************************/
/************************* KAAPI *************************/
/*********************************************************/

#ifdef __FFLASFFPACK_USE_KAAPI // KAAPI

#define SPAWN(f,N) CONCATENATE_ARGS(ka::Spawn<Task ## f, N)
#define CONCATENATE_ARGS(f, N) f ## N

// TASK definition
#define TASK(r, w, rw, f, Args...) CONCATENATE_ARGS(spawner,f)(Args)

// WAIT do nothing in kaapi
#define WAIT

// BARRIER means synchronization in kaapi (waits for the execution of all tasks)
#define BARRIER do{				\
    ka::Sync();					\
}while(0)

#define PAR_BLOCK
#define PARFOR1D for

#  define THREAD_INDEX kaapi_get_thread_num()

// Number of threads
#  define NUM_THREADS kaapi_getconcurrency_cpu()
#  define SET_THREADS(numthreads) {}
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

#define SYNCH_GROUP(Args...) {{Args};}



#endif // KAAPI macros



/*********************************************************/
/********************* common macros *********************/
/*********************************************************/

#define COMMA ,
#define MODE(...)  __VA_ARGS__
#define RETURNPARAM(f, P1, Args...) P1=f(Args)

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

#define NOSPLIT() FFLAS::ParSeqHelper::Sequential()

// overload of SPLITTER
#define splitting_0() FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads>()
#define splitting_1(a) FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads>(a)
#define splitting_2(a,c) FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,c>(a)
#define splitting_3(a,b,c) FFLAS::ParSeqHelper::Parallel<b,c>(a)

#define splitt(_1,_2,_3, NAME,...) NAME

#define SPLITTER(...) splitt(__VA_ARGS__, splitting_3, splitting_2, splitting_1, splitting_0)(__VA_ARGS__)


#include "fflas-ffpack/paladin/blockcuts.inl"

#endif //__FFLASFFPACK_fflas_parallel_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
