

#ifdef __FFLASFFPACK_USE_OPENMP
#include "omp.h"
#endif

#ifdef __FFLASFFPACK_USE_KAAPI
#include "kaapi++"
//#include "kaapi-routines.h"
#endif


//OMP
#define comma ,
#define read(Args...) comma Args
#define write(Args...) comma Args
#define readwrite(Args...) comma Args
#define noread(void) 
#define nowrite(void) 
#define noreadwrite(void) 


// KAAPI

//#define DSLREADARRAY(A,pointer, range) A,pointer,range 

#define GlobalShared(a, Args...) shared(Args)
#ifdef __FFLASFFPACK_USE_OPENMP
    #define TASK(r, w, rw, f, Args...) \
    PRAGMA_OMP_TASK_IMPL( omp task GlobalShared(x  r w rw) )	\
    f(Args)
#else 
   #ifdef __FFLASFFPACK_USE_KAAPI
      #define TASK(r, w, rw, f,Args...)  ka::Spawn<Task ## f <Field> >()(Args)
   #else
      #define Task(r,w,rw,f,Args...) f(Args)
   #endif
#endif

#define SHARED(a,b,c,d)

#define STR_IMPL(a) #a
#define STR(a) STR_IMPL(a)

#define PRAGMA_OMP_TASK(F,A,B,C, ARGS...) PRAGMA_OMP_TASK_IMPL( omp task shared( F, A, B, C ) )
#define PRAGMA_OMP_TASK_IMPL( F ) _Pragma( #F )

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
