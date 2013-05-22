#define RBLOCKSIZE 512
#define CBLOCKSIZE 100000
#include <omp.h>
namespace FFLAS {
template<class Field>
inline typename Field::Element*
pfgemm( const Field& F,
               const FFLAS_TRANSPOSE ta,
               const FFLAS_TRANSPOSE tb,
               const size_t m,
               const size_t n,
               const size_t k,
               const typename Field::Element alpha,
               const typename Field::Element* A, const size_t lda,
               const typename Field::Element* B, const size_t ldb,
               const typename Field::Element beta,
               typename Field::Element* C, const size_t ldc,
               const size_t w){
    
    size_t NrowBlocks = m/RBLOCKSIZE;
    size_t LastrowBlockSize = m % RBLOCKSIZE;
    if (LastrowBlockSize) 
        NrowBlocks++;
    else 
        LastrowBlockSize = RBLOCKSIZE;
    size_t NcolBlocks = n/CBLOCKSIZE;
    size_t LastcolBlockSize = n % CBLOCKSIZE;
    if (LastcolBlockSize) 
        NcolBlocks++;
    else 
        LastcolBlockSize = CBLOCKSIZE;

#pragma omp parallel for default (none) shared (A, B, C, F, NcolBlocks, NrowBlocks, LastcolBlockSize, LastrowBlockSize)
    for (size_t t = 0; t < NrowBlocks*NcolBlocks; ++t){
        size_t i = t / NcolBlocks;
        size_t j = t % NcolBlocks;
        size_t BlockRowDim = RBLOCKSIZE;
        if (i == NrowBlocks-1)
            BlockRowDim = LastrowBlockSize;
        size_t BlockColDim = CBLOCKSIZE;
        if (j == NcolBlocks-1)
            BlockColDim = LastcolBlockSize;

        fgemm( F, ta, tb, BlockRowDim, BlockColDim, k, alpha, A + RBLOCKSIZE * i*lda, lda, B + CBLOCKSIZE * j, ldb, beta, C+ RBLOCKSIZE*i*ldc+j*CBLOCKSIZE, ldc, w);
    }
    return C;
}
} // FFLAS
