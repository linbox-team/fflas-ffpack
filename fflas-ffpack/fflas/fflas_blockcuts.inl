/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* fflas/fflas_bounds.inl
 * Copyright (C) 2013 Jean-Guillaume Dumas
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
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

#ifndef __FFLASFFPACK_fflas_blockcuts_INL
#define __FFLASFFPACK_fflas_blockcuts_INL


#define __FFLASFFPACK_MINBLOCKCUTS 256

namespace FFLAS {

    template<CuttingStrategy Method>
    void BlockCuts(size_t& RBLOCKSIZE, size_t& CBLOCKSIZE,
                   const size_t m, const size_t n,
                   const size_t numthreads);

    template<>
    void BlockCuts<ROW_FIXED>(size_t& RBLOCKSIZE,
                              size_t& CBLOCKSIZE,
                              const size_t m, const size_t n,
                              const size_t numthreads) {
        RBLOCKSIZE = __FFLASFFPACK_MINBLOCKCUTS;
        CBLOCKSIZE = n;
    }


    template<>
    void BlockCuts<COLUMN_FIXED>(size_t& RBLOCKSIZE,
                                 size_t& CBLOCKSIZE,
                                 const size_t m, const size_t n,
                                 const size_t numthreads) {
        RBLOCKSIZE = m;
        CBLOCKSIZE =__FFLASFFPACK_MINBLOCKCUTS;
    }


    template<>
    void BlockCuts<BLOCK_FIXED>(size_t& RBLOCKSIZE,
                                size_t& CBLOCKSIZE,
                                const size_t m, const size_t n,
                                const size_t numthreads) {
        RBLOCKSIZE =__FFLASFFPACK_MINBLOCKCUTS;
        CBLOCKSIZE =__FFLASFFPACK_MINBLOCKCUTS;
    }

    template<>
    void BlockCuts<ROW_THREADS>(size_t& RBLOCKSIZE,
                                size_t& CBLOCKSIZE,
                                const size_t m, const size_t n,
                                const size_t numthreads) {
        RBLOCKSIZE = std::max(m/numthreads,(size_t)1);
        CBLOCKSIZE = n;
    }


    template<>
    void BlockCuts<COLUMN_THREADS>(size_t& RBLOCKSIZE,
                                   size_t& CBLOCKSIZE,
                                   const size_t m, const size_t n,
                                   const size_t numthreads) {
        RBLOCKSIZE = m;
        CBLOCKSIZE = std::max(n/numthreads,(size_t)1);
    }

    template<>
    void BlockCuts<BLOCK_THREADS>(size_t& RBLOCKSIZE,
                                  size_t& CBLOCKSIZE,
                                  const size_t m, const size_t n,
                                  const size_t numthreads) {
        if (numthreads<65) {
                //CP: Let's not compute these values all the time
            const short maxtr[64] = {1,1,1,2,1,2,1,2,3,2,1,3,1,2,3,4,1,3,1,4,3,2,1,4,5,2,3,4,1,5,1,4,3,2,5,6,1,2,3,5,1,6,1,4,5,2,1,6,7,5,3,4,1,6,5,7,3,2,1,6,1,2,7,8};
            const short maxtc[64] = {1,2,3,2,5,3,7,4,3,5,11,4,13,7,5,4,17,6,19,5,7,11,23,6,5,13,9,7,29,6,31,8,11,17,7,6,37,19,13,8,41,7,43,11,9,23,47,8,7,10,17,13,53,9,11,8,19,29,59,10,61,31,9,8};
            
            RBLOCKSIZE=std::max(m/maxtr[numthreads-1],(size_t)1);
            CBLOCKSIZE=std::max(n/maxtc[numthreads-1],(size_t)1);
        } else {
            const size_t maxt = (size_t)sqrt((double)numthreads);
            size_t maxtr=maxt,maxtc=maxt;
            for(size_t i=maxt; i>=1; --i) {
                size_t j=maxt;
                size_t newpr = i*j;
                for( ; newpr < numthreads; ++j, newpr+=i ) {}
                if (newpr == numthreads) {
                    maxtr = i;
                    maxtc = j;
                    break;
                }
            }
            RBLOCKSIZE=std::max(m/maxtr,(size_t)1);
            CBLOCKSIZE=std::max(n/maxtc,(size_t)1);
        }
    }
    
    void BlockCuts(size_t& r, size_t& c,
                   size_t m, size_t n,
                   const CuttingStrategy method,
                   const size_t t) {
        switch(method) {
            case BLOCK_THREADS: BlockCuts<BLOCK_THREADS>(r,c,m,n,t); break;
            case ROW_THREADS: BlockCuts<ROW_THREADS>(r,c,m,n,t); break;
            case ROW_FIXED: BlockCuts<ROW_FIXED>(r,c,m,n,t); break;
            case BLOCK_FIXED: BlockCuts<BLOCK_FIXED>(r,c,m,n,t); break;
            case COLUMN_THREADS: BlockCuts<COLUMN_THREADS>(r,c,m,n,t); break;
            case COLUMN_FIXED: BlockCuts<COLUMN_FIXED>(r,c,m,n,t); break;
            default: BlockCuts<BLOCK_THREADS>(r,c,m,n,t);
        };

    }



    void BlockCuts(size_t& rowBlockSize, size_t& colBlockSize,
                   size_t& lastRBS, size_t& lastCBS,
                   size_t& numRowBlock, size_t& numColBlock,
                   size_t m, size_t n,
                   const CuttingStrategy method,
                   const size_t numthreads) {
        BlockCuts(rowBlockSize, colBlockSize, m, n, method, numthreads);
        numRowBlock = m/rowBlockSize;
        numColBlock = n/colBlockSize;
        lastRBS = m % rowBlockSize;
        lastCBS = n % colBlockSize;

            // Better preserve numRowBlock and numColBlock
        if (lastRBS) {
            lastRBS = m-rowBlockSize*numRowBlock;
            ++rowBlockSize;
        } else lastRBS = rowBlockSize;
        if (lastCBS) {
            lastCBS = n-colBlockSize*numColBlock;
            ++colBlockSize;
        } else lastCBS = colBlockSize;

//             // Better preserve rowBlockSize and colBlockSize
//         if (lastRBS) ++numRowBlock; else lastRBS = rowBlockSize;
//         if (lastCBS) ++numColBlock; else lastCBS = colBlockSize;
    }


}



namespace FFLAS {

    struct ForStrategy1D {
        ForStrategy1D(const size_t n, const CuttingStrategy method, const size_t numthreads) {
//             std::cout<<"FS1D n : "<<n<<std::endl;
//             std::cout<<"FS1D method    : "<<method<<std::endl;
//             std::cout<<"FS1D numthreads : "<<numthreads<<std::endl;

            if ( method == BLOCK_THREADS || method == ROW_THREADS || method == COLUMN_THREADS) {
                numBlock = std::max(numthreads,(size_t)1);
            } else {
                numBlock = std::max(n/__FFLASFFPACK_MINBLOCKCUTS,(size_t)1);
            }
            firstBlockSize = n/numBlock;
            if (firstBlockSize<1) {
                firstBlockSize = (size_t)1;
                numBlock = n;
            }
            changeBS = n - numBlock*firstBlockSize;
            lastBlockSize = firstBlockSize;
            if (changeBS) ++firstBlockSize;

//             std::cout<<"FS1D 1BLOCKSIZE : "<<firstBlockSize<<std::endl;
//             std::cout<<"FS1D 2BLOCKSIZE : "<<lastBlockSize<<std::endl;
//             std::cout<<"FS1D changeBS : "<<changeBS<<std::endl;
//             std::cout<<"FS1D NBlocks : "<<numBlock<<std::endl;
        }
        
        size_t begin() { 
            ibeg = 0; iend = firstBlockSize; 
//             std::cout << "FS1D 0   : " << 0 << std::endl;
//             std::cout << "FS1D ibeg: " << ibeg << std::endl;
//             std::cout << "FS1D iend: " << iend << std::endl;
            
            return current = 0; 
        }
        bool end() const { return current == numBlock; }
        size_t operator++() { 
            ibeg = iend;
            iend += ++current<changeBS?firstBlockSize:lastBlockSize;

//             std::cout << "FS1D i   : " << current << std::endl;
//             std::cout << "FS1D ibeg: " << ibeg << std::endl;
//             std::cout << "FS1D iend: " << iend << std::endl;
            

            return current;
        }
        
        size_t ibeg, iend;

    protected:
        size_t current;
        size_t firstBlockSize,lastBlockSize; 
        size_t changeBS;
        size_t numBlock;

    };
    

    struct ForStrategy2D {
        ForStrategy2D(const size_t m, const size_t n, const CuttingStrategy method, const size_t numthreads) {
            BlockCuts(rowBlockSize, colBlockSize, 
                      lastRBS, lastCBS,
                      numRowBlock, numColBlock,
                      m, n, method, numthreads);

            BLOCKS = numRowBlock * numColBlock;
//             std::cout<<"RBLOCKSIZE : "<<rowBlockSize<<std::endl;
//             std::cout<<"CBLOCKSIZE : "<<colBlockSize<<std::endl;
//             std::cout<<"lastRBS    : "<<lastRBS<<std::endl;
//             std::cout<<"lastCBS    : "<<lastCBS<<std::endl;
//             std::cout<<"NrowBlocks : "<<numRowBlock<<std::endl;
//             std::cout<<"NcolBlocks : "<<numColBlock<<std::endl;
        }
        
        size_t begin() { current = 0; return setCurrentBlock(); }
        bool end() const { return current == BLOCKS; }
        size_t setCurrentBlock() { 
			ibeg = current / numColBlock;
			jbeg = current % numColBlock;
			size_t BlockRowDim = rowBlockSize;
			if (ibeg == numRowBlock-1)
				BlockRowDim = lastRBS;
			size_t BlockColDim = colBlockSize;
			if (jbeg == numColBlock-1)
				BlockColDim = lastCBS;
            ibeg *= rowBlockSize;
            jbeg *= colBlockSize;
            iend = ibeg + BlockRowDim;
            jend = jbeg + BlockColDim;
            return current;
        }
        
        size_t operator++() { ++current; return setCurrentBlock(); }
        
        size_t ibeg, iend, jbeg, jend;

    protected:
        size_t current;
        size_t rowBlockSize; size_t colBlockSize;
        size_t lastRBS; size_t lastCBS;
        size_t numRowBlock; size_t numColBlock;
        size_t BLOCKS;
        
   };
    
}



#endif

