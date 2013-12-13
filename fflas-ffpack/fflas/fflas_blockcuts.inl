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
        RBLOCKSIZE = MAX(m/numthreads,1);
        CBLOCKSIZE = n;
    }


    template<>
    void BlockCuts<COLUMN_THREADS>(size_t& RBLOCKSIZE,
                                   size_t& CBLOCKSIZE,
                                   const size_t m, const size_t n,
                                   const size_t numthreads) {
        RBLOCKSIZE = m;
        CBLOCKSIZE = MAX(n/numthreads,1);
    }

    template<>
    void BlockCuts<BLOCK_THREADS>(size_t& RBLOCKSIZE,
                                  size_t& CBLOCKSIZE,
                                  const size_t m, const size_t n,
                                  const size_t numthreads) {
	    //CP: Let's not compute these values all the time
            const short maxtr[64] = {1,1,1,2,1,2,1,2,3,2,1,3,1,2,3,4,1,3,1,4,3,2,1,4,5,2,3,4,1,5,1,4,3,2,5,6,1,2,3,5,1,6,1,4,5,2,1,6,7,5,3,4,1,6,5,7,3,2,1,6,1,2,7,8};
            const short maxtc[64] = {1,2,3,2,5,3,7,4,3,5,11,4,13,7,5,4,17,6,19,5,7,11,23,6,5,13,9,7,29,6,31,8,11,17,7,6,37,19,13,8,41,7,43,11,9,23,47,8,7,10,17,13,53,9,11,8,19,29,59,10,61,31,9,8};

       // const size_t maxt = (size_t)sqrt((double)numthreads);
       // 	size_t maxtr=maxt,maxtc=maxt;
       // 	for(size_t i=maxt; i>=1; --i) {
       // 		size_t j=maxt;
       //          size_t newpr = i*j;
       // 		for( ; newpr < numthreads; ++j, newpr+=i ) {}
       //          if (newpr == numthreads) {
       //              maxtr = i;
       //              maxtc = j;
       //              break;
       //          }
       //  }

	    RBLOCKSIZE=MAX(m/maxtr[numthreads-1],1);
	    CBLOCKSIZE=MAX(n/maxtc[numthreads-1],1);
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
        if (lastRBS) ++numRowBlock; else lastRBS = rowBlockSize;
        lastCBS = n % colBlockSize;
        if (lastCBS) ++numColBlock; else lastCBS = colBlockSize;
    }
    

}



#endif

