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

#include <fflas-ffpack/fflas/fflas_enum.h>
#include <math.h>
#include <cassert>

#define __FFLASFFPACK_MINBLOCKCUTS ((size_t)256)

namespace FFLAS {
    namespace CuttingStrategy{
        struct Single{};
        struct Row{};
        struct Column{};
        struct Block{};
        struct Recursive{};
        typedef Row RNSModulus;
    }

    namespace StrategyParameter{
        struct Fixed{};
        struct Threads{};
        struct Grain{};
        struct TwoD{};
        struct TwoDAdaptive{};
        struct ThreeD{};
        struct ThreeDInPlace{};
        struct ThreeDAdaptive{};
    }

    /*! ParSeqHelper for both fgemm and ftrsm
    */
    /*! ParSeqHelper for both fgemm and ftrsm
    */
    namespace ParSeqHelper {
        template <typename C=CuttingStrategy::Block, typename P=StrategyParameter::Threads>
        struct Parallel{
            typedef C Cut;
            typedef P Param;

            Parallel(size_t n=NUM_THREADS):_numthreads(n){}

            friend std::ostream& operator<<(std::ostream& out, const Parallel& p) {
                return out << "Parallel: " << p.numthreads();
            }
            size_t numthreads() const { return _numthreads; }
            size_t& set_numthreads(size_t n) { return _numthreads=n; }
            // CuttingStrategy method() const { return _method; }
            // StrategyParameter strategy() const { return _param; }
        private:
            size_t _numthreads;
            // CuttingStrategy _method;
            // StrategyParameter _param;

        };
        struct Sequential{
            Sequential() {}
            Sequential(size_t nth) { assert (nth == 1); }
            template<class Cut,class Param>
            Sequential(Parallel<Cut,Param>& ) {}
            friend std::ostream& operator<<(std::ostream& out, const Sequential&) {
                return out << "Sequential";
            }
            size_t numthreads() const { return 1; }
            // 	CuttingStrategy method() const { return SINGLE; }
            // // numthreads==1 ==> a single block
            // 	StrategyParameter strategy() const { return THREADS; }
        };
        template <typename H1=Sequential, typename H2=Sequential>
        struct Compose{

            Compose() : _comp1 (), _comp2 () {}
            Compose(const Compose & other) : _comp1 (other.first_component()), _comp2 (other.second_component()) {}
            Compose(const Sequential & S) : _comp1 (1), _comp2 (1) {}
            Compose(size_t th1, size_t th2) : _comp1 (th1), _comp2 (th2) {}
            Compose(const H1 & o1, const H2 & o2) : _comp1 (o1), _comp2 (o2) {}

            H1 first_component () const { return _comp1; }
            H2 second_component () const { return _comp2; }

            friend std::ostream& operator<<(std::ostream& o, const Compose& c) {
                return o << "Compose: (" << c.first_component() << ", "
                << c.second_component() << ")";
            }

        private:
            H1 _comp1;
            H2 _comp2;
        };
    }


    template<class Cut=CuttingStrategy::Block, class Strat=StrategyParameter::Threads>
    inline void BlockCuts(size_t& RBLOCKSIZE, size_t& CBLOCKSIZE,
                          const size_t m, const size_t n,
                          const size_t numthreads);

    template<>
    inline void BlockCuts<CuttingStrategy::Single,StrategyParameter::Threads>(size_t& RBLOCKSIZE,
                                                                              size_t& CBLOCKSIZE,
                                                                              const size_t m, const size_t n,
                                                                              const size_t numthreads) {
        assert(numthreads==1);
        RBLOCKSIZE = std::max(m,(size_t)1);
        CBLOCKSIZE = std::max(n,(size_t)1);
    }


    template<>
    inline void BlockCuts<CuttingStrategy::Row,StrategyParameter::Fixed>(size_t& RBLOCKSIZE,
                                                                         size_t& CBLOCKSIZE,
                                                                         const size_t m, const size_t n,
                                                                         const size_t numthreads) {
        RBLOCKSIZE = std::max(std::min(m,__FFLASFFPACK_MINBLOCKCUTS),(size_t)1);
        CBLOCKSIZE = std::max(n,(size_t)1);
    }


    template<>
    inline void BlockCuts<CuttingStrategy::Row,StrategyParameter::Grain>(size_t& RBLOCKSIZE,
                                                                         size_t& CBLOCKSIZE,
                                                                         const size_t m, const size_t n,
                                                                         const size_t grainsize) {
        RBLOCKSIZE = std::max(std::min(m,grainsize),(size_t)1);
        CBLOCKSIZE = std::max(n,(size_t)1);
    }

    template<>
    inline void BlockCuts<CuttingStrategy::Block,StrategyParameter::Grain>(size_t& RBLOCKSIZE,
                                                                           size_t& CBLOCKSIZE,
                                                                           const size_t m, const size_t n,
                                                                           const size_t grainsize) {
        RBLOCKSIZE = std::max(std::min(m,grainsize),(size_t)1);
        CBLOCKSIZE = std::max(std::min(n,grainsize),(size_t)1);
    }


    template<>
    inline void BlockCuts<CuttingStrategy::Column,StrategyParameter::Fixed>(size_t& RBLOCKSIZE,
                                                                            size_t& CBLOCKSIZE,
                                                                            const size_t m, const size_t n,
                                                                            const size_t numthreads) {
        RBLOCKSIZE = std::max(m,(size_t)1);
        CBLOCKSIZE = std::max(std::min(n,__FFLASFFPACK_MINBLOCKCUTS),(size_t)1);
    }


    template<>
    inline void BlockCuts<CuttingStrategy::Column,StrategyParameter::Grain>(size_t& RBLOCKSIZE,
                                                                            size_t& CBLOCKSIZE,
                                                                            const size_t m, const size_t n,
                                                                            const size_t grainsize) {
        RBLOCKSIZE = std::max(m,(size_t)1);
        CBLOCKSIZE = std::max(std::min(n,grainsize),(size_t)1);
    }


    template<>
    inline void BlockCuts<CuttingStrategy::Block,StrategyParameter::Fixed>(size_t& RBLOCKSIZE,
                                                                           size_t& CBLOCKSIZE,
                                                                           const size_t m, const size_t n,
                                                                           const size_t numthreads) {
        RBLOCKSIZE = std::max(std::min(m,__FFLASFFPACK_MINBLOCKCUTS),(size_t)1);
        CBLOCKSIZE = std::max(std::min(n,__FFLASFFPACK_MINBLOCKCUTS),(size_t)1);
    }

    template<>
    inline void BlockCuts<CuttingStrategy::Row,StrategyParameter::Threads>(size_t& RBLOCKSIZE,
                                                                           size_t& CBLOCKSIZE,
                                                                           const size_t m, const size_t n,
                                                                           const size_t numthreads) {
        RBLOCKSIZE = std::max(m/numthreads,(size_t)1);
        CBLOCKSIZE = std::max(n,(size_t)1);
    }


    template<>
    inline void BlockCuts<CuttingStrategy::Column,StrategyParameter::Threads>(size_t& RBLOCKSIZE,
                                                                              size_t& CBLOCKSIZE,
                                                                              const size_t m, const size_t n,
                                                                              const size_t numthreads) {
        RBLOCKSIZE = std::max(m,(size_t)1);
        CBLOCKSIZE = std::max(n/numthreads,(size_t)1);
    }

    template<>
    inline void BlockCuts<CuttingStrategy::Block,StrategyParameter::Threads>(size_t& RBLOCKSIZE,
                                                                             size_t& CBLOCKSIZE,
                                                                             const size_t m, const size_t n,
                                                                             const size_t numthreads) {
        if (numthreads<65) {
            //CP: Let's not compute these values all the time
            const short maxtc[64] = {1,2,3,2,5,3,7,4,3,5,11,4,13,7,5,4,17,6,19,5,7,11,23,6,5,13,9,7,29,6,31,8,11,17,7,6,37,19,13,8,41,7,43,11,9,23,47,8,7,10,17,13,53,9,11,8,19,29,59,10,61,31,9,8};
            const short maxtr[64] = {1,1,1,2,1,2,1,2,3,2,1,3,1,2,3,4,1,3,1,4,3,2,1,4,5,2,3,4,1,5,1,4,3,2,5,6,1,2,3,5,1,6,1,4,5,2,1,6,7,5,3,4,1,6,5,7,3,2,1,6,1,2,7,8};

            RBLOCKSIZE=std::max(m/(size_t)maxtr[numthreads-1],(size_t)1);
            CBLOCKSIZE=std::max(n/(size_t)maxtc[numthreads-1],(size_t)1);
        } else {
            const size_t maxt = (size_t)sqrt((double)numthreads);
            size_t maxtr=maxt,maxtc=maxt;
            for(size_t i=maxt; i>=1; --i) {
                size_t j=maxt;
                size_t newpr = i*j;
                for( ; newpr < numthreads; ++j, newpr+=i ) {}
                if (newpr == numthreads) {
                    maxtc = j;
                    maxtr = i;
                    break;
                }
            }
            RBLOCKSIZE=std::max(m/maxtr,(size_t)1);
            CBLOCKSIZE=std::max(n/maxtc,(size_t)1);
        }
    }

    // inline void BlockCuts(size_t& r, size_t& c,
    //                size_t m, size_t n,
    //                const CuttingStrategy method,
    //                const StrategyParameter strategy,
    //                const size_t t) {
    //     switch(method) {
    //         case CuttingStrategy::Block:
    // 		    switch(strategy) {
    // 			case StrategyParameter::Threads: BlockCuts<CuttingStrategy::Block,StrategyParameter::Threads>(r,c,m,n,t); break;
    // 			case StrategyParameter::Grain: BlockCuts<CuttingStrategy::Block,StrategyParameter::Grain>(r,c,m,n,t); break;
    // 			case StrategyParameter::Fixed: BlockCuts<CuttingStrategy::Block,StrategyParameter::Fixed>(r,c,m,n,t); break;
    // 			default: BlockCuts<CuttingStrategy::Block,StrategyParameter::Threads>(r,c,m,n,t);
    // 		    }
    // 		    break;
    //         case CuttingStrategy::Row:
    // 		    switch(strategy) {
    // 			case StrategyParameter::Threads: BlockCuts<CuttingStrategy::Row,StrategyParameter::Threads>(r,c,m,n,t); break;
    // 			case StrategyParameter::Grain: BlockCuts<CuttingStrategy::Row,StrategyParameter::Grain>(r,c,m,n,t); break;
    // 			case StrategyParameter::Fixed: BlockCuts<CuttingStrategy::Row,StrategyParameter::Fixed>(r,c,m,n,t); break;
    // 			default: BlockCuts<CuttingStrategy::Row,StrategyParameter::Threads>(r,c,m,n,t);
    // 		    }
    // 		    break;
    //         case CuttingStrategy::Column:
    // 		    switch(strategy) {
    // 			case StrategyParameter::Threads: BlockCuts<CuttingStrategy::Column,StrategyParameter::Threads>(r,c,m,n,t); break;
    // 			case StrategyParameter::Grain: BlockCuts<CuttingStrategy::Column,StrategyParameter::Grain>(r,c,m,n,t); break;
    // 			case StrategyParameter::Fixed: BlockCuts<CuttingStrategy::Column,StrategyParameter::Fixed>(r,c,m,n,t); break;
    // 			default: BlockCuts<CuttingStrategy::Column,StrategyParameter::Threads>(r,c,m,n,t);
    // 		    }
    // 		    break;
    //         default: BlockCuts<CuttingStrategy::Block,StrategyParameter::Threads>(r,c,m,n,t);
    //     };
    // }


    template<class Cut=CuttingStrategy::Block, class Param=StrategyParameter::Threads>
    inline void BlockCuts(size_t& rowBlockSize, size_t& colBlockSize,
                          size_t& lastRBS, size_t& lastCBS,
                          size_t& changeRBS, size_t& changeCBS,
                          size_t& numRowBlock, size_t& numColBlock,
                          size_t m, size_t n,
                          //                 const CuttingStrategy method,
                          //                 const StrategyParameter strategy,
                          const size_t numthreads) {
        BlockCuts<Cut,Param>(rowBlockSize, colBlockSize, m, n, numthreads);
        numRowBlock = m/rowBlockSize;
        numColBlock = n/colBlockSize;

        changeRBS = m-rowBlockSize*numRowBlock;
        lastRBS = rowBlockSize;
        if (changeRBS) ++rowBlockSize;

        changeCBS = n-colBlockSize*numColBlock;
        lastCBS = colBlockSize;
        if (changeCBS) ++colBlockSize;


        /*
        // Better preserve numRowBlock and numColBlock
        if (lastRBS) {
        lastRBS = m-rowBlockSize*numRowBlock;
        ++rowBlockSize;
        } else lastRBS = rowBlockSize;
        if (lastCBS) {
        lastCBS = n-colBlockSize*numColBlock;
        ++colBlockSize;
        } else lastCBS = colBlockSize;
        */



        //             // Better preserve rowBlockSize and colBlockSize
        //         lastRBS = m % rowBlockSize;
        //         lastCBS = n % colBlockSize;
        //          if (lastRBS) ++numRowBlock; else lastRBS = rowBlockSize;
        //          if (lastCBS) ++numColBlock; else lastCBS = colBlockSize;
    }


}



namespace FFLAS {
    template <typename blocksize_t=size_t, typename Cut=CuttingStrategy::Block, typename Param=StrategyParameter::Threads>
    struct ForStrategy1D {
        ForStrategy1D(const blocksize_t n, const ParSeqHelper::Parallel<Cut,Param> H) {
            build(n,H);
        }
        ForStrategy1D(const blocksize_t b, const blocksize_t e, const ParSeqHelper::Parallel<Cut,Param> H) {
            build(e-b,H);
        }

        void build(const blocksize_t n, const ParSeqHelper::Parallel<Cut,Param> H) {

            if ( Protected::AreEqual<Param, StrategyParameter::Threads>::value ) {
                numBlock = std::max((blocksize_t)(H.numthreads()),(blocksize_t)1);
            } else if ( Protected::AreEqual<Param,StrategyParameter::Grain>::value ) {
                numBlock = std::max(n/ (blocksize_t)(H.numthreads()), (blocksize_t)1);
            } else {
                numBlock = std::max(n/(blocksize_t)(__FFLASFFPACK_MINBLOCKCUTS),(blocksize_t)1);
            }
            firstBlockSize = n/numBlock;
            if (firstBlockSize<1) {
                firstBlockSize = (blocksize_t)1;
                numBlock = n;
            }
            changeBS = n - numBlock*firstBlockSize;
            lastBlockSize = firstBlockSize;
            if (changeBS) ++firstBlockSize;
        }

        blocksize_t initialize() {
            ibeg = 0; iend = firstBlockSize;
            return current = 0;
        }
        bool isTerminated() const { return current == numBlock; }

        blocksize_t begin() const { return ibeg; }
        blocksize_t end() const { return iend; }

        blocksize_t numblocks() const { return numBlock; }
        blocksize_t blockindex() const { return current; }


        blocksize_t operator++() {
            ibeg = iend;
            iend += (++current<changeBS?firstBlockSize:lastBlockSize);
            return current;
        }

    protected:
        blocksize_t ibeg, iend;

        blocksize_t current;
        blocksize_t firstBlockSize,lastBlockSize;
        blocksize_t changeBS;
        blocksize_t numBlock;

    };

    template <typename blocksize_t=size_t, typename Cut=CuttingStrategy::Block, typename Param=StrategyParameter::Threads>
    struct ForStrategy2D {
        ForStrategy2D(const blocksize_t m, const blocksize_t n, const ParSeqHelper::Parallel<Cut,Param> H)
        {
            BlockCuts<Cut,Param>(rowBlockSize, colBlockSize,
                                 lastRBS, lastCBS,
                                 changeRBS, changeCBS,
                                 numRowBlock, numColBlock,
                                 m, n,
                                 //						     H.method(), H.strategy(),
                                 H.numthreads());

            BLOCKS = numRowBlock * numColBlock;
        }


        blocksize_t initialize() {
            _ibeg = 0; _iend = rowBlockSize;
            _jbeg = 0; _jend = colBlockSize;
            return current = 0;
        }
        bool isTerminated() const { return current == BLOCKS; }

        blocksize_t ibegin() const { return _ibeg; }
        blocksize_t jbegin() const { return _jbeg; }
        blocksize_t iend() const { return _iend; }
        blocksize_t jend() const { return _jend; }


        blocksize_t operator++() {
            ++current;
            blocksize_t icurr = current/numColBlock;
            blocksize_t jcurr = current%numColBlock;
            if (jcurr) {
                _jbeg = _jend;
                _jend += (jcurr<changeCBS?colBlockSize:lastCBS);
            } else {
                _ibeg = _iend;
                _iend += (icurr<changeRBS?rowBlockSize:lastRBS);
                _jbeg = 0;
                _jend = colBlockSize;
            }
            return current;
        }

        friend std::ostream& operator<<(std::ostream& out, const ForStrategy2D& FS2D) {
            out<<"RBLOCKSIZE: "<<FS2D.rowBlockSize<<std::endl;
            out<<"CBLOCKSIZE: "<<FS2D.colBlockSize<<std::endl;
            out<<"changeRBS : "<<FS2D.changeRBS<<std::endl;
            out<<"changeCBS : "<<FS2D.changeCBS<<std::endl;
            out<<"lastRBS   : "<<FS2D.lastRBS<<std::endl;
            out<<"lastCBS   : "<<FS2D.lastCBS<<std::endl;
            out<<"NrowBlocks: "<<FS2D.numRowBlock<<std::endl;
            out<<"NcolBlocks: "<<FS2D.numColBlock<<std::endl;
            out<<"curr: " << FS2D.current << '/' << FS2D.BLOCKS << std::endl;
            out<<"_ibeg: " << FS2D._ibeg << std::endl;
            out<<"_iend: " << FS2D._iend << std::endl;
            out<<"_jbeg: " << FS2D._jbeg << std::endl;
            out<<"_jend: " << FS2D._jend << std::endl;
            return out;
        }

        blocksize_t rownumblocks() const { return numRowBlock; }
        blocksize_t colnumblocks() const { return numColBlock; }
        blocksize_t blockindex() const { return current; }
        blocksize_t rowblockindex() const { return blockindex()/colnumblocks(); }
        blocksize_t colblockindex() const { return blockindex()%colnumblocks(); }


    protected:
        blocksize_t _ibeg, _iend, _jbeg, _jend;
        blocksize_t rowBlockSize, colBlockSize;

        blocksize_t current;
        blocksize_t lastRBS; blocksize_t lastCBS;
        blocksize_t changeRBS; blocksize_t changeCBS;
        blocksize_t numRowBlock; blocksize_t numColBlock;
        blocksize_t BLOCKS;

    };

}



#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
