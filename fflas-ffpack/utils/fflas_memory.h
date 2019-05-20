/* fflas/fflas_memory.h
 * Copyright (C) 2014 fflas-ffpack group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * The cache size detection has been copied from the Eigen library,
 * a lightweight C++ template library for linear algebra, licensed under
 * the Mozilla
 * Public License v. 2.0. If a copy of the MPL was not distributed
 * with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
 * Copyright (C) 2008-2010 Gael Guennebaud <gael.guennebaud@inria.fr>
 * Copyright (C) 2008-2009 Benoit Jacob <jacob.benoit.1@gmail.com>
 * Copyright (C) 2009 Kenneth Riddile <kfriddile@yahoo.com>
 * Copyright (C) 2010 Hauke Heibel <hauke.heibel@gmail.com>
 * Copyright (C) 2010 Thomas Capricelli <orzel@freehackers.org>
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

#ifndef __FFLASFFPACK_memory_H
#define __FFLASFFPACK_memory_H

#include "fflas-ffpack/utils/align-allocator.h"
#include <givaro/givinteger.h>

namespace FFLAS{

    template<class Element>
    inline bool alignable() {
        return true ;
    }

    // BB : segfault in Givaro::Integer::logcpy otherwise
    template<>
    inline bool alignable<Givaro::Integer*>() {
        return false;
    }

    template<class Field>
    inline typename Field::Element_ptr fflas_new (const Field& F, const size_t m, const Alignment align = Alignment::DEFAULT)
    {
        if (alignable<typename Field::Element_ptr>() ) {
            return malloc_align<typename Field::Element>(m, align);
        }
        else {
            return new typename Field::Element[m];
        }
    }

    template<class Field>
    inline typename Field::Element_ptr fflas_new (const Field& F, const size_t m, const size_t n, const Alignment align = Alignment::DEFAULT)
    {
        return fflas_new(F, m*n, align);
    }

    template<class Element >
    inline Element* fflas_new (const size_t m, const Alignment align = Alignment::DEFAULT)
    {
        if (alignable<Element*>() ) {
            return malloc_align<Element>(m, align);
        }
        else {
            return new Element[m];
        }

    }

    template<class Element_ptr>
    inline void fflas_delete(Element_ptr A)
    {
        if (alignable<Element_ptr>() )
            free(A);
        else
            delete[] A;
    }

    template<class Ptr, class ...Args>
    inline void fflas_delete(Ptr p, Args ... args){
        fflas_delete(p);
        fflas_delete(std::forward<Args>(args)...);
    }

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
    inline void prefetch(const int64_t* addr) { _mm_prefetch((const char*)(addr), _MM_HINT_T0); }
#else
    inline void prefetch(const int64_t*) {}
#endif

#if ( defined(__i386__) || defined(__x86_64__) )
#define __CPUID(abcd,func,id) \
    __asm__ __volatile__ ("cpuid": "=a" (abcd[0]), "=b" (abcd[1]), "=c" (abcd[2]), "=d" (abcd[3]) : "a" (func), "c" (id) );

    inline void getCacheSize(int& l1, int& l2, int& l3)
    {
        int abcd[4];
        l1 = l2 = l3 = 0;
        int cache_id = 0;
        int cache_type = 0;
        do {
            abcd[0] = abcd[1] = abcd[2] = abcd[3] = 0;
            __CPUID(abcd,0x4,cache_id);
            cache_type  = (abcd[0] & 0x0F) >> 0;
            if(cache_type==1||cache_type==3) // data or unified cache
            {
                int cache_level = (abcd[0] & 0xE0) >> 5;  // A[7:5]
                int ways        = (abcd[1] & 0xFFC00000) >> 22; // B[31:22]
                int partitions  = (abcd[1] & 0x003FF000) >> 12; // B[21:12]
                int line_size   = (abcd[1] & 0x00000FFF) >>  0; // B[11:0]
                int sets        = (abcd[2]);                    // C[31:0]
                int cache_size = (ways+1) * (partitions+1) * (line_size+1) * (sets+1);
                switch(cache_level)
                {
                case 1: l1 = cache_size; break;
                case 2: l2 = cache_size; break;
                case 3: l3 = cache_size; break;
                default: break;
                }
            }
            cache_id++;
        } while(cache_type>0 && cache_id<16);
    }

    inline void getTLBSize(int& tlb)
    {
        int abcd[4]={};
        int sTLB=0;
        int lTLB;
        __CPUID(abcd,0x2,0);
        unsigned char * bytes = reinterpret_cast<unsigned char *>(abcd)+2;
        for(int i=0; i<14; ++i)
            switch(bytes[i]){
            case 0x03: sTLB=64; break;
            case 0x04: lTLB=8;  break;
            case 0x05: lTLB=32; break;
            case 0x56: lTLB=16; break;
            case 0x57: sTLB=16; break;
            case 0x59: sTLB=16; break;
            case 0x5A: lTLB=32; break;
            case 0x5B: sTLB=lTLB=64;  break;
            case 0x5C: sTLB=lTLB=128; break;
            case 0x5D: sTLB=lTLB=256; break;
            case 0xB3: sTLB=128; break;
            case 0xB4: sTLB=256; break;
            case 0xBA: sTLB=64;  break;
            case 0xC0: sTLB=lTLB=8; break;
            case 0xCA: sTLB=512; break;
            default: break;
            }
        //cout<<"small TLB: "<<sTLB<<endl;
        //cout<<"large TLB: "<<lTLB<<endl;
        tlb=sTLB*4096;
    }
#else
    inline void getTLBSize(int& tlb){tlb = 0;} // not implemented in non x86 archs
#endif
    //---------- Cache sizes ----------

#if !defined(EIGEN_NO_CPUID)
#  if defined(__GNUC__) && ( defined(__i386__) || defined(__x86_64__) )
#    if defined(__PIC__) && defined(__i386__)
    // Case for x86 with PIC
#      define EIGEN_CPUID(abcd,func,id) \
    __asm__ __volatile__ ("xchgl %%ebx, %%esi;cpuid; xchgl %%ebx,%%esi": "=a" (abcd[0]), "=S" (abcd[1]), "=c" (abcd[2]), "=d" (abcd[3]) : "a" (func), "c" (id));
#    else
    // Case for x86_64 or x86 w/o PIC
#      define EIGEN_CPUID(abcd,func,id) \
    __asm__ __volatile__ ("cpuid": "=a" (abcd[0]), "=b" (abcd[1]), "=c" (abcd[2]), "=d" (abcd[3]) : "a" (func), "c" (id) );
#    endif
#  elif defined(_MSC_VER)
#    if (_MSC_VER > 1500) && ( defined(_M_IX86) || defined(_M_X64) )
#      define EIGEN_CPUID(abcd,func,id) __cpuidex((int*)abcd,func,id)
#    endif
#  endif
#endif


#ifdef EIGEN_CPUID

    inline bool cpuid_is_vendor(int abcd[4], const char* vendor)
    {
        return abcd[1]==(reinterpret_cast<const int*>(vendor))[0] && abcd[3]==(reinterpret_cast<const int*>(vendor))[1] && abcd[2]==(reinterpret_cast<const int*>(vendor))[2];
    }

    inline void queryCacheSizes_intel_direct(int& l1, int& l2, int& l3)
    {
        int abcd[4];
        l1 = l2 = l3 = 0;
        int cache_id = 0;
        int cache_type = 0;
        do {
            abcd[0] = abcd[1] = abcd[2] = abcd[3] = 0;
            EIGEN_CPUID(abcd,0x4,cache_id);
            cache_type  = (abcd[0] & 0x0F) >> 0;
            if(cache_type==1||cache_type==3) // data or unified cache
            {
                int cache_level = (abcd[0] & 0xE0) >> 5;  // A[7:5]
                int ways        = (abcd[1] & 0xFFC00000) >> 22; // B[31:22]
                int partitions  = (abcd[1] & 0x003FF000) >> 12; // B[21:12]
                int line_size   = (abcd[1] & 0x00000FFF) >>  0; // B[11:0]
                int sets        = (abcd[2]);                    // C[31:0]

                int cache_size = (ways+1) * (partitions+1) * (line_size+1) * (sets+1);

                switch(cache_level)
                {
                case 1: l1 = cache_size; break;
                case 2: l2 = cache_size; break;
                case 3: l3 = cache_size; break;
                default: break;
                }
            }
            cache_id++;
        } while(cache_type>0 && cache_id<16);
    }

    inline void queryCacheSizes_intel_codes(int& l1, int& l2, int& l3)
    {
        int abcd[4];
        abcd[0] = abcd[1] = abcd[2] = abcd[3] = 0;
        l1 = l2 = l3 = 0;
        EIGEN_CPUID(abcd,0x00000002,0);
        unsigned char * bytes = reinterpret_cast<unsigned char *>(abcd)+2;
        bool check_for_p2_core2 = false;
        for(int i=0; i<14; ++i)
        {
            switch(bytes[i])
            {
            case 0x0A: l1 = 8; break;   // 0Ah   data L1 cache, 8 KB, 2 ways, 32 byte lines
            case 0x0C: l1 = 16; break;  // 0Ch   data L1 cache, 16 KB, 4 ways, 32 byte lines
            case 0x0E: l1 = 24; break;  // 0Eh   data L1 cache, 24 KB, 6 ways, 64 byte lines
            case 0x10: l1 = 16; break;  // 10h   data L1 cache, 16 KB, 4 ways, 32 byte lines (IA-64)
            case 0x15: l1 = 16; break;  // 15h   code L1 cache, 16 KB, 4 ways, 32 byte lines (IA-64)
            case 0x2C: l1 = 32; break;  // 2Ch   data L1 cache, 32 KB, 8 ways, 64 byte lines
            case 0x30: l1 = 32; break;  // 30h   code L1 cache, 32 KB, 8 ways, 64 byte lines
            case 0x60: l1 = 16; break;  // 60h   data L1 cache, 16 KB, 8 ways, 64 byte lines, sectored
            case 0x66: l1 = 8; break;   // 66h   data L1 cache, 8 KB, 4 ways, 64 byte lines, sectored
            case 0x67: l1 = 16; break;  // 67h   data L1 cache, 16 KB, 4 ways, 64 byte lines, sectored
            case 0x68: l1 = 32; break;  // 68h   data L1 cache, 32 KB, 4 ways, 64 byte lines, sectored
            case 0x1A: l2 = 96; break;   // code and data L2 cache, 96 KB, 6 ways, 64 byte lines (IA-64)
            case 0x22: l3 = 512; break;   // code and data L3 cache, 512 KB, 4 ways (!), 64 byte lines, dual-sectored
            case 0x23: l3 = 1024; break;   // code and data L3 cache, 1024 KB, 8 ways, 64 byte lines, dual-sectored
            case 0x25: l3 = 2048; break;   // code and data L3 cache, 2048 KB, 8 ways, 64 byte lines, dual-sectored
            case 0x29: l3 = 4096; break;   // code and data L3 cache, 4096 KB, 8 ways, 64 byte lines, dual-sectored
            case 0x39: l2 = 128; break;   // code and data L2 cache, 128 KB, 4 ways, 64 byte lines, sectored
            case 0x3A: l2 = 192; break;   // code and data L2 cache, 192 KB, 6 ways, 64 byte lines, sectored
            case 0x3B: l2 = 128; break;   // code and data L2 cache, 128 KB, 2 ways, 64 byte lines, sectored
            case 0x3C: l2 = 256; break;   // code and data L2 cache, 256 KB, 4 ways, 64 byte lines, sectored
            case 0x3D: l2 = 384; break;   // code and data L2 cache, 384 KB, 6 ways, 64 byte lines, sectored
            case 0x3E: l2 = 512; break;   // code and data L2 cache, 512 KB, 4 ways, 64 byte lines, sectored
            case 0x40: l2 = 0; break;   // no integrated L2 cache (P6 core) or L3 cache (P4 core)
            case 0x41: l2 = 128; break;   // code and data L2 cache, 128 KB, 4 ways, 32 byte lines
            case 0x42: l2 = 256; break;   // code and data L2 cache, 256 KB, 4 ways, 32 byte lines
            case 0x43: l2 = 512; break;   // code and data L2 cache, 512 KB, 4 ways, 32 byte lines
            case 0x44: l2 = 1024; break;   // code and data L2 cache, 1024 KB, 4 ways, 32 byte lines
            case 0x45: l2 = 2048; break;   // code and data L2 cache, 2048 KB, 4 ways, 32 byte lines
            case 0x46: l3 = 4096; break;   // code and data L3 cache, 4096 KB, 4 ways, 64 byte lines
            case 0x47: l3 = 8192; break;   // code and data L3 cache, 8192 KB, 8 ways, 64 byte lines
            case 0x48: l2 = 3072; break;   // code and data L2 cache, 3072 KB, 12 ways, 64 byte lines
            case 0x49: if(l2!=0) l3 = 4096; else {check_for_p2_core2=true; l3 = l2 = 4096;} break;// code and data L3 cache, 4096 KB, 16 ways, 64 byte lines (P4) or L2 for core2
            case 0x4A: l3 = 6144; break;   // code and data L3 cache, 6144 KB, 12 ways, 64 byte lines
            case 0x4B: l3 = 8192; break;   // code and data L3 cache, 8192 KB, 16 ways, 64 byte lines
            case 0x4C: l3 = 12288; break;   // code and data L3 cache, 12288 KB, 12 ways, 64 byte lines
            case 0x4D: l3 = 16384; break;   // code and data L3 cache, 16384 KB, 16 ways, 64 byte lines
            case 0x4E: l2 = 6144; break;   // code and data L2 cache, 6144 KB, 24 ways, 64 byte lines
            case 0x78: l2 = 1024; break;   // code and data L2 cache, 1024 KB, 4 ways, 64 byte lines
            case 0x79: l2 = 128; break;   // code and data L2 cache, 128 KB, 8 ways, 64 byte lines, dual-sectored
            case 0x7A: l2 = 256; break;   // code and data L2 cache, 256 KB, 8 ways, 64 byte lines, dual-sectored
            case 0x7B: l2 = 512; break;   // code and data L2 cache, 512 KB, 8 ways, 64 byte lines, dual-sectored
            case 0x7C: l2 = 1024; break;   // code and data L2 cache, 1024 KB, 8 ways, 64 byte lines, dual-sectored
            case 0x7D: l2 = 2048; break;   // code and data L2 cache, 2048 KB, 8 ways, 64 byte lines
            case 0x7E: l2 = 256; break;   // code and data L2 cache, 256 KB, 8 ways, 128 byte lines, sect. (IA-64)
            case 0x7F: l2 = 512; break;   // code and data L2 cache, 512 KB, 2 ways, 64 byte lines
            case 0x80: l2 = 512; break;   // code and data L2 cache, 512 KB, 8 ways, 64 byte lines
            case 0x81: l2 = 128; break;   // code and data L2 cache, 128 KB, 8 ways, 32 byte lines
            case 0x82: l2 = 256; break;   // code and data L2 cache, 256 KB, 8 ways, 32 byte lines
            case 0x83: l2 = 512; break;   // code and data L2 cache, 512 KB, 8 ways, 32 byte lines
            case 0x84: l2 = 1024; break;   // code and data L2 cache, 1024 KB, 8 ways, 32 byte lines
            case 0x85: l2 = 2048; break;   // code and data L2 cache, 2048 KB, 8 ways, 32 byte lines
            case 0x86: l2 = 512; break;   // code and data L2 cache, 512 KB, 4 ways, 64 byte lines
            case 0x87: l2 = 1024; break;   // code and data L2 cache, 1024 KB, 8 ways, 64 byte lines
            case 0x88: l3 = 2048; break;   // code and data L3 cache, 2048 KB, 4 ways, 64 byte lines (IA-64)
            case 0x89: l3 = 4096; break;   // code and data L3 cache, 4096 KB, 4 ways, 64 byte lines (IA-64)
            case 0x8A: l3 = 8192; break;   // code and data L3 cache, 8192 KB, 4 ways, 64 byte lines (IA-64)
            case 0x8D: l3 = 3072; break;   // code and data L3 cache, 3072 KB, 12 ways, 128 byte lines (IA-64)

            default: break;
            }
        }
        if(check_for_p2_core2 && l2 == l3)
            l3 = 0;
        l1 *= 1024;
        l2 *= 1024;
        l3 *= 1024;
    }

    inline void queryCacheSizes_intel(int& l1, int& l2, int& l3, int max_std_funcs)
    {
        if(max_std_funcs>=4)
            queryCacheSizes_intel_direct(l1,l2,l3);
        else
            queryCacheSizes_intel_codes(l1,l2,l3);
    }

    inline void queryCacheSizes_amd(int& l1, int& l2, int& l3)
    {
        int abcd[4];
        abcd[0] = abcd[1] = abcd[2] = abcd[3] = 0;
        EIGEN_CPUID(abcd,0x80000005,0);
        l1 = (abcd[2] >> 24) * 1024; // C[31:24] = L1 size in KB
        abcd[0] = abcd[1] = abcd[2] = abcd[3] = 0;
        EIGEN_CPUID(abcd,0x80000006,0);
        l2 = (abcd[2] >> 16) * 1024; // C[31;16] = l2 cache size in KB
        l3 = ((abcd[3] & 0xFFFC000) >> 18) * 512 * 1024; // D[31;18] = l3 cache size in 512KB
    }
#endif

    /** \internal
     * Queries and returns the cache sizes in Bytes of the L1, L2, and L3 data caches respectively */
    inline void queryCacheSizes(int& l1, int& l2, int& l3)
    {
#ifdef EIGEN_CPUID
        int abcd[4];

        // identify the CPU vendor
        EIGEN_CPUID(abcd,0x0,0);
        int max_std_funcs = abcd[1];
        if(cpuid_is_vendor(abcd,"GenuineIntel"))
            queryCacheSizes_intel(l1,l2,l3,max_std_funcs);
        else if(cpuid_is_vendor(abcd,"AuthenticAMD") || cpuid_is_vendor(abcd,"AMDisbetter!")
          || cpuid_is_vendor(abcd, "HygonGenuine"))
            queryCacheSizes_amd(l1,l2,l3);
        else
            // by default let's use Intel's API
            queryCacheSizes_intel(l1,l2,l3,max_std_funcs);

        // here is the list of other vendors:
        //   ||cpuid_is_vendor(abcd,"VIA VIA VIA ")
        //   ||cpuid_is_vendor(abcd,"CyrixInstead")
        //   ||cpuid_is_vendor(abcd,"CentaurHauls")
        //   ||cpuid_is_vendor(abcd,"GenuineTMx86")
        //   ||cpuid_is_vendor(abcd,"TransmetaCPU")
        //   ||cpuid_is_vendor(abcd,"RiseRiseRise")
        //   ||cpuid_is_vendor(abcd,"Geode by NSC")
        //   ||cpuid_is_vendor(abcd,"SiS SiS SiS ")
        //   ||cpuid_is_vendor(abcd,"UMC UMC UMC ")
        //   ||cpuid_is_vendor(abcd,"NexGenDriven")
#else
        l1 = l2 = l3 = -1;
#endif
    }

    /** \internal
     * \returns the size in Bytes of the L1 data cache */
    inline int queryL1CacheSize()
    {
        int l1(-1), l2, l3;
        queryCacheSizes(l1,l2,l3);
        return l1;
    }

    /** \internal
     * \returns the size in Bytes of the L2 or L3 cache if this later is present */
    inline int queryTopLevelCacheSize()
    {
        int l1, l2(-1), l3(-1);
        queryCacheSizes(l1,l2,l3);
        return (std::max)(l2,l3);
    }

} // namespace FFLAS
#endif // __FFLASFFPACK_memory_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
