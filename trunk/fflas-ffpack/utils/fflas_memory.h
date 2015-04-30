/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_memory.h
 * Copyright (C) 2014 fflas-ffpack group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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
    inline typename Field::Element_ptr fflas_new (const Field& F, const size_t m, const size_t n, const Alignment align = Alignment::DEFAULT)
    {
	    if (alignable<typename Field::Element_ptr>() ) {
	    	    return malloc_align<typename Field::Element>(m*n, align);
	    }
	    else {
	    	    return new typename Field::Element[m*n];
	    }
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

#ifdef __FFLASFFPACK_USE_SIMD
    inline void prefetch(const int64_t* addr) { _mm_prefetch((const char*)(addr), _MM_HINT_T0); }
#endif


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
	    int sTLB, lTLB;
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

} // namespace FFLAS
#endif // __FFLASFFPACK_memory_H
