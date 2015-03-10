/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/*
 * Copyright (C) 2013  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 * the code is inspired and adapted from the Eigen library
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

#include <immintrin.h>
#include <stdlib.h>

#include "fflas-ffpack/fflas/fflas_igebp.h"


template<size_t k>
void pack_lhs(int64_t* XX, const int64_t* X, size_t ldx, size_t rows, size_t cols);

template<size_t k>
void pack_rhs(int64_t* XX, const int64_t* X, size_t ldx, size_t rows, size_t cols);

void gebp(size_t rows, size_t cols, size_t depth,int64_t* C, size_t ldc, const int64_t* blockA, size_t lda,
	  const int64_t* BlockB, size_t ldb, int64_t* BlockW);

void BlockingFactor(size_t& m, size_t& n, size_t& k);


// Assume matrices A,B,C are stored in column major order
void igemm_colmajor(size_t rows, size_t cols, size_t depth, int64_t* C, size_t ldc, int64_t* A, size_t lda, int64_t* B, size_t ldb){

	size_t mc,kc,nc;
	mc=rows;
	nc=cols;
	kc=depth;
	BlockingFactor(mc,nc,kc);
	size_t sizeA = mc*kc;
	size_t sizeB = kc*cols;
	size_t sizeW = _vect_s*kc*_nr; // store data duplicated by the number of elements fitting in vector register

	// these data must be (_vect_align) byte aligned
	int64_t *blockA, *blockB, *blockW;
	posix_memalign(reinterpret_cast<void**>(&blockA),_vect_align, sizeof(int64_t)*sizeA);;
	posix_memalign(reinterpret_cast<void**>(&blockB),_vect_align, sizeof(int64_t)*sizeB);;
	posix_memalign(reinterpret_cast<void**>(&blockW),_vect_align, sizeof(int64_t)*sizeW);;


	// For each horizontal panel of B, and corresponding vertical panel of A
	for(size_t k2=0; k2<depth; k2+=kc){

		const size_t actual_kc = std::min(k2+kc,depth)-k2;

		// pack horizontal panel of B into sequential memory (L2 cache)
		pack_rhs<_nr>(blockB, B+k2, ldc, actual_kc, cols);

		// For each mc x kc block of the lhs's vertical panel...
		for(size_t i2=0; i2<rows; i2+=mc){

			const size_t actual_mc = std::min(i2+mc,rows)-i2;

			//cout<<"mc= "<<actual_mc<<" kc= "<<actual_mc<<endl;

			// pack a chunk of the vertical panel of A into a sequential memory (L1 cache)
			pack_lhs<_mr>(blockA, A+i2+k2*lda, lda, actual_mc, actual_kc);
			// call block*panel kernel
			igebp(actual_mc, cols, actual_kc, C+i2, ldc, blockA, actual_kc, blockB, actual_kc, blockW);
		}
	}
	free(blockA);
	free(blockB);
	free(blockW);
}


#define SSE_LOADU(A,B)  A=_mm_loadu_si128(reinterpret_cast<const __m128i*>(&B));
#define SSE_STORE(A,B)    _mm_store_si128(reinterpret_cast<__m128i*>(&A),B);


// store each rows x k submatrices of Rhs in row major mode
// if k does not divide cols, the remaining column are not packed
template<size_t k>
void pack_rhs(int64_t* XX, const int64_t* X, size_t ldx, size_t rows, size_t cols){
	size_t cols_by_k=(cols/k)*k;
	size_t p=0;
	// pack by k columns
	for(size_t j=0;j<cols_by_k;j+=k){
		for(size_t i=0;i<rows;i++)
			for (size_t l=0;l<k;l++,p++)
				XX[p]=X[i+(j+l)*ldx];
	}
	// the remaining columns are not packed
	for(size_t j=cols_by_k;j<cols;j++)
		for(size_t i=0;i<rows;i++,p++)
			XX[p]=X[i+j*ldx];
}


// store each k x cols submatrices of Lhs in column major mode
// if k does not divide rows, the remaining rows are not packed
template<size_t k>
void pack_lhs(int64_t* XX, const int64_t* X, size_t ldx, size_t rows, size_t cols){
	size_t p=0;
	size_t rows_by_k=(rows/k)*k;
	// pack rows by group of k
	for(size_t i=0;i<rows_by_k;i+=k)
		for(size_t j=0;j<cols;j++)
			//for (size_t l=0;l<k;l++,p++) XX[p]=X[i+l+j*ldx];
				for (size_t l=0;l<k;l+=4){
				__m128i T0,T1;
				SSE_LOADU(T0,X[i+l+  j*ldx]);
				SSE_LOADU(T1,X[i+l+2+j*ldx]);
				SSE_STORE(XX[p],T0);p+=2;
				SSE_STORE(XX[p],T1);p+=2;
			}
	// the remaining rows are packed by group of StepA (if possible)
	if (rows-rows_by_k>=StepA){
		for(size_t j=0;j<cols;j++)
			for (size_t l=0;l<StepA;l+=2){
				__m128i T0;
				SSE_LOADU(T0,X[rows_by_k+l+j*ldx]);
				SSE_STORE(XX[p],T0);p+=2;
			}
		rows_by_k+=StepA;
	}
	for(size_t i=rows_by_k;i<rows;i++)
		for(size_t j=0;j<cols;j++,p++){
			XX[p]=X[i+j*ldx];
		}

}

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

void BlockingFactor(size_t& m, size_t& n, size_t& k)
{
	int l1, l2, l3, tlb;
	getCacheSize(l1,l2,l3);
	getTLBSize(tlb);
	/*
	cout<<"Cache size: ";
	cout<<"L1 ("<<l1<<") ";
	cout<<"L2 ("<<l2<<") ";
	cout<<"L3 ("<<l3<<") ";
	cout<<"TLB ("<<tlb<<") ";
	cout<<endl;
	*/
	l2=std::max(l2,l3);
	l2=std::min(l2,tlb);
	size_t kc,mc;
	// kc * 2*(_mr+_nr) must fit in L1 cache
	// kc * (n+mc) must fit in L2 cache and in TLB
	size_t kdiv= 2*(_nr+_mr)*sizeof(int64_t);
	kc = std::min(k, l1/kdiv);
	mc = std::min(m, l2/(sizeof(int64_t) * kc));
	k=kc;
       	m=mc;
	//cout<<"mc="<<m<<endl;
	//cout<<"kc="<<k<<endl;
}


void igemm(size_t rows, size_t cols, size_t depth, int64_t* C, size_t ldc, int64_t* A, size_t lda, int64_t* B, size_t ldb){
	igemm_colmajor(cols, rows, depth, C, ldc, B, ldb, A, lda);
}

