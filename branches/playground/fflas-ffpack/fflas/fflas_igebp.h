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

#ifdef __AVX2__
#define _vect_t __m256i
#define _vect_align 32
#define _vect_s 4
#define _nr 4
#define _mr 8
#define StepA 4
#define StepB 4
#define VEC_MADD_32(C,A,B,X)  X= _mm256_mul_epi32(A,B); C= _mm256_add_epi64(C,X);
#define VEC_ADDIN_64(A,B)     A= _mm256_add_epi64(A,B);
#define VEC_LOAD(C,X)         C= _mm256_load_si256(reinterpret_cast<const __m256i*>(X));
#define VEC_LOADU(C,X)        C= _mm256_loadu_si256(reinterpret_cast<const __m256i*>(X));
#define VEC_STORE(C,X)        _mm256_store_si256(reinterpret_cast<__m256i*>(C),X);
#define VEC_STOREU(C,X)       _mm256_storeu_si256(reinterpret_cast<__m256i*>(C),X);
#define VEC_SET_64(C,X)	      C= _mm256_set1_epi64x(X);
#else
#define _vect_t __m128i
#define _vect_align 16
#define _vect_s 2
#define _nr 4
#define _mr 4
#define StepA 2
#define StepB 2
#define VEC_MADD_32(C,A,B,X)  X= _mm_mul_epi32(A,B); C= _mm_add_epi64(C,X);
#define VEC_ADDIN_64(A,B)     A= _mm_add_epi64(A,B);
#define VEC_LOAD(C,X)         C= _mm_load_si128(reinterpret_cast<const __m128i*>(X));
#define VEC_LOADU(C,X)        C= _mm_loadu_si128(reinterpret_cast<const __m128i*>(X));
#define VEC_STORE(C,X)        _mm_store_si128(reinterpret_cast<__m128i*>(C),X);
#define VEC_STOREU(C,X)       _mm_storeu_si128(reinterpret_cast<__m128i*>(C),X);
#define VEC_SET_64(C,X)	      C= _mm_set1_epi64x(X);

#endif



// I ASSUME _vect_s <= 4 (meaning that we store at most 4 element per vector register)
// -> this is to explicitely treat the dimension that are not divisible by _vect_s


inline void igebb44(size_t i, size_t j, size_t depth, size_t pdeth, int64_t* C, size_t ldc, const int64_t *blA, const int64_t* blB);
inline void igebb24(size_t i, size_t j, size_t depth, size_t pdeth, int64_t* C, size_t ldc, const int64_t *blA, const int64_t* blB);
inline void igebb14(size_t i, size_t j, size_t depth, size_t pdeth, int64_t* C, size_t ldc, const int64_t *blA, const int64_t* blB);
inline void igebb41(size_t i, size_t j, size_t depth, size_t pdeth, int64_t* C, size_t ldc, const int64_t *blA, const int64_t* blB);
inline void igebb21(size_t i, size_t j, size_t depth, size_t pdeth, int64_t* C, size_t ldc, const int64_t *blA, const int64_t* blB);
inline void igebb11(size_t i, size_t j, size_t depth, size_t pdeth, int64_t* C, size_t ldc, const int64_t *blA, const int64_t* blB);

void prefetch(const int64_t* addr) { _mm_prefetch((const char*)(addr), _MM_HINT_T0); }

// duplicate each entry into vector register
template<size_t N>
inline void duplicate_vect (int64_t* XX, const int64_t* X, size_t n){}

template<>
inline void duplicate_vect<2>(int64_t* XX, const int64_t* X, size_t n){
	int64_t *p=XX;
	for(size_t i=0;i<n;i++){
		p[0]=X[i];
		p[1]=X[i];
		p+=2;
	}
}

template<>
inline void duplicate_vect<4>(int64_t* XX, const int64_t* X, size_t n){
	int64_t *p=XX;
	for(size_t i=0;i<n;i++){
		p[0]=X[i];
		p[1]=X[i];
		p[2]=X[i];
		p[3]=X[i];
		p+=4;
	}
}


/*************************
 *  MAIN GEBP OPERATION  *
 ************************/

void igebp( size_t rows, size_t cols, size_t depth, int64_t* C, size_t ldc, const int64_t* blockA, size_t lda,
	    const int64_t* blockB, size_t ldb, int64_t* blockW){

	size_t i,j,k;
	size_t prows,pcols,pdepth;
	prows=(rows/_mr)*_mr;
	pcols=(cols/_nr)*_nr;
	pdepth=(depth/4)*4;
	// process columns by pack of _nr
	for(j=0;j<pcols;j+=_nr){
		duplicate_vect<_vect_s>(blockW, blockB+j*ldb,depth*_nr);
		prefetch(blockW);
		// process rows by pack of _mr
		for (i=0;i<prows;i+=_mr){
			const int64_t* blA = blockA+i*lda;
			prefetch(blA);
			igebb44(i, j, depth, pdepth, C, ldc, blA, blockW);
		}
		i=prows;
		// process the (rows%_mr) remainings rows
		int rem=rows-prows;
		while (rem >0) {
			if (rem>=_vect_s){
				igebb24(i  ,j,depth, pdepth, C, ldc, blockA+i*lda, blockW);
				i+=_vect_s;
				rem-=_vect_s;
			}
			else{	// use blockB since no vectorization
				igebb14(i,j,depth, pdepth, C, ldc, blockA+i*lda, blockB+j*ldb);
				i++;
				rem--;
			}
		}
	}
	// process the (columns%_nr) remaining columns one by one
	for (;j<cols;j++){
		duplicate_vect<_vect_s>(blockW, blockB+j*ldb,depth);
		prefetch(blockW);
		// process rows by pack of _mr
		for (i=0;i<prows;i+=_mr){
			const int64_t* blA = blockA+i*lda;
			prefetch(blA);
			igebb41(i, j, depth, pdepth, C, ldc, blA, blockW);
		}
		i=prows;
		// process the (rows%_mr) remainings rows
		int rem=rows-prows;
		while (rem >0) {
			if (rem>=_vect_s){
				igebb21(i  ,j,depth, pdepth, C, ldc, blockA+i*lda, blockW);
				i+=_vect_s;
				rem-=_vect_s;
			}
			else{   // use blockB since no vectorization
				igebb11(i,j,depth, pdepth, C, ldc, blockA+i*lda, blockB+j*ldb);
				i++;
				rem--;
			}
		}
	}
}


/********************************************************
 * KERNEL FOR MATMUL USING SIMD OPERATION AND REGISTERS *
 ********************************************************/

inline void igebb44(size_t i, size_t j, size_t depth, size_t pdepth, int64_t* C, size_t ldc, const int64_t *blA, const int64_t* blB){
	size_t k;
	_vect_t C0,C1,C2,C3,C4,C5,C6,C7;
	VEC_SET_64(C0,0L);
	VEC_SET_64(C1,0L);
	VEC_SET_64(C2,0L);
	VEC_SET_64(C3,0L);
	VEC_SET_64(C4,0L);
	VEC_SET_64(C5,0L);
	VEC_SET_64(C6,0L);
	VEC_SET_64(C7,0L);
	int64_t *r0 = C+j*ldc+i;
	int64_t *r1 = r0+ldc;
	int64_t *r2 = r1+ldc;
	int64_t *r3 = r2+ldc;
	prefetch(r0+_vect_s);
	prefetch(r1+_vect_s);
	prefetch(r2+_vect_s);
	prefetch(r3+_vect_s);
	// process the loop by (_mrx4) by (4x4) matrix mul
	for (k=0;k<pdepth;k+=4){
		_vect_t A0,A1;
		_vect_t B0,B1,B2,B3;
		_vect_t Tmp;
		VEC_LOAD(A0, blA+0*StepA);
		VEC_LOAD(A1, blA+1*StepA);
		VEC_LOAD(B0, blB+0*StepB);
		VEC_LOAD(B1, blB+1*StepB);
		VEC_MADD_32(C0,A0,B0,Tmp);
		VEC_LOAD(B2, blB+2*StepB);
		VEC_MADD_32(C4,A1,B0,B0);
		VEC_LOAD(B3, blB+3*StepB);
		VEC_LOAD(B0, blB+4*StepB);
		VEC_MADD_32(C1,A0,B1,Tmp);
		VEC_MADD_32(C5,A1,B1,B1);
		VEC_LOAD(B1, blB+5*StepB);
		VEC_MADD_32(C2,A0,B2,Tmp);
		VEC_MADD_32(C6,A1,B2,B2);
		VEC_LOAD(B2, blB+6*StepB);
		VEC_MADD_32(C3,A0,B3,Tmp);
		VEC_LOAD(A0, blA+2*StepA);
		VEC_MADD_32(C7,A1,B3,B3);
		VEC_LOAD(A1, blA+3*StepA);
		VEC_LOAD(B3, blB+7*StepB);
		VEC_MADD_32(C0,A0,B0,Tmp);
		VEC_MADD_32(C4,A1,B0,B0);
		VEC_LOAD(B0, blB+8*StepB);
		VEC_MADD_32(C1,A0,B1,Tmp);
		VEC_MADD_32(C5,A1,B1,B1);
		VEC_LOAD(B1, blB+9*StepB);
		VEC_MADD_32(C2,A0,B2,Tmp);
		VEC_MADD_32(C6,A1,B2,B2);
		VEC_LOAD(B2, blB+10*StepB);
		VEC_MADD_32(C3,A0,B3,Tmp);
		VEC_LOAD(A0, blA+4*StepA);
		VEC_MADD_32(C7,A1,B3,B3);
		VEC_LOAD(A1, blA+5*StepA);
		VEC_LOAD(B3, blB+11*StepB);
		VEC_MADD_32(C0,A0,B0,Tmp);
		VEC_MADD_32(C4,A1,B0,B0);
		VEC_LOAD(B0, blB+12*StepB);
		VEC_MADD_32(C1,A0,B1,Tmp);
		VEC_MADD_32(C5,A1,B1,B1);
		VEC_LOAD(B1, blB+13*StepB);
		VEC_MADD_32(C2,A0,B2,Tmp);
		VEC_MADD_32(C6,A1,B2,B2);
		VEC_LOAD(B2, blB+14*StepB);
		VEC_MADD_32(C3,A0,B3,Tmp);
		VEC_LOAD(A0, blA+6*StepA);
		VEC_MADD_32(C7,A1,B3,B3);
		VEC_LOAD(A1, blA+7*StepA);
		VEC_LOAD(B3, blB+15*StepB);
		VEC_MADD_32(C0,A0,B0,Tmp);
		VEC_MADD_32(C4,A1,B0,B0);
		VEC_MADD_32(C1,A0,B1,Tmp);
		VEC_MADD_32(C5,A1,B1,B1);
		VEC_MADD_32(C2,A0,B2,Tmp);
		VEC_MADD_32(C6,A1,B2,B2);
		VEC_MADD_32(C3,A0,B3,Tmp);
		VEC_MADD_32(C7,A1,B3,B3);
		blA+= 8*StepA;
		blB+=16*StepB;
	}
	// process (depth mod 4) remaining entries by  (_mrx1) by (1x4) matrix mul
	for(;k<depth;k++){
		_vect_t A0,A1;
		_vect_t B0,B1,B2,B3;
		_vect_t Tmp;
		VEC_LOAD(A0, blA+0*StepA);
		VEC_LOAD(A1, blA+1*StepA);
		VEC_LOAD(B0, blB+0*StepB);
		VEC_LOAD(B1, blB+1*StepB);
		VEC_MADD_32(C0,A0,B0,Tmp);
		VEC_LOAD(B2, blB+2*StepB);
		VEC_MADD_32(C4,A1,B0,B0);
		VEC_LOAD(B3, blB+3*StepB);
		VEC_MADD_32(C1,A0,B1,Tmp);
		VEC_MADD_32(C5,A1,B1,B1);
		VEC_MADD_32(C2,A0,B2,Tmp);
		VEC_MADD_32(C6,A1,B2,B2);
		VEC_MADD_32(C3,A0,B3,Tmp);
		VEC_MADD_32(C7,A1,B3,B3);
		blA+=2*StepA;
		blB+=4*StepB;
	}
	_vect_t R0, R1, R2, R3, R4, R5, R6;
	VEC_LOADU(R0, r0);
	VEC_LOADU(R1, r1);
	VEC_LOADU(R2, r2);
	VEC_LOADU(R3, r3);
	VEC_LOADU(R4, r0+_vect_s);
	VEC_LOADU(R5, r1+_vect_s);
	VEC_LOADU(R6, r2+_vect_s);
	VEC_ADDIN_64(R0,C0);
	VEC_STOREU(r0,R0);
	VEC_LOADU(R0, r3+_vect_s);
	VEC_ADDIN_64(R1,C1);
	VEC_ADDIN_64(R2,C2);
	VEC_ADDIN_64(R3,C3);
	VEC_ADDIN_64(R4,C4);
	VEC_ADDIN_64(R5,C5);
	VEC_ADDIN_64(R6,C6);
	VEC_ADDIN_64(R0,C7);
	VEC_STOREU(r1,R1);
	VEC_STOREU(r2,R2);
	VEC_STOREU(r3,R3);
	VEC_STOREU(r0+_vect_s,R4);
	VEC_STOREU(r1+_vect_s,R5);
	VEC_STOREU(r2+_vect_s,R6);
	VEC_STOREU(r3+_vect_s,R0);

}


inline void igebb24(size_t i, size_t j, size_t depth, size_t pdepth, int64_t* C, size_t ldc, const int64_t *blA, const int64_t* blB){
//cout<<"aligned 32:"<< int64_t( blA)% 32 <<endl;
	size_t k;
	_vect_t C0,C1,C2,C3;
	VEC_SET_64(C0,0L);
	VEC_SET_64(C1,0L);
	VEC_SET_64(C2,0L);
	VEC_SET_64(C3,0L);
	int64_t *r0 = C+j*ldc+i;
	int64_t *r1 = r0+ldc;
	int64_t *r2 = r1+ldc;
	int64_t *r3 = r2+ldc;
	// process the loop by (1/2_mrx4) by (4x4) matrix mul
	for (k=0;k<pdepth;k+=4){
		_vect_t A0;
		_vect_t B0,B1,B2,B3;
		_vect_t Tmp;
		VEC_LOAD(A0, blA+0*StepA);
		VEC_LOAD(B0, blB+0*StepB);
		VEC_LOAD(B1, blB+1*StepB);
		VEC_MADD_32(C0,A0,B0,Tmp);
		VEC_LOAD(B2, blB+2*StepB);
		VEC_LOAD(B3, blB+3*StepB);
		VEC_LOAD(B0, blB+4*StepB);
		VEC_MADD_32(C1,A0,B1,Tmp);
		VEC_LOAD(B1, blB+5*StepB);
		VEC_MADD_32(C2,A0,B2,Tmp);
		VEC_LOAD(B2, blB+6*StepB);
		VEC_MADD_32(C3,A0,B3,Tmp);
		VEC_LOAD(A0, blA+1*StepA);
		VEC_LOAD(B3, blB+7*StepB);
		VEC_MADD_32(C0,A0,B0,Tmp);
		VEC_LOAD(B0, blB+8*StepB);
		VEC_MADD_32(C1,A0,B1,Tmp);
		VEC_LOAD(B1, blB+9*StepB);
		VEC_MADD_32(C2,A0,B2,Tmp);
		VEC_LOAD(B2, blB+10*StepB);
		VEC_MADD_32(C3,A0,B3,Tmp);
		VEC_LOAD(A0, blA+2*StepA);
		VEC_LOAD(B3, blB+11*StepB);
		VEC_MADD_32(C0,A0,B0,Tmp);
		VEC_LOAD(B0, blB+12*StepB);
		VEC_MADD_32(C1,A0,B1,Tmp);
		VEC_LOAD(B1, blB+13*StepB);
		VEC_MADD_32(C2,A0,B2,Tmp);
		VEC_LOAD(B2, blB+14*StepB);
		VEC_MADD_32(C3,A0,B3,Tmp);
		VEC_LOAD(A0, blA+3*StepA);
		VEC_LOAD(B3, blB+15*StepB);
		VEC_MADD_32(C0,A0,B0,Tmp);
		VEC_MADD_32(C1,A0,B1,Tmp);
		VEC_MADD_32(C2,A0,B2,Tmp);
		VEC_MADD_32(C3,A0,B3,Tmp);
		blA+= 4*StepA;
		blB+=16*StepB;
	}
	// process (depth mod 4) remaining entries by  (1/2_mrx1) by (1x4) matrix mul
	for(;k<depth;k++){
		_vect_t A0;
		_vect_t B0,B1,B2,B3;
		_vect_t Tmp;
		VEC_LOAD(A0, blA+0*StepA);
		VEC_LOAD(B0, blB+0*StepB);
		VEC_LOAD(B1, blB+1*StepB);
		VEC_MADD_32(C0,A0,B0,Tmp);
		VEC_LOAD(B2, blB+2*StepB);
		VEC_LOAD(B3, blB+3*StepB);
		VEC_MADD_32(C1,A0,B1,Tmp);
		VEC_MADD_32(C2,A0,B2,Tmp);
		VEC_MADD_32(C3,A0,B3,Tmp);
		blA+=StepA;
		blB+=4*StepB;
	}
	_vect_t R0, R1, R2, R3;
	VEC_LOADU(R0, r0);
	VEC_LOADU(R1, r1);
	VEC_LOADU(R2, r2);
	VEC_LOADU(R3, r3);
	VEC_ADDIN_64(R0,C0);
	VEC_ADDIN_64(R1,C1);
	VEC_ADDIN_64(R2,C2);
	VEC_ADDIN_64(R3,C3);
	VEC_STOREU(r0,R0);
	VEC_STOREU(r1,R1);
	VEC_STOREU(r2,R2);
	VEC_STOREU(r3,R3);

}


inline void igebb14(size_t i, size_t j, size_t depth, size_t pdepth, int64_t* C, size_t ldc, const int64_t *blA, const int64_t* blB){
	size_t k;
	int64_t *r0 = C+j*ldc+i;
	int64_t *r1 = r0+ldc;
	int64_t *r2 = r1+ldc;
	int64_t *r3 = r2+ldc;
	for(k=0;k<depth;k++){
		r0[0]+=blA[0]*blB[0];
		r1[0]+=blA[0]*blB[1];
		r2[0]+=blA[0]*blB[2];
		r3[0]+=blA[0]*blB[3];
		blA++;
		blB+=4;
	}
}


inline void igebb41(size_t i, size_t j, size_t depth, size_t pdepth, int64_t* C, size_t ldc, const int64_t *blA, const int64_t* blB){
	size_t k;
	_vect_t C0,C4;
	VEC_SET_64(C0,0L);
	VEC_SET_64(C4,0L);
	int64_t *r0 = C+j*ldc+i;
	int64_t *r4 = r0+_vect_s;

	// process the loop by (_mrx1) by (1x1) matrix mul
	for (k=0;k<depth;k++){
		_vect_t A0,A1;
		_vect_t B0;
		_vect_t Tmp;
		VEC_LOAD(A0, blA+0*StepA);
		VEC_LOAD(A1, blA+1*StepA);
		VEC_LOAD(B0, blB+0*StepB);
		VEC_MADD_32(C0,A0,B0,Tmp);
		VEC_MADD_32(C4,A1,B0,B0);
		blA+= 2*StepA;
		blB+= 1*StepB;
	}
	_vect_t R0, R4;
	VEC_LOADU(R0, r0);
	VEC_LOADU(R4, r4);
	VEC_ADDIN_64(R0,C0);
	VEC_ADDIN_64(R4,C4);
	VEC_STOREU(r0,R0);
	VEC_STOREU(r4,R4);
}


inline void igebb21(size_t i, size_t j, size_t depth, size_t pdepth, int64_t* C, size_t ldc, const int64_t *blA, const int64_t* blB){
	size_t k;
	_vect_t C0;
	VEC_SET_64(C0,0L);
	int64_t *r0 = C+j*ldc+i;

	// process the loop by (1/2_mrx1) by (1x1) matrix mul
	for (k=0;k<depth;k++){
		_vect_t A0;
		_vect_t B0;
		_vect_t Tmp;
		VEC_LOAD(A0, blA+0*StepA);
		VEC_LOAD(B0, blB+0*StepB);
		VEC_MADD_32(C0,A0,B0,Tmp);
		blA+= 1*StepA;
		blB+= 1*StepB;
	}
	_vect_t R0;
	VEC_LOADU(R0, r0);
	VEC_ADDIN_64(R0,C0);
	VEC_STOREU(r0,R0);
}


inline void igebb11(size_t i, size_t j, size_t depth, size_t pdepth, int64_t* C, size_t ldc, const int64_t *blA, const int64_t* blB){
	size_t k;
	int64_t *r0 = C+j*ldc+i;
	for(k=0;k<depth;k++){
		r0[0]+=blA[k]*blB[k];
	}
}


#include <iostream>
ostream& operator<<(ostream &os, _vect_t& x){
	int64_t T[_vect_s];
	VEC_STORE(T,x);
	os <<"[";
	for(size_t i=0;i<_vect_s-1;i++)
		os<<T[i]<<" ";
	os<<T[_vect_s-1]<<"]";
	return os;
}

