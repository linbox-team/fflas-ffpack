/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2013,2014  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 * the code is inspired and adapted from the Eigen library
 * modified by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

#ifndef __FFLASFFPACK_fflas_igemm_igemm_kernels_INL
#define __FFLASFFPACK_fflas_igemm_igemm_kernels_INL


#ifdef __AVX2__
#define _nr 4
#define _mr 8
#define StepA 4
#define StepB 4
#elif defined(__SSE4_1__) or defined(__AVX__)
#define _nr 4
#define _mr 4
#define StepA 2
#define StepB 2
#else
#error "kernels not supported"
#endif

#include "fflas-ffpack/utils/fflas_memory.h"
#include "igemm_tools.h"

/********************************************************
 * KERNEL FOR MATMUL USING SIMD OPERATION AND REGISTERS *
 ********************************************************/

namespace FFLAS { namespace details { /*  kernels */

	template<enum number_kind K>
	inline void igebb44(size_t i, size_t j, size_t depth, size_t pdepth
			    , const int64_t alpha
			    , const int64_t *blA, const int64_t* blB
			    , int64_t* C, size_t ldc
			   )
	{
		using simd = Simd<int64_t>;
		using vect_t =  typename simd::vect_t;
		size_t k;
		vect_t C0,C1,C2,C3,C4,C5,C6,C7;
		C0 = simd::zero();
		C1 = simd::zero();
		C2 = simd::zero();
		C3 = simd::zero();
		C4 = simd::zero();
		C5 = simd::zero();
		C6 = simd::zero();
		C7 = simd::zero();
		int64_t *r0 = C+j*ldc+i;
		int64_t *r1 = r0+ldc;
		int64_t *r2 = r1+ldc;
		int64_t *r3 = r2+ldc;
		prefetch(r0+simd::vect_size);
		prefetch(r1+simd::vect_size);
		prefetch(r2+simd::vect_size);
		prefetch(r3+simd::vect_size);
		// process the loop by (_mrx4) by (4x4) matrix mul
		for (k=0;k<pdepth;k+=4){
			vect_t A0,A1;
			vect_t B0,B1,B2,B3;
			A0 = simd::load( blA+0*StepA);
			A1 = simd::load( blA+1*StepA);
			B0 = simd::load( blB+0*StepB);
			B1 = simd::load( blB+1*StepB);
			simd::fmaddxin(C0,A0,B0);
			B2 = simd::load( blB+2*StepB);
			simd::fmaddxin(C4,A1,B0); // B0
			B3 = simd::load( blB+3*StepB);
			B0 = simd::load( blB+4*StepB);
			simd::fmaddxin(C1,A0,B1);
			simd::fmaddxin(C5,A1,B1); // B1
			B1 = simd::load( blB+5*StepB);
			simd::fmaddxin(C2,A0,B2);
			simd::fmaddxin(C6,A1,B2); // B2
			B2 = simd::load( blB+6*StepB);
			simd::fmaddxin(C3,A0,B3);
			A0 = simd::load( blA+2*StepA);
			simd::fmaddxin(C7,A1,B3); // B3
			A1 = simd::load( blA+3*StepA);
			B3 = simd::load( blB+7*StepB);
			simd::fmaddxin(C0,A0,B0);
			simd::fmaddxin(C4,A1,B0); // B0
			B0 = simd::load( blB+8*StepB);
			simd::fmaddxin(C1,A0,B1);
			simd::fmaddxin(C5,A1,B1); // B1
			B1 = simd::load( blB+9*StepB);
			simd::fmaddxin(C2,A0,B2);
			simd::fmaddxin(C6,A1,B2); // B2
			B2 = simd::load( blB+10*StepB);
			simd::fmaddxin(C3,A0,B3);
			A0 = simd::load( blA+4*StepA);
			simd::fmaddxin(C7,A1,B3); // B3
			A1 = simd::load( blA+5*StepA);
			B3 = simd::load( blB+11*StepB);
			simd::fmaddxin(C0,A0,B0);
			simd::fmaddxin(C4,A1,B0); // B0
			B0 = simd::load( blB+12*StepB);
			simd::fmaddxin(C1,A0,B1);
			simd::fmaddxin(C5,A1,B1); // B1
			B1 = simd::load( blB+13*StepB);
			simd::fmaddxin(C2,A0,B2);
			simd::fmaddxin(C6,A1,B2); // B2
			B2 = simd::load( blB+14*StepB);
			simd::fmaddxin(C3,A0,B3);
			A0 = simd::load( blA+6*StepA);
			simd::fmaddxin(C7,A1,B3); // B3
			A1 = simd::load( blA+7*StepA);
			B3 = simd::load( blB+15*StepB);
			simd::fmaddxin(C0,A0,B0);
			simd::fmaddxin(C4,A1,B0); // B0
			simd::fmaddxin(C1,A0,B1);
			simd::fmaddxin(C5,A1,B1); // B1
			simd::fmaddxin(C2,A0,B2);
			simd::fmaddxin(C6,A1,B2); // B2
			simd::fmaddxin(C3,A0,B3);
			simd::fmaddxin(C7,A1,B3); // B3
			blA+= 8*StepA;
			blB+=16*StepB;
		}
		// process (depth mod 4) remaining entries by  (_mrx1) by (1x4) matrix mul
		for(;k<depth;k++){
			vect_t A0,A1;
			vect_t B0,B1,B2,B3;
			A0 = simd::load( blA+0*StepA);
			A1 = simd::load( blA+1*StepA);
			B0 = simd::load( blB+0*StepB);
			B1 = simd::load( blB+1*StepB);
			simd::fmaddxin(C0,A0,B0);
			B2 = simd::load( blB+2*StepB);
			simd::fmaddxin(C4,A1,B0); // B0
			B3 = simd::load( blB+3*StepB);
			simd::fmaddxin(C1,A0,B1);
			simd::fmaddxin(C5,A1,B1);  // B1
			simd::fmaddxin(C2,A0,B2);
			simd::fmaddxin(C6,A1,B2); // B2
			simd::fmaddxin(C3,A0,B3);
			simd::fmaddxin(C7,A1,B3); // B3
			blA+=2*StepA;
			blB+=4*StepB;
		}
		vect_t R0, R1, R2, R3, R4, R5, R6;
		vect_t A0 ;
		A0 = simd::set1(alpha);
		R0 = simd::loadu( r0);
		R1 = simd::loadu( r1);
		R2 = simd::loadu( r2);
		R3 = simd::loadu( r3);
		R4 = simd::loadu( r0+simd::vect_size);
		R5 = simd::loadu( r1+simd::vect_size);
		R6 = simd::loadu( r2+simd::vect_size);
		if (K == number_kind::one) {
			simd::addin(R0,C0);
		}
		if (K == number_kind::mone) {
			simd::subin(R0,C0);
		}
		if (K == number_kind::other) {
			simd::fmaddxin(R0,A0,C0);
		}
		simd::storeu(r0,R0);
		R0 = simd::loadu( r3+simd::vect_size);
		if (K == number_kind::one) {
			simd::addin(R1,C1);
			simd::addin(R2,C2);
			simd::addin(R3,C3);
			simd::addin(R4,C4);
			simd::addin(R5,C5);
			simd::addin(R6,C6);
			simd::addin(R0,C7);
		}
		if (K == number_kind::mone) {
			simd::subin(R1,C1);
			simd::subin(R2,C2);
			simd::subin(R3,C3);
			simd::subin(R4,C4);
			simd::subin(R5,C5);
			simd::subin(R6,C6);
			simd::subin(R0,C7);
		}
		if (K == number_kind::other) {
			simd::fmaddxin(R1,A0,C1);
			simd::fmaddxin(R2,A0,C2);
			simd::fmaddxin(R3,A0,C3);
			simd::fmaddxin(R4,A0,C4);
			simd::fmaddxin(R5,A0,C5);
			simd::fmaddxin(R6,A0,C6);
			simd::fmaddxin(R0,A0,C7);
		}
		simd::storeu(r1,R1);
		simd::storeu(r2,R2);
		simd::storeu(r3,R3);
		simd::storeu(r0+simd::vect_size,R4);
		simd::storeu(r1+simd::vect_size,R5);
		simd::storeu(r2+simd::vect_size,R6);
		simd::storeu(r3+simd::vect_size,R0);

	}


	template<enum number_kind K>
	inline void igebb24(size_t i, size_t j, size_t depth, size_t pdepth
			    , const int64_t alpha
			    , const int64_t *blA, const int64_t* blB
			    , int64_t* C, size_t ldc
			   )
	{
		using simd = Simd<int64_t>;
		using vect_t =  typename simd::vect_t;

		//cout<<"aligned 32:"<< int64_t( blA)% 32 <<endl;
		size_t k;
		vect_t C0,C1,C2,C3;
		C0 = simd::zero();
		C1 = simd::zero();
		C2 = simd::zero();
		C3 = simd::zero();
		int64_t *r0 = C+j*ldc+i;
		int64_t *r1 = r0+ldc;
		int64_t *r2 = r1+ldc;
		int64_t *r3 = r2+ldc;
		// process the loop by (1/2_mrx4) by (4x4) matrix mul
		for (k=0;k<pdepth;k+=4){
			vect_t A0;
			vect_t B0,B1,B2,B3;
			A0 = simd::load( blA+0*StepA);
			B0 = simd::load( blB+0*StepB);
			B1 = simd::load( blB+1*StepB);
			simd::fmaddxin(C0,A0,B0);
			B2 = simd::load( blB+2*StepB);
			B3 = simd::load( blB+3*StepB);
			B0 = simd::load( blB+4*StepB);
			simd::fmaddxin(C1,A0,B1);
			B1 = simd::load( blB+5*StepB);
			simd::fmaddxin(C2,A0,B2);
			B2 = simd::load( blB+6*StepB);
			simd::fmaddxin(C3,A0,B3);
			A0 = simd::load( blA+1*StepA);
			B3 = simd::load( blB+7*StepB);
			simd::fmaddxin(C0,A0,B0);
			B0 = simd::load( blB+8*StepB);
			simd::fmaddxin(C1,A0,B1);
			B1 = simd::load( blB+9*StepB);
			simd::fmaddxin(C2,A0,B2);
			B2 = simd::load( blB+10*StepB);
			simd::fmaddxin(C3,A0,B3);
			A0 = simd::load( blA+2*StepA);
			B3 = simd::load( blB+11*StepB);
			simd::fmaddxin(C0,A0,B0);
			B0 = simd::load( blB+12*StepB);
			simd::fmaddxin(C1,A0,B1);
			B1 = simd::load( blB+13*StepB);
			simd::fmaddxin(C2,A0,B2);
			B2 = simd::load( blB+14*StepB);
			simd::fmaddxin(C3,A0,B3);
			A0 = simd::load( blA+3*StepA);
			B3 = simd::load( blB+15*StepB);
			simd::fmaddxin(C0,A0,B0);
			simd::fmaddxin(C1,A0,B1);
			simd::fmaddxin(C2,A0,B2);
			simd::fmaddxin(C3,A0,B3);
			blA+= 4*StepA;
			blB+=16*StepB;
		}
		// process (depth mod 4) remaining entries by  (1/2_mrx1) by (1x4) matrix mul
		for(;k<depth;k++){
			vect_t A0;
			vect_t B0,B1,B2,B3;
			A0 = simd::load( blA+0*StepA);
			B0 = simd::load( blB+0*StepB);
			B1 = simd::load( blB+1*StepB);
			simd::fmaddxin(C0,A0,B0);
			B2 = simd::load( blB+2*StepB);
			B3 = simd::load( blB+3*StepB);
			simd::fmaddxin(C1,A0,B1);
			simd::fmaddxin(C2,A0,B2);
			simd::fmaddxin(C3,A0,B3);
			blA+=StepA;
			blB+=4*StepB;
		}
		vect_t R0, R1, R2, R3;
		vect_t A0 ;
		A0 = simd::set1(alpha);
		R0 = simd::loadu( r0);
		R1 = simd::loadu( r1);
		R2 = simd::loadu( r2);
		R3 = simd::loadu( r3);
		if ( K == number_kind::one) {
			simd::addin(R0,C0);
			simd::addin(R1,C1);
			simd::addin(R2,C2);
			simd::addin(R3,C3);
		}
		if ( K == number_kind::mone) {
			simd::subin(R0,C0);
			simd::subin(R1,C1);
			simd::subin(R2,C2);
			simd::subin(R3,C3);
		}
		if ( K == number_kind::other) {
			simd::fmaddxin(R0,A0,C0);
			simd::fmaddxin(R1,A0,C1);
			simd::fmaddxin(R2,A0,C2);
			simd::fmaddxin(R3,A0,C3);
		}
		simd::storeu(r0,R0);
		simd::storeu(r1,R1);
		simd::storeu(r2,R2);
		simd::storeu(r3,R3);

	}


	template<enum number_kind K>
	inline void igebb14(size_t i, size_t j, size_t depth, size_t pdepth
			    , const int64_t alpha
			    , const int64_t *blA, const int64_t* blB
			    , int64_t* C, size_t ldc
			   )
	{
		// using simd = Simd<int64_t>;
		// using vect_t =  typename simd::vect_t;

		size_t k;
		int64_t *r0 = C+j*ldc+i;
		int64_t *r1 = r0+ldc;
		int64_t *r2 = r1+ldc;
		int64_t *r3 = r2+ldc;
		for(k=0;k<depth;k++){
			if (K == number_kind::one) {
				r0[0]+=blA[0]*blB[0];
				r1[0]+=blA[0]*blB[1];
				r2[0]+=blA[0]*blB[2];
				r3[0]+=blA[0]*blB[3];
			}
			if (K == number_kind::mone) {
				r0[0]-=blA[0]*blB[0];
				r1[0]-=blA[0]*blB[1];
				r2[0]-=blA[0]*blB[2];
				r3[0]-=blA[0]*blB[3];
			}
			if (K == number_kind::other) {
				int64_t abla = alpha*blA[0];
				r0[0]+=abla*blB[0];
				r1[0]+=abla*blB[1];
				r2[0]+=abla*blB[2];
				r3[0]+=abla*blB[3];
			}

			blA++;
			blB+=4;
		}
	}


	template<enum number_kind K>
	inline void igebb41(size_t i, size_t j, size_t depth, size_t pdepth
			    , const int64_t alpha
			    , const int64_t *blA, const int64_t* blB
			    , int64_t* C, size_t ldc
			   )
	{
		using simd = Simd<int64_t>;
		using vect_t =  typename simd::vect_t;

		size_t k;
		vect_t C0,C4;
		C0 = simd::zero();
		C4 = simd::zero();
		int64_t *r0 = C+j*ldc+i;
		int64_t *r4 = r0+simd::vect_size;

		// process the loop by (_mrx1) by (1x1) matrix mul
		for (k=0;k<depth;k++){
			vect_t A0,A1;
			vect_t B0;
			A0 = simd::load( blA+0*StepA);
			A1 = simd::load( blA+1*StepA);
			B0 = simd::load( blB+0*StepB);
			simd::fmaddxin(C0,A0,B0);
			simd::fmaddxin(C4,A1,B0); //! bug ,B0 dans VEC_MADD_32 ?
			blA+= 2*StepA;
			blB+= 1*StepB;
		}
		vect_t R0, R4;
		R0 = simd::loadu( r0);
		R4 = simd::loadu( r4);
		vect_t A0 ;
		A0 = simd::set1(alpha);
		if (K == number_kind::one) {
			simd::addin(R0,C0);
			simd::addin(R4,C4);
		}
		if (K == number_kind::mone) {
			simd::subin(R0,C0);
			simd::subin(R4,C4);
		}
		if (K == number_kind::other) {
			simd::fmaddxin(R0,A0,C0);
			simd::fmaddxin(R4,A0,C4);
		}
		simd::storeu(r0,R0);
		simd::storeu(r4,R4);
	}


	template<enum number_kind K>
	inline void igebb21(size_t i, size_t j, size_t depth, size_t pdepth
			    , const int64_t alpha
			    , const int64_t *blA, const int64_t* blB
			    , int64_t* C, size_t ldc
			   )
	{
		using simd = Simd<int64_t>;
		using vect_t =  typename simd::vect_t;

		size_t k;
		vect_t C0;
		C0 = simd::zero();
		int64_t *r0 = C+j*ldc+i;

		// process the loop by (1/2_mrx1) by (1x1) matrix mul
		for (k=0;k<depth;k++){
			vect_t A0;
			vect_t B0;
			A0 = simd::load( blA+0*StepA);
			B0 = simd::load( blB+0*StepB);
			simd::fmaddxin(C0,A0,B0);
			blA+= 1*StepA;
			blB+= 1*StepB;
		}
		vect_t R0;
		vect_t A0 ;
		A0 = simd::set1(alpha);

		R0 = simd::loadu( r0);
		if ( K == number_kind::one)
			simd::addin(R0,C0);
		if ( K == number_kind::mone)
			simd::subin(R0,C0);
		if ( K == number_kind::other)
			simd::fmaddxin(R0,A0,C0);
		simd::storeu(r0,R0);
	}


	template<enum number_kind K>
	inline void igebb11(size_t i, size_t j, size_t depth, size_t pdepth
			    , const int64_t alpha
			    , const int64_t *blA, const int64_t* blB
			    , int64_t* C, size_t ldc
			   )
	{
		// using simd = Simd<int64_t>;
		// using vect_t =  typename simd::vect_t;
		size_t k;
		int64_t *r0 = C+j*ldc+i;
		for(k=0;k<depth;k++){
			if (K == number_kind::one)
				r0[0]+=blA[k]*blB[k];
			if ( K == number_kind::mone)
				r0[0]-=blA[k]*blB[k];
			if ( K == number_kind::other)
				r0[0]+=alpha*blA[k]*blB[k];
		}
	}


} // details
} // FFLAS


/*************************
 *  MAIN GEBP OPERATION  *
 ************************/

namespace FFLAS { namespace details { /*  main */

	template<enum number_kind K>
	void igebp( size_t rows, size_t cols, size_t depth
		    , const int64_t alpha
		    , const int64_t* blockA, size_t lda,
		    const int64_t* blockB, size_t ldb,
		    int64_t* C, size_t ldc,
		    int64_t* blockW)
	{
		using simd = Simd<int64_t>;
		// using vect_t =  typename simd::vect_t;
		size_t i,j;
		size_t prows,pcols,pdepth;
		prows=(rows/_mr)*_mr;
		pcols=(cols/_nr)*_nr;
		pdepth=(depth/4)*4;
		// process columns by pack of _nr
		for(j=0;j<pcols;j+=_nr){
			duplicate_vect<simd::vect_size>(blockW, blockB+j*ldb,depth*_nr);
			prefetch(blockW);
			// process rows by pack of _mr
			for (i=0;i<prows;i+=_mr){
				const int64_t* blA = blockA+i*lda;
				prefetch(blA);
				igebb44<K>(i, j, depth, pdepth, alpha, blA, blockW, C, ldc);
			}
			i=prows;
			// process the (rows%_mr) remainings rows
			int rem=(int)(rows-prows);
			while (rem >0) {
				if (rem>=(int)simd::vect_size){
					igebb24<K>(i  ,j,depth, pdepth, alpha , blockA+i*lda, blockW, C, ldc);
					i+=simd::vect_size;
					rem-=(int)simd::vect_size;
				}
				else{	// use blockB since no vectorization
					igebb14<K>(i,j,depth, pdepth, alpha, blockA+i*lda, blockB+j*ldb, C, ldc);
					i++;
					rem--;
				}
			}
		}
		// process the (columns%_nr) remaining columns one by one
		for (;j<cols;j++){
                        duplicate_vect<simd::vect_size>(blockW, blockB+j*ldb,depth);
			prefetch(blockW);
			// process rows by pack of _mr
			for (i=0;i<prows;i+=_mr){
				const int64_t* blA = blockA+i*lda;
				prefetch(blA);
				igebb41<K>(i, j, depth, pdepth, alpha, blA, blockW, C, ldc);
			}
			i=prows;
			// process the (rows%_mr) remainings rows
			int rem=(int)(rows-prows);
			while (rem >0) {
				if (rem>=(int)simd::vect_size){
					igebb21<K>(i  ,j,depth, pdepth, alpha, blockA+i*lda, blockW, C, ldc);
					i+=simd::vect_size;
					rem-=(int)(simd::vect_size);
				}
				else{   // use blockB since no vectorization
					igebb11<K>(i,j,depth, pdepth, alpha, blockA+i*lda, blockB+j*ldb, C, ldc);
					i++;
					rem--;
				}
			}
		}
	}


} // details
} // FFLAS

#endif // __FFLASFFPACK_fflas_igemm_igemm_kernels_INL
