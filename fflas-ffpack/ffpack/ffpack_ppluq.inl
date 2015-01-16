/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* ffpack/ffpack_ppluq.inl
 * Copyright (C) 2014 Ziad Sultan
 *
 * Written by Ziad.Sultan@imag.fr
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
 * License along with this library; if not, WRITE to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __FFLASFFPACK_ffpack_ppluq_INL
#define __FFLASFFPACK_ffpack_ppluq_INL



#ifdef __FFLASFFPACK_USE_OPENMP

#define __FFLAS__TRSM_READONLY

#define PBASECASE_K 256
//#define PBASECASE_K 5

namespace FFPACK {

	template<class Field>
	void  threads_fgemm(const size_t m, const size_t n, const size_t r, int nbthreads, size_t * W1, size_t * W2, size_t * W3, size_t gamma)
	{
		size_t H1, H2, H3;
		size_t M2 = m>>1;
		size_t N2 = n>>1;

		H1 = ((m-N2)*r*(N2-r))<<1;
		H2 = ((M2-r)*r*(n-N2))<<1;
		H3 = ((m-M2)*r*(n-N2))<<1;

		// if we take into account 2 concurrent pluq calls....
		size_t h;
		size_t z1= h*((m-M2)*(N2-r)*(N2-r)-(N2-r)*(N2-r)*(N2-r)/3);
		size_t z2= h*((n-N2)*(M2-r)*(M2-r)-(M2-r)*(M2-r)*(M2-r)/3);

		H1+= z1;
		H2+= z2;

		// compute number of threads for each fgemm call
		*W1=std::max(H1*nbthreads/(H1+H2+H3),(size_t)1);
		*W2=std::max(H2*nbthreads/(H1+H2+H3),(size_t)1);
		*W3=std::max(nbthreads-*W1-*W2,(size_t)1);

		// add gamma factor to change number of threads for pluq calls
		W1-= gamma*z1/(z1+z2);
		W2-= gamma*(1-z1/(z1+z2));
		W3+= gamma;

	}

	template<class Field>
	void threads_ftrsm(const size_t m, const size_t n, int nbthreads, size_t * t1, size_t * t2)
	{
		*t1 = nbthreads*m/(m+n);
		*t2 = nbthreads-(int)*t1;
	}


	template<class Field>
	inline size_t
	pPLUQ(const Field& Fi, const FFLAS::FFLAS_DIAG Diag,
	      const size_t M, const size_t N,
	      typename Field::Element_ptr A, const size_t lda,
	      size_t* P, size_t* Q, int nt)
	  {


    for (size_t i=0; i<M; ++i) P[i] = i;
    for (size_t i=0; i<N; ++i) Q[i] = i;
    if (std::min(M,N) == 0) return 0;
    if (std::max (M,N) == 1) return (Fi.isZero(*A))? 0 : 1;
    if (M == 1){
      size_t piv = 0;
      while ((piv < N) && Fi.isZero (A[piv])) piv++;
      if (piv == N)
	return 0;
      if (piv){
	Q[0] = piv;
	Fi.assign (*A, A[piv]);
	Fi.assign (A[piv], Fi.zero);
      }
      if (Diag== FFLAS::FflasUnit){
	typename Field::Element invpivot;
	Fi.inv(invpivot, *A);
	// for (size_t i=piv+1; i<N; ++i)
	  // Fi.mulin (A[i], invpivot);
	FFLAS::fscalin(Fi,N-piv-1,invpivot,A+piv+1,1);
      }
      return 1;
    }
    if (N == 1){
      size_t piv = 0;
      while ((piv < M) && Fi.isZero (A[piv*lda])) piv++;
      if (piv == M)
	return 0;
      if (piv){
	P[0] = piv;
	Fi.assign (*A, *(A+piv*lda));
	Fi.assign (*(A+piv*lda), Fi.zero);
      }
      if (Diag== FFLAS::FflasNonUnit){
	typename Field::Element invpivot;
	Fi.inv(invpivot, *A);
	// for (size_t i=piv+1; i<M; ++i)
	  // Fi.mulin (*(A+i*lda), invpivot);
	FFLAS::fscalin(Fi,M-piv-1,invpivot,A+(piv+1)*lda,lda);
      }
      return 1;
    }

    #ifdef PBASECASE_K
    if (std::min(M,N) < PBASECASE_K)
      return PLUQ_basecaseCrout (Fi, Diag, M, N, A, lda, P, Q);
    #endif
    FFLAS::FFLAS_DIAG OppDiag = (Diag == FFLAS::FflasUnit)? FFLAS::FflasNonUnit : FFLAS::FflasUnit;

    size_t M2 = M >> 1;
    size_t N2 = N >> 1;
    size_t * P1 = new size_t [M2];
    size_t * Q1 = new size_t [N2];
    size_t R1,R2,R3,R4;
	//if (M>8000) std::cerr<<"PLUQ 1"<<std::endl;
    // A1 = P1 [ L1 ] [ U1 V1 ] Q1
    //        [ M1 ]
    R1 = pPLUQ (Fi, Diag, M2, N2, A, lda, P1, Q1,nt);

    typename Field::Element * A2 = A + N2;
    typename Field::Element * A3 = A + M2*lda;
    typename Field::Element * A4 = A3 + N2;
    typename Field::Element * F = A2 + R1*lda;
    typename Field::Element * G = A3 + R1;

//    const FFLAS::CuttingStrategy meth = FFLAS::BLOCK_THREADS;
    const FFLAS::CuttingStrategy meth = FFLAS::TWO_D_ADAPT;
//    const FFLAS::CuttingStrategy meth = FFLAS::THREE_D_ADAPT;
    /*
    FFLAS::MMHelper<Field,FFLAS::MMHelperAlgo::Winograd,
		    typename FFLAS::FieldTraits<Field>::value,
		    FFLAS::ParSeqHelper::Parallel>
	    pWH (Fi, -1, FFLAS::ParSeqHelper::Parallel(MAX_THREADS, method));
    */
    typename FFLAS::ParSeqHelper::Parallel pWH (nt, meth);
    typename FFLAS::ParSeqHelper::Parallel PH (std::max(nt,1), meth);

#ifdef __FFLASFFPACK_USE_OPENMP4
    //#pragma omp parallel firstprivate(M2, N2) shared(P1, Q1, A2, A3, Fi)

      // [ B1 ] <- P1^T A2
      // [ B2 ]
#pragma omp task shared(P1, A2, Fi) depend(inout:A2) depend(in:P1)
      //#pragma omp task depend(in:A2, Fi)
      papplyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N-N2, 0, M2, A2, lda, P1);
      // [ C1 C2 ] <- A3 Q1^T
#pragma omp task  shared(Q1, A3, Fi) depend(inout:A3) depend(in:Q1)
      papplyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M-M2, 0, N2, A3, lda, Q1);

      //#pragma omp taskwait

      // D <- L1^-1 B1
//    if (M>8000) std::cerr<<"TRSMs D,E"<<std::endl;
#pragma omp task  shared( A2, A, Fi, N2, R1, PH) depend(in:A) depend(inout:A2)
      ftrsm (Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R1, N-N2, Fi.one, A, lda, A2, lda, PH);

      // E <- C1 U1^-1
#pragma omp task  shared(A, A3, Fi, M2, R1, PH) depend(in:A) depend(inout:A3)
      ftrsm (Fi, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag, M-M2, R1, Fi.one, A, lda, A3, lda, PH);

//#pragma omp taskwait


    typename Field::Element * AR1lda = A + R1*lda;
    typename Field::Element * AR1 = A + R1;
    //    typename Field::Element * A3R1 = A3 + R1;

//    if (M>8000) std::cerr<<"GEMMs F,G"<<std::endl;
    // F <- B2 - M1 D
#pragma omp task  shared(A, A2, Fi, F, pWH) depend(in: A2)  depend(inout:F) //depend(inout:A2[R1*lda:M2-R1])
    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M2-R1, N-N2, R1, Fi.mOne, AR1lda, lda, A2, lda, Fi.one, F, lda, pWH);
    // G <- C2 - E V1
#pragma omp task shared(A, A3, Fi, G, pWH) depend(in:A3) depend(inout:G)
    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2, N2-R1, R1, Fi.mOne, A3, lda, AR1, lda, Fi.one, G, lda, pWH);

//    if (M>8000) std::cerr<<"GEMM H"<<std::endl;
    // H <- A4 - ED
#pragma omp task shared(A3, A2, A4, Fi, pWH) depend(in:A3,A2) depend(inout:A4)
    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2, N-N2, R1, Fi.mOne, A3, lda, A2, lda, Fi.one, A4, lda, pWH);

//#pragma omp taskwait
//    if (M>8000) std::cerr<<"PLUQ 2"<<std::endl;
    // F = P2 [ L2 ] [ U2 V2 ] Q2
    //        [ M2 ]
    size_t * P2 = new size_t [M2-R1];
    size_t * Q2 = new size_t [N-N2];
    typename Field::Element * A4R2 = 0;
#pragma omp task shared(P2, Q2, F, Fi, R2, A4R2) depend(inout:F) depend(out: P2,Q2,R2, A4R2) // P2 Q2 should be out only
    {
	    R2 = pPLUQ (Fi, Diag, M2-R1, N-N2, F, lda, P2, Q2,nt/2);
	    A4R2 = A4 + R2;
    }
    // G = P3 [ L3 ] [ U3 V3 ] Q3
    //        [ M3 ]
    size_t * P3 = new size_t [M-M2];
    size_t * Q3 = new size_t [N2-R1];
//    if (M>8000) std::cerr<<"PLUQ 3"<<std::endl;

#pragma omp task shared(P3, Q3, G, Fi, R3) depend(inout:P3,Q3,G) depend(out: R3)
    R3 = pPLUQ (Fi, Diag, M-M2, N2-R1, G, lda, P3, Q3,nt/2);


    // [ H1 H2 ] <- P3^T H Q2^T
    // [ H3 H4 ]
    #pragma omp task shared(Fi, A4, Q2, P3) depend(inout:A4) depend(in:P3, Q2)
    {
    papplyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M-M2, 0, N-N2, A4, lda, Q2);
    papplyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N-N2, 0, M-M2, A4, lda, P3);
    }
    #pragma omp taskwait
    // [ E1 ] <- P3^T E
    // [ E2 ]
#pragma omp task shared( P3, A3, Fi) depend(in:P3) depend(inout:A3)
    papplyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, R1, 0, M-M2, A3, lda, P3);

    // [ M11 ] <- P2^T M1
    // [ M12 ]
#pragma omp task shared( P2, A, Fi) depend(inout:AR1lda) depend(in:P2)
    papplyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, R1, 0, M2-R1, AR1lda, lda, P2); // pas de dependence ?????

    // [ D1 D2 ] <- D Q2^T
#pragma omp task  shared( Q2, A2, Fi) depend(in:Q2) depend(inout:A2)
    papplyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, R1, 0, N-N2, A2, lda, Q2);

    // [ V1 V2 ] <- V1 Q3^T
#pragma omp task shared(Q3, A, Fi) depend(in:Q3) depend(inout:AR1)
    papplyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, R1, 0, N2-R1, AR1, lda, Q3);


      //#pragma omp task firstprivate(M2, R2, lda) shared(F, A4, Fi)
    // I <- H U2^-1
    // K <- H3 U2^-1
//    if (M>8000) std::cerr<<"TRSM I, K"<<std::endl;

#pragma omp task shared(Fi, A4, F, R2, PH) depend(in:F, R2) depend(inout:A4)
    ftrsm (Fi, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag, M-M2, R2, Fi.one, F, lda, A4, lda, PH);
    //    #pragma omp taskwait
    //#pragma omp task shared(temp, A4R2) depend(in:R2, R3) depend(out:temp, A4R2)


    typename Field::Element * temp = 0;
	/*
    for (size_t i=0; i<R3; ++i)
      fassign (Fi, R2, temp + i*R2, 1, A4 + i*lda, 1);
    */
#pragma omp task  shared(temp, Fi, R3, R2, A4) depend(in:A4, R2, R3) depend(out:temp)
    {
	    temp = new typename Field::Element [R3*R2];
	    FFLAS::fassign (Fi, R3, R2, A4 , lda, temp , R2);
    }
//    if (M>8000) std::cerr<<"TRSMs J, N"<<std::endl;
    // J <- L3^-1 I (in a temp)
#pragma omp task  shared(G, temp, Fi, R3, R2, A4) depend(in:A4, G, R2, R3) depend(inout:temp)
    ftrsm (Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R3, R2, Fi.one, G, lda, temp, R2, PH);
    // N <- L3^-1 H2
#pragma omp task  shared(G, A4R2, Fi, R3, R2) depend(in:G, R3, R2) depend(inout:A4R2)
    ftrsm (Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R3, N-N2-R2, Fi.one, G, lda, A4R2, lda, PH);
	// #pragma omp taskwait

    //    typename Field::Element * FR2 = F + R2;
    // O <- N - J V2
//    if (M>8000) std::cerr<<"GEMM O"<<std::endl;

#pragma omp task shared(F,R2, A4R2, Fi,temp, R3) depend(inout:A4R2,temp) depend(in:F,R2, R3)
    {
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, R3, N-N2-R2, R2, Fi.mOne, temp, R2, F+R2, lda, Fi.one, A4R2, lda, pWH);
	    delete[] temp;
	    temp=0;
    }
//    if (M>8000)  std::cerr<<"done GEMM O"<<std::endl;
    // R <- H4 - K V2 - M3 O
    typename Field::Element * R = 0;


#pragma omp task shared(A4R2, F, R2, R, Fi, R3) depend(inout:R) depend(in:F,R2, A4R2,R3)
    {
	    R = A4 + R2 + R3*lda;
	    typename Field::Element * A4R3lda = A4 + R3*lda;
	    typename Field::Element * GR3lda = G+ R3*lda;
//	    if (M>8000) std::cerr<<"GEMMs R"<<std::endl;

	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2-R3, N-N2-R2, R2, Fi.mOne, A4R3lda, lda, F+R2, lda, Fi.one, R, lda, pWH);
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2-R3, N-N2-R2, R3, Fi.mOne, GR3lda, lda, A4R2, lda, Fi.one, R, lda, pWH);
    }
#pragma omp taskwait // R2, R, P4 and Q4 are computed
    // H4 = P4 [ L4 ] [ U4 V4 ] Q4
    //         [ M4 ]
    size_t * P4 = new size_t [M-M2-R3];
    size_t * Q4 = new size_t [N-N2-R2];
    typename Field::Element * A3R3lda = A3 + R3*lda;
    typename Field::Element * A2R2 = A2 + R2;

//#pragma omp task shared(Fi, P4, Q4, R) depend(inout:R) depend(out:R4, P4, Q4)
//    if (M>8000) std::cerr<<"PLUQ 4"<<std::endl;
    R4 = pPLUQ (Fi, Diag, M-M2-R3, N-N2-R2, R, lda, P4, Q4,nt);

    // [ E21 M31 0 K1 ] <- P4^T [ E2 M3 0 K ]
    // [ E22 M32 0 K2 ]
#pragma omp task shared(A3R3lda, P4, Fi) depend(in:P4) depend(inout:A3R3lda)
    papplyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N2+R2, 0, M-M2-R3, A3R3lda, lda, P4);
    // [ D21 D22 ]     [ D2 ]
    // [ V21 V22 ]  <- [ V2 ] Q4^T
    // [  0   0  ]     [  0 ]
    // [ O1   O2 ]     [  O ]
#pragma omp task shared(A2R2, Q4, Fi) depend(in:Q4) depend(inout:A2R2)
    papplyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M2+R3, 0, N-N2-R2, A2R2, lda, Q4);
#pragma omp taskwait // To be removed?

    // P <- Diag (P1 [ I_R1    ] , P3 [ I_R3    ])
    //               [      P2 ]      [      P4 ]
    size_t* MathP = new size_t[M];
    composePermutationsP (MathP, P1, P2, R1, M2);
    composePermutationsP (MathP+M2, P3, P4, R3, M-M2);

    for (size_t i=M2; i<M; ++i)
      MathP[i] += M2;
    if (R1+R2 < M2){
      // P <- P S
#pragma omp task shared(MathP) depend(inout:MathP)
      PermApplyS (MathP, 1,1, M2, R1, R2, R3, R4);

      // A <-  S^T A
#pragma omp task shared(Fi, A) depend(inout:A)
      pMatrixApplyS(Fi, A, lda, N, M2, R1, R2, R3, R4);
      //      applyS (A, lda, N, M2, R1, R2, R3, R4);
    }

    // Q<- Diag ( [ I_R1    ] Q1,  [ I_R2    ] Q2 )
    //            [      Q3 ]      [      P4 ]
    size_t * MathQ = new size_t [N];
    composePermutationsQ (MathQ, Q1, Q3, R1, N2);
    composePermutationsQ (MathQ+N2, Q2, Q4, R2, N-N2);
    for (size_t i=N2; i<N; ++i)
      MathQ[i] += N2;

    if (R1 < N2){
//#pragma omp task shared(MathQ) depend(inout:MathQ)
		// Q <- T Q
	    PermApplyT (MathQ, 1,1,N2, R1, R2, R3, R4);
      // A <-   A T^T
#pragma omp task shared (Fi, A) depend(inout:A)
	    pMatrixApplyT( Fi, A, lda, M, N2, R1, R2, R3, R4);

		//      applyT (A, lda, M, N2, R1, R2, R3, R4);
    }
#pragma omp taskwait
    MathPerm2LAPACKPerm (Q, MathQ, N);
    MathPerm2LAPACKPerm (P, MathP, M);

    delete[] Q1;
    delete[] Q2;
    delete[] Q3;
    delete[] Q4;
    delete[] MathP;
    delete[] MathQ;
    delete[] P1;
    delete[] P2;
    delete[] P3;
    delete[] P4;
    return R1+R2+R3+R4;


#else // use __FFLASFFPACK_USE_OPENMP4

		  // [ B1 ] <- P1^T A2
		  // [ B2 ]
    TASK(MODE(READ(P1) REFERENCE(Fi, P1, A2) READWRITE(A2)), papplyP( Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N-N2, 0, M2, A2, lda, P1));
		  //applyP( Fi, FflasLeft, FflasNoTrans, N-N2, 0, M2, A2, lda, P1);

		  // [ C1 C2 ] <- A3 Q1^T
    TASK(MODE(READ( Q1) REFERENCE(Fi, Q1, A3) READWRITE(A3)), papplyP( Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M-M2, 0, N2, A3, lda, Q1));
		  //papplyP( Fi, FflasRight, FflasTrans, M-M2, 0, N2, A3, lda, Q1);

		  WAIT;
#if(DEBUG!=0)
struct timespec tsi;
clock_gettime(CLOCK_REALTIME, &tsi);
std::cerr << "R1 : " << tsi.tv_nsec << std::endl;
#endif
		  // D <- L1^-1 B1
 TASK(MODE(READ(A, R1, PH) REFERENCE(Fi, A, R1, PH, A2) READWRITE(A2)), ftrsm( Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R1, N-N2, Fi.one, A, lda, A2, lda, PH));
		  //    pftrsm( Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R1, N-N2, Fi.one, A, lda, A2 , lda,  method, NUM);
		  //ftrsm( Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R1, N-N2, Fi.one, A, lda, A2 , lda);

		  // E <- C1 U1^-1
 TASK(MODE(READ(R1, A, PH) REFERENCE(Fi, R1, A, PH, A3) READWRITE(A3)), ftrsm(Fi, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag, M-M2, R1, Fi.one, A, lda, A3, lda,  PH));
		  //pftrsm(Fi, FflasRight, FflasUpper, FflasNoTrans, Diag, M-M2, R1, Fi.one, A, lda, A3, lda,  method, NUM);
		  //ftrsm(Fi, FflasRight, FflasUpper, FflasNoTrans, Diag, M-M2, R1, Fi.one, A, lda, A3, lda);

		  WAIT;
#if(DEBUG!=0)
clock_gettime(CLOCK_REALTIME, &tsi);
std::cerr << "F  : " << tsi.tv_nsec << std::endl;
#endif
		  // F <- B2 - M1 D
 TASK(MODE(READ(A, pWH) REFERENCE(Fi, A, pWH, F, A2) WRITE(F) READWRITE(A2)),  fgemm( Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M2-R1, N-N2, R1, Fi.mOne, A + R1*lda, lda, A2, lda, Fi.one, A2+R1*lda, lda, pWH));
		  //fgemm( Fi, FflasNoTrans, FflasNoTrans, M2-R1, N-N2, R1, Fi.mOne, A + R1*lda, lda, A2, lda, Fi.one, A2+R1*lda, lda);

		  // G <- C2 - E V1
 TASK(MODE(READ(R1, A, pWH) REFERENCE(Fi, R1, G, A3, A, pWH) WRITE(G) READWRITE(A3)), fgemm( Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2, N2-R1, R1, Fi.mOne, A3, lda, A+R1, lda, Fi.one, A3+R1, lda, pWH));
		  //fgemm( Fi, FflasNoTrans, FflasNoTrans, M-M2, N2-R1, R1, Fi.mOne, A3, lda, A+R1, lda, Fi.one, A3+R1, lda);

		  WAIT;

#if(DEBUG!=0)
clock_gettime(CLOCK_REALTIME, &tsi);
std::cerr << "R2 : " << tsi.tv_nsec << std::endl;
#endif
		  size_t * P2 = FFLAS::fflas_new<size_t>(M2-R1);
		  size_t * Q2 = FFLAS::fflas_new<size_t>(N-N2);
		  // F = P2 [ L2 ] [ U2 V2 ] Q2
		  //        [ M2 ]
		  TASK(MODE(REFERENCE(Fi, R2, F, P2, Q2) WRITE(R2) READWRITE(F, P2, Q2)),
		       R2 = pPLUQ(Fi, Diag, M2-R1, N-N2, F, lda, P2, Q2,nt/2);
		       );
		  //R2 = PLUQ (Fi, Diag, M2-R1, N-N2, F, lda, P2, Q2);

		  size_t * P3 = FFLAS::fflas_new<size_t>(M-M2);
		  size_t * Q3 = FFLAS::fflas_new<size_t>(N2-R1);
		  // G = P3 [ L3 ] [ U3 V3 ] Q3
		  //        [ M3 ]
		  TASK(MODE(REFERENCE(Fi, R3, G, P3, Q3) WRITE(R3) READWRITE(G, P3, Q3)),
		       R3 = pPLUQ(Fi, Diag, M-M2, N2-R1, G, lda, P3, Q3,nt/2);
		       );
		  //		  R3 = PLUQ (Fi, Diag, M-M2, N2-R1, G, lda, P3, Q3);



		  // H <- A4 - ED
		  TASK(MODE(READ(M2, N2, R1, A3, A2) REFERENCE(Fi, M2, N2, R1, A3, A2, A4)  READWRITE(A4)), fgemm( Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2, N-N2, R1, Fi.mOne, A3, lda, A2, lda, Fi.one, A4, lda, pWH));
		  //fgemm( Fi, FflasNoTrans, FflasNoTrans, M-M2, N-N2, R1, Fi.mOne, A3, lda, A2, lda, Fi.one, A4, lda);
		  //		  std::cout<<"NUM "<<NUM_THREADS<<std::endl;

		  WAIT;
#if(DEBUG!=0)
clock_gettime(CLOCK_REALTIME, &tsi);
std::cerr << "aq2: " << tsi.tv_nsec << std::endl;
#endif

		  // [ H1 H2 ] <- P3^T H Q2^T
		  // [ H3 H4 ]
 TASK(MODE(REFERENCE(Fi, Q2, A4) READ(Q2) READWRITE(A4)), papplyP( Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M-M2, 0, N-N2, A4, lda, Q2));
		  //applyP( Fi, FflasRight, FflasTrans, M-M2, 0, N-N2, A4, lda, Q2);
		  WAIT;

		  TASK(MODE(READ(P3) REFERENCE(Fi, P3, A4) READWRITE(A4)), papplyP( Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N-N2, 0, M-M2, A4, lda, P3));
		  //applyP( Fi, FflasLeft, FflasNoTrans, N-N2, 0, M-M2, A4, lda, P3);
		  WAIT;

#if(DEBUG!=0)
clock_gettime(CLOCK_REALTIME, &tsi);
std::cerr << "ap3: " << tsi.tv_nsec << std::endl;
#endif
      // [ E1 ] <- P3^T E
      // [ E2 ]
 TASK(MODE(READ(P3) READWRITE(A3) REFERENCE(Fi, P3, A3)), papplyP(Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, R1, 0, M-M2, A3, lda, P3));
		  //applyP(  Fi, FflasLeft, FflasNoTrans, R1, 0, M-M2, A3, lda, P3);

      // [ M11 ] <- P2^T M1
      // [ M12 ]
 TASK(MODE(READ(P2) REFERENCE(Fi, P2, A) READWRITE(A)), papplyP(Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, R1, 0, M2-R1, A+R1*lda, lda, P2));
		  //applyP(Fi, FflasLeft, FflasNoTrans, R1, 0, M2-R1, A+R1*lda, lda, P2);

      // [ D1 D2 ] <- D Q2^T
 TASK(MODE(READ(Q2) REFERENCE(Fi, Q2, A2) READWRITE(A2)), papplyP(Fi, FFLAS::FflasRight, FFLAS::FflasTrans, R1, 0, N-N2, A2, lda, Q2));
		  //papplyP( Fi, FflasRight, FflasTrans, R1, 0, N-N2, A2, lda, Q2);

      // [ V1 V2 ] <- V1 Q3^T
 TASK(MODE(READ(Q3) REFERENCE(Fi, Q3, A) READWRITE(A)), papplyP( Fi, FFLAS::FflasRight, FFLAS::FflasTrans, R1, 0, N2-R1, A+R1, lda, Q3));
		  //applyP( Fi, FflasRight, FflasTrans, R1, 0, N2-R1, A+R1, lda, Q3);

      // I <- H1 U2^-1
      // K <- H3 U2^-1
 TASK(MODE(READ(R2, F) REFERENCE(Fi, R2, F, A4) READWRITE(A4)), ftrsm( Fi, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag, M-M2, R2, Fi.one, F, lda, A4, lda, PH));
		  //pftrsm( Fi, FflasRight, FflasUpper, FflasNoTrans, Diag, M-M2, R2, Fi.one, F, lda, A4, lda,  method, NUM);
		  //ftrsm( Fi, FflasRight, FflasUpper, FflasNoTrans, Diag, M-M2, R2, Fi.one, F, lda, A4, lda);
		  WAIT;

#if(DEBUG!=0)
clock_gettime(CLOCK_REALTIME, &tsi);
std::cerr << "R3 : " << tsi.tv_nsec << std::endl;
#endif
		  typename Field::Element_ptr temp = FFLAS::fflas_new (Fi, R3, R2);
		  /*    for (size_t i=0; i<R3; ++i)
			fassign (Fi, R2, temp + i*R2, 1, A4 + i*lda, 1);
		  */
		  FFLAS::fassign (Fi, R3, R2, A4 , lda, temp , R2);

    // J <- L3^-1 I (in a temp)
		  TASK(MODE(READ(R2, G) REFERENCE(Fi, R2, G, temp) READWRITE(temp)), ftrsm(  Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R3, R2, Fi.one, G, lda, temp, R2, PH));
		  //pftrsm( Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R3, R2, Fi.one, G, lda, temp, R2,  method, NUM);
		  //ftrsm( Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R3, R2, Fi.one, G, lda, temp, R2);

   // N <- L3^-1 H2
		  TASK(MODE(READ(R3, G) REFERENCE(Fi, R3, G, A4) READWRITE(A4)), ftrsm(Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R3, N-N2-R2, Fi.one, G, lda, A4+R2, lda, PH));
		  //    pftrsm(Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R3, N-N2-R2, Fi.one, G, lda, A4+R2, lda,  method, NUM);
		  //ftrsm(Fi, FflasLeft, FflasLower, FflasNoTrans, OppDiag, R3, N-N2-R2, Fi.one, G, lda, A4+R2, lda);
		  WAIT;

#if(DEBUG!=0)
clock_gettime(CLOCK_REALTIME, &tsi);
std::cerr << "mm3: " << tsi.tv_nsec << std::endl;
#endif
   // O <- N - J V2
 TASK(MODE(READ(R2, temp, F) REFERENCE(Fi, F, temp, R2, A4) READWRITE(A4)), fgemm( Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, R3, N-N2-R2, R2, Fi.mOne, temp, R2, F+R2, lda, Fi.one, A4+R2, lda, pWH));
		  //fgemm( Fi, FflasNoTrans, FflasNoTrans, R3, N-N2-R2, R2, Fi.mOne, temp, R2, F+R2, lda, Fi.one, A4+R2, lda);

		  typename Field::Element_ptr R = A4 + R2 + R3*lda;
  // R <- H4 - K V2
		  TASK(MODE(READ(R2, A4, F) REFERENCE(Fi, R2, A4, F, R) READWRITE(R)), fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2-R3, N-N2-R2, R2, Fi.mOne, A4+R3*lda, lda, F+R2, lda, Fi.one, R, lda, pWH));
		  //fgemm( Fi, FflasNoTrans, FflasNoTrans, M-M2-R3, N-N2-R2, R2, Fi.mOne, A4+R3*lda, lda, F+R2, lda, Fi.one, R, lda);
		  WAIT;

#if(DEBUG!=0)
clock_gettime(CLOCK_REALTIME, &tsi);
std::cerr << "mm4: " << tsi.tv_nsec << std::endl;
#endif
		  FFLAS::fflas_delete (temp);
    // R <- R - M3 O
		  TASK(MODE(READ(R3, A4, G) REFERENCE(Fi, R3, A4, G, R) READWRITE(R)), fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2-R3, N-N2-R2, R3, Fi.mOne, G+R3*lda, lda, A4+R2, lda, Fi.one, R, lda, pWH));
		  //fgemm( Fi, FflasNoTrans, FflasNoTrans, M-M2-R3, N-N2-R2, R3, Fi.mOne, G+R3*lda, lda, A4+R2, lda, Fi.one, R, lda);

		  WAIT;
#if(DEBUG!=0)
clock_gettime(CLOCK_REALTIME, &tsi);
std::cerr << "P4 : " << tsi.tv_nsec << std::endl;
#endif

		  size_t * P4 = FFLAS::fflas_new<size_t>(M-M2-R3);
		  size_t * Q4 = FFLAS::fflas_new<size_t>(N-N2-R2);

    // H4 = P4 [ L4 ] [ U4 V4 ] Q4
    //         [ M4 ]
		  //TASK(READ(Fi), NOWRITE(R4), READWRITE(R, P4, Q4), PPLUQ, R4, Fi, Diag, M-M2-R3, N-N2-R2, R, lda, P4, Q4);
		  R4 = pPLUQ (Fi, Diag, M-M2-R3, N-N2-R2, R, lda, P4, Q4,nt);

#if(DEBUG!=0)
clock_gettime(CLOCK_REALTIME, &tsi);
std::cerr << "R4 : " << tsi.tv_nsec << std::endl;
#endif

    // [ E21 M31 0 K1 ] <- P4^T [ E2 M3 0 K ]
    // [ E22 M32 0 K2 ]
 TASK(MODE(READ(P4) REFERENCE(Fi, P4, A3) READWRITE(A3)), papplyP( Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N2+R2, 0, M-M2-R3, A3+R3*lda, lda, P4));
		  //applyP( Fi, FflasLeft, FflasNoTrans, N2+R2, 0, M-M2-R3, A3+R3*lda, lda, P4);

    // [ D21 D22 ]     [ D2 ]
    // [ V21 V22 ]  <- [ V2 ] Q4^T
    // [  0   0  ]     [  0 ]
    // [ O1   O2 ]     [  O ]
 TASK(MODE(READ(Q4) REFERENCE(Fi, A2, Q4) READWRITE(A2)), papplyP( Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M2+R3, 0, N-N2-R2, A2+R2, lda, Q4));
		  //applyP( Fi, FflasRight, FflasTrans, M2+R3, 0, N-N2-R2, A2+R2, lda, Q4);

		  size_t* MathP = new size_t[M];
    // P <- Diag (P1 [ I_R1    ] , P3 [ I_R3    ])
    //               [      P2 ]      [      P4 ]
		  composePermutationsP (MathP, P1, P2, R1, M2);
		  composePermutationsP (MathP+M2, P3, P4, R3, M-M2);
		  for (size_t i=M2; i<M; ++i)
			  MathP[i] += M2;

		  WAIT;
#if(DEBUG!=0)
clock_gettime(CLOCK_REALTIME, &tsi);
std::cerr << "M2 : " << tsi.tv_nsec << std::endl;
#endif
		  if (R1+R2 < M2){
			  // A <-  S^T A
			  TASK(MODE(READ(R1, R2, R3, R4) REFERENCE(Fi, R1, R2, R3, R4, A) READWRITE(A)), pMatrixApplyS( Fi, A, lda, N, M2, R1, R2, R3, R4));
			  //MatrixApplyS(Fi, A, lda, N, M2, R1, R2, R3, R4);

			  // P <- P S
			  PermApplyS( MathP, 1,1, M2, R1, R2, R3, R4);
		  }

		  // Q<- Diag ( [ I_R1    ] Q1,  [ I_R2    ] Q2 )
		  //            [      Q3 ]      [      P4 ]
		  size_t * MathQ = FFLAS::fflas_new<size_t>(N);
		  composePermutationsQ (MathQ, Q1, Q3, R1, N2);
		  composePermutationsQ (MathQ+N2, Q2, Q4, R2, N-N2);
		  for (size_t i=N2; i<N; ++i)
			  MathQ[i] += N2;

		  if (R1 < N2){
			  // Q <- T Q
			  PermApplyT (MathQ, 1,1,N2, R1, R2, R3, R4);
			  WAIT;
#if(DEBUG!=0)
clock_gettime(CLOCK_REALTIME, &tsi);
std::cerr << "pMatAppT: " << tsi.tv_nsec << std::endl;
#endif
			  // A <-   A T^T
 TASK(MODE(READ(R1, R2, R3, R4) REFERENCE(Fi, A, R1, R2, R3, R4) READWRITE(A)), pMatrixApplyT( Fi, A, lda, M, N2, R1, R2, R3, R4));
			  //			  MatrixApplyT(Fi, A, lda, M, N2, R1, R2, R3, R4);
		  }
		  MathPerm2LAPACKPerm (Q, MathQ, N);
		  MathPerm2LAPACKPerm (P, MathP, M);
		  WAIT;
#if(DEBUG!=0)
clock_gettime(CLOCK_REALTIME, &tsi);
std::cerr << "de : " << tsi.tv_nsec << std::endl;
#endif

    FFLAS::fflas_delete( MathQ);
    FFLAS::fflas_delete( MathP);
    FFLAS::fflas_delete( P1);
    FFLAS::fflas_delete( P2);
    FFLAS::fflas_delete( P3);
    FFLAS::fflas_delete( P4);
    FFLAS::fflas_delete( Q1);
    FFLAS::fflas_delete( Q2);
    FFLAS::fflas_delete( Q3);
    FFLAS::fflas_delete( Q4);
    return R1+R2+R3+R4;
#endif
	  }

}// namespace FFPACK

#endif
#endif // __FFLASFFPACK_ffpack_ppluq_INL
