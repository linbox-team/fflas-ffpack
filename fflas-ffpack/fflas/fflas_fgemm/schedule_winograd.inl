/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2014 the LinBox group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *            Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *            Ziad Sultan <ziad.sultan@imag.fr>
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

/** @file fflas/fflas_fgemm/winograd.inl
 * @ingroup MMalgos
 * @brief Winograd implementation
 * @bib ISSAC09 Scheduling
 */

#ifndef __FFLASFFPACK_fgemm_winograd_INL
#define __FFLASFFPACK_fgemm_winograd_INL

namespace FFLAS { namespace BLAS3 {

        template < class Field, class FieldTrait, class Strat, class Param >
        inline typename Field::Element_ptr
        WinoPar (const Field& F,
				 const FFLAS_TRANSPOSE ta,
				 const FFLAS_TRANSPOSE tb,
				 const size_t mr, const size_t nr, const size_t kr,
				 const typename Field::Element alpha,
				 typename Field::ConstElement_ptr A,const size_t lda,
				 typename Field::ConstElement_ptr B,const size_t ldb,
				 const typename Field::Element  beta,
				 typename Field::Element_ptr C, const size_t ldc,
					 // const size_t kmax, const size_t w, const FFLAS_BASE base
                 MMHelper<Field, MMHelperAlgo::WinogradPar, FieldTrait, ParSeqHelper::Parallel<Strat,Param> > & WH
				 )
        {
            FFLASFFPACK_check(F.isZero(beta));

				//          typedef MMHelper<Field, MMHelperAlgo::WinogradPar, FieldTrait > MMH_t;
            typedef MMHelper<Field, MMHelperAlgo::WinogradPar, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Recursive,StrategyParameter::TwoDAdaptive> > MMH_t;
            const typename MMH_t::DelayedField & DF = WH.delayedField;
            typedef typename  MMH_t::DelayedField::Element DFElt;

            size_t lb, cb, la, ca, ldX2;
                // size_t x3rd = std::max(mr,kr);
            typename Field::ConstElement_ptr A11=A, A12, A21, A22;
            typename Field::ConstElement_ptr B11=B, B12, B21, B22;
            typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;

            size_t x1rd = std::max(nr,kr);
            size_t ldX1;
            if (ta == FflasTrans) {
                A21 = A + mr;
                A12 = A + kr*lda;
                A22 = A12 + mr;
                la = kr;
                ca = mr;
                ldX1 = mr;
            } else {
                A12 = A + kr;
                A21 = A + mr*lda;
                A22 = A21 + kr;
                la = mr;
                ca = kr;
                ldX1  = x1rd;
            }
            if (tb == FflasTrans) {
                B21 = B + kr;
                B12 = B + nr*ldb;
                B22 = B12 + kr;
                lb = nr;
                cb = kr;
                ldX2 = kr;
            } else {
                B12 = B + nr;
                B21 = B + kr*ldb;
                B22 = B21 + nr;
                lb = kr;
                ldX2 = cb = nr;
            }

                // 11 temporary submatrices are required
            typename Field::Element_ptr X21 = fflas_new (F, kr, nr);
            typename Field::Element_ptr X11 = fflas_new (F,mr,x1rd);
        
            typename Field::Element_ptr X22 = fflas_new (F, kr, nr);
            typename Field::Element_ptr X12 = fflas_new (F,mr,x1rd);

            typename Field::Element_ptr X23 = fflas_new (F, kr, nr);
            typename Field::Element_ptr X13 = fflas_new (F,mr,x1rd);

            typename Field::Element_ptr X24 = fflas_new (F, kr, nr);
            typename Field::Element_ptr X14 = fflas_new (F,mr,x1rd);
            typename Field::Element_ptr X15 = fflas_new (F,mr,x1rd);

            typename Field::Element_ptr C_11 = fflas_new (F,mr,nr);
            typename Field::Element_ptr CC_11 = fflas_new (F,mr,nr);
            SYNCH_GROUP(
                        
					// T3 = B22 - B12 in X21  and S3 = A11 - A21 in X11
				TASK(MODE(READ(B22, B12) WRITE(X21) CONSTREFERENCE(DF)),
					 pfsub(DF,lb,cb,B22,ldb,B12,ldb,X21,ldX2, NUM_THREADS););
				TASK(MODE(READ(A11, A21) WRITE(X11) CONSTREFERENCE(DF)),
					 pfsub(DF,la,ca,A11,lda,A21,lda,X11,ldX1, NUM_THREADS););

					// T1 = B12 - B11 in X22 and  S1 = A21 + A22 in X12
				TASK(MODE(READ(B11, B12) WRITE(X22) CONSTREFERENCE(DF)),
					 pfsub(DF,lb,cb,B12,ldb,B11,ldb,X22,ldX2, NUM_THREADS););
				TASK(MODE(READ(A12, A22) WRITE(X12) CONSTREFERENCE(DF)),
					 pfadd(DF,la,ca,A21,lda,A22,lda,X12,ldX1, NUM_THREADS););

				CHECK_DEPENDENCIES;

					// T2 = B22 - T1 in X23 and  S2 = S1 - A11 in X13
				TASK(MODE(READ(B22, X22) READWRITE(X23) CONSTREFERENCE(DF)),
					 pfsub(DF,lb,cb,B22,ldb,X22,ldX2,X23,ldX2, NUM_THREADS););
				TASK(MODE(READ(A11, X12) READWRITE(X13) CONSTREFERENCE(DF)),
						 //          fsub(DF,la,ca,A11,lda,X12,ldX1,X13,ldX1););
					 pfsub(DF,la,ca,X12,ldX1,A11,lda,X13,ldX1, NUM_THREADS););
					/*
					  fsub(DF,lb,cb,B22,ldb,X2,ldX2,X2,ldX2);
					  fsubin(DF,la,ca,A11,lda,X1,ldX1););
					*/
				CHECK_DEPENDENCIES;

					// T4 = T2 - B21 in X2 and S4 = A12 -S2 in X1
				TASK(MODE(READ(B21, X23) READWRITE(X24) CONSTREFERENCE(DF)),
						 //          fsub(DF,lb,cb,B21,ldb,X23,ldX2,X24,ldX2);
					 pfsub(DF,lb,cb,X23,ldX2,B21,ldb,X24,ldX2, NUM_THREADS););
				TASK(MODE(READ(A12, X13) READWRITE(X14) CONSTREFERENCE(DF)),
					 pfsub(DF,la,ca,A12,lda,X13,ldX1,X14,ldX1, NUM_THREADS););

					/*
					  fsubin(DF,lb,cb,B21,ldb,X2,ldX2);
					  fsub(DF,la,ca,A12,lda,X1,ldX1,X1,ldX1););
					*/
				CHECK_DEPENDENCIES;

					// P1 = alpha . A11 * B11 in X1

				MMH_t H1(F, WH.recLevel-1, WH.Amin, WH.Amax, WH.Bmin, WH.Bmax, 0, 0);
				MMH_t H7(F, WH.recLevel-1, -(WH.Amax-WH.Amin), WH.Amax-WH.Amin, -(WH.Bmax-WH.Bmin), WH.Bmax-WH.Bmin, 0,0);
				MMH_t H5(F, WH.recLevel-1, 2*WH.Amin, 2*WH.Amax, -(WH.Bmax-WH.Bmin), WH.Bmax-WH.Bmin, 0, 0);
				MMH_t H6(F, WH.recLevel-1, 2*WH.Amin-WH.Amax, 2*WH.Amax-WH.Amin, 2*WH.Bmin-WH.Bmax, 2*WH.Bmax-WH.Bmin, 0, 0);
				MMH_t H3(F, WH.recLevel-1, 2*WH.Amin-2*WH.Amax, 2*WH.Amax-2*WH.Amin, WH.Bmin, WH.Bmax, 0, 0);
				MMH_t H4(F, WH.recLevel-1, WH.Amin, WH.Amax, 2*WH.Bmin-2*WH.Bmax, 2*WH.Bmax-2*WH.Bmin, 0, 0);
				MMH_t H2(F, WH.recLevel-1, WH.Amin, WH.Amax, WH.Bmin, WH.Bmax, 0, 0);

				size_t nt = WH.parseq.numthreads();
				size_t nt_rec = nt/7;
				size_t nt_mod = nt % 7 ;
				H1.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
				H2.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
				H3.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
				H4.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
				H5.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
				H6.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
				H7.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));

				TASK(MODE(READ(A11, B11) WRITE(X15) CONSTREFERENCE(F,H1)),
					 fgemm (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X15, x1rd, H1););
					// P7 = alpha . S3 * T3  in C21
				TASK(MODE(READ(X11, X21) WRITE(C21) CONSTREFERENCE(F,H7)),
					 fgemm (F, ta, tb, mr, nr, kr, alpha, X11, ldX1, X21, ldX2, F.zero, C21, ldc, H7););

					// P5 = alpha . S1*T1 in C22
				TASK(MODE(READ(X12, X22) WRITE(C22) CONSTREFERENCE(F,H5)),
					 fgemm (F, ta, tb, mr, nr, kr, alpha, X12, ldX1, X22, ldX2, F.zero, C22, ldc, H5););

					// P6 = alpha . S2 * T2 in C12
				TASK(MODE(READ(X13, X23) WRITE(C12) CONSTREFERENCE(F,H6)),
					 fgemm (F, ta, tb, mr, nr, kr, alpha, X13, ldX1, X23, ldX2, F.zero, C12, ldc, H6););

					// P3 = alpha . S4*B22 in CC_11
				TASK(MODE(READ(X14, B22) WRITE(CC_11) CONSTREFERENCE(F,H3)),
					 fgemm (F, ta, tb, mr, nr, kr, alpha, X14, ldX1, B22, ldb, F.zero, CC_11, nr, H3););

					// P4 = alpha . A22 * T4 in C_11
				TASK(MODE(READ(A22) WRITE(C_11) READWRITE(X24, X22, X23, X21) CONSTREFERENCE(F,H4)),
					 fgemm (F, ta, tb, mr, nr, kr, alpha, A22, lda, X24, ldX2, F.zero, C_11, nr, H4);
					 );

					// P2 = alpha . A12 * B21  in C11
				TASK(MODE(READ(A12, B21) WRITE(C11) CONSTREFERENCE(F,H2)),
					 fgemm (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.zero, C11, ldc, H2););
				CHECK_DEPENDENCIES;

				DFElt U2Min, U2Max;
				DFElt U3Min, U3Max;
				DFElt U4Min, U4Max;
				DFElt U7Min, U7Max;
				DFElt U5Min, U5Max;
					// U2 = P1 + P6 in C12  and
					// U3 = P7 + U2 in C21  and
					// U4 = P5 + U2 in C12    and
					// U7 = P5 + U3 in C22    and
					// U5 = P3 + U4 in C12
					// BIG TASK with 5 Addin function calls
//      TASK(MODE(READWRITE(X15, C12) CONSTREFERENCE(F, DF, WH, U2Min, U2Max, H1.Outmin, H1.Outmax, H6.Outmin, H6.Outmax)),
				if (Protected::NeedPreAddReduction(U2Min, U2Max, H1.Outmin, H1.Outmax, H6.Outmin, H6.Outmax, WH)){
					TASK(MODE(READWRITE(X15) CONSTREFERENCE(F)),
						 pfreduce (F, mr, x1rd, X15, x1rd, NUM_THREADS);
						 );
					TASK(MODE(READWRITE(C12) CONSTREFERENCE(F)),
						 pfreduce (F, mr, nr, C12, ldc, NUM_THREADS);
						 );
					CHECK_DEPENDENCIES;
				}
				TASK(MODE(READWRITE(X15, C12) CONSTREFERENCE(DF)),
					 pfaddin(DF,mr,nr,X15,x1rd,C12,ldc, NUM_THREADS);
					 );
				CHECK_DEPENDENCIES;
//      TASK(MODE(READWRITE(C12, C21) CONSTREFERENCE(F, DF, WH, U3Min, U3Max, U2Min, U2Max)),
				if (Protected::NeedPreAddReduction(U3Min, U3Max, U2Min, U2Max, H7.Outmin, H7.Outmax, WH)){
					TASK(MODE(READWRITE(C12) CONSTREFERENCE(F)),
						 pfreduce (F, mr, nr, C12, ldc, NUM_THREADS);
						 );
					TASK(MODE(READWRITE(C21) CONSTREFERENCE(F)),
						 pfreduce (F, mr, nr, C21, ldc, NUM_THREADS);
						 );
					CHECK_DEPENDENCIES;
				}
				TASK(MODE(READWRITE(C12, C21) CONSTREFERENCE(DF)),
					 pfaddin(DF,mr,nr,C12,ldc,C21,ldc, NUM_THREADS);
					 );
				CHECK_DEPENDENCIES;
//      TASK(MODE(READWRITE(C12, C22) CONSTREFERENCE(F, DF, WH) VALUE(U4Min, U4Max, U2Min, U2Max)),
				if (Protected::NeedPreAddReduction(U4Min, U4Max, U2Min, U2Max, H5.Outmin, H5.Outmax, WH)){
					TASK(MODE(READWRITE(C22) CONSTREFERENCE(F)),
						 pfreduce (F, mr, nr, C22, ldc, NUM_THREADS);
						 );
					TASK(MODE(READWRITE(C12) CONSTREFERENCE(F)),
						 pfreduce (F, mr, nr, C12, ldc, NUM_THREADS);
						 );
					CHECK_DEPENDENCIES;
				}
				TASK(MODE(READWRITE(C12, C22) CONSTREFERENCE(DF, WH)),
					 pfaddin(DF,mr,nr,C22,ldc,C12,ldc, NUM_THREADS);
					 );
				CHECK_DEPENDENCIES;
//      TASK(MODE(READWRITE(C22, C21) CONSTREFERENCE(F, DF, WH) VALUE(U3Min, U3Max, U7Min, U7Max)),
				if (Protected::NeedPreAddReduction (U7Min,U7Max, U3Min, U3Max, H5.Outmin,H5.Outmax, WH) ){
					TASK(MODE(READWRITE(C21) CONSTREFERENCE(F)),
						 pfreduce (F, mr, nr, C21, ldc, NUM_THREADS);
						 );
					TASK(MODE(READWRITE(C22) CONSTREFERENCE(F)),
						 pfreduce (F, mr, nr, C22, ldc, NUM_THREADS);
						 );
					CHECK_DEPENDENCIES;
				}
				TASK(MODE(READWRITE(C22, C21) CONSTREFERENCE(DF, WH)),
					 pfaddin(DF,mr,nr,C21,ldc,C22,ldc, NUM_THREADS);
					 );
//      TASK(MODE(READWRITE(C12, CC_11) CONSTREFERENCE(F, DF, WH) VALUE(U5Min, U5Max, U4Min, U4Max)),
				if (Protected::NeedPreAddReduction (U5Min,U5Max, U4Min, U4Max, H3.Outmin, H3.Outmax, WH) ){
					TASK(MODE(READWRITE(C12) CONSTREFERENCE(F)),
						 pfreduce (F, mr, nr, C12, ldc, NUM_THREADS);
						 );
					TASK(MODE(READWRITE(CC_11) CONSTREFERENCE(F)),
						 pfreduce (F, mr, nr, CC_11, nr, NUM_THREADS);
						 );
					CHECK_DEPENDENCIES;
				}
				TASK(MODE(READWRITE(C12, CC_11) CONSTREFERENCE(DF, WH)),
					 pfaddin(DF,mr,nr,CC_11,nr,C12,ldc, NUM_THREADS);
					 );
				CHECK_DEPENDENCIES;

					// U6 = U3 - P4 in C21
				DFElt U6Min, U6Max;
//      TASK(MODE(READWRITE(C_11, C21) CONSTREFERENCE(F, DF, WH) VALUE(U6Min, U6Max, U3Min, U3Max)),
				if (Protected::NeedPreSubReduction (U6Min,U6Max, U3Min, U3Max, H4.Outmin,H4.Outmax, WH) ){
					TASK(MODE(READWRITE(CC_11) CONSTREFERENCE(F)),
						 pfreduce (F, mr, nr, C_11, nr, NUM_THREADS);
						 );
					TASK(MODE(READWRITE(C21) CONSTREFERENCE(F)),
						 pfreduce (F, mr, nr, C21, ldc, NUM_THREADS);
						 );
					CHECK_DEPENDENCIES
						}
				TASK(MODE(READWRITE(C_11, C21) CONSTREFERENCE(DF, WH) ),
					 pfsubin(DF,mr,nr,C_11,nr,C21,ldc, NUM_THREADS);
					 );
        
					//CHECK_DEPENDENCIES;

					//  U1 = P2 + P1 in C11
				DFElt U1Min, U1Max;
//      TASK(MODE(READWRITE(C11, X15/*, X14, X13, X12, X11*/) CONSTREFERENCE(F, DF, WH) VALUE(U1Min, U1Max)),
				if (Protected::NeedPreAddReduction (U1Min, U1Max, H1.Outmin, H1.Outmax, H2.Outmin,H2.Outmax, WH) ){
					TASK(MODE(READWRITE(X15) CONSTREFERENCE(F)),
						 pfreduce (F, mr, nr, X15, x1rd, NUM_THREADS);
						 );
					TASK(MODE(READWRITE(C11) CONSTREFERENCE(F)),
						 pfreduce (F, mr, nr, C11, ldc, NUM_THREADS);
						 );
					CHECK_DEPENDENCIES
						}
				TASK(MODE(READWRITE(C11, X15) CONSTREFERENCE(DF, WH)),
					 pfaddin(DF,mr,nr,X15,x1rd,C11,ldc, NUM_THREADS);
					 );

				WH.Outmin = std::min (U1Min, std::min (U5Min, std::min (U6Min, U7Min)));
				WH.Outmax = std::max (U1Max, std::max (U5Max, std::max (U6Max, U7Max)));

                        );
//          WAIT;
            
        
            fflas_delete (CC_11);
            fflas_delete (C_11);
            fflas_delete (X15);
            fflas_delete (X14);
            fflas_delete (X24);
            fflas_delete (X13);
            fflas_delete (X23);
            fflas_delete (X12);
            fflas_delete (X22);
            fflas_delete (X11);
            fflas_delete (X21);
             
            return C;
        } //wino parallel


        template < class Field, class FieldTrait >
        inline void Winograd (const Field& F,
                              const FFLAS_TRANSPOSE ta,
                              const FFLAS_TRANSPOSE tb,
                              const size_t mr, const size_t nr, const size_t kr,
                              const typename Field::Element alpha,
                              typename Field::ConstElement_ptr A,const size_t lda,
                              typename Field::ConstElement_ptr B,const size_t ldb,
                              const typename Field::Element  beta,
                              typename Field::Element_ptr C, const size_t ldc,
                                  // const size_t kmax, const size_t w, const FFLAS_BASE base
                              MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait> & WH
                              )
        {

            FFLASFFPACK_check(F.isZero(beta));

            typedef MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > MMH_t;
            //     //typedef typename MMgr::DelayedField::Element_ptr DFEptr;
            //     //typedef typename MMgr::DelayedField::ConstElement_ptr DFCEptr;
            // typedef typename MMgr::DFElt DFElt;
            // typedef typename MMgr::DFEptr DFEptr;

            // const typename MMgr::DelayedField & DF = WH.ModeManager.delayedField;
            typename MMH_t::ModeMgr_t& WHMM = WH.ModeManager;
            size_t lb, cb, la, ca, ldX2;
                // size_t x3rd = std::max(mr,kr);
            typename Field::ConstElement_ptr A11=A, A12, A21, A22;
            typename Field::ConstElement_ptr B11=B, B12, B21, B22;
            typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;

            size_t x1rd = std::max(nr,kr);
            size_t ldX1;
            if (ta == FflasTrans) {
                A21 = A + mr;
                A12 = A + kr*lda;
                A22 = A12 + mr;
                la = kr;
                ca = mr;
                ldX1 = mr;
            } else {
                A12 = A + kr;
                A21 = A + mr*lda;
                A22 = A21 + kr;
                la = mr;
                ca = kr;
                ldX1  = x1rd;
            }
            if (tb == FflasTrans) {
                B21 = B + kr;
                B12 = B + nr*ldb;
                B22 = B12 + kr;
                lb = nr;
                cb = kr;
                ldX2 = kr;
            } else {
                B12 = B + nr;
                B21 = B + kr*ldb;
                B22 = B21 + nr;
                lb = kr;
                ldX2 = cb = nr;
            }
                // Two temporary submatrices are required

                // T3 = B22 - B12 in X2
            typename Field::Element_ptr X2 = fflas_new (F, kr, nr);
            AddSubHelper<Field, ModeCategories::LazyTag> T3H(F, WHMM.B, WHMM.B);
            fsub (F, lb, cb,  B22, ldb,  B12, ldb, X2, ldX2, T3H);

                // S3 = A11 - A21 in X1
            typename Field::Element_ptr X1 = fflas_new (F,mr,x1rd);
            AddSubHelper<Field, ModeCategories::LazyTag> S3H(F, WHMM.A, WHMM.A);
            fsub (F, la, ca, A11, lda, A21, lda, X1, ldX1, S3H);

                // P7 = alpha . S3 * T3  in C21
            MMH_t H7(F, WH.AlgoManager.recLevel-1, S3H.Out, T3H.Out);
            fgemm (F, ta, tb, mr, nr, kr, alpha, X1, ldX1, X2, ldX2, F.zero, C21, ldc, H7);

                // T1 = B12 - B11 in X2
            AddSubHelper<Field, ModeCategories::LazyTag> T1H(F, WHMM.B, WHMM.B);
            fsub (F, lb, cb, B12, ldb, B11, ldb, X2, ldX2, T1H);

                // S1 = A21 + A22 in X1
            AddSubHelper<Field, ModeCategories::LazyTag> S1H(F, WHMM.A, WHMM.A);
            fadd (F, la, ca, A21, lda, A22, lda, X1, ldX1, S1H);

                // P5 = alpha . S1*T1 in C22
            MMH_t H5(F, WH.AlgoManager.recLevel-1, S1H.Out, T1H.Out);
            fgemm (F, ta, tb, mr, nr, kr, alpha, X1, ldX1, X2, ldX2, F.zero, C22, ldc, H5);

                // T2 = B22 - T1 in X2
            AddSubHelper<Field, ModeCategories::LazyTag> T2H(F, WHMM.B, T1H.Out);
            fsub (F, lb, cb, B22, ldb, X2, ldX2, X2, ldX2, T2H);

                // S2 = S1 - A11 in X1
            AddSubHelper<Field, ModeCategories::LazyTag> S2H(F, S1H.Out, WHMM.A);
            fsubin(F, la, ca, A11, lda, X1, ldX1, S2H);

                // P6 = alpha . S2 * T2 in C12
            MMH_t H6(F, WH.AlgoManager.recLevel-1, S2H.Out, T2H.Out);
            fgemm (F, ta, tb, mr, nr, kr, alpha, X1, ldX1, X2, ldX2, F.zero, C12, ldc, H6);

                // S4 = A12 -S2 in X1
            AddSubHelper<Field, ModeCategories::LazyTag> S4H(F, WHMM.A, S2H.Out);
            fsub (F, la, ca, A12, lda, X1, ldX1, X1, ldX1, S4H);

                // P3 = alpha . S4*B22 in C11
            MMH_t H3 (F, WH.AlgoManager.recLevel-1, S4H.Out, WHMM.B);
            fgemm (F, ta, tb, mr, nr, kr, alpha, X1, ldX1, B22, ldb, F.zero, C11, ldc, H3);

                // P1 = alpha . A11 * B11 in X1
            MMH_t H1(F, WH.AlgoManager.recLevel-1, WHMM.A, WHMM.B);
            fgemm (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X1, nr, H1);

                // U2 = P1 + P6 in C12  and
            // DFElt U2Min, U2Max;
            //     // This test will be optimized out
            // if (Protected::NeedPreAddReduction(U2Min, U2Max, H1.ModeManager.Outmin, H1.ModeManager.Outmax, H6.ModeManager.Outmin, H6.ModeManager.Outmax, WH)){
            //     freduce (F, mr, nr, X1, nr);
            //     freduce (F, mr, nr, C12, ldc);
            // }
            AddSubHelper<Field,ModeCategories::LazyTag> U2H (F, H1.ModeManager.Out, H6.ModeManager.Out);
            faddin (F, mr, nr, X1, nr, C12, ldc, U2H);

                // U3 = P7 + U2 in C21  and
            // DFElt U3Min, U3Max;
            //     // This test will be optimized out
            // if (Protected::NeedPreAddReduction(U3Min, U3Max, U2Min, U2Max, H7.ModeManager.Outmin, H7.ModeManager.Outmax, WH)){
            //     freduce (F, mr, nr, C12, ldc);
            //     freduce (F, mr, nr, C21, ldc);
            // }
            AddSubHelper<Field,ModeCategories::LazyTag> U3H (F, H7.ModeManager.Out, U2H.Out);
            faddin (F, mr, nr, C12, ldc, C21, ldc, U3H);


                // U4 = P5 + U2 in C12    and
            // DFElt U4Min, U4Max;
            //     // This test will be optimized out
            // if (Protected::NeedPreAddReduction(U4Min, U4Max, U2Min, U2Max, H5.ModeManager.Outmin, H5.ModeManager.Outmax, WH)){
            //     freduce (F, mr, nr, C22, ldc);
            //     freduce (F, mr, nr, C12, ldc);
            // }
            AddSubHelper<Field,ModeCategories::LazyTag> U4H (F, H5.ModeManager.Out, U2H.Out);
            faddin (F, mr, nr, C22, ldc, C12, ldc, U4H);

                // U7 = P5 + U3 in C22    and
            // DFElt U7Min, U7Max;
            //     // This test will be optimized out
            // if (Protected::NeedPreAddReduction (U7Min,U7Max, U3Min, U3Max, H5.ModeManager.Outmin,H5.ModeManager.Outmax, WH) ){
            //     freduce (F, mr, nr, C21, ldc);
            //     freduce (F, mr, nr, C22, ldc);
            // }
            AddSubHelper<Field,ModeCategories::LazyTag> U7H (F, H5.ModeManager.Out, U3H.Out);
            faddin (F, mr, nr, C21, ldc, C22, ldc, U7H);

                // U5 = P3 + U4 in C12
            // DFElt U5Min, U5Max;
            //     // This test will be optimized out
            // if (Protected::NeedPreAddReduction (U5Min,U5Max, U4Min, U4Max, H3.ModeManager.Outmin, H3.ModeManager.Outmax, WH) ){
            //     freduce (F, mr, nr, C12, ldc);
            //     freduce (F, mr, nr, C11, ldc);
            // }
            AddSubHelper<Field,ModeCategories::LazyTag> U5H (F, H3.ModeManager.Out, U4H.Out);
            faddin (F, mr, nr, C11, ldc, C12, ldc, U5H);

                // T4 = T2 - B21 in X2
            AddSubHelper<Field,ModeCategories::LazyTag> T4H (F, T2H.Out, WHMM.B);
            fsubin (F, lb, cb, B21, ldb, X2, ldX2, T4H);

                // P4 = alpha . A22 * T4 in C11
            MMH_t H4(F, WH.AlgoManager.recLevel-1, WHMM.A,  T4H.Out);
            fgemm (F, ta, tb, mr, nr, kr, alpha, A22, lda, X2, ldX2, F.zero, C11, ldc, H4);

            fflas_delete (X2);

                // U6 = U3 - P4 in C21
            // DFElt U6Min, U6Max;
            //     // This test will be optimized out
            // if (Protected::NeedPreSubReduction (U6Min,U6Max, U3Min, U3Max, H4.ModeManager.Outmin,H4.ModeManager.Outmax, WH) ){
            //     freduce (F, mr, nr, C11, ldc);
            //     freduce (F, mr, nr, C21, ldc);
            // }
            AddSubHelper<Field,ModeCategories::LazyTag> U6H (F, U3H.Out, H4.ModeManager.Out);
            fsubin (F, mr, nr, C11, ldc, C21, ldc, U6H);

                // P2 = alpha . A12 * B21  in C11
            MMH_t H2 (F, WH.AlgoManager.recLevel-1, WHMM.A, WHMM.B);
            fgemm (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.zero, C11, ldc, H2);

                //  U1 = P2 + P1 in C11
            // DFElt U1Min, U1Max;
            //     // This test will be optimized out
            // if (Protected::NeedPreAddReduction (U1Min, U1Max, H1.ModeManager.Outmin, H1.ModeManager.Outmax, H2.ModeManager.Outmin,H2.ModeManager.Outmax, WH) ){
            //     freduce (F, mr, nr, X1, nr);
            //     freduce (F, mr, nr, C11, ldc);
            // }
            AddSubHelper<Field,ModeCategories::LazyTag> U1H (F, H2.ModeManager.Out, H1.ModeManager.Out);
            faddin (F, mr, nr, X1, nr, C11, ldc, U1H);

            fflas_delete (X1);

            WH.ModeManager.setOutBoundsMM(U1H, U5H, U6H, U7H);
            
        } // Winograd

    } // BLAS3


} // FFLAS

#endif // __FFLASFFPACK_fgemm_winograd_INL

