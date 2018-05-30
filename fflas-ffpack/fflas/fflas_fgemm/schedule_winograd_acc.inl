/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* Copyright (C) 2014 the LinBox group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

/** @file fflas/fflas_fgemm/winograd_acc.inl
 * @ingroup MMalgos
 * @brief Winograd implementation
 * @bib ISSAC09 Scheduling
 */

#ifndef __FFLASFFPACK_fgemm_winograd_acc_INL
#define __FFLASFFPACK_fgemm_winograd_acc_INL

namespace FFLAS { namespace BLAS3 {


	// 3 temps and 23 ops
        // TODO: Add check for modular reductions before final additions
	template < class Field,class FieldTrait >
	inline void WinogradAcc_3_23 (const Field& F,
				      const FFLAS_TRANSPOSE ta,
				      const FFLAS_TRANSPOSE tb,
				      const size_t mr, const size_t nr, const size_t kr,
				      const typename Field::Element alpha,
				      typename Field::ConstElement_ptr A,const size_t lda,
				      typename Field::ConstElement_ptr B,const size_t ldb,
				      const typename Field::Element  beta,
				      typename Field::Element_ptr C, const size_t ldc,
				      MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > & WH
				     )
	{
		MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > H = WH ;
		H.recLevel = H.recLevel - 1 ;

		FFLASFFPACK_check(!F.isZero(beta));

		typename Field::Element mbeta  ;
		F.neg(mbeta,beta);

		size_t lb, cb, la, ca;
		size_t x3rd = std::max(mr,kr);
		typename Field::ConstElement_ptr A11=A, A12, A21, A22;
		typename Field::ConstElement_ptr B11=B, B12, B21, B22;
		typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;



		size_t ldX3;
		// Three temporary submatrices are required

		if (ta == FflasTrans) {
			A21 = A + mr;
			A12 = A + kr*lda;
			A22 = A12 + mr;
			la = kr;
			ca = mr;
		}
		else { // ta == FflasNoTrans
			A12 = A + kr;
			A21 = A + mr*lda;
			A22 = A21 + kr;
			la = mr;
			ca = kr;
		}
		if (tb == FflasTrans) {
			B21 = B + kr;
			B12 = B + nr*ldb;
			B22 = B12 + kr;
			lb = nr;
			cb = kr;
			ldX3 = x3rd;
		}
		else { // ta == FflasNoTrans
			B12 = B + nr;
			B21 = B + kr*ldb;
			B22 = B21 + nr;
			lb = kr;
			ldX3 = cb = nr;
		}

		// P2 = alpha . A12 * B21 + beta . C11  in C11
		fgemm (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, beta, C11, ldc, H);

		typename Field::Element_ptr X3 = fflas_new (F, x3rd, nr);

		// T3 = B22 - B12 in X3
		fsub(F,lb,cb,B22,ldb,B12,ldb,X3,ldX3);

		typename Field::Element_ptr X2 = fflas_new (F, mr, kr);

		// S3 = A11 - A21 in X2
		fsub(F,la,ca,A11,lda,A21,lda,X2,ca);

		// C22 = C22 - C12 if beta != 0
		fsubin(F,mr,nr,C12,ldc,C22,ldc);

		// C21 = C21 - C22
		fsubin(F,mr,nr,C22,ldc,C21,ldc);

		// P7 = alpha . S3 * T3 + beta . C22 in C22
		fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, X3, ldX3, beta, C22, ldc, H);

		// T1 = B12 - B11 in X3
		fsub(F,lb,cb,B12,ldb,B11,ldb,X3,ldX3);

		// S1 = A21 + A22 in X2
		fadd(F,la,ca,A21,lda,A22,lda,X2,ca);

		// P5 = alpha . S1*T1 + beta . C12 in C12
		fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, X3, ldX3, beta, C12, ldc, H);

		// T2 = B22 - T1 in X3
		fsub(F,lb,cb,B22,ldb,X3,ldX3,X3,ldX3);

		// S2 = S1 - A11 in X2
		fsubin(F,la,ca,A11,lda,X2,ca);

		typename Field::Element_ptr X1 = fflas_new (F, mr, nr);

		// P6 = alpha . S2 * T2 in X1
		fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, X3, ldX3, F.zero, X1, nr, H);

		// T4 = T2 - B21 in X3
		fsubin(F,lb,cb,B21,ldb,X3,ldX3);

		// S4 = A12 -S2 in X2
		fsub(F,la,ca,A12,lda,X2,ca,X2,ca);

		// P4 = alpha . A22 * T4 - beta . C21 in C21
		fgemm (F, ta, tb, mr, nr, kr, alpha, A22, lda, X3, ldX3, mbeta, C21, ldc, H);

		// P1 = alpha . A11 * B11 in X3
		fgemm (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X3, nr, H);

		//  U1 = P2 + P1 in C11
		faddin(F,mr,nr,X3,nr,C11,ldc);

		// U2 = P1 + P6 in tmpU2/X1  and
		faddin(F, mr, nr, X3, nr, X1, nr);

		// U3 = P7 + U2 in tmpU3/X3  and
		fadd(F, mr, nr, X1, nr, C22, ldc, X3, nr);

		// U7 = P5 + U3 in C22    and
		fadd(F, mr, nr, C12, ldc, X3, nr, C22, ldc);

		// U4 = P5 + U2 in C12    and
		faddin(F, mr, nr, X1, nr, C12, ldc);

		fflas_delete (X1);

		// U6 = U3 - P4 in C21    and
		fsub(F, mr, nr, X3, nr, C21, ldc, C21, ldc);

		fflas_delete (X3);

		// P3 = alpha . S4*B22 in X1
		fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, B22, ldb, F.one, C12, ldc, H);

		fflas_delete (X2);

	} // WinogradAccOld

	// 3 temps and 21 ops
	template < class Field, class ModeTrait>
	inline void WinogradAcc_3_21 (const Field& F,
				      const FFLAS_TRANSPOSE ta,
				      const FFLAS_TRANSPOSE tb,
				      const size_t mr, const size_t nr, const size_t kr,
				      const typename Field::Element alpha,
				      typename Field::ConstElement_ptr A,const size_t lda,
				      typename Field::ConstElement_ptr B,const size_t ldb,
				      const typename Field::Element  beta,
				      typename Field::Element_ptr C, const size_t ldc,
				      MMHelper<Field, MMHelperAlgo::Winograd, ModeTrait > & WH
				     )
	{
		typedef MMHelper<Field, MMHelperAlgo::Winograd, ModeTrait > MMH_t;
		// typedef typename  MMH_t::ModeMgr_t::DelayedField::Element_ptr DFEptr;
		// typedef typename  MMH_t::ModeMgr_t::DelayedField::ConstElement_ptr DFCEptr;
		// typedef typename  MMH_t::ModeMgr_t::DelayedField::Element DFElt;

		// const typename MMH_t::ModeMgr_t::DelayedField & DF = WH.ModeManager.delayedField;

		FFLASFFPACK_check(!F.isZero(beta));

		typename MMH_t::ModeMgr_t& WHMM = WH.ModeManager;

		size_t lb, cb, la, ca;
		size_t x3rd = std::max(mr,kr);
		typename Field::ConstElement_ptr A11=A, A12, A21, A22;
		typename Field::ConstElement_ptr B11=B, B12, B21, B22;
		typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;

		typename Field::Element mbeta;
		F.neg(mbeta,beta);
		// DFElt betadf;
		// if (F.isMOne(beta))
		// 	DF.assign(betadf,DF.mOne);
		// else
		// 	DF.init(betadf, beta);

		size_t ldX3;

		if (ta == FflasTrans) {
			A21 = A + mr;
			A12 = A + kr*lda;
			A22 = A12 + mr;
			la = kr;
			ca = mr;
		} else { // ta == FflasNoTrans
			A12 = A + kr;
			A21 = A + mr*lda;
			A22 = A21 + kr;
			la = mr;
			ca = kr;
		}
		if (tb == FflasTrans) {
			B21 = B + kr;
			B12 = B + nr*ldb;
			B22 = B12 + kr;
			lb = nr;
			cb = kr;
			ldX3 = x3rd;
		} else { // ta == FflasNoTrans
			B12 = B + nr;
			B21 = B + kr*ldb;
			B22 = B21 + nr;
			lb = kr;
			ldX3 = cb = nr;
		}

		// Three temporary submatrices are required

		// T1 = B12 - B11 in X3
		typename Field::Element_ptr X3 = fflas_new (F, x3rd, nr);
		AddSubHelper<Field, typename TryLazy<Field>::value> T1H(F, WHMM.B, WHMM.B);
		fsub (F, lb, cb, B12, ldb, B11, ldb, X3, ldX3, T1H);

		// S1 = A21 + A22 in X2
		typename Field::Element_ptr X2 = fflas_new(F,mr,kr);
		AddSubHelper<Field, typename TryLazy<Field>::value> S1H(F, WHMM.A, WHMM.A);
		fadd (F, la, ca, A21, lda, A22, lda, X2, ca, S1H);

                // P5 = alpha . S1*T1  in X1
		typename Field::Element_ptr X1 = fflas_new(F,mr,nr);
		MMH_t H5(F, WH.AlgoManager.recLevel-1, S1H.ModeManager.Out, T1H.ModeManager.Out);
		fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, X3, ldX3, F.zero, X1, nr, H5);

		// DFElt C22Min, C22Max;
		// DFElt C12Min, C12Max;
		// // This test will be optimized out
		// if (Protected::NeedDoublePreAddReduction (C12Min, C12Max, H5.ModeManager.Outmin, H5.ModeManager.Outmax, WH.ModeManager.Cmin, WH.ModeManager.Cmax, betadf, WH)){
		// 	freduce(F,mr,nr,X1,nr);
		// 	H5.ModeManager.initOut();
		// }
		// C22Min = C12Min; C22Max = C12Max;


		// Local copy of WHMM.C: if C12 or C22 gets reduded, then WHMM.C is reset
		// but C21 or C11 are not. WHMM.C tracks C11 only.
		Operand<Field, ModeTrait> Cx2H(WHMM.C);
		Operand<Field, ModeTrait> C21H(WHMM.C);
                // C22 = P5 + beta C22 in C22
		AddSubHelper<Field, typename TryLazy<Field>::value> P5H(F, H5.ModeManager.Out, Cx2H);
		fadd (F, mr, nr, X1, nr, beta, C22, ldc, C22, ldc, P5H);

		// C12 = P5 + beta C12 in C12
		AddSubHelper<Field, typename TryLazy<Field>::value> P5H2(F, H5.ModeManager.Out, Cx2H);
		std::cerr<<"Avant P5H2 = "<<P5H2<<std::endl;
		fadd  (F, mr, nr, X1, nr, beta, C12, ldc, C12, ldc, P5H2);
		std::cerr<<"Avant P5H2 = "<<P5H2<<std::endl;

		// P1 = alpha . A11 * B11 in X1
		MMH_t H1(F, WH.AlgoManager.recLevel-1, WHMM.A, WHMM.B);
		fgemm (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X1, nr, H1);

		// P2 = alpha . A12 * B21 + beta . C11  in C11
		MMH_t H2(F, WH.AlgoManager.recLevel-1, WHMM.A, WHMM.B, WHMM.C);
		fgemm (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, beta, C11, ldc, H2);

		//  U1 = P2 + P1 in C11
		// DFElt U1Min, U1Max;
		// if (Protected::NeedPreAddReduction (U1Min,U1Max, H1.ModeManager.Outmin, H1.ModeManager.Outmax, H2.ModeManager.Outmin,H2.ModeManager.Outmax, WH) ){
		// 	freduce(F,mr,nr,X1,nr);
		// 	freduce(F,mr,nr,C11,ldc);
		// }
		AddSubHelper<Field,typename TryLazy<Field>::value> U1H (F, H1.ModeManager.Out, H1.ModeManager.Out);
		faddin (F, mr, nr, X1, nr, C11, ldc, U1H);

		// T2 = B22 - T1 in X3
		AddSubHelper<Field, typename TryLazy<Field>::value> T2H(F, WHMM.B, T1H.ModeManager.Out);
		fsub (F, lb, cb, B22, ldb, X3, ldX3, X3, ldX3, T2H);

		// S2 = S1 - A11 in X2
		AddSubHelper<Field, typename TryLazy<Field>::value> S2H(F, S1H.ModeManager.Out, WHMM.A);
		fsubin (F, la, ca, A11, lda, X2, ca, S2H);

		// U2 = P6 + P1 = alpha . S2 * T2 + P1 in X1
		MMH_t H6(F, WH.AlgoManager.recLevel-1, S2H.ModeManager.Out, T2H.ModeManager.Out, H1.ModeManager.Out);
			 // 2*WH.ModeManager.Amin-WH.ModeManager.Amax, 2*WH.ModeManager.Amax-WH.ModeManager.Amin,
			 // 2*WH.ModeManager.Bmin-WH.ModeManager.Bmax, 2*WH.ModeManager.Bmax-WH.ModeManager.Bmin,
			 // H1.ModeManager.Outmin, H1.ModeManager.Outmax);
		fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, X3, ldX3, F.one, X1, nr, H6);

		// U4 = U2 + C12 in C12
		// DFElt U4Min, U4Max;
		// if (Protected::NeedPreAddReduction (U4Min, U4Max, H6.ModeManager.Outmin, H6.ModeManager.Outmax, C12Min, C12Max, WH)){
		// 	freduce(F,mr,nr,C12,ldc);
		// 	freduce(F,mr,nr,X1,nr);
		// }
		AddSubHelper<Field, typename TryLazy<Field>::value> U4H(F, H6.ModeManager.Out, P5H2.ModeManager.Out);
		std::cerr<<"Avant U4H = "<<U4H<<std::endl;
		WriteMatrix(std::cerr<<"U2 = "<<std::endl,F,mr,nr,X1,nr);
		WriteMatrix(std::cerr<<"C12 = "<<std::endl,F,mr,nr,C12,ldc);
		faddin (F, mr, nr, X1, nr, C12, ldc, U4H);
		std::cerr<<"Apres U4H = "<<U4H<<std::endl;
		WriteMatrix(std::cerr<<"C12 = "<<std::endl,F,mr,nr,C12,ldc);

		// T4 = T2 - B21 in X3
		AddSubHelper<Field, typename TryLazy<Field>::value> T4H(F, T2H.ModeManager.Out, WHMM.B);
		fsubin (F, lb, cb, B21, ldb, X3, ldX3, T4H);

		// S4 = A12 -S2 in X2
		AddSubHelper<Field, typename TryLazy<Field>::value> S4H(F, WHMM.A, S2H.ModeManager.Out);
		fsub (F, la, ca, A12, lda, X2, ca, X2, ca, S4H);

		// P4 = alpha . A22 * T4 - beta . C21 in C21
		MMH_t H4(F, WH.AlgoManager.recLevel-1, S4H.ModeManager.Out, T4H.ModeManager.Out, C21H);
		fgemm (F, ta, tb, mr, nr, kr, alpha, A22, lda, X3, ldX3, mbeta, C21, ldc, H4);

		// U5 = P3 + U4 = alpha . S4*B22 + U4 in C12
		MMH_t H3(F, WH.AlgoManager.recLevel-1, S4H.ModeManager.Out, WHMM.B, U4H.ModeManager.Out);
			 // 2*WH.ModeManager.Amin-2*WH.ModeManager.Amax, 2*WH.ModeManager.Amax-2*WH.ModeManager.Amin,
			 // WH.ModeManager.Bmin, WH.ModeManager.Bmax,
			 // U4Min, U4Max);
		std::cerr<<"Avant H3 = "<<H3<<std::endl;
		WriteMatrix(std::cerr<<"S4 = "<<std::endl,F,mr,kr,X2,ca);
		WriteMatrix(std::cerr<<"B22 = "<<std::endl,F,kr,nr,B22,ldb);
		WriteMatrix(std::cerr<<"U4 = "<<std::endl,F,mr,nr,C12,ldc);
		fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, B22, ldb, F.one, C12, ldc, H3);
		std::cerr<<"Apres H3 = "<<H3<<std::endl;
		WriteMatrix(std::cerr<<"U5 = "<<std::endl,F,mr,nr,C12,ldc);

                // T3 = B22 - B12 in X3
		AddSubHelper<Field, typename TryLazy<Field>::value> T3H (F, WHMM.B, WHMM.B);
		fsub (F, lb, cb, B22, ldb, B12, ldb, X3, ldX3, T3H);

		// S3 = A11 - A21 in X2
		AddSubHelper<Field, typename TryLazy<Field>::value> S3H (F, WHMM.A, WHMM.A);
		fsub (F, la, ca, A11, lda, A21, lda, X2, ca, S3H);

		// U3 = P7 + U2  = alpha . S3 * T3 + U2 in X1
		MMH_t H7(F, WH.AlgoManager.recLevel-1, S3H.ModeManager.Out, T3H.ModeManager.Out, H6.ModeManager.Out);
			 // WH.ModeManager.Amin-WH.ModeManager.Amax, WH.ModeManager.Amax-WH.ModeManager.Amin,
			 // WH.ModeManager.Bmin-WH.ModeManager.Bmax, WH.ModeManager.Bmax-WH.ModeManager.Bmin,
			 // H6.ModeManager.Outmin, H6.ModeManager.Outmax);
		fgemm (F, ta, tb, mr, nr, kr, alpha, X2, ca, X3, ldX3, F.one, X1, nr, H7);

		fflas_delete (X2);
		fflas_delete (X3);

		// U7 =  U3 + C22 in C22
		// DFElt U7Min, U7Max;
		// if (Protected::NeedPreAddReduction (U7Min, U7Max, H7.ModeManager.Outmin, H7.ModeManager.Outmax, C22Min, C22Max, WH)){
		// 	freduce(F,mr,nr,X1,nr);
		// 	freduce(F,mr,nr,C22,ldc);
		// }
		AddSubHelper<Field,typename TryLazy<Field>::value> U7H (F, H7.ModeManager.Out, P5H.ModeManager.Out);
		faddin (F, mr, nr, X1, nr, C22, ldc, U7H);

		// U6 = U3 - P4 in C21
		// DFElt U6Min, U6Max;
		// if (Protected::NeedPreSubReduction(U6Min, U6Max, H7.ModeManager.Outmin, H7.ModeManager.Outmax, H4.ModeManager.Outmin, H4.ModeManager.Outmax, WH)){
		// 	freduce(F,mr,nr,X1,nr);
		// 	freduce(F,mr,nr,C21,ldc);
		// }
		AddSubHelper<Field,typename TryLazy<Field>::value> U6H (F, H7.ModeManager.Out, H4.ModeManager.Out);
		fsub (F,mr,nr,X1,nr,C21,ldc,C21,ldc, U6H);

		fflas_delete (X1);

		// Updating WH with Outmin, Outmax of the result
		WH.ModeManager.setOutBoundsMM(U1H, H3, U6H, U7H);
	} // WinogradAcc


	// 2 temps and 24 ops
        // TODO: Add check for modular reductions before final additions
	template < class Field, class FieldTrait >
	inline void WinogradAcc_2_24 (const Field& F,
				      const FFLAS_TRANSPOSE ta,
				      const FFLAS_TRANSPOSE tb,
				      const size_t mr, const size_t nr, const size_t kr,
				      const typename Field::Element alpha,
				      const typename Field::Element_ptr A,const size_t lda,
				      const typename Field::Element_ptr B,const size_t ldb,
				      const typename Field::Element  beta,
				      typename Field::Element_ptr C, const size_t ldc,
				      MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > & WH
				     )
	{
		MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > H = WH ;
		H.recLevel = H.recLevel - 1 ;

		FFLASFFPACK_check(!F.isZero(beta));

		typename Field::Element malpha ;
		F.neg(malpha,alpha);

		// A, B and c submatrices
		const typename Field::Element_ptr A11=A, A12, A21, A22;
		const typename Field::Element_ptr B11=B, B12, B21, B22;
		typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;



		size_t la, ca, lb, cb; // lines and columns in A,B sub matrices

		// Three temporary submatrices are required

		if (ta == FflasTrans) {
			A21 = A + mr;
			A12 = A + kr*lda;
			A22 = A12 + mr;
			la = kr ;
			ca = mr ;
		}
		else { // ta == FflasNoTrans
			A12 = A + kr;
			A21 = A + mr*lda;
			A22 = A21 + kr;
			la = mr ;
			ca = kr ;

		}
		if (tb == FflasTrans) {
			B21 = B + kr;
			B12 = B + nr*ldb;
			B22 = B12 + kr;
			lb = nr ;
			cb = kr ;

		}
		else { // ta == FflasNoTrans
			B12 = B + nr;
			B21 = B + kr*ldb;
			B22 = B21 + nr;
			lb = kr ;
			cb = nr ;
		}

		// Z1 = C22 - C12         in C22
		fsubin(F,mr,nr,C12,ldc,C22,ldc);
		// Z3 = C12-C21           in C12
		fsubin(F,mr,nr,C21,ldc,C12,ldc);
		// S1 = A21 + A22         in X
		typename Field::Element_ptr X = fflas_new(F,mr,std::max(nr,kr));
		fadd(F,la,ca,A21,lda,A22,lda,X,ca);
		// T1 = B12 - B11         in Y
		typename Field::Element_ptr Y = fflas_new(F,nr,kr);
		fsub(F,lb,cb,B12,ldb,B11,ldb,Y,cb);
		// P5 = a S1 T1 + b Z3    in C12
		fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, Y, cb, beta, C12, ldc, H);
		// S2 = S1 - A11          in X
		fsubin(F,la,ca,A11,lda,X,ca);
		// T2 = B22 - T1          in Y
		fsub(F,lb,cb,B22,ldb,Y,cb,Y,cb);
		// P6 = a S2 T2 + b C21   in C21
		fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, Y, cb, beta, C21, ldc, H);
		// S4 = A12 - S2          in X
		fsub(F,la,ca,A12,lda,X,ca,X,ca);
		// W1 = P5 + beta Z1      in C22
		fadd(F,mr,nr,C12,ldc,beta,C22,ldc,C22,ldc);
		// P3 = a S4 B22 + P5     in C12
		fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, B22, ldb, F.one, C12, ldc, H);
		// P1 = a A11 B11         in X
		fgemm (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X, nr, H);
		// U2 = P6 + P1           in C21
		faddin(F,mr,nr,X,nr,C21,ldc);
		// P2 = a A12 B21 + b C11 in C11
		fgemm (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, beta, C11, ldc, H);
		// U1 = P1 + P2           in C11
		faddin(F,mr,nr,X,nr,C11,ldc);
		// U5 = U2 + P3           in C12
		faddin(F,mr,nr,C21,ldc,C12,ldc);
		// S3 =  A11 - A21        in X ;
		fsub(F,la,ca,A11,lda,A21,lda,X,ca);
		// T3 = B22 - B12         in Y
		fsub(F,lb,cb,B22,ldb,B12,ldb,Y,cb);
		// U3 = a S3 T3 + U2      in C21
		fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, Y, cb, F.one, C21, ldc, H);
		fflas_delete (X);
		// U7 = U3 + W1           in C22
		faddin(F,mr,nr,C21,ldc,C22,ldc);
		// T1_ = B12 - B11        in Y
		fsub(F,lb,cb,B12,ldb,B11,ldb,Y,cb);
		// T2_ = B22 - T1_        in Y
		fsub(F,lb,cb,B22,ldb,Y,cb,Y,cb);
		// T4 = T2_ - B21         in Y
		fsub(F,lb,cb,Y,cb,B21,ldb,Y,cb);
		// U6 = -a A22 T4 + U3    in C21;
		fgemm (F, ta, tb, mr, nr, kr, malpha, A22, lda, Y, cb, F.one, C21, ldc, H);
		fflas_delete (Y);


	} // WinogradAccOld

	// 2 temps and 27 ops
        // TODO: Add check for modular reductions before final additions
	template < class Field, class FieldTrait >
	inline void WinogradAcc_2_27 (const Field& F,
				      const FFLAS_TRANSPOSE ta,
				      const FFLAS_TRANSPOSE tb,
				      const size_t mr, const size_t nr, const size_t kr,
				      const typename Field::Element alpha,
				      const typename Field::Element_ptr A,const size_t lda,
				      const typename Field::Element_ptr B,const size_t ldb,
				      const typename Field::Element  beta,
				      typename Field::Element_ptr C, const size_t ldc,
				      MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > & WH)
	{
		MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > H = WH ;
		H.recLevel = H.recLevel - 1 ;

		FFLASFFPACK_check(!F.isZero(beta));

		typename Field::Element malpha ;
		F.neg(malpha,alpha);

		// A, B and c submatrices
		const typename Field::Element_ptr A11=A, A12, A21, A22;
		const typename Field::Element_ptr B11=B, B12, B21, B22;
		typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;



		size_t la, ca, lb, cb; // lines and columns in A,B sub matrices

		// Three temporary submatrices are required

		if (ta == FflasTrans) {
			A21 = A + mr;
			A12 = A + kr*lda;
			A22 = A12 + mr;
			la = kr ;
			ca = mr ;
		}
		else { // ta == FflasNoTrans
			A12 = A + kr;
			A21 = A + mr*lda;
			A22 = A21 + kr;
			la = mr ;
			ca = kr ;

		}
		if (tb == FflasTrans) {
			B21 = B + kr;
			B12 = B + nr*ldb;
			B22 = B12 + kr;
			lb = nr ;
			cb = kr ;

		}
		else { // ta == FflasNoTrans
			B12 = B + nr;
			B21 = B + kr*ldb;
			B22 = B21 + nr;
			lb = kr ;
			cb = nr ;
		}

		// Z1 = C22 - C12         in C22
		fsubin(F,mr,nr,C12,ldc,C22,ldc);
		// Z3 = C12-C21           in C12
		fsubin(F,mr,nr,C21,ldc,C12,ldc);
		// S1 = A21 + A22         in X
		typename Field::Element_ptr X = fflas_new(F,mr,std::max(nr,kr));
		fadd(F,la,ca,A21,lda,A22,lda,X,ca);
		// T1 = B12 - B11         in Y
		typename Field::Element_ptr Y = fflas_new(F,nr,std::max(kr,mr));
		fsub(F,lb,cb,B12,ldb,B11,ldb,Y,cb);
		// P5 = a S1 T1 + b Z3    in C12
		fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, Y, cb, beta, C12, ldc, H);
		// S2 = S1 - A11          in X
		fsubin(F,la,ca,A11,lda,X,ca);
		// T2 = B22 - T1          in Y
		fsub(F,lb,cb,B22,ldb,Y,cb,Y,cb);
		// P6 = a S2 T2 + b C21   in C21
		fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, Y, cb, beta, C21, ldc, H);
		// S4 = A12 - S2          in X
		fsub(F,la,ca,A12,lda,X,ca,X,ca);
		// W1 = P5 + beta Z1      in C22
		fadd(F,mr,nr,C12,ldc,beta,C22,ldc,C22,ldc);
		// P3 = a S4 B22 + P5     in C12
		fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, B22, ldb, F.zero, Y, nr, H);
		fadd(F,mr,nr,Y,nr,C12,ldc,C12,ldc);
		// P1 = a A11 B11         in X
		fgemm (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X, nr, H);
		// U2 = P6 + P1           in C21
		faddin(F,mr,nr,X,nr,C21,ldc);
		// P2 = a A12 B21 + b C11 in C11
		fgemm (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.zero, Y, nr, H);
		fadd(F,mr,nr,Y,nr,beta,C11,ldc,C11,ldc);
		// U1 = P1 + P2           in C11
		faddin(F,mr,nr,X,nr,C11,ldc);
		// U5 = U2 + P3           in C12
		faddin(F,mr,nr,C21,ldc,C12,ldc);
		// S3 =  A11 - A21        in X ;
		fsub(F,la,ca,A11,lda,A21,lda,X,ca);
		// T3 = B22 - B12         in Y
		fsub(F,lb,cb,B22,ldb,B12,ldb,Y,cb);
		// U3 = a S3 T3 + U2      in C21
		fgemm (F, ta, tb, mr, nr, kr, alpha, X, ca, Y, cb, F.one, C21, ldc, H);
		// U7 = U3 + W1           in C22
		faddin(F,mr,nr,C21,ldc,C22,ldc);
		// T1_ = B12 - B11        in Y
		fsub(F,lb,cb,B12,ldb,B11,ldb,Y,cb);
		// T2_ = B22 - T1_        in Y
		fsub(F,lb,cb,B22,ldb,Y,cb,Y,cb);
		// T4 = T2_ - B21         in Y
		fsub(F,lb,cb,Y,cb,B21,ldb,Y,cb);
		// U6 = -a A22 T4 + U3    in C21;
		fgemm (F, ta, tb, mr, nr, kr, alpha, A22, lda, Y, cb, F.zero, X, nr, H);
		fflas_delete (Y);
		fsub(F,mr,nr,C21,ldc,X,nr,C21,ldc);
		fflas_delete (X);


	} // WinogradAcc3


} // BLAS3

} // FFLAS

#endif // __FFLASFFPACK_fgemm_winograd_acc_INL

