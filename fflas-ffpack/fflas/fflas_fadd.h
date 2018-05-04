/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2014 FFLAS-FFPACK group
 *
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

#ifndef __FFLASFFPACK_fadd_H
#define __FFLASFFPACK_fadd_H

namespace FFLAS {

	template<class T>
	struct support_simd_add  : public std::false_type {} ;

// #ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
	template<>
	struct support_simd_add<float> : public std::true_type {} ;
	template<>
	struct support_simd_add<double> : public std::true_type {} ;
 #ifdef SIMD_INT
	template<>
	struct support_simd_add<int64_t> : public std::true_type {} ;
	template<>
	struct support_simd_add<int32_t> : public std::true_type {} ;

 #endif  // SIMD_INT

// #endif // __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

} // FFLAS

#include "fflas_fadd.inl"

namespace FFLAS {

	/***************************/
	/*         LEVEL 1         */
	/***************************/

	template <class Field>
	void
	fadd (const Field & F,  const size_t N,
	      typename Field::ConstElement_ptr A, const size_t inca,
	      typename Field::ConstElement_ptr B, const size_t incb,
	      typename Field::Element_ptr C, const size_t incc)
	{
		details::fadd<Field, true>(F,N,A,inca,B,incb,C,incc
					   , typename FieldTraits<Field>::category() );
	}



	template <class Field>
	void
	faddin (const Field& F,  const size_t N,
		typename Field::ConstElement_ptr B, const size_t incb,
		typename Field::Element_ptr C, const size_t incc)
	{
		fadd(F,N,B,incb,C,incc,C,incc);
		return;
	}

	template <class Field>
	void
	fsub(const Field & F,  const size_t N,
	      typename Field::ConstElement_ptr A, const size_t inca,
	      typename Field::ConstElement_ptr B, const size_t incb,
	      typename Field::Element_ptr C, const size_t incc)
	{
		details::fadd<Field, false>(F,N,A,inca,B,incb,C,incc, typename FieldTraits<Field>::category() );
	}



	template <class Field>
	void
	fsubin (const Field& F,  const size_t N,
		typename Field::ConstElement_ptr B, const size_t incb,
		typename Field::Element_ptr C, const size_t incc)
	{
		fsub(F,N,C,incc,B,incb,C,incc);
		return;
	}

	// C = A + a B
	template <class Field>
	void
	fadd (const Field& F, const size_t N,
	      typename Field::ConstElement_ptr A, const size_t inca,
	      const typename Field::Element alpha,
	      typename Field::ConstElement_ptr B, const size_t incb,
	      typename Field::Element_ptr C, const size_t incc)
	{
		if (C == A && inca == incc)
			return faxpy(F,N,alpha,B,incb,C,incc);
		if (F.isOne(alpha))
			return fadd(F,N,A,inca,B,incb,C,incc);
		if (F.isMOne(alpha)){
			return fsub(F,N,A,inca,B,incb,C,incc);
		}
		if (F.isZero(alpha))
			return fassign(F,N,A,inca,C,incc);

		if (inca == 1 && incb == 1 && incc == 1) {
			for (size_t i = 0 ; i < N ; ++i) {
				//!@todo optimise here
				F.mul(C[i],alpha,B[i]);
				F.addin(C[i],A[i]);
			}
			return;
		}

		typename Field::ConstElement_ptr Ai = A, Bi = B;
		typename Field::Element_ptr Ci = C;
		for (; Ai < A+N*inca; Ai+=inca, Bi+=incb, Ci+=incc) {
			F.mul(*Ci,alpha,*Bi);
			F.addin (*Ci, *Ai);
		}
	}


	/***************************/
	/*         LEVEL 2         */
	/***************************/


	template <class Field>
        void
        pfadd (const Field & F,  const size_t M, const size_t N,
               typename Field::ConstElement_ptr A, const size_t lda,
               typename Field::ConstElement_ptr B, const size_t ldb,
               typename Field::Element_ptr C, const size_t ldc, const size_t numths){
            SYNCH_GROUP(
              FORBLOCK1D(iter, M, SPLITTER(numths),
			    size_t rowsize= iter.end()-iter.begin();
                TASK(MODE(CONSTREFERENCE(F) READWRITE(C[iter.begin()*ldc]) READ(A[iter.begin()*lda], B[iter.begin()*ldb])),
                fadd(F, rowsize, N, A+iter.begin()*lda, lda, B+iter.begin()*ldb, ldb, C+iter.begin()*ldc, ldc);
                     );
                         );
                        );
        }

        template <class Field>
        void
        pfsub (const Field & F,  const size_t M, const size_t N,
               typename Field::ConstElement_ptr A, const size_t lda,
               typename Field::ConstElement_ptr B, const size_t ldb,
               typename Field::Element_ptr C, const size_t ldc, const size_t numths){
            SYNCH_GROUP(
              FORBLOCK1D(iter, M, SPLITTER(numths),
                size_t rowsize= iter.end()-iter.begin();
                TASK(MODE(CONSTREFERENCE(F) READWRITE(C[iter.begin()*ldc]) READ(A[iter.begin()*lda], B[iter.begin()*ldb])),
                fsub(F, rowsize, N, A+iter.begin()*lda, lda, B+iter.begin()*ldb, ldb, C+iter.begin()*ldc, ldc);
                     );
                         );
                        );
        }


        template <class Field>
        void
        pfaddin (const Field& F, const size_t M, const size_t N,
                typename Field::ConstElement_ptr B, const size_t ldb,
                 typename Field::Element_ptr C, const size_t ldc, size_t numths){

            SYNCH_GROUP(
              FORBLOCK1D(iter, M, SPLITTER(numths),
                size_t rowsize= iter.end()-iter.begin();
                TASK(MODE(CONSTREFERENCE(F) READWRITE(C[iter.begin()*ldc]) READ(B[iter.begin()*ldb])),
                     faddin(F, rowsize, N, B+iter.begin()*ldb, ldb, C+iter.begin()*ldc, ldc);
                     );
                         );
                        );
        }

        template <class Field>
        void
        pfsubin (const Field& F, const size_t M, const size_t N,
                typename Field::ConstElement_ptr B, const size_t ldb,
                 typename Field::Element_ptr C, const size_t ldc, size_t numths){
            SYNCH_GROUP(
              FORBLOCK1D(iter, M, SPLITTER(numths),
              size_t rowsize= iter.end()-iter.begin();
              TASK(MODE(CONSTREFERENCE(F) READWRITE(C[iter.begin()*ldc]) READ(B[iter.begin()*ldb])),
              fsubin(F, rowsize, N, B+iter.begin()*ldb, ldb, C+iter.begin()*ldc, ldc);
                   );
                         );
                        );
        }

	template <class Field>
	void fadd (const Field& F, const size_t M, const size_t N,
		   typename Field::ConstElement_ptr A, const size_t lda,
		   typename Field::ConstElement_ptr B, const size_t ldb,
		   typename Field::Element_ptr C, const size_t ldc) {
		AddSubHelper<Field> ASH(F);
		fadd (F, M, N, A, lda, B, ldb, C, ldc, ASH);
	}
	template <class Field,class ModeT>
	void fadd (const Field& F, const size_t M, const size_t N,
		   typename Field::ConstElement_ptr A, const size_t lda,
		   typename Field::ConstElement_ptr B, const size_t ldb,
		   typename Field::Element_ptr C, const size_t ldc,
		   AddSubHelper<Field,ModeT,ParSeqHelper::Sequential>& H) {
		H.ModeManager.setOutBoundsAdd();
		std::cerr<<"in general fadd H="<<H<<std::endl;
		if (N == lda && N == ldb && N == ldc)
			return fadd(F,M*N,A,1,B,1,C,1);
		typename Field::ConstElement_ptr Ai = A, Bi = B;
		typename Field::Element_ptr Ci = C;
		for (; Ai < A+M*lda; Ai+=lda, Bi+=ldb, Ci+=ldc)
			fadd(F,N,Ai,1,Bi,1,Ci,1);

	}
	
	template <class Field>
	void fadd (const Field& F, const size_t M, const size_t N,
		   typename Field::ConstElement_ptr A, const size_t lda,
		   typename Field::ConstElement_ptr B, const size_t ldb,
		   typename Field::Element_ptr C, const size_t ldc,
		   AddSubHelper<Field,ModeCategories::LazyTag,ParSeqHelper::Sequential>& H) {

		std::cerr<<"fadd lazy"<<std::endl;
		typename AddSubHelper<Field,ModeCategories::LazyTag,ParSeqHelper::Sequential>::ModeMgr_t& MM = H.ModeManager;
		if (MM.MaxStorableValue - MM.A.max < -MM.B.max ||
		    MM.MaxStorableValue - MM.A.min < -MM.B.min){
			MM.initA();
			MM.initB();
			    // TODO merge these freduce with the fsub (to reduce cache misses)
			freduce (F, M, N, const_cast<typename Field::Element_ptr>(A), lda);
			freduce (F, M, N, const_cast<typename Field::Element_ptr>(B), lda);
		}
		AddSubHelper<typename associatedDelayedField<const Field>::field, ModeCategories::DefaultBoundedTag> Hdf(H);
		fadd (MM.delayedField, M, N, A, lda, B, ldb, C, ldc, Hdf);
		std::cerr<<"MM = "<<MM<<std::endl;
		MM.Out = Hdf.ModeManager.Out;
	}

	template <class Field>
	void fsub (const Field& F, const size_t M, const size_t N,
		   typename Field::ConstElement_ptr A, const size_t lda,
		   typename Field::ConstElement_ptr B, const size_t ldb,
		   typename Field::Element_ptr C, const size_t ldc) {
		AddSubHelper<Field> ASH(F);
		fsub(F, M, N, A, lda, B, ldb, C, ldc, ASH);
	}

	template <class Field, class ModeT>
	void fsub (const Field& F, const size_t M, const size_t N,
		   typename Field::ConstElement_ptr A, const size_t lda,
		   typename Field::ConstElement_ptr B, const size_t ldb,
		   typename Field::Element_ptr C, const size_t ldc,
		   AddSubHelper<Field,ModeT,ParSeqHelper::Sequential>& H) {
		H.ModeManager.setOutBoundsSub();
		if (N == lda && N == ldb && N == ldc)
			return fsub(F,M*N,A,1,B,1,C,1);
		typename Field::ConstElement_ptr Ai = A, Bi = B;
		typename Field::Element_ptr Ci = C;
		for (; Ai < A+M*lda; Ai+=lda, Bi+=ldb, Ci+=ldc)
			fsub(F,N,Ai,1,Bi,1,Ci,1);
	}

	template <class Field>
	void fsub (const Field& F, const size_t M, const size_t N,
		   typename Field::ConstElement_ptr A, const size_t lda,
		   typename Field::ConstElement_ptr B, const size_t ldb,
		   typename Field::Element_ptr C, const size_t ldc,
		   AddSubHelper<Field,ModeCategories::LazyTag,ParSeqHelper::Sequential>& H) {

		typename AddSubHelper<Field,ModeCategories::LazyTag,ParSeqHelper::Sequential>::ModeMgr_t& MM = H.ModeManager;
		if (MM.MaxStorableValue - MM.A.max < -MM.B.min ||
		    MM.MaxStorableValue - MM.B.max < -MM.A.min){
			MM.initA();
			MM.initB();
			    // TODO merge these freduce with the fsub (to reduce cache misses)
			freduce (F, M, N, const_cast<typename Field::Element_ptr>(A), lda);
			freduce (F, M, N, const_cast<typename Field::Element_ptr>(B), lda);
		}
		AddSubHelper<typename associatedDelayedField<const Field>::field,
			     ModeCategories::DefaultBoundedTag> Hdf(H);
		fsub (MM.delayedField, M, N, A, lda, B, ldb, C, ldc, Hdf);
		MM.Out = Hdf.ModeManager.Out;
	}

	template <class Field>
	void faddin (const Field& F, const size_t M, const size_t N,
		     typename Field::ConstElement_ptr B, const size_t ldb,
		     typename Field::Element_ptr C, const size_t ldc) {
		AddSubHelper<Field> ASH(F);
		faddin(F, M, N, B, ldb, C, ldc, ASH);
	}

	template <class Field, class ModeT>
	void faddin (const Field& F, const size_t M, const size_t N,
		typename Field::ConstElement_ptr B, const size_t ldb,
		typename Field::Element_ptr C, const size_t ldc,
		AddSubHelper<Field,ModeT,ParSeqHelper::Sequential>& H){
		H.ModeManager.setOutBoundsAdd();
		if (N == ldb && N == ldc)
			return faddin(F,M*N,B,1,C,1);
		typename Field::ConstElement_ptr Bi = B;
		typename Field::Element_ptr Ci = C;
		for (; Bi < B+M*ldb;  Bi+=ldb, Ci+=ldc)
			faddin(F,N,Bi,1,Ci,1);

	}
	template <class Field>
	void faddin (const Field& F, const size_t M, const size_t N,
		     typename Field::ConstElement_ptr B, const size_t ldb,
		     typename Field::Element_ptr C, const size_t ldc,
		     AddSubHelper<Field,ModeCategories::LazyTag,ParSeqHelper::Sequential>& H) {
	
		typename AddSubHelper<Field,ModeCategories::LazyTag,ParSeqHelper::Sequential>::ModeMgr_t& MM = H.ModeManager;
		if (MM.MaxStorableValue - MM.A.max < -MM.B.max ||
		    MM.MaxStorableValue - MM.A.min < -MM.B.min){
			MM.initA();
			MM.initB();
			    // TODO merge these freduce with the fsub (to reduce cache misses)
			freduce (F, M, N, const_cast<typename Field::Element_ptr>(B), ldb);
			freduce (F, M, N, C, ldc);
		}
		AddSubHelper<typename associatedDelayedField<const Field>::field, ModeCategories::DefaultBoundedTag> Hdf(H);
		faddin (MM.delayedField, M, N, B, ldb, C, ldc, Hdf);
		MM.Out = Hdf.ModeManager.Out;
	}


	template <class Field>
	void fsubin (const Field& F, const size_t M, const size_t N,
		     typename Field::ConstElement_ptr B, const size_t ldb,
		     typename Field::Element_ptr C, const size_t ldc) {
		AddSubHelper<Field> ASH(F);
		fsubin(F, M, N, B, ldb, C, ldc, ASH);
	}
	template <class Field, class ModeT>
	void fsubin (const Field& F, const size_t M, const size_t N,
		     typename Field::ConstElement_ptr B, const size_t ldb,
		     typename Field::Element_ptr C, const size_t ldc,
		     AddSubHelper<Field,ModeT,ParSeqHelper::Sequential>& H){
		H.ModeManager.setOutBoundsSub();
		if (N == ldb && N == ldc)
			return fsubin(F,M*N,B,1,C,1);
		typename Field::ConstElement_ptr Bi = B;
		typename Field::Element_ptr Ci = C;
		for (; Bi < B+M*ldb;  Bi+=ldb, Ci+=ldc)
			fsubin(F,N,Bi,1,Ci,1);

	}

	template <class Field>
	void fsubin (const Field& F, const size_t M, const size_t N,
		     typename Field::ConstElement_ptr B, const size_t ldb,
		     typename Field::Element_ptr C, const size_t ldc,
		     AddSubHelper<Field,ModeCategories::LazyTag,ParSeqHelper::Sequential>& H) {

		typename AddSubHelper<Field,ModeCategories::LazyTag,ParSeqHelper::Sequential>::ModeMgr_t& MM = H.ModeManager;
		if (MM.MaxStorableValue - MM.A.max < -MM.B.min ||
		    MM.MaxStorableValue - MM.B.max < -MM.A.min){
			MM.initA();
			MM.initB();
			    // TODO merge these freduce with the fsub (to reduce cache misses)
			freduce (F, M, N, const_cast<typename Field::Element_ptr>(B), ldb);
			freduce (F, M, N, C, ldc);
		}
		AddSubHelper<typename associatedDelayedField<const Field>::field, ModeCategories::DefaultBoundedTag> Hdf(H);
		fsubin (MM.delayedField, M, N, B, ldb, C, ldc, Hdf);
		MM.Out = Hdf.ModeManager.Out;
	}
	// C = A + a B
	template <class Field>
	void fadd (const Field& F, const size_t M, const size_t N,
		   typename Field::ConstElement_ptr A, const size_t lda,
		   const typename Field::Element alpha,
		   typename Field::ConstElement_ptr B, const size_t ldb,
		   typename Field::Element_ptr C, const size_t ldc){
		AddSubHelper<Field> ASH(F);
		return fadd(F, M, N, A, lda, alpha, B, ldb, C, ldc, ASH);
	}

	template <class Field, class ModeT>
	void fadd (const Field& F, const size_t M, const size_t N,
		   typename Field::ConstElement_ptr A, const size_t lda,
		   const typename Field::Element alpha,
		   typename Field::ConstElement_ptr B, const size_t ldb,
		   typename Field::Element_ptr C, const size_t ldc,
		   AddSubHelper<Field,ModeT,ParSeqHelper::Sequential>& H){
		if (C == A && lda == ldc)
			return faxpy(F,M,N,alpha,B,ldb,C,ldc);
		if (F.isOne(alpha))
			return fadd(F,M,N,A,lda,B,ldb,C,ldc);
		if (F.isMOne(alpha))
			return fsub(F,M,N,A,lda,B,ldb,C,ldc);
		if (F.isZero(alpha))
			return fassign(F,M,N,A,lda,C,ldc);

		if (N == lda && N == ldb && N == ldc)
			return fadd(F,M*N,A,1,alpha,B,1,C,1);

		typename Field::ConstElement_ptr Ai = A, Bi = B;
		typename Field::Element_ptr Ci = C;
		for (; Ai < A+M*lda; Ai+=lda, Bi+=ldb, Ci+=ldc)
			fadd (F, N, Ai, 1, alpha, Bi, 1, Ci, 1);
		    //TODO update Helper.ModeManager.setOutBounds
	}

	template <class Field>
	void fadd (const Field& F, const size_t M, const size_t N,
		   typename Field::ConstElement_ptr A, const size_t lda,
		   const typename Field::Element alpha,
		   typename Field::ConstElement_ptr B, const size_t ldb,
		   typename Field::Element_ptr C, const size_t ldc,
		   AddSubHelper<Field,ModeCategories::LazyTag,ParSeqHelper::Sequential>& H){
		AddSubHelper<typename associatedDelayedField<const Field>::field, ModeCategories::DefaultBoundedTag> Hdf(H);
		fadd (H.ModeManager.delayedField, M, N, A, lda, alpha, B, ldb, C, ldc, Hdf);
		H.ModeManager.Out = Hdf.ModeManager.Out;
	}
} // FFLAS


#endif // __FFLASFFPACK_fadd_H

