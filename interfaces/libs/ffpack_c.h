/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2015 FFLAS-FFPACK
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

/** @file ffpack-c.h
 * @author  Brice Boyer
 * @brief C functions calls for FFPACK
 * @see ffpack/ffpack.h
 */

#ifndef __FFLASFFPACK_interfaces_libs_ffpack_c_H
#define __FFLASFFPACK_interfaces_libs_ffpack_c_H

#include <stdbool.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif



#ifndef __FFLASFFPACK_interfaces_libs_fflas_c_H

enum FFLAS_order {
	FflasRowMajor=101,
	FflasColMajor=102
};
enum FFLAS_transpose {
	FflasNoTrans = 111,
	FflasTrans   = 112
};
enum FFLAS_uplo {
	FflasUpper = 121,
	FflasLower = 122
};
enum FFLAS_diag {
	FflasNonUnit = 131,
	FflasUnit    = 132
};
enum FFLAS_side {
	FflasLeft  = 141,
	FflasRight = 142
};
enum FFLAS_base {
	FflasDouble  = 151,
	FflasFloat   = 152,
	FflasGeneric = 153
};

#endif // __FFLASFFPACK_interfaces_libs_fflas_c_H

enum FFPACK_lu_tag
{
	FfpackSlabRecursive = 1,
	FfpackTileRecursive = 2,
	FfpackSingular = 3
};
enum FFPACK_charpoly_tag
{
	FfpackLUK=1,
	FfpackKG=2,
	FfpackHybrid=3,
	FfpackKGFast=4,
	FfpackDanilevski=5,
	FfpackArithProg=6,
	FfpackKGFastG=7
};
enum FFPACK_minpoly_tag
{
	FfpackDense=1,
	FfpackKGF=2
};

/*****************/
/* PERMUTATIONS  */
/*****************/


void LAPACKPerm2MathPerm (size_t * MathP, const size_t * LapackP,
			  const size_t N);

void MathPerm2LAPACKPerm (size_t * LapackP, const size_t * MathP,
			  const size_t N);

void MatrixApplyS_modular_double (const double p, double * A, const size_t lda, const size_t width,
				  const size_t M2,
				  const size_t R1, const size_t R2,
				  const size_t R3, const size_t R4
				  , bool positive );

void PermApplyS_double (double * A, const size_t lda, const size_t width,
			const size_t M2,
			const size_t R1, const size_t R2,
			const size_t R3, const size_t R4);


void MatrixApplyT_modular_double (const double p, double * A, const size_t lda, const size_t width,
				  const size_t N2,
				  const size_t R1, const size_t R2,
				  const size_t R3, const size_t R4
				  , bool positive );

void PermApplyT_double (double * A, const size_t lda, const size_t width,
			const size_t N2,
			const size_t R1, const size_t R2,
			const size_t R3, const size_t R4);

void composePermutationsP (size_t * MathP,
			   const size_t * P1,
			   const size_t * P2,
			   const size_t R, const size_t N);

void composePermutationsQ (size_t * MathP,
			   const size_t * Q1,
			   const size_t * Q2,
			   const size_t R, const size_t N);

void cyclic_shift_mathPerm (size_t * P,  const size_t s);

#if 0
template<typename Base_t>
void cyclic_shift_row_col(Base_t * A, size_t m, size_t n, size_t lda);
#endif


void cyclic_shift_row_modular_double(const double p, double * A, size_t m, size_t n, size_t lda
				     , bool positive );


void cyclic_shift_col_modular_double(const double p, double * A, size_t m, size_t n, size_t lda
				     , bool positive );



void
applyP_modular_double( const double p,
		       const enum FFLAS_side Side,
		       const enum FFLAS_transpose Trans,
		       const size_t M, const size_t ibeg, const size_t iend,
		       double * A, const size_t lda, const size_t * P
		       , bool positive  );





/* fgetrs, fgesv */

void
fgetrsin_modular_double (const double p,
		       const enum FFLAS_side Side,
		       const size_t M, const size_t N, const size_t R,
		       double * A, const size_t lda,
		       const size_t *P, const size_t *Q,
		       double * B, const size_t ldb,
		       int * info
		       , bool positive );


double *
fgetrs_modular_double (const double p,
		       const enum FFLAS_side Side,
		       const size_t M, const size_t N, const size_t NRHS, const size_t R,
		       double * A, const size_t lda,
		       const size_t *P, const size_t *Q,
		       double * X, const size_t ldx,
		       const double * B, const size_t ldb,
		       int * info
		       , bool positive );


size_t
fgesvin_modular_double (const double p,
		      const enum FFLAS_side Side,
		      const size_t M, const size_t N,
		      double * A, const size_t lda,
		      double * B, const size_t ldb,
		      int * info
		      , bool positive );


size_t
fgesv_modular_double (const double p,
		      const enum FFLAS_side Side,
		      const size_t M, const size_t N, const size_t NRHS,
		      double * A, const size_t lda,
		      double * X, const size_t ldx,
		      const double * B, const size_t ldb,
		      int * info);

/* ftrtr */


void
ftrtri_modular_double (const double p, const enum FFLAS_uplo Uplo, const enum FFLAS_diag Diag,
		       const size_t N, double * A, const size_t lda
		       , bool positive );


void trinv_left_modular_double( const double p, const size_t N, const double * L, const size_t ldl,
				double * X, const size_t ldx
				, bool positive  );


void
ftrtrm_modular_double (const double p, const enum FFLAS_diag diag, const size_t N,
		       double * A, const size_t lda
		       , bool positive );



/* PLUQ */


size_t
PLUQ_modular_double (const double p, const enum FFLAS_diag Diag,
		     const size_t M, const size_t N,
		     double * A, const size_t lda,
		     size_t*P, size_t *Q
		     , bool positive );



size_t
LUdivine_modular_double (const double p, const enum FFLAS_diag Diag,  const enum FFLAS_transpose trans,
			 const size_t M, const size_t N,
			 double * A, const size_t lda,
			 size_t* P, size_t* Qt,
			 const enum FFPACK_lu_tag LuTag,
			 const size_t cutoff
			 , bool positive );


size_t
LUdivine_small_modular_double (const double p, const enum FFLAS_diag Diag,  const enum FFLAS_transpose trans,
			       const size_t M, const size_t N,
			       double * A, const size_t lda,
			       size_t* P, size_t* Q,
			       const enum FFPACK_lu_tag LuTag
			       , bool positive );


size_t
LUdivine_gauss_modular_double (const double p, const enum FFLAS_diag Diag,
			       const size_t M, const size_t N,
			       double * A, const size_t lda,
			       size_t* P, size_t* Q,
			       const enum FFPACK_lu_tag LuTag
			       , bool positive );



/*****************/
/* ECHELON FORMS */
/*****************/


size_t
ColumnEchelonForm_modular_double (const double p, const size_t M, const size_t N,
				  double * A, const size_t lda,
				  size_t* P, size_t* Qt, bool transform,
				  const enum FFPACK_lu_tag LuTag
				  , bool positive );


size_t
RowEchelonForm_modular_double (const double p, const size_t M, const size_t N,
			       double * A, const size_t lda,
			       size_t* P, size_t* Qt, const bool transform,
			       const enum FFPACK_lu_tag LuTag
			       , bool positive );


size_t
ReducedColumnEchelonForm_modular_double (const double p, const size_t M, const size_t N,
					 double * A, const size_t lda,
					 size_t* P, size_t* Qt, const bool transform,
					 const enum FFPACK_lu_tag LuTag
					 , bool positive );


size_t
ReducedRowEchelonForm_modular_double (const double p, const size_t M, const size_t N,
				      double * A, const size_t lda,
				      size_t* P, size_t* Qt, const bool transform,
				      const enum FFPACK_lu_tag LuTag
				      , bool positive );


size_t
ReducedRowEchelonForm2_modular_double (const double p, const size_t M, const size_t N,
				       double * A, const size_t lda,
				       size_t* P, size_t* Qt, const bool transform
				       , bool positive );


size_t
REF_modular_double (const double p, const size_t M, const size_t N,
		    double * A, const size_t lda,
		    const size_t colbeg, const size_t rowbeg, const size_t colsize,
		    size_t* Qt, size_t* P
		    , bool positive );


/*****************/
/*   INVERSION   */
/*****************/


double *
Invertin_modular_double (const double p, const size_t M,
			 double * A, const size_t lda,
			 int * nullity
			 , bool positive );


double *
Invert_modular_double (const double p, const size_t M,
		       const double * A, const size_t lda,
		       double * X, const size_t ldx,
		       int* nullity
		       , bool positive );


double *
Invert2_modular_double( const double p, const size_t M,
			double * A, const size_t lda,
			double * X, const size_t ldx,
			int* nullity
			, bool positive );

/*****************************/
/* CHARACTERISTIC POLYNOMIAL */
/*****************************/


#if 0 /*  pas pour le moment */
template <class Polynomial>
std::list<Polynomial>&
CharPoly( const double p, std::list<Polynomial>& charp, const size_t N,
	  double * A, const size_t lda,
	  const enum FFPACK_CHARPOLY_TAG CharpTag= FfpackArithProg);

template<class Polynomial>
Polynomial & mulpoly_modular_double(const double p, Polynomial &res, const Polynomial & P1, const Polynomial & P2);

template <class Polynomial>
Polynomial&
CharPoly_modular_double( const double p, Polynomial& charp, const size_t N,
			 double * A, const size_t lda,
			 const enum FFPACK_CHARPOLY_TAG CharpTag= FfpackArithProg);



template <class Polynomial>
std::list<Polynomial>&
CharpolyArithProg_modular_double (const double p, std::list<Polynomial>& frobeniusForm,
				  const size_t N, double * A, const size_t lda, const size_t c);
#endif



/**********************/
/* MINIMAL POLYNOMIAL */
/**********************/

#if 0 /*  pas pour le moment */
template <class Polynomial>
Polynomial&
MinPoly_modular_double( const double p, Polynomial& minP, const size_t N,
			const double * A, const size_t lda,
			double * X, const size_t ldx, size_t* P,
			const enum FFPACK::FFPACK_MINPOLY_TAG MinTag= FFPACK::FfpackDense,
			const size_t kg_mc=0, const size_t kg_mb=0, const size_t kg_j=0 );
#endif


/* Krylov Elim */


size_t KrylovElim_modular_double( const double p, const size_t M, const size_t N,
				  double * A, const size_t lda, size_t*P,
				  size_t *Q, const size_t deg, size_t *iterates, size_t * inviterates, const size_t maxit,size_t virt
				  , bool positive );


size_t  SpecRankProfile_modular_double (const double p, const size_t M, const size_t N,
					double * A, const size_t lda, const size_t deg, size_t *rankProfile
					, bool positive );


/********/
/* RANK */
/********/


size_t
Rank_modular_double( const double p, const size_t M, const size_t N,
		     double * A, const size_t lda
		     , bool positive ) ;

/********/
/* DET  */
/********/


bool
IsSingular_modular_double( const double p, const size_t M, const size_t N,
			   double * A, const size_t lda
			   , bool positive );


double
Det_modular_double( const double p, const size_t M, const size_t N,
		    double * A, const size_t lda
		    , bool positive );

/*********/
/* SOLVE */
/*********/



double *
Solve_modular_double( const double p, const size_t M,
		      double * A, const size_t lda,
		      double * x, const int incx,
		      const double * b, const int incb
		      , bool positive  );



void
solveLB_modular_double( const double p, const enum FFLAS_side Side,
			const size_t M, const size_t N, const size_t R,
			double * L, const size_t ldl,
			const size_t * Q,
			double * B, const size_t ldb );


void
solveLB2_modular_double( const double p, const enum FFLAS_side Side,
			 const size_t M, const size_t N, const size_t R,
			 double * L, const size_t ldl,
			 const size_t * Q,
			 double * B, const size_t ldb
			 , bool positive  );


/*************/
/* NULLSPACE */
/*************/


void RandomNullSpaceVector_modular_double (const double p, const enum FFLAS_side Side,
					   const size_t M, const size_t N,
					   double * A, const size_t lda,
					   double * X, const size_t incX
					   , bool positive );


size_t NullSpaceBasis_modular_double (const double p, const enum FFLAS_side Side,
				      const size_t M, const size_t N,
				      double * A, const size_t lda,
				      double ** NS, size_t* ldn,
				      size_t * NSdim
				      , bool positive );

/*****************/
/* RANK PROFILES */
/*****************/


size_t RowRankProfile_modular_double (const double p, const size_t M, const size_t N,
				      double * A, const size_t lda,
				      size_t ** rkprofile,
				      const enum FFPACK_lu_tag LuTag
				      , bool positive );



size_t ColumnRankProfile_modular_double (const double p, const size_t M, const size_t N,
					 double * A, const size_t lda,
					 size_t ** rkprofile,
					 const enum FFPACK_lu_tag LuTag
					 , bool positive );

void RankProfileFromLU (const size_t* P, const size_t N, const size_t R,
			size_t* rkprofile, const enum FFPACK_lu_tag LuTag);

size_t LeadingSubmatrixRankProfiles (const size_t M, const size_t N, const size_t R,
				     const size_t LSm, const size_t LSn,
				     const size_t* P, const size_t* Q,
				     size_t* RRP, size_t* CRP);


size_t RowRankProfileSubmatrixIndices_modular_double (const double p,
						      const size_t M, const size_t N,
						      double * A,
						      const size_t lda,
						      size_t ** rowindices,
						      size_t ** colindices,
						      size_t * R
						      , bool positive );


size_t ColRankProfileSubmatrixIndices_modular_double (const double p,
						      const size_t M, const size_t N,
						      double * A,
						      const size_t lda,
						      size_t** rowindices,
						      size_t** colindices,
						      size_t* R
						      , bool positive );


size_t RowRankProfileSubmatrix_modular_double (const double p,
					       const size_t M, const size_t N,
					       double * A,
					       const size_t lda,
					       double ** X, size_t* R
					       , bool positive );


size_t ColRankProfileSubmatrix_modular_double (const double p, const size_t M, const size_t N,
					       double * A, const size_t lda,
					       double ** X, size_t* R
					       , bool positive );

/*********************************************/
/* Accessors to Triangular and Echelon forms */
/*********************************************/


void
getTriangular_modular_double (const double p, const enum FFLAS_uplo Uplo,
			      const enum FFLAS_diag diag,
			      const size_t M, const size_t N, const size_t R,
			      const double * A, const size_t lda,
			      double * T, const size_t ldt,
			      const bool OnlyNonZeroVectors
			      , bool positive );


void
getTriangularin_modular_double (const double p, const enum FFLAS_uplo Uplo,
				const enum FFLAS_diag diag,
				const size_t M, const size_t N, const size_t R,
				double * A, const size_t lda
				, bool positive );


void
getEchelonForm_modular_double (const double p, const enum FFLAS_uplo Uplo,
			       const enum FFLAS_diag diag,
			       const size_t M, const size_t N, const size_t R, const size_t* P,
			       const double * A, const size_t lda,
			       double * T, const size_t ldt,
			       const bool OnlyNonZeroVectors,
			       const enum FFPACK_lu_tag LuTag
			       , bool positive );


void
getEchelonFormin_modular_double (const double p, const enum FFLAS_uplo Uplo,
				 const enum FFLAS_diag diag,
				 const size_t M, const size_t N, const size_t R, const size_t* P,
				 double * A, const size_t lda,
				 const enum FFPACK_lu_tag LuTag
				 , bool positive );

void
getEchelonTransform_modular_double (const double p, const enum FFLAS_uplo Uplo,
				    const enum FFLAS_diag diag,
				    const size_t M, const size_t N, const size_t R, const size_t* P, const size_t* Q,
				    const double * A, const size_t lda,
				    double * T, const size_t ldt,
				    const enum FFPACK_lu_tag LuTag
				    , bool positive );


void
getReducedEchelonForm_modular_double (const double p, const enum FFLAS_uplo Uplo,
				      const size_t M, const size_t N, const size_t R, const size_t* P,
				      const double * A, const size_t lda,
				      double * T, const size_t ldt,
				      const bool OnlyNonZeroVectors,
				      const enum FFPACK_lu_tag LuTag
				      , bool positive );


void
getReducedEchelonFormin_modular_double (const double p, const enum FFLAS_uplo Uplo,
					const size_t M, const size_t N, const size_t R, const size_t* P,
					double * A, const size_t lda,
					const enum FFPACK_lu_tag LuTag
					, bool positive );


void
getReducedEchelonTransform_modular_double (const double p, const enum FFLAS_uplo Uplo,
					   const size_t M, const size_t N, const size_t R, const size_t* P, const size_t* Q,
					   const double * A, const size_t lda,
					   double * T, const size_t ldt,
					   const enum FFPACK_lu_tag LuTag
					   , bool positive );

void
PLUQtoEchelonPermutation (const size_t N, const size_t R, const size_t * P, size_t * outPerm);

#ifdef __cplusplus
}
#endif

#endif // __FFLASFFPACK_interfaces_libs_ffpack_c_H

