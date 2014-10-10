/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2014 FFLAS-FFPACK group
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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

#ifndef __FFLASFFPACK_fscal_mp_INL
#define __FFLASFFPACK_fscal_mp_INL

// activate only if FFLAS-FFPACK haves multiprecision integer
#ifdef __FFLASFFPACK_HAVE_INTEGER
#include "fflas-ffpack/field/rns-integer.h"
#include "fflas-ffpack/field/rns-integer-mod.h"

namespace FFLAS {

	// specialization of the level1 fscalin function for the field RNSInteger<rns_double>
	template<>
	void fscalin(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F,  const size_t n,
		     const FFPACK::rns_double::Element alpha,
		     FFPACK::rns_double::Element_ptr A, const size_t inc) 
	{
#ifdef BENCH_PERF_SCAL_MP
		FFLAS::Timer chrono;chrono.start();
#endif
		for (size_t i=0;i<F.size();i++)
			fscalin(F.rns()._field_rns[i], n, alpha._ptr[i*alpha._stride], A._ptr+i*A._stride,inc);
#ifdef BENCH_PERF_SCAL_MP
		chrono.stop();F.t_scal+=chrono.usertime();
#endif
		finit(F,n,A,inc);
	}
	template<>
	void fscal(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F,  const size_t n,
		   const FFPACK::rns_double::Element alpha,
		   FFPACK::rns_double::ConstElement_ptr A, const size_t Ainc,
		   FFPACK::rns_double::Element_ptr B, const size_t Binc) 
	{
#ifdef BENCH_PERF_SCAL_MP
		FFLAS::Timer chrono;chrono.start();
#endif
		for (size_t i=0;i<F.size();i++)
			fscal(F.rns()._field_rns[i], n, alpha._ptr[i*alpha._stride], A._ptr+i*A._stride,Ainc, B._ptr+i*B._stride,Binc);
#ifdef BENCH_PERF_SCAL_MP
		chrono.stop();F.t_scal+=chrono.usertime();
#endif
		finit(F,n,B,Binc);
		
	}
	
	// specialization of the level2 fscalin function for the field RNSInteger<rns_double>
	template<>
	void fscalin(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F,  const size_t m, const size_t n,
		     const FFPACK::rns_double::Element alpha,
		     FFPACK::rns_double::Element_ptr A, const size_t lda) {
#ifdef BENCH_PERF_SCAL_MP
		FFLAS::Timer chrono;chrono.start();
#endif
		for (size_t i=0;i<F.size();i++)
			fscalin(F.rns()._field_rns[i], m, n, alpha._ptr[i*alpha._stride], A._ptr+i*A._stride,lda);
#ifdef BENCH_PERF_SCAL_MP
				chrono.stop();F.t_scal+=chrono.usertime();
#endif
		finit(F,m,n,A,lda);
	}
	template<>
	void fscal(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F, const size_t m, const size_t n,
		   const FFPACK::rns_double::Element alpha,
		   FFPACK::rns_double::ConstElement_ptr A, const size_t lda,
		   FFPACK::rns_double::Element_ptr B, const size_t ldb) {
#ifdef BENCH_PERF_SCAL_MP
		FFLAS::Timer chrono;chrono.start();
#endif
		for (size_t i=0;i<F.size();i++)
			fscal(F.rns()._field_rns[i], m, n, alpha._ptr[i*alpha._stride], A._ptr+i*A._stride, lda, B._ptr+i*B._stride, ldb);
#ifdef BENCH_PERF_SCAL_MP
		chrono.stop();F.t_scal+=chrono.usertime();
#endif
		finit(F,n,B,ldb);
	}

} //end of namespace FFLAS

#endif // __FFLASFFPACK_HAVE_INTEGER


#endif 
