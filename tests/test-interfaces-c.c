/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) FFLAS-FFPACK
 * This file is Free Software and part of FFLAS-FFPACK.
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
#include <fflas-ffpack/interfaces/libs/fflas_c.h>
#include <fflas-ffpack/interfaces/libs/ffpack_c.h>

#include <stdlib.h>
#include <stdio.h>

int main() {
	double * A = (double*)malloc(4*sizeof(double));
	A[0] = A[2] = 1 ;
	A[1] = A[3] = 0 ;
	size_t * P = (size_t*) malloc(2*sizeof(size_t));
	size_t * Qt = (size_t*) malloc(2*sizeof(size_t));
            
        size_t r = RowEchelonForm_modular_double(101,2,2,A,2,P,Qt,false,FfpackSlabRecursive,true);
        freducein_2_modular_double(101,2,2,A,2,false);
	freducein_1_modular_double(101,4,A,1,false);
	fsquare_3_modular_double(101,FflasNoTrans,2,1,A,2,1,A,1,true);
	return !(r==1);
}

