/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) the FFLAS-FFPACK group 2017
 * Written by Clément Pernet
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


#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>

#include "fflas-ffpack/utils/fflas_io.h"

#include "fflas-ffpack/utils/Matio.h"
using namespace std;
using namespace FFLAS;
using Givaro::Modular;

int main(){

	typedef Modular<double> Field;
	Field F(101);
	string file_dense = "data/mat.dense";
	string file_sms = "data/mat.sms";
	string outfile_dense = "data/out.dense";

	Field::Element_ptr A=NULL;
	size_t m,n;
	ReadMatrix (file_dense, F, m, n, A);

	WriteMatrix (std::cout<<"A = "<<std::endl,F,m,n,A,n,FflasSMS);

	WriteMatrix (outfile_dense,F,m,n,A,n);

	fflas_delete(A);	

	return 0;
}
