/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2014 FFLAS-FFPACK
 * Written by :
 *        BB <brice.boyer@lip6.fr>
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

#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas/fflas_fspmv.h"
#include "fflas-ffpack/utils/args-parser.h"

using namespace FFLAS ;
using namespace FFPACK ;

/* 9 x 13 ; 20                                */
/*     I           V              X           */
/*   _______________________________________  */
/* 1 | 0  0  0  0  0  0  0  0  0  0  0  0  0  */
/*   | 0  1  0  0  2  0  0  0  3  0  0  0  0  */
/*   | 0  0  0  0  0  0  0  0  0  0  0  0  0  */
/*   | 0  0 -1  0  0  0 10  0  0  0  0  0  0  */
/* 5 | 1  2  3  4  0  0  7  8  9  0 11 12 13  */
/*   | 0  0  0 -2 -3  4  5  0  0  0  0  0  0  */
/*   | 0  0  0  0  0  0  0  0  0  0  0  0  0  */
/*   | 0  0  0  0  0  0  0  0  0  0  0  0  1  */
/*   | 0  0  0  0  0  0  0  0  0  0  0  0  0  */
/*   | 0  0  0  0  0  0  0  0  0  0  0  8  0  */

int test1_coo()
{
	Modular<double> F(1051);

	COO<double> Mat ;
	size_t nbnz = 21 ;
	Mat.m = 10 ;
	Mat.n = 13 ;
	Mat.row = FFLAS::fflas_new<index_t>(nbnz) ;
	Mat.col = FFLAS::fflas_new<index_t>(nbnz);
	Mat.dat = FFLAS::fflas_new<double >(nbnz);

	Mat.dat[0] = 1 ;
	Mat.dat[1] = 2 ;
	Mat.dat[2] = 3 ;
	Mat.dat[3] = -1 ;
	Mat.dat[4] = 10 ;
	Mat.dat[5] = 1 ;
	Mat.dat[6] = 2 ;
	Mat.dat[7] = 3 ;
	Mat.dat[8] = 4 ;
	Mat.dat[9] = 7 ;
	Mat.dat[10] = 8 ;
	Mat.dat[11] = 9 ;
	Mat.dat[12] = 11 ;
	Mat.dat[13] = 12 ;
	Mat.dat[14] = 13 ;
	Mat.dat[15] = -2 ;
	Mat.dat[16] = -3 ;
	Mat.dat[17] = 4 ;
	Mat.dat[18] = 5 ;
	Mat.dat[19] = 1 ;
	Mat.dat[20] = 8 ;

	Mat.row[0 ] = 1 ;
	Mat.row[1 ] = 1 ;
	Mat.row[2 ] = 1 ;
	Mat.row[3 ] = 3 ;
	Mat.row[4 ] = 3 ;
	Mat.row[5 ] = 4 ;
	Mat.row[6 ] = 4 ;
	Mat.row[7 ] = 4 ;
	Mat.row[8 ] = 4 ;
	Mat.row[9 ] = 4 ;
	Mat.row[10] = 4 ;
	Mat.row[11] = 4 ;
	Mat.row[12] = 4 ;
	Mat.row[13] = 4 ;
	Mat.row[14] = 4 ;
	Mat.row[15] = 5 ;
	Mat.row[16] = 5 ;
	Mat.row[17] = 5 ;
	Mat.row[18] = 5 ;
	Mat.row[19] = 7 ;
	Mat.row[20] = 9 ;

	Mat.col[0 ] = 1 ;
	Mat.col[1 ] = 4 ;
	Mat.col[2 ] = 8 ;
	Mat.col[3 ] = 2 ;
	Mat.col[4 ] = 6 ;
	Mat.col[5 ] = 0 ;
	Mat.col[6 ] = 1 ;
	Mat.col[7 ] = 2 ;
	Mat.col[8 ] = 3 ;
	Mat.col[9 ] = 6 ;
	Mat.col[10] = 7 ;
	Mat.col[11] = 8 ;
	Mat.col[12] = 10 ;
	Mat.col[13] = 11 ;
	Mat.col[14] = 12 ;
	Mat.col[15] = 3 ;
	Mat.col[16] = 4 ;
	Mat.col[17] = 5 ;
	Mat.col[18] = 6 ;
	Mat.col[19] = 12 ;
	Mat.col[20] = 11 ;



	Mat.z =  nbnz;


	VECT<double> x,y ;
	x.m = Mat.n ;
	x.inc = 1 ;
	x.dat = FFLAS::fflas_new<double>(x.m);

	y.m = Mat.m ;
	y.inc = 1 ;
	y.dat = FFLAS::fflas_new<double>(y.m);

	for (size_t i = 0 ; i < x.m ; ++i) {
		F.init(x.dat[i],i+1) ;
	}

	for (size_t i = 0 ; i < y.m ; ++i) {
		F.init(y.dat[i] ,i+1) ;
	}

	// y
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax + 0 y
	sp_fgemv(F,Mat,x,0,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;


	// y = Ax + y ( y = 2 Ax)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax - y (y = - Ax)
	sp_fgemv(F,Mat,x,-1,y);

	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + 3 y (y = - Ax)
	sp_fgemv(F,Mat,x,2,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + y (y = 0)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;



	return 0 ;


}

int test2_coo()
{
	Modular<double> F(1051);

	COO_sub<double> Mat ;
	size_t nbnz = 21 ;
	Mat.m = 10 ;
	Mat.n = 13 ;
	Mat.row = FFLAS::fflas_new<index_t>(nbnz) ;
	Mat.col = FFLAS::fflas_new<index_t>(nbnz);
	Mat.dat = FFLAS::fflas_new<double >(nbnz);

	Mat.dat[0] = 1 ;
	Mat.dat[1] = 2 ;
	Mat.dat[2] = 3 ;
	Mat.dat[3] = -1 ;
	Mat.dat[4] = 10 ;
	Mat.dat[5] = 1 ;
	Mat.dat[6] = 2 ;
	Mat.dat[7] = 3 ;
	Mat.dat[8] = 4 ;
	Mat.dat[9] = 7 ;
	Mat.dat[10] = 8 ;
	Mat.dat[11] = 9 ;
	Mat.dat[12] = 11 ;
	Mat.dat[13] = 12 ;
	Mat.dat[14] = 13 ;
	Mat.dat[15] = -2 ;
	Mat.dat[16] = -3 ;
	Mat.dat[17] = 4 ;
	Mat.dat[18] = 5 ;
	Mat.dat[19] = 1 ;
	Mat.dat[20] = 8 ;

	Mat.row[0 ] = 1 ;
	Mat.row[1 ] = 1 ;
	Mat.row[2 ] = 1 ;
	Mat.row[3 ] = 3 ;
	Mat.row[4 ] = 3 ;
	Mat.row[5 ] = 4 ;
	Mat.row[6 ] = 4 ;
	Mat.row[7 ] = 4 ;
	Mat.row[8 ] = 4 ;
	Mat.row[9 ] = 4 ;
	Mat.row[10] = 4 ;
	Mat.row[11] = 4 ;
	Mat.row[12] = 4 ;
	Mat.row[13] = 4 ;
	Mat.row[14] = 4 ;
	Mat.row[15] = 5 ;
	Mat.row[16] = 5 ;
	Mat.row[17] = 5 ;
	Mat.row[18] = 5 ;
	Mat.row[19] = 7 ;
	Mat.row[20] = 9 ;

	Mat.col[0 ] = 1 ;
	Mat.col[1 ] = 4 ;
	Mat.col[2 ] = 8 ;
	Mat.col[3 ] = 2 ;
	Mat.col[4 ] = 6 ;
	Mat.col[5 ] = 0 ;
	Mat.col[6 ] = 1 ;
	Mat.col[7 ] = 2 ;
	Mat.col[8 ] = 3 ;
	Mat.col[9 ] = 6 ;
	Mat.col[10] = 7 ;
	Mat.col[11] = 8 ;
	Mat.col[12] = 10 ;
	Mat.col[13] = 11 ;
	Mat.col[14] = 12 ;
	Mat.col[15] = 3 ;
	Mat.col[16] = 4 ;
	Mat.col[17] = 5 ;
	Mat.col[18] = 6 ;
	Mat.col[19] = 12 ;
	Mat.col[20] = 11 ;



	Mat.z =  nbnz;


	VECT<double> x,y ;
	x.m = Mat.n ;
	x.inc = 1 ;
	x.dat = FFLAS::fflas_new<double>(x.m);

	y.m = Mat.m ;
	y.inc = 1 ;
	y.dat = FFLAS::fflas_new<double>(y.m);

	for (size_t i = 0 ; i < x.m ; ++i) {
		F.init(x.dat[i],i+1) ;
	}

	for (size_t i = 0 ; i < y.m ; ++i) {
		F.init(y.dat[i] ,i+1) ;
	}

	// y
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax + 0 y
	sp_fgemv(F,Mat,x,0,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;


	// y = Ax + y ( y = 2 Ax)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax - y (y = - Ax)
	sp_fgemv(F,Mat,x,-1,y);

	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + 3 y (y = - Ax)
	sp_fgemv(F,Mat,x,2,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + y (y = 0)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;


	return 0 ;


}

int test3_coo(int CST)
{

	{
		Modular<double> F(1051);

		COO<double> Mat ;
		size_t nbnz = 21 ;
		Mat.m = 10 ;
		Mat.n = 13 ;
		Mat.row = FFLAS::fflas_new<index_t>(nbnz) ;
		Mat.col = FFLAS::fflas_new<index_t>(nbnz);
		Mat.dat = FFLAS::fflas_new<double >(nbnz);

		Mat.dat[0] = CST ;
		Mat.dat[1] = CST ;
		Mat.dat[2] = CST ;
		Mat.dat[3] = CST ;
		Mat.dat[4] = CST ;
		Mat.dat[5] = CST ;
		Mat.dat[6] = CST ;
		Mat.dat[7] = CST ;
		Mat.dat[8] = CST ;
		Mat.dat[9] = CST ;
		Mat.dat[10] = CST ;
		Mat.dat[11] = CST ;
		Mat.dat[12] = CST ;
		Mat.dat[13] = CST ;
		Mat.dat[14] = CST ;
		Mat.dat[15] = CST ;
		Mat.dat[16] = CST ;
		Mat.dat[17] = CST ;
		Mat.dat[18] = CST ;
		Mat.dat[19] = CST ;
		Mat.dat[20] = CST ;

		Mat.row[0 ] = 1 ;
		Mat.row[1 ] = 1 ;
		Mat.row[2 ] = 1 ;
		Mat.row[3 ] = 3 ;
		Mat.row[4 ] = 3 ;
		Mat.row[5 ] = 4 ;
		Mat.row[6 ] = 4 ;
		Mat.row[7 ] = 4 ;
		Mat.row[8 ] = 4 ;
		Mat.row[9 ] = 4 ;
		Mat.row[10] = 4 ;
		Mat.row[11] = 4 ;
		Mat.row[12] = 4 ;
		Mat.row[13] = 4 ;
		Mat.row[14] = 4 ;
		Mat.row[15] = 5 ;
		Mat.row[16] = 5 ;
		Mat.row[17] = 5 ;
		Mat.row[18] = 5 ;
		Mat.row[19] = 7 ;
		Mat.row[20] = 9 ;

		Mat.col[0 ] = 1 ;
		Mat.col[1 ] = 4 ;
		Mat.col[2 ] = 8 ;
		Mat.col[3 ] = 2 ;
		Mat.col[4 ] = 6 ;
		Mat.col[5 ] = 0 ;
		Mat.col[6 ] = 1 ;
		Mat.col[7 ] = 2 ;
		Mat.col[8 ] = 3 ;
		Mat.col[9 ] = 6 ;
		Mat.col[10] = 7 ;
		Mat.col[11] = 8 ;
		Mat.col[12] = 10 ;
		Mat.col[13] = 11 ;
		Mat.col[14] = 12 ;
		Mat.col[15] = 3 ;
		Mat.col[16] = 4 ;
		Mat.col[17] = 5 ;
		Mat.col[18] = 6 ;
		Mat.col[19] = 12 ;
		Mat.col[20] = 11 ;



		Mat.z =  nbnz;


		VECT<double> x,y ;
		x.m = Mat.n ;
		x.inc = 1 ;
		x.dat = FFLAS::fflas_new<double>(x.m);

		y.m = Mat.m ;
		y.inc = 1 ;
		y.dat = FFLAS::fflas_new<double>(y.m);

		for (size_t i = 0 ; i < x.m ; ++i) {
			F.init(x.dat[i],i+1) ;
		}

		for (size_t i = 0 ; i < y.m ; ++i) {
			F.init(y.dat[i] ,i+1) ;
		}

		// y
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;

		// y = Ax + 0 y
		sp_fgemv(F,Mat,x,0,y);
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;


		// y = Ax + y ( y = 2 Ax)
		// sp_fgemv(F,Mat,x,1,y);
		// for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		// std::cout << std::endl;

		// y = Ax - y (y = - Ax)
		// sp_fgemv(F,Mat,x,-1,y);

		// for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		// std::cout << std::endl;

		// y += Ax + 3 y (y = - Ax)
		// sp_fgemv(F,Mat,x,2,y);
		// for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		// std::cout << std::endl;

		// y += Ax + y (y = 0)
		// sp_fgemv(F,Mat,x,1,y);
		// for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		// std::cout << std::endl;


	}

	{
		Modular<double> F(1051);
		COO_ZO<double> Mat ;

		size_t nbnz = 21 ;
		Mat.m = 10 ;
		Mat.n = 13 ;
		Mat.row = FFLAS::fflas_new<index_t>(nbnz) ;
		Mat.col = FFLAS::fflas_new<index_t>(nbnz);
		Mat.cst = CST ;

		Mat.row[0 ] = 1 ;
		Mat.row[1 ] = 1 ;
		Mat.row[2 ] = 1 ;
		Mat.row[3 ] = 3 ;
		Mat.row[4 ] = 3 ;
		Mat.row[5 ] = 4 ;
		Mat.row[6 ] = 4 ;
		Mat.row[7 ] = 4 ;
		Mat.row[8 ] = 4 ;
		Mat.row[9 ] = 4 ;
		Mat.row[10] = 4 ;
		Mat.row[11] = 4 ;
		Mat.row[12] = 4 ;
		Mat.row[13] = 4 ;
		Mat.row[14] = 4 ;
		Mat.row[15] = 5 ;
		Mat.row[16] = 5 ;
		Mat.row[17] = 5 ;
		Mat.row[18] = 5 ;
		Mat.row[19] = 7 ;
		Mat.row[20] = 9 ;

		Mat.col[0 ] = 1 ;
		Mat.col[1 ] = 4 ;
		Mat.col[2 ] = 8 ;
		Mat.col[3 ] = 2 ;
		Mat.col[4 ] = 6 ;
		Mat.col[5 ] = 0 ;
		Mat.col[6 ] = 1 ;
		Mat.col[7 ] = 2 ;
		Mat.col[8 ] = 3 ;
		Mat.col[9 ] = 6 ;
		Mat.col[10] = 7 ;
		Mat.col[11] = 8 ;
		Mat.col[12] = 10 ;
		Mat.col[13] = 11 ;
		Mat.col[14] = 12 ;
		Mat.col[15] = 3 ;
		Mat.col[16] = 4 ;
		Mat.col[17] = 5 ;
		Mat.col[18] = 6 ;
		Mat.col[19] = 12 ;
		Mat.col[20] = 11 ;



		Mat.z =  nbnz;


		VECT<double> x,y ;
		x.m = Mat.n ;
		x.inc = 1 ;
		x.dat = FFLAS::fflas_new<double>(x.m);

		y.m = Mat.m ;
		y.inc = 1 ;
		y.dat = FFLAS::fflas_new<double>(y.m);

		for (size_t i = 0 ; i < x.m ; ++i) {
			F.init(x.dat[i],i+1) ;
		}

		for (size_t i = 0 ; i < y.m ; ++i) {
			F.init(y.dat[i] ,i+1) ;
		}

		// y
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;

		// y = Ax + 0 y
		sp_fgemv(F,Mat,x,0,y);
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;


		// y = Ax + y ( y = 2 Ax)
		sp_fgemv(F,Mat,x,1,y);
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;

		// y = Ax - y (y = - Ax)
		sp_fgemv(F,Mat,x,-1,y);

		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;

		// y += Ax + 3 y (y = - Ax)
		sp_fgemv(F,Mat,x,2,y);
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;

		// y += Ax + y (y = 0)
		sp_fgemv(F,Mat,x,1,y);
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;




	}

	return 0 ;

}

int test1_csr()
{
	Modular<double> F(1051);

	CSR<double> Mat ;
	size_t nbnz = 21 ;
	Mat.m = 10 ;
	Mat.n = 13 ;
	Mat.st = FFLAS::fflas_new<index_t>(Mat.m+1) ;
	Mat.col = FFLAS::fflas_new<index_t>(nbnz);
	Mat.dat = FFLAS::fflas_new<double >(nbnz);


	Mat.dat[0] = 1 ;
	Mat.dat[1] = 2 ;
	Mat.dat[2] = 3 ;
	Mat.dat[3] = -1 ;
	Mat.dat[4] = 10 ;
	Mat.dat[5] = 1 ;
	Mat.dat[6] = 2 ;
	Mat.dat[7] = 3 ;
	Mat.dat[8] = 4 ;
	Mat.dat[9] = 7 ;
	Mat.dat[10] = 8 ;
	Mat.dat[11] = 9 ;
	Mat.dat[12] = 11 ;
	Mat.dat[13] = 12 ;
	Mat.dat[14] = 13 ;
	Mat.dat[15] = -2 ;
	Mat.dat[16] = -3 ;
	Mat.dat[17] = 4 ;
	Mat.dat[18] = 5 ;
	Mat.dat[19] = 1 ;
	Mat.dat[20] = 8 ;

	Mat.st[0 ] = 0 ;
	Mat.st[1 ] = 0 ;
	Mat.st[2 ] = 3 ;
	Mat.st[3 ] = 3 ;
	Mat.st[4 ] = 5 ;
	Mat.st[5 ] = 15 ;
	Mat.st[6 ] = 15 ;
	Mat.st[7 ] = 19 ;
	Mat.st[8 ] = 20 ;
	Mat.st[9 ] = 20 ;
	Mat.st[10] = 21 ;

	Mat.col[0 ] = 1 ;
	Mat.col[1 ] = 4 ;
	Mat.col[2 ] = 8 ;
	Mat.col[3 ] = 2 ;
	Mat.col[4 ] = 6 ;
	Mat.col[5 ] = 0 ;
	Mat.col[6 ] = 1 ;
	Mat.col[7 ] = 2 ;
	Mat.col[8 ] = 3 ;
	Mat.col[9 ] = 6 ;
	Mat.col[10] = 7 ;
	Mat.col[11] = 8 ;
	Mat.col[12] = 10 ;
	Mat.col[13] = 11 ;
	Mat.col[14] = 12 ;
	Mat.col[15] = 3 ;
	Mat.col[16] = 4 ;
	Mat.col[17] = 5 ;
	Mat.col[18] = 6 ;
	Mat.col[19] = 12 ;
	Mat.col[20] = 11 ;



	VECT<double> x,y ;
	x.m = Mat.n ;
	x.inc = 1 ;
	x.dat = FFLAS::fflas_new<double>(x.m);

	y.m = Mat.m ;
	y.inc = 1 ;
	y.dat = FFLAS::fflas_new<double>(y.m);

	for (size_t i = 0 ; i < x.m ; ++i) {
		F.init(x.dat[i],i+1) ;
	}

	for (size_t i = 0 ; i < y.m ; ++i) {
		F.init(y.dat[i] ,i+1) ;
	}

	sp_fgemv(F,Mat,x,0,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;


	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	sp_fgemv(F,Mat,x,-1,y);

	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	sp_fgemv(F,Mat,x,2,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;



	return 0 ;


}

int test2_csr()
{
	Modular<double> F(1051);

	CSR_sub<double> Mat ;
	size_t nbnz = 21 ;
	Mat.m = 10 ;
	Mat.n = 13 ;
	Mat.st = FFLAS::fflas_new<index_t>(Mat.m+1) ;
	Mat.col = FFLAS::fflas_new<index_t>(nbnz);
	Mat.dat = FFLAS::fflas_new<double >(nbnz);


	Mat.dat[0] = 1 ;
	Mat.dat[1] = 2 ;
	Mat.dat[2] = 3 ;
	Mat.dat[3] = -1 ;
	Mat.dat[4] = 10 ;
	Mat.dat[5] = 1 ;
	Mat.dat[6] = 2 ;
	Mat.dat[7] = 3 ;
	Mat.dat[8] = 4 ;
	Mat.dat[9] = 7 ;
	Mat.dat[10] = 8 ;
	Mat.dat[11] = 9 ;
	Mat.dat[12] = 11 ;
	Mat.dat[13] = 12 ;
	Mat.dat[14] = 13 ;
	Mat.dat[15] = -2 ;
	Mat.dat[16] = -3 ;
	Mat.dat[17] = 4 ;
	Mat.dat[18] = 5 ;
	Mat.dat[19] = 1 ;
	Mat.dat[20] = 8 ;

	Mat.st[0 ] = 0 ;
	Mat.st[1 ] = 0 ;
	Mat.st[2 ] = 3 ;
	Mat.st[3 ] = 3 ;
	Mat.st[4 ] = 5 ;
	Mat.st[5 ] = 15 ;
	Mat.st[6 ] = 15 ;
	Mat.st[7 ] = 19 ;
	Mat.st[8 ] = 20 ;
	Mat.st[9 ] = 20 ;
	Mat.st[10] = 21 ;

	Mat.col[0 ] = 1 ;
	Mat.col[1 ] = 4 ;
	Mat.col[2 ] = 8 ;
	Mat.col[3 ] = 2 ;
	Mat.col[4 ] = 6 ;
	Mat.col[5 ] = 0 ;
	Mat.col[6 ] = 1 ;
	Mat.col[7 ] = 2 ;
	Mat.col[8 ] = 3 ;
	Mat.col[9 ] = 6 ;
	Mat.col[10] = 7 ;
	Mat.col[11] = 8 ;
	Mat.col[12] = 10 ;
	Mat.col[13] = 11 ;
	Mat.col[14] = 12 ;
	Mat.col[15] = 3 ;
	Mat.col[16] = 4 ;
	Mat.col[17] = 5 ;
	Mat.col[18] = 6 ;
	Mat.col[19] = 12 ;
	Mat.col[20] = 11 ;



	VECT<double> x,y ;
	x.m = Mat.n ;
	x.inc = 1 ;
	x.dat = FFLAS::fflas_new<double>(x.m);

	y.m = Mat.m ;
	y.inc = 1 ;
	y.dat = FFLAS::fflas_new<double>(y.m);

	for (size_t i = 0 ; i < x.m ; ++i) {
		F.init(x.dat[i],i+1) ;
	}

	for (size_t i = 0 ; i < y.m ; ++i) {
		F.init(y.dat[i] ,i+1) ;
	}

	sp_fgemv(F,Mat,x,0,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;


	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	sp_fgemv(F,Mat,x,-1,y);

	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	sp_fgemv(F,Mat,x,2,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;



	return 0 ;


}

int test3_csr(int CST)
{

	{
		Modular<double> F(1051);

		CSR<double> Mat ;
		size_t nbnz = 21 ;
		Mat.m = 10 ;
		Mat.n = 13 ;
		Mat.st = FFLAS::fflas_new<index_t>(Mat.m+1) ;
		Mat.col = FFLAS::fflas_new<index_t>(nbnz);
		Mat.dat = FFLAS::fflas_new<double >(nbnz);

		Mat.dat[0] = CST ;
		Mat.dat[1] = CST ;
		Mat.dat[2] = CST ;
		Mat.dat[3] = CST ;
		Mat.dat[4] = CST ;
		Mat.dat[5] = CST ;
		Mat.dat[6] = CST ;
		Mat.dat[7] = CST ;
		Mat.dat[8] = CST ;
		Mat.dat[9] = CST ;
		Mat.dat[10] = CST ;
		Mat.dat[11] = CST ;
		Mat.dat[12] = CST ;
		Mat.dat[13] = CST ;
		Mat.dat[14] = CST ;
		Mat.dat[15] = CST ;
		Mat.dat[16] = CST ;
		Mat.dat[17] = CST ;
		Mat.dat[18] = CST ;
		Mat.dat[19] = CST ;
		Mat.dat[20] = CST ;

		Mat.st[0 ] = 0 ;
		Mat.st[1 ] = 0 ;
		Mat.st[2 ] = 3 ;
		Mat.st[3 ] = 3 ;
		Mat.st[4 ] = 5 ;
		Mat.st[5 ] = 15 ;
		Mat.st[6 ] = 15 ;
		Mat.st[7 ] = 19 ;
		Mat.st[8 ] = 20 ;
		Mat.st[9 ] = 20 ;
		Mat.st[10] = 21 ;


		Mat.col[0 ] = 1 ;
		Mat.col[1 ] = 4 ;
		Mat.col[2 ] = 8 ;
		Mat.col[3 ] = 2 ;
		Mat.col[4 ] = 6 ;
		Mat.col[5 ] = 0 ;
		Mat.col[6 ] = 1 ;
		Mat.col[7 ] = 2 ;
		Mat.col[8 ] = 3 ;
		Mat.col[9 ] = 6 ;
		Mat.col[10] = 7 ;
		Mat.col[11] = 8 ;
		Mat.col[12] = 10 ;
		Mat.col[13] = 11 ;
		Mat.col[14] = 12 ;
		Mat.col[15] = 3 ;
		Mat.col[16] = 4 ;
		Mat.col[17] = 5 ;
		Mat.col[18] = 6 ;
		Mat.col[19] = 12 ;
		Mat.col[20] = 11 ;



		// Mat.z =  nbnz;


		VECT<double> x,y ;
		x.m = Mat.n ;
		x.inc = 1 ;
		x.dat = FFLAS::fflas_new<double>(x.m);

		y.m = Mat.m ;
		y.inc = 1 ;
		y.dat = FFLAS::fflas_new<double>(y.m);

		for (size_t i = 0 ; i < x.m ; ++i) {
			F.init(x.dat[i],i+1) ;
		}

		for (size_t i = 0 ; i < y.m ; ++i) {
			F.init(y.dat[i] ,i+1) ;
		}

		// y
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;

		// y = Ax + 0 y
		sp_fgemv(F,Mat,x,0,y);
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;


		// y = Ax + y ( y = 2 Ax)
		// sp_fgemv(F,Mat,x,1,y);
		// for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		// std::cout << std::endl;

		// y = Ax - y (y = - Ax)
		// sp_fgemv(F,Mat,x,-1,y);

		// for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		// std::cout << std::endl;

		// y += Ax + 3 y (y = - Ax)
		// sp_fgemv(F,Mat,x,2,y);
		// for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		// std::cout << std::endl;

		// y += Ax + y (y = 0)
		// sp_fgemv(F,Mat,x,1,y);
		// for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		// std::cout << std::endl;


	}

	{
		Modular<double> F(1051);
		CSR_ZO<double> Mat ;

		size_t nbnz = 21 ;
		Mat.m = 10 ;
		Mat.n = 13 ;
		Mat.st = FFLAS::fflas_new<index_t>(Mat.m+1) ;
		Mat.col = FFLAS::fflas_new<index_t>(nbnz);
		Mat.cst = CST ;

		Mat.col[0 ] = 1 ;
		Mat.col[1 ] = 4 ;
		Mat.col[2 ] = 8 ;
		Mat.col[3 ] = 2 ;
		Mat.col[4 ] = 6 ;
		Mat.col[5 ] = 0 ;
		Mat.col[6 ] = 1 ;
		Mat.col[7 ] = 2 ;
		Mat.col[8 ] = 3 ;
		Mat.col[9 ] = 6 ;
		Mat.col[10] = 7 ;
		Mat.col[11] = 8 ;
		Mat.col[12] = 10 ;
		Mat.col[13] = 11 ;
		Mat.col[14] = 12 ;
		Mat.col[15] = 3 ;
		Mat.col[16] = 4 ;
		Mat.col[17] = 5 ;
		Mat.col[18] = 6 ;
		Mat.col[19] = 12 ;
		Mat.col[20] = 11 ;

		Mat.st[0 ] = 0 ;
		Mat.st[1 ] = 0 ;
		Mat.st[2 ] = 3 ;
		Mat.st[3 ] = 3 ;
		Mat.st[4 ] = 5 ;
		Mat.st[5 ] = 15 ;
		Mat.st[6 ] = 15 ;
		Mat.st[7 ] = 19 ;
		Mat.st[8 ] = 20 ;
		Mat.st[9 ] = 20 ;
		Mat.st[10] = 21 ;



		// Mat.z =  nbnz;


		VECT<double> x,y ;
		x.m = Mat.n ;
		x.inc = 1 ;
		x.dat = FFLAS::fflas_new<double>(x.m);

		y.m = Mat.m ;
		y.inc = 1 ;
		y.dat = FFLAS::fflas_new<double>(y.m);

		for (size_t i = 0 ; i < x.m ; ++i) {
			F.init(x.dat[i],i+1) ;
		}

		for (size_t i = 0 ; i < y.m ; ++i) {
			F.init(y.dat[i] ,i+1) ;
		}

		// y
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;

		// y = Ax + 0 y
		sp_fgemv(F,Mat,x,0,y);
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;


		// y = Ax + y ( y = 2 Ax)
		sp_fgemv(F,Mat,x,1,y);
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;

		// y = Ax - y (y = - Ax)
		sp_fgemv(F,Mat,x,-1,y);

		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;

		// y += Ax + 3 y (y = - Ax)
		sp_fgemv(F,Mat,x,2,y);
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;

		// y += Ax + y (y = 0)
		sp_fgemv(F,Mat,x,1,y);
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;




	}

	return 0 ;

}

int test1_ell()
{
	Modular<double> F(1051);

	ELL<double,false> Mat ;
	// size_t nbnz = 20 ;
	Mat.m = 10 ;
	Mat.n = 13 ;
	size_t ld = 10 ;
	Mat.ld = ld ;
	Mat.col = FFLAS::fflas_new<index_t>(ld*Mat.m);
	Mat.dat = FFLAS::fflas_new<double >(ld*Mat.m);
	for (size_t i = 0 ; i < ld*Mat.m ; ++i) Mat.col[i]  = 0 ;
	for (size_t i = 0 ; i < ld*Mat.m ; ++i) Mat.dat[i]  = 0 ;

	size_t j = 0 ;
	++j ;
	Mat.dat[j*ld+0] = 1 ;
	Mat.dat[j*ld+1] = 2 ;
	Mat.dat[j*ld+2] = 3 ;
	++j;
	++j;
	Mat.dat[j*ld+0] = -1 ;
	Mat.dat[j*ld+1] = 10 ;
	++j;
	Mat.dat[j*ld+0] = 1 ;
	Mat.dat[j*ld+1] = 2 ;
	Mat.dat[j*ld+2] = 3 ;
	Mat.dat[j*ld+3] = 4 ;
	Mat.dat[j*ld+4] = 7 ;
	Mat.dat[j*ld+5] = 8 ;
	Mat.dat[j*ld+6] = 9 ;
	Mat.dat[j*ld+7] = 11 ;
	Mat.dat[j*ld+8] = 12 ;
	Mat.dat[j*ld+9] = 13 ;
	++j;
	Mat.dat[j*ld+0] = -2 ;
	Mat.dat[j*ld+1] = -3 ;
	Mat.dat[j*ld+2] = 4 ;
	Mat.dat[j*ld+3] = 5 ;
	++j;
	++j;
	Mat.dat[j*ld+0] = 1 ;
	++j;
	++j;
	Mat.dat[j*ld+0] = 8 ;

	j = 0 ;
	++j ;
	Mat.col[j*ld+0] = 1 ;
	Mat.col[j*ld+1] = 4 ;
	Mat.col[j*ld+2] = 8 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 2 ;
	Mat.col[j*ld+1] = 6 ;
	++j;
	Mat.col[j*ld+0] = 0 ;
	Mat.col[j*ld+1] = 1 ;
	Mat.col[j*ld+2] = 2 ;
	Mat.col[j*ld+3] = 3 ;
	Mat.col[j*ld+4] = 6 ;
	Mat.col[j*ld+5] = 7 ;
	Mat.col[j*ld+6] = 8 ;
	Mat.col[j*ld+7] = 10 ;
	Mat.col[j*ld+8] = 11 ;
	Mat.col[j*ld+9] = 12 ;
	++j;
	Mat.col[j*ld+0] = 3 ;
	Mat.col[j*ld+1] = 4 ;
	Mat.col[j*ld+2] = 5 ;
	Mat.col[j*ld+3] = 6 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 12 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 11 ;



	VECT<double> x,y ;
	x.m = Mat.n ;
	x.inc = 1 ;
	x.dat = FFLAS::fflas_new<double>(x.m);

	y.m = Mat.m ;
	y.inc = 1 ;
	y.dat = FFLAS::fflas_new<double>(y.m);

	for (size_t i = 0 ; i < x.m ; ++i) {
		F.init(x.dat[i],i+1) ;
	}

	for (size_t i = 0 ; i < y.m ; ++i) {
		F.init(y.dat[i] ,i+1) ;
	}

	// y
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax + 0 y
	sp_fgemv(F,Mat,x,0,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;


	// y = Ax + y ( y = 2 Ax)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax - y (y = - Ax)
	sp_fgemv(F,Mat,x,-1,y);

	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + 3 y (y = - Ax)
	sp_fgemv(F,Mat,x,2,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + y (y = 0)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;



	return 0 ;


}

int test2_ell()
{
	Modular<double> F(1051);

	ELL_sub<double,false> Mat ;
	// size_t nbnz = 20 ;
	Mat.m = 10 ;
	Mat.n = 13 ;
	size_t ld = 10 ;
	Mat.ld = ld ;
	Mat.col = FFLAS::fflas_new<index_t>(ld*Mat.m);
	Mat.dat = FFLAS::fflas_new<double >(ld*Mat.m);
	for (size_t i = 0 ; i < ld*Mat.m ; ++i) Mat.col[i]  = 0 ;
	for (size_t i = 0 ; i < ld*Mat.m ; ++i) Mat.dat[i]  = 0 ;

	size_t j = 0 ;
	++j ;
	Mat.dat[j*ld+0] = 1 ;
	Mat.dat[j*ld+1] = 2 ;
	Mat.dat[j*ld+2] = 3 ;
	++j;
	++j;
	Mat.dat[j*ld+0] = -1 ;
	Mat.dat[j*ld+1] = 10 ;
	++j;
	Mat.dat[j*ld+0] = 1 ;
	Mat.dat[j*ld+1] = 2 ;
	Mat.dat[j*ld+2] = 3 ;
	Mat.dat[j*ld+3] = 4 ;
	Mat.dat[j*ld+4] = 7 ;
	Mat.dat[j*ld+5] = 8 ;
	Mat.dat[j*ld+6] = 9 ;
	Mat.dat[j*ld+7] = 11 ;
	Mat.dat[j*ld+8] = 12 ;
	Mat.dat[j*ld+9] = 13 ;
	++j;
	Mat.dat[j*ld+0] = -2 ;
	Mat.dat[j*ld+1] = -3 ;
	Mat.dat[j*ld+2] = 4 ;
	Mat.dat[j*ld+3] = 5 ;
	++j;
	++j;
	Mat.dat[j*ld+0] = 1 ;
	++j;
	++j;
	Mat.dat[j*ld+0] = 8 ;

	j = 0 ;
	++j ;
	Mat.col[j*ld+0] = 1 ;
	Mat.col[j*ld+1] = 4 ;
	Mat.col[j*ld+2] = 8 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 2 ;
	Mat.col[j*ld+1] = 6 ;
	++j;
	Mat.col[j*ld+0] = 0 ;
	Mat.col[j*ld+1] = 1 ;
	Mat.col[j*ld+2] = 2 ;
	Mat.col[j*ld+3] = 3 ;
	Mat.col[j*ld+4] = 6 ;
	Mat.col[j*ld+5] = 7 ;
	Mat.col[j*ld+6] = 8 ;
	Mat.col[j*ld+7] = 10 ;
	Mat.col[j*ld+8] = 11 ;
	Mat.col[j*ld+9] = 12 ;
	++j;
	Mat.col[j*ld+0] = 3 ;
	Mat.col[j*ld+1] = 4 ;
	Mat.col[j*ld+2] = 5 ;
	Mat.col[j*ld+3] = 6 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 12 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 11 ;



	VECT<double> x,y ;
	x.m = Mat.n ;
	x.inc = 1 ;
	x.dat = FFLAS::fflas_new<double>(x.m);

	y.m = Mat.m ;
	y.inc = 1 ;
	y.dat = FFLAS::fflas_new<double>(y.m);

	for (size_t i = 0 ; i < x.m ; ++i) {
		F.init(x.dat[i],i+1) ;
	}

	for (size_t i = 0 ; i < y.m ; ++i) {
		F.init(y.dat[i] ,i+1) ;
	}

	// y
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax + 0 y
	sp_fgemv(F,Mat,x,0,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;


	// y = Ax + y ( y = 2 Ax)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax - y (y = - Ax)
	sp_fgemv(F,Mat,x,-1,y);

	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + 3 y (y = - Ax)
	sp_fgemv(F,Mat,x,2,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + y (y = 0)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;



	return 0 ;


}

#ifdef __FFLASFFPACK_USE_SIMD
#define projette(a,b) \
	(a)/(simd::vect_size)*ld*simd::vect_size+(b)*simd::vect_size+(a)%(simd::vect_size)

int test1_ell_simd()
{
	Modular<double> F(1051);
	using simd = Simd<double>;
	// using vect_t = typename simd::vect_t;

	ELL<double,true> Mat ;
	// size_t nbnz = 20 ;
	Mat.m = 10 ;
	Mat.n = 13 ;
	size_t ld = 10 ;

	{

		Mat.ld = ld ;
		Mat.col = FFLAS::fflas_new<index_t>(ld*Mat.m);
		Mat.dat = FFLAS::fflas_new<double >(ld*Mat.m);
		for (size_t i = 0 ; i < ld*Mat.m ; ++i) Mat.col[i]  = 0 ;
		for (size_t i = 0 ; i < ld*Mat.m ; ++i) Mat.dat[i]  = 0 ;

		size_t i = 0 ;
		++i ;
		Mat.dat[projette(i,0)] = 1 ;
		Mat.dat[projette(i,1)] = 2 ;
		Mat.dat[projette(i,2)] = 3 ;
		++i;
		++i;
		Mat.dat[projette(i,0)] = -1 ;
		Mat.dat[projette(i,1)] = 10 ;
		++i;
		Mat.dat[projette(i,0)] = 1 ;
		Mat.dat[projette(i,1)] = 2 ;
		Mat.dat[projette(i,2)] = 3 ;
		Mat.dat[projette(i,3)] = 4 ;
		Mat.dat[projette(i,4)] = 7 ;
		Mat.dat[projette(i,5)] = 8 ;
		Mat.dat[projette(i,6)] = 9 ;
		Mat.dat[projette(i,7)] = 11 ;
		Mat.dat[projette(i,8)] = 12 ;
		Mat.dat[projette(i,9)] = 13 ;
		++i;
		Mat.dat[projette(i,0)] = -2 ;
		Mat.dat[projette(i,1)] = -3 ;
		Mat.dat[projette(i,2)] = 4 ;
		Mat.dat[projette(i,3)] = 5 ;
		++i;
		++i;
		Mat.dat[projette(i,0)] = 1 ;
		++i;
		++i;
		Mat.dat[projette(i,0)] = 8 ;

		i = 0 ;
		++i ;
		Mat.col[projette(i,0)] = 1 ;
		Mat.col[projette(i,1)] = 4 ;
		Mat.col[projette(i,2)] = 8 ;
		++i;
		++i;
		Mat.col[projette(i,0)] = 2 ;
		Mat.col[projette(i,1)] = 6 ;
		++i;
		Mat.col[projette(i,0)] = 0 ;
		Mat.col[projette(i,1)] = 1 ;
		Mat.col[projette(i,2)] = 2 ;
		Mat.col[projette(i,3)] = 3 ;
		Mat.col[projette(i,4)] = 6 ;
		Mat.col[projette(i,5)] = 7 ;
		Mat.col[projette(i,6)] = 8 ;
		Mat.col[projette(i,7)] = 10 ;
		Mat.col[projette(i,8)] = 11 ;
		Mat.col[projette(i,9)] = 12 ;
		++i;
		Mat.col[projette(i,0)] = 3 ;
		Mat.col[projette(i,1)] = 4 ;
		Mat.col[projette(i,2)] = 5 ;
		Mat.col[projette(i,3)] = 6 ;
		++i;
		++i;
		Mat.col[projette(i,0)] = 12 ;
		++i;
		++i;
		Mat.col[projette(i,0)] = 11 ;

#if 0
		for (size_t ii = 0 ; ii < Mat.m/4 ; ++ ii) {
			for (size_t jj = 0 ; jj < ld ; ++jj) {
				for (size_t kk = 0 ; kk < simd::vect_size ; ++kk) {
					std::cout << Mat.dat[ii*ld*simd::vect_size+jj*simd::vect_size+kk] << ' ' ;
				}
				std::cout << '|' ;
			}
			std::cout << std::endl;
		}
			std::cout << std::endl;
#endif

	}

	VECT<double> x,y ;
	x.m = Mat.n ;
	x.inc = 1 ;
	x.dat = FFLAS::fflas_new<double>(x.m);

	y.m = Mat.m ;
	y.inc = 1 ;
	y.dat = FFLAS::fflas_new<double>(y.m);

	for (size_t i = 0 ; i < x.m ; ++i) {
		F.init(x.dat[i],i+1) ;
	}

	for (size_t i = 0 ; i < y.m ; ++i) {
		F.init(y.dat[i] ,i+1) ;
	}

	// y
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax + 0 y
	sp_fgemv(F,Mat,x,0,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;


	// y = Ax + y ( y = 2 Ax)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax - y (y = - Ax)
	sp_fgemv(F,Mat,x,-1,y);

	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + 3 y (y = - Ax)
	sp_fgemv(F,Mat,x,2,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + y (y = 0)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;



	return 0 ;


}

int test2_ell_simd()
{
	Modular<double> F(1051);
	using simd = Simd<double>;
	// using vect_t = typename simd::vect_t;

	ELL_sub<double,true> Mat ;
	// size_t nbnz = 20 ;
	Mat.m = 10 ;
	Mat.n = 13 ;
	size_t ld = 10 ;

	{

		Mat.ld = ld ;
		Mat.col = FFLAS::fflas_new<index_t>(ld*Mat.m);
		Mat.dat = FFLAS::fflas_new<double>(ld*Mat.m);
		for (size_t i = 0 ; i < ld*Mat.m ; ++i) Mat.col[i]  = 0 ;
		for (size_t i = 0 ; i < ld*Mat.m ; ++i) Mat.dat[i]  = 0 ;

		size_t i = 0 ;
		++i ;
		Mat.dat[projette(i,0)] = 1 ;
		Mat.dat[projette(i,1)] = 2 ;
		Mat.dat[projette(i,2)] = 3 ;
		++i;
		++i;
		Mat.dat[projette(i,0)] = -1 ;
		Mat.dat[projette(i,1)] = 10 ;
		++i;
		Mat.dat[projette(i,0)] = 1 ;
		Mat.dat[projette(i,1)] = 2 ;
		Mat.dat[projette(i,2)] = 3 ;
		Mat.dat[projette(i,3)] = 4 ;
		Mat.dat[projette(i,4)] = 7 ;
		Mat.dat[projette(i,5)] = 8 ;
		Mat.dat[projette(i,6)] = 9 ;
		Mat.dat[projette(i,7)] = 11 ;
		Mat.dat[projette(i,8)] = 12 ;
		Mat.dat[projette(i,9)] = 13 ;
		++i;
		Mat.dat[projette(i,0)] = -2 ;
		Mat.dat[projette(i,1)] = -3 ;
		Mat.dat[projette(i,2)] = 4 ;
		Mat.dat[projette(i,3)] = 5 ;
		++i;
		++i;
		Mat.dat[projette(i,0)] = 1 ;
		++i;
		++i;
		Mat.dat[projette(i,0)] = 8 ;

		i = 0 ;
		++i ;
		Mat.col[projette(i,0)] = 1 ;
		Mat.col[projette(i,1)] = 4 ;
		Mat.col[projette(i,2)] = 8 ;
		++i;
		++i;
		Mat.col[projette(i,0)] = 2 ;
		Mat.col[projette(i,1)] = 6 ;
		++i;
		Mat.col[projette(i,0)] = 0 ;
		Mat.col[projette(i,1)] = 1 ;
		Mat.col[projette(i,2)] = 2 ;
		Mat.col[projette(i,3)] = 3 ;
		Mat.col[projette(i,4)] = 6 ;
		Mat.col[projette(i,5)] = 7 ;
		Mat.col[projette(i,6)] = 8 ;
		Mat.col[projette(i,7)] = 10 ;
		Mat.col[projette(i,8)] = 11 ;
		Mat.col[projette(i,9)] = 12 ;
		++i;
		Mat.col[projette(i,0)] = 3 ;
		Mat.col[projette(i,1)] = 4 ;
		Mat.col[projette(i,2)] = 5 ;
		Mat.col[projette(i,3)] = 6 ;
		++i;
		++i;
		Mat.col[projette(i,0)] = 12 ;
		++i;
		++i;
		Mat.col[projette(i,0)] = 11 ;

#if 0
		for (size_t ii = 0 ; ii < Mat.m/4 ; ++ ii) {
			for (size_t jj = 0 ; jj < ld ; ++jj) {
				for (size_t kk = 0 ; kk < simd::vect_size ; ++kk) {
					std::cout << Mat.dat[ii*ld*simd::vect_size+jj*simd::vect_size+kk] << ' ' ;
				}
				std::cout << '|' ;
			}
			std::cout << std::endl;
		}
			std::cout << std::endl;
#endif

	}

	VECT<double> x,y ;
	x.m = Mat.n ;
	x.inc = 1 ;
	x.dat = FFLAS::fflas_new<double>(x.m);

	y.m = Mat.m ;
	y.inc = 1 ;
	y.dat = FFLAS::fflas_new<double>(y.m);

	for (size_t i = 0 ; i < x.m ; ++i) {
		F.init(x.dat[i],i+1) ;
	}

	for (size_t i = 0 ; i < y.m ; ++i) {
		F.init(y.dat[i] ,i+1) ;
	}

	// y
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax + 0 y
	sp_fgemv(F,Mat,x,0,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;


	// y = Ax + y ( y = 2 Ax)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax - y (y = - Ax)
	sp_fgemv(F,Mat,x,-1,y);

	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + 3 y (y = - Ax)
	sp_fgemv(F,Mat,x,2,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + y (y = 0)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;



	return 0 ;


}

#endif

int test1_ellr()
{
	Modular<double> F(1051);

	ELLR<double> Mat ;
	// size_t nbnz = 20 ;
	Mat.m = 10 ;
	Mat.n = 13 ;
	size_t ld = 10 ;
	Mat.ld = ld ;
	Mat.col = FFLAS::fflas_new<index_t>(ld*Mat.m);
	Mat.row = FFLAS::fflas_new<index_t>(Mat.m);
	Mat.dat = FFLAS::fflas_new<double >(ld*Mat.m);
	for (size_t i = 0 ; i < ld*Mat.m ; ++i) Mat.col[i]  = 0 ;
	for (size_t i = 0 ; i < ld*Mat.m ; ++i) Mat.dat[i]  = 0 ;
	for (size_t i = 0 ; i < Mat.m ; ++i) Mat.row[i]  = 0 ;

	size_t j = 0 ;
	++j ;
	Mat.dat[j*ld+0] = 1 ;
	Mat.dat[j*ld+1] = 2 ;
	Mat.dat[j*ld+2] = 3 ;
	++j;
	++j;
	Mat.dat[j*ld+0] = -1 ;
	Mat.dat[j*ld+1] = 10 ;
	++j;
	Mat.dat[j*ld+0] = 1 ;
	Mat.dat[j*ld+1] = 2 ;
	Mat.dat[j*ld+2] = 3 ;
	Mat.dat[j*ld+3] = 4 ;
	Mat.dat[j*ld+4] = 7 ;
	Mat.dat[j*ld+5] = 8 ;
	Mat.dat[j*ld+6] = 9 ;
	Mat.dat[j*ld+7] = 11 ;
	Mat.dat[j*ld+8] = 12 ;
	Mat.dat[j*ld+9] = 13 ;
	++j;
	Mat.dat[j*ld+0] = -2 ;
	Mat.dat[j*ld+1] = -3 ;
	Mat.dat[j*ld+2] = 4 ;
	Mat.dat[j*ld+3] = 5 ;
	++j;
	++j;
	Mat.dat[j*ld+0] = 1 ;
	++j;
	++j;
	Mat.dat[j*ld+0] = 8 ;

	j = 0 ;
	++j ;
	Mat.col[j*ld+0] = 1 ;
	Mat.col[j*ld+1] = 4 ;
	Mat.col[j*ld+2] = 8 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 2 ;
	Mat.col[j*ld+1] = 6 ;
	++j;
	Mat.col[j*ld+0] = 0 ;
	Mat.col[j*ld+1] = 1 ;
	Mat.col[j*ld+2] = 2 ;
	Mat.col[j*ld+3] = 3 ;
	Mat.col[j*ld+4] = 6 ;
	Mat.col[j*ld+5] = 7 ;
	Mat.col[j*ld+6] = 8 ;
	Mat.col[j*ld+7] = 10 ;
	Mat.col[j*ld+8] = 11 ;
	Mat.col[j*ld+9] = 12 ;
	++j;
	Mat.col[j*ld+0] = 3 ;
	Mat.col[j*ld+1] = 4 ;
	Mat.col[j*ld+2] = 5 ;
	Mat.col[j*ld+3] = 6 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 12 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 11 ;

	Mat.row[0] = 0 ;
	Mat.row[1] = 3 ;
	Mat.row[2] = 0 ;
	Mat.row[3] = 2 ;
	Mat.row[4] = 10 ;
	Mat.row[5] = 4 ;
	Mat.row[6] = 0 ;
	Mat.row[7] = 1 ;
	Mat.row[8] = 0 ;
	Mat.row[9] = 1 ;


	VECT<double> x,y ;
	x.m = Mat.n ;
	x.inc = 1 ;
	x.dat = FFLAS::fflas_new<double>(x.m);

	y.m = Mat.m ;
	y.inc = 1 ;
	y.dat = FFLAS::fflas_new<double>(y.m);

	for (size_t i = 0 ; i < x.m ; ++i) {
		F.init(x.dat[i],i+1) ;
	}

	for (size_t i = 0 ; i < y.m ; ++i) {
		F.init(y.dat[i] ,i+1) ;
	}

	// y
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax + 0 y
	sp_fgemv(F,Mat,x,0,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;


	// y = Ax + y ( y = 2 Ax)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax - y (y = - Ax)
	sp_fgemv(F,Mat,x,-1,y);

	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + 3 y (y = - Ax)
	sp_fgemv(F,Mat,x,2,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + y (y = 0)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;



	return 0 ;


}

int test2_ellr()
{
	Modular<double> F(1051);

	ELLR_sub<double> Mat ;
	// size_t nbnz = 20 ;
	Mat.m = 10 ;
	Mat.n = 13 ;
	size_t ld = 10 ;
	Mat.ld = ld ;
	Mat.col = FFLAS::fflas_new<index_t>(ld*Mat.m);
	Mat.dat = FFLAS::fflas_new<double >(ld*Mat.m);
	Mat.row = FFLAS::fflas_new<index_t>(Mat.m);
	for (size_t i = 0 ; i < ld*Mat.m ; ++i) Mat.col[i]  = 0 ;
	for (size_t i = 0 ; i < ld*Mat.m ; ++i) Mat.dat[i]  = 0 ;
	for (size_t i = 0 ; i < Mat.m ; ++i) Mat.row[i]  = 0 ;

	size_t j = 0 ;
	++j ;
	Mat.dat[j*ld+0] = 1 ;
	Mat.dat[j*ld+1] = 2 ;
	Mat.dat[j*ld+2] = 3 ;
	++j;
	++j;
	Mat.dat[j*ld+0] = -1 ;
	Mat.dat[j*ld+1] = 10 ;
	++j;
	Mat.dat[j*ld+0] = 1 ;
	Mat.dat[j*ld+1] = 2 ;
	Mat.dat[j*ld+2] = 3 ;
	Mat.dat[j*ld+3] = 4 ;
	Mat.dat[j*ld+4] = 7 ;
	Mat.dat[j*ld+5] = 8 ;
	Mat.dat[j*ld+6] = 9 ;
	Mat.dat[j*ld+7] = 11 ;
	Mat.dat[j*ld+8] = 12 ;
	Mat.dat[j*ld+9] = 13 ;
	++j;
	Mat.dat[j*ld+0] = -2 ;
	Mat.dat[j*ld+1] = -3 ;
	Mat.dat[j*ld+2] = 4 ;
	Mat.dat[j*ld+3] = 5 ;
	++j;
	++j;
	Mat.dat[j*ld+0] = 1 ;
	++j;
	++j;
	Mat.dat[j*ld+0] = 8 ;

	j = 0 ;
	++j ;
	Mat.col[j*ld+0] = 1 ;
	Mat.col[j*ld+1] = 4 ;
	Mat.col[j*ld+2] = 8 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 2 ;
	Mat.col[j*ld+1] = 6 ;
	++j;
	Mat.col[j*ld+0] = 0 ;
	Mat.col[j*ld+1] = 1 ;
	Mat.col[j*ld+2] = 2 ;
	Mat.col[j*ld+3] = 3 ;
	Mat.col[j*ld+4] = 6 ;
	Mat.col[j*ld+5] = 7 ;
	Mat.col[j*ld+6] = 8 ;
	Mat.col[j*ld+7] = 10 ;
	Mat.col[j*ld+8] = 11 ;
	Mat.col[j*ld+9] = 12 ;
	++j;
	Mat.col[j*ld+0] = 3 ;
	Mat.col[j*ld+1] = 4 ;
	Mat.col[j*ld+2] = 5 ;
	Mat.col[j*ld+3] = 6 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 12 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 11 ;

	Mat.row[0] = 0 ;
	Mat.row[1] = 3 ;
	Mat.row[2] = 0 ;
	Mat.row[3] = 2 ;
	Mat.row[4] = 10 ;
	Mat.row[5] = 4 ;
	Mat.row[6] = 0 ;
	Mat.row[7] = 1 ;
	Mat.row[8] = 0 ;
	Mat.row[9] = 1 ;




	VECT<double> x,y ;
	x.m = Mat.n ;
	x.inc = 1 ;
	x.dat = FFLAS::fflas_new<double>(x.m);

	y.m = Mat.m ;
	y.inc = 1 ;
	y.dat = FFLAS::fflas_new<double>(y.m);

	for (size_t i = 0 ; i < x.m ; ++i) {
		F.init(x.dat[i],i+1) ;
	}

	for (size_t i = 0 ; i < y.m ; ++i) {
		F.init(y.dat[i] ,i+1) ;
	}

	// y
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax + 0 y
	sp_fgemv(F,Mat,x,0,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;


	// y = Ax + y ( y = 2 Ax)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y = Ax - y (y = - Ax)
	sp_fgemv(F,Mat,x,-1,y);

	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + 3 y (y = - Ax)
	sp_fgemv(F,Mat,x,2,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;

	// y += Ax + y (y = 0)
	sp_fgemv(F,Mat,x,1,y);
	for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
	std::cout << std::endl;



	return 0 ;


}

int test3_ellr(int CST)
{

	{
		Modular<double> F(1051);

		ELLR<double> Mat ;
	Mat.m = 10 ;
	Mat.n = 13 ;
	size_t ld = 10 ;
	Mat.ld = ld ;
	Mat.col = FFLAS::fflas_new<index_t>(ld*Mat.m);
	Mat.dat = FFLAS::fflas_new<double >(ld*Mat.m);
	Mat.row = FFLAS::fflas_new<index_t>(Mat.m);
	for (size_t i = 0 ; i < ld*Mat.m ; ++i) Mat.col[i]  = 0 ;
	for (size_t i = 0 ; i < ld*Mat.m ; ++i) Mat.dat[i]  = 0 ;
	for (size_t i = 0 ; i < Mat.m ; ++i) Mat.row[i]  = 0 ;

	size_t j = 0 ;
	++j ;
	Mat.dat[j*ld+0] = CST ;
	Mat.dat[j*ld+1] = CST ;
	Mat.dat[j*ld+2] = CST ;
	++j;
	++j;
	Mat.dat[j*ld+0] = CST ;
	Mat.dat[j*ld+1] = CST ;
	++j;
	Mat.dat[j*ld+0] = CST ;
	Mat.dat[j*ld+1] = CST ;
	Mat.dat[j*ld+2] = CST ;
	Mat.dat[j*ld+3] = CST ;
	Mat.dat[j*ld+4] = CST ;
	Mat.dat[j*ld+5] = CST ;
	Mat.dat[j*ld+6] = CST ;
	Mat.dat[j*ld+7] = CST ;
	Mat.dat[j*ld+8] = CST ;
	Mat.dat[j*ld+9] = CST ;
	++j;
	Mat.dat[j*ld+0] = CST ;
	Mat.dat[j*ld+1] = CST ;
	Mat.dat[j*ld+2] = CST ;
	Mat.dat[j*ld+3] = CST ;
	++j;
	++j;
	Mat.dat[j*ld+0] = CST ;
	++j;
	++j;
	Mat.dat[j*ld+0] = CST ;

	j = 0 ;
	++j ;
	Mat.col[j*ld+0] = 1 ;
	Mat.col[j*ld+1] = 4 ;
	Mat.col[j*ld+2] = 8 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 2 ;
	Mat.col[j*ld+1] = 6 ;
	++j;
	Mat.col[j*ld+0] = 0 ;
	Mat.col[j*ld+1] = 1 ;
	Mat.col[j*ld+2] = 2 ;
	Mat.col[j*ld+3] = 3 ;
	Mat.col[j*ld+4] = 6 ;
	Mat.col[j*ld+5] = 7 ;
	Mat.col[j*ld+6] = 8 ;
	Mat.col[j*ld+7] = 10 ;
	Mat.col[j*ld+8] = 11 ;
	Mat.col[j*ld+9] = 12 ;
	++j;
	Mat.col[j*ld+0] = 3 ;
	Mat.col[j*ld+1] = 4 ;
	Mat.col[j*ld+2] = 5 ;
	Mat.col[j*ld+3] = 6 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 12 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 11 ;

	Mat.row[0] = 0 ;
	Mat.row[1] = 3 ;
	Mat.row[2] = 0 ;
	Mat.row[3] = 2 ;
	Mat.row[4] = 10 ;
	Mat.row[5] = 4 ;
	Mat.row[6] = 0 ;
	Mat.row[7] = 1 ;
	Mat.row[8] = 0 ;
	Mat.row[9] = 1 ;



		// Mat.z =  nbnz;


		VECT<double> x,y ;
		x.m = Mat.n ;
		x.inc = 1 ;
		x.dat = FFLAS::fflas_new<double>(x.m);

		y.m = Mat.m ;
		y.inc = 1 ;
		y.dat = FFLAS::fflas_new<double>(y.m);

		for (size_t i = 0 ; i < x.m ; ++i) {
			F.init(x.dat[i],i+1) ;
		}

		for (size_t i = 0 ; i < y.m ; ++i) {
			F.init(y.dat[i] ,i+1) ;
		}

		// y
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;

		// y = Ax + 0 y
		sp_fgemv(F,Mat,x,0,y);
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;


		// y = Ax + y ( y = 2 Ax)
		// sp_fgemv(F,Mat,x,1,y);
		// for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		// std::cout << std::endl;

		// y = Ax - y (y = - Ax)
		// sp_fgemv(F,Mat,x,-1,y);

		// for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		// std::cout << std::endl;

		// y += Ax + 3 y (y = - Ax)
		// sp_fgemv(F,Mat,x,2,y);
		// for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		// std::cout << std::endl;

		// y += Ax + y (y = 0)
		// sp_fgemv(F,Mat,x,1,y);
		// for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		// std::cout << std::endl;


	}

	{
		Modular<double> F(1051);
		ELLR_ZO<double> Mat ;

	Mat.m = 10 ;
	Mat.n = 13 ;
	size_t ld = 10 ;
	Mat.ld = ld ;
	Mat.col = FFLAS::fflas_new<index_t>(ld*Mat.m);
	Mat.row = FFLAS::fflas_new<index_t>(Mat.m);
	for (size_t i = 0 ; i < ld*Mat.m ; ++i) Mat.col[i]  = 0 ;
	for (size_t i = 0 ; i < Mat.m ; ++i) Mat.row[i]  = 0 ;

	size_t j = 0 ;

	++j ;
	Mat.col[j*ld+0] = 1 ;
	Mat.col[j*ld+1] = 4 ;
	Mat.col[j*ld+2] = 8 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 2 ;
	Mat.col[j*ld+1] = 6 ;
	++j;
	Mat.col[j*ld+0] = 0 ;
	Mat.col[j*ld+1] = 1 ;
	Mat.col[j*ld+2] = 2 ;
	Mat.col[j*ld+3] = 3 ;
	Mat.col[j*ld+4] = 6 ;
	Mat.col[j*ld+5] = 7 ;
	Mat.col[j*ld+6] = 8 ;
	Mat.col[j*ld+7] = 10 ;
	Mat.col[j*ld+8] = 11 ;
	Mat.col[j*ld+9] = 12 ;
	++j;
	Mat.col[j*ld+0] = 3 ;
	Mat.col[j*ld+1] = 4 ;
	Mat.col[j*ld+2] = 5 ;
	Mat.col[j*ld+3] = 6 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 12 ;
	++j;
	++j;
	Mat.col[j*ld+0] = 11 ;

	Mat.row[0] = 0 ;
	Mat.row[1] = 3 ;
	Mat.row[2] = 0 ;
	Mat.row[3] = 2 ;
	Mat.row[4] = 10 ;
	Mat.row[5] = 4 ;
	Mat.row[6] = 0 ;
	Mat.row[7] = 1 ;
	Mat.row[8] = 0 ;
	Mat.row[9] = 1 ;

	Mat.cst = CST ;



		// Mat.z =  nbnz;


		VECT<double> x,y ;
		x.m = Mat.n ;
		x.inc = 1 ;
		x.dat = FFLAS::fflas_new<double>(x.m);

		y.m = Mat.m ;
		y.inc = 1 ;
		y.dat = FFLAS::fflas_new<double>(y.m);

		for (size_t i = 0 ; i < x.m ; ++i) {
			F.init(x.dat[i],i+1) ;
		}

		for (size_t i = 0 ; i < y.m ; ++i) {
			F.init(y.dat[i] ,i+1) ;
		}

		// y
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;

		// y = Ax + 0 y
		sp_fgemv(F,Mat,x,0,y);
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;


		// y = Ax + y ( y = 2 Ax)
		sp_fgemv(F,Mat,x,1,y);
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;

		// y = Ax - y (y = - Ax)
		sp_fgemv(F,Mat,x,-1,y);

		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;

		// y += Ax + 3 y (y = - Ax)
		sp_fgemv(F,Mat,x,2,y);
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;

		// y += Ax + y (y = 0)
		sp_fgemv(F,Mat,x,1,y);
		for (size_t i = 0 ; i < y.m ; ++i) std::cout << y.dat[i] << ' '  ;
		std::cout << std::endl;




	}

	return 0 ;

}



int main()
{
	std::cout << "test (COO)         " << std::endl;
	test1_coo();

	std::cout << "test (COO) SUB     " << std::endl;
	test2_coo();

	std::cout << "test (COO) ZO      " << std::endl;
	test3_coo(1);
	test3_coo(2);

	std::cout << "test (CSR)         " << std::endl;
	test1_csr();

	std::cout << "test (CSR) SUB     " << std::endl;
	test2_csr();

	std::cout << "test (CSR) ZO      " << std::endl;
	test3_csr(1);
	test3_csr(2);

	std::cout << "test (ELL)         " << std::endl;
	test1_ell();

	std::cout << "test (ELL) SUB     " << std::endl;
	test2_ell();

#ifdef __FFLASFFPACK_USE_SIMD
	std::cout << "test (ELL)     SIMD" << std::endl;
	test1_ell_simd();

	std::cout << "test (ELL) SUB SIMD" << std::endl;
	test2_ell_simd();

#endif

	std::cout << "test (ELLR)        " << std::endl;
	test1_ellr();

	std::cout << "test (ELLR) SUB    " << std::endl;
	test2_ellr();

	std::cout << "test (ELLR) ZO     " << std::endl;
	test3_ellr(1);
	test3_ellr(2);




}
