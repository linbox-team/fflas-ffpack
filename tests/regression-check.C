/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* tests/regression-check.C
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by all reporters of bugs (see ffpack-devel@googlegroups.com)
 *
 * ------------------------------------
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
#include "fflas-ffpack/fflas-ffpack.h"

/*  #1  */
bool check1 () ;

/*  #2  */
bool check2()
{
	FFPACK::Modular<double> F(2);
	FFPACK::Modular<double>::RandIter R(F);

	size_t ok = 0 ;
	size_t tot = 500 ;
	for (size_t i = 0 ; i < tot ; ++i) {
		double elt ;
		R.random(elt);
		if (elt == 1) ++ok ;
	}
	double f = (double) ok / (double) tot ;
	if (f < 0.3 or f > 0.7) return false ;

	return true ;

}

/*  #3  */
bool check3()
{
	FFPACK::Modular<double> F(2);
	double * A = NULL ;
	double d = FFPACK::Det(F,0,0,A,0);
	return F.areEqual(d,F.one);

}

/*  #4  */
bool check4()
{
        typedef int32_t Element;
	FFPACK::Modular<Element> F(2);
	Element * A = NULL ;
	Element * X = NULL ;
	int nul;
	FFPACK::Invert2(F,0,A,0,X,0,nul);
	return true ;
}


int main() {
	bool pass = true ;
	pass &= check2();
	pass &= check3();
	pass &= check4();
	return !pass;
}

