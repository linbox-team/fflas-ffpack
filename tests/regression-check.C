/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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
 
#include "fflas-ffpack/fflas-ffpack-config.h"
#include <givaro/modular.h>
#include "fflas-ffpack/fflas-ffpack.h"
using namespace Givaro;
using namespace FFLAS;
using namespace FFPACK;
/*  #1  */
bool check1 () ;

/*  #2  */
bool check2()
{
	Modular<double> F(2);
	Modular<double>::RandIter R(F);

	size_t ok = 0 ;
	size_t tot = 500 ;
	for (size_t i = 0 ; i < tot ; ++i) {
		double elt ;
		R.random(elt);
		if (elt == 1) ++ok ;
	}
	double f = (double) ok / (double) tot ;
	if (f < 0.3 or f > 0.7) return false ;

	return true;

}

/*  #3  */
bool check3()
{
	Modular<double> F(2);
	double * A = NULL ;
	double d = Det(F,0,0,A,0);
	return F.areEqual(d,F.one);

}

/*  #4  */
bool check4()
{
    typedef int32_t Element;
	Modular<Element> F(2);
	Element * A = NULL ;
	Element * X = NULL ;
	int nul;
	Invert2(F,0,A,0,X,0,nul);
	return true ;
}


bool checkZeroDimCharpoly(){
	Modular<double> F(101);
	double * A = fflas_new(F,0,0);
	Poly1Dom<Modular<double> > PR (F);
	Poly1Dom<Modular<double> >::Element charp;
	CharPoly(PR, charp, 0, A, 0);
	return PR.isOne(charp);
}
bool checkZeroDimMinPoly(){
	Modular<double> F(101);
	double * A = fflas_new(F,0,0);
	Poly1Dom<Modular<double> > PR (F);
	Poly1Dom<Modular<double> >::Element minp;
	MinPoly(F, minp, 0, A, 0);
	return PR.isOne(minp);
}
int main() {
	bool pass = true ;
	pass = pass && check2();
	pass = pass && check3();
	pass = pass && check4();
	pass = pass && checkZeroDimCharpoly();
	pass = pass && checkZeroDimMinPoly();
	return !pass;
}

