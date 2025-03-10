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
    double d;
    Det(F,d,0,A,0);
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
    fflas_delete(A);
    return PR.isOne(charp);
}
bool checkZeroDimMinPoly(){
    Modular<double> F(101);
    double * A = fflas_new(F,0,0);
    Poly1Dom<Modular<double> > PR (F);
    Poly1Dom<Modular<double> >::Element minp;
    MinPoly(F, minp, 0, A, 0);
    fflas_delete(A);
    return PR.isOne(minp);
}

bool gf2ModularBalanced(){

    typedef Givaro::Modular<int64_t> Field;

    Field F(2);

    int h_[] = {1,0,1,0,1,0,1, 0,1,1,0,0,1,1, 0,0,0,1,1,1,1};

    Field::Element_ptr G, H;
    H = fflas_new(F,3,7);

    for (unsigned i=0; i<21; i++)
        F.init(H[i], h_[i]);

    size_t NSdim;
    size_t ldn;
    NullSpaceBasis (F, FflasLeft, 3, 7, H, 7, G, ldn, NSdim);

    fflas_delete(G,H);
    return true;
}

bool charpolyArithProg(){

    Modular<double> F(37);
    Poly1Dom<Modular<double> >PolRing(F,'X');

    // Reading the matrix from a file
    double* A;
    size_t m, n;
    std::string file("data/regression_charpoly.sms");
    ReadMatrix(file.c_str(), F, m, n, A);
    // here m=n=35
    Poly1Dom<Modular<double> >::Element charp(n+1);

    CharPoly(PolRing, charp, n, A, n);

    fflas_delete(A);
    bool pass = F.isOne(charp[n]);
    for (size_t i = 0; i<n; i++)
            pass = pass && (F.isZero(charp[i]));
    if (!pass) std::cerr<<"charpolyArithProg FAILED"<<std::endl;
    return pass;
}

int main() {
    bool pass = true ;
    pass = pass && check2();
    pass = pass && check3();
    pass = pass && check4();
    pass = pass && checkZeroDimCharpoly();
    pass = pass && checkZeroDimMinPoly();
    pass = pass && gf2ModularBalanced();
    pass = pass && charpolyArithProg();
    return !pass;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
