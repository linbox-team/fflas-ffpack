/*
 * Copyright (C) FFLAS-FFPACK
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


//--------------------------------------------------------------------------
//          Test for the krylov-elimination
//--------------------------------------------------------------------------
// usage: test-krylov-elim p A, to compute the rank profile of the (n+m)xn matrix B
// formed by the n identity vectors and the mxn matrix A over Z/pZ
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <list>
#include <vector>
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/timer.h"
using namespace std;
#include "givaro/modular.h"
#include "fflas-ffpack/ffpack/ffpack.h"


using namespace FFPACK;
typedef Givaro::Modular<double> Field;

template<class T>
std::ostream& printvect(std::ostream& o, vector<T>& vect){
    for(size_t i=0; i < vect.size()-1; ++i)
        o << vect[i] << " " ;
    return o << vect[vect.size()-1] << std::endl;
}

int main(int argc, char** argv){

    size_t m,n;
    cout<<setprecision(20);

    if (argc!=4){
        cerr<<"usage : test-frobenius <p> <A> <c>"<<endl
        <<"         to compute the frobenius normal form of the matrix A over Z/pZ, with conditonning parameter c"
        <<endl;
        exit(-1);
    }
    Field F( atoi(argv[1]) );
    Field::Element* A;
    FFLAS::ReadMatrix (argv[2],F,m,n,A);
    size_t c = atoi(argv[3]);

    std::list<vector<Field::Element> > frobForm;
    FFLAS::Timer tim;
    tim.clear();
    tim.start();
    FFPACK::CharpolyArithProg (F, frobForm, n, A, n, c);
    tim.stop();
    std::list<vector<Field::Element> >::iterator it = frobForm.begin();
    while(it != frobForm.end()){
        printvect (cout, *(it++));
    }
    cerr<<c<<" "<<tim.usertime()<<" "<<4.55*n*n/1000000.0*n/tim.usertime()<<endl;
    FFLAS::fflas_delete( A);
    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
