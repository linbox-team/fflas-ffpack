/* Copyright (c) FFLAS-FFPACK
* Written by  

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
*/


/*    ===============================================
      |    Finit Field Linar Algebra Subprograms    |
      ===============================================
The puspose of this file is to be able to handle the basics for
using fflas.
 */

// ==========Includes=========

// fflas-ffpack modules
//#include <fflas-ffpack/fflas-ffpack-config.h>
#include <fflas-ffpack/fflas/fflas.h>

// Givaro provides some finite fields such as modular rings.
// Here we include two types of modular rings.
//    1. Normal   modular rings with elements numbered from 0 to p
//    2. Balanced modular rings with elemetns numbered form -(p+1)/2 to (p+1)/2
// Typically p is a prime number but Givaro works with any
// interger n as well.
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>



#include <iostream>

// All of the functions used come from the FFLAS namespace.
using namespace FFLAS;

int main(int argc, char** argv) {
  
  // ========== Defining finite fields  =========

  // Defining rings using Givaro.
  typedef Givaro::Modular<float> Float_Ring;
  Float_Ring F1(11); 
  // F1 is a ring of floats defined modulo 11.
  
  typedef Givaro::Modular<double> Double_Balanced_Ring;
  Double_Balanced_Ring F2(15);
  // F2 is a ring of doubles defined modulo 15;
  
  // For other ways of defining fields using Givaro please refer to Givaro's documentation.

  // ===== Elements in finite fields ======
  
  // Let a, b and c be three elements of F1
  // and d, e and f be three elements of F2.

  Float_Ring::Element a,b;
  Double_Balanced_Ring::Element d,e,f;
  
  // The elements need to be assigned.
  
  F1.init(a,2);
  F1.init(b,25);
  // a and b are initialised with the elements of F1 that
  // correspond to 2 and 25.
  // Since F1 is made of floats modulo 11, these are 2.0 and 3.0.
  
  // This means the second parameter will be converted into an
  // actual element of the field.

  // A field has two elements that are neutral with respect to
  // the addition and multiplication laws, ie : zero and one.
  // An element of the field can be assigned to an element variable.
  
  F2.assign(d, F2.zero);
  F2.assign(e, F2.one);
  F2.assign(f, e);
  
  // Operations on elements
  Float_Ring::Element c;
  Double_Balanced_Ring::Element g;

  F1.add(c,a,b);       // c = a + b
  F2.mul(g,e,5);       // g = e * 5
  
} //end main
