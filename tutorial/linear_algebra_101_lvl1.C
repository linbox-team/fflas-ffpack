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



/*    =====================================================
      |    Finit Field Linar Algebra Subprograms Level 1  |
      =====================================================
The puspose of this file is to be able to handle the first level of fflas.
Here 
 */

// For information on the inclusions check the file 101-fflas.C
#include <fflas-ffpack/fflas/fflas.h>
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include <iostream>

using namespace FFLAS;

int main(int argc, char** argv) {
  std::cout << "" << std::endl;
  
  typedef Givaro::Modular<float> Float_Ring;
  Float_Ring F(120);
  Float_Ring::Element alpha,beta;    // scalars
  Float_Ring::Element * X, * Y, * Z; // vectors
  
  X = fflas_new(F,2,1);
  Y = fflas_new(F,2,1);
  Z = fflas_new(F,2,1);

  // ===== manual initialisation =====
  
  F.init(*(X+0),  52);               // X[0] = 52
  F.init(  X[1], 183);               // X[1] = 63

  std::cout << "\nInitialisation" << std::endl;
  write_field(F, std::cout << "X:=", X, 2, 1, 1,true) << std::endl;
  write_field(F, std::cout << "Y:=", Y, 2, 1, 1,true) << std::endl;
  

  // ===== initialisation from a vector ===== ???
  /*  finit(F,2,Y,1,X,1);

  write_field(F, std::cout << "X:=", X, 2, 1, 1,true) << std::endl;
  write_field(F, std::cout << "Y:=", Y, 2, 1, 1,true) << std::endl;
  */

  
  

  // ===== modular reduction =====
  std::cout << "\n=== Modular reductions ===" << std::endl;
  
  typedef Givaro::Modular<float> Ring;
  Ring F2(2);
  freduce(F2,2,Y,1);      // reduce Y modulo F2 inplace
  freduce(F2,2,X,1,Y,1);  // reduce X modulo F2 and store in Y

  
  std::cout << "\nY mod 2 inplace" << std::endl;
  write_field(F, std::cout << "Y:=", Y, 2, 1, 1,true) << std::endl;
  std::cout << "\nY := X mod 2" << std::endl;
  write_field(F, std::cout << "Y:=", Y, 2, 1, 1,true) << std::endl;


  // ===== negation =====
  fnegin(F,2,X,1);        // X := -X
  fneg(F,2,X,1,Y,1);      // Y := -X

  std::cout << "\n=== Negations ===" << std::endl;
  std::cout << "\nX := -X" << std::endl;
  write_field(F, std::cout << "X:=", X, 2, 1, 1,true) << std::endl;
  std::cout << "\nY := -X" << std::endl;
  write_field(F, std::cout << "Y:=", Y, 2, 1, 1,true) << std::endl;
  

  // ===== addition / substraction =====
  // both of these exist inplace with the suffix "in"
  F.init(X[0],3); F.init(X[1],4);       // X := (3;4)
  F.init(Y[0],2); F.assign(Y[1],F.one); // Y := (2;1)

  std::cout << "\n=== Additions / Substractions ===" << std::endl;
  write_field(F, std::cout << "X:=", X, 2, 1, 1,true) << std::endl;
  write_field(F, std::cout << "Y:=", Y, 2, 1, 1,true) << std::endl;
  write_field(F, std::cout << "Z:=", Z, 2, 1, 1,true) << std::endl;

  fadd(F,2,X,1,Y,1,Z,1);  // Z := X + Y
  fsub(F,2,Z,1,Y,1,X,1);  // X := Z - Y

  std::cout << "\nZ := X + Y"  << std::endl;
  write_field(F, std::cout << "Z:=", Z, 2, 1, 1,true) << std::endl;
  std::cout << "\nX := Z - Y"  << std::endl;
  write_field(F, std::cout << "X:=", X, 2, 1, 1,true) << std::endl;

  // ===== scalar times vector =====
  F.init(alpha,42); // the scalar is scalar with respect to F
  fscal(F,2,alpha,X,1,Y,1); // Y = alpha * X

  std::cout << "\n=== Scalar times vector ===" << std::endl;
  write_field(F, std::cout << "X:=", X, 2, 1, 1,true) << std::endl;
  write_field(F, std::cout << "Y:=", Y, 2, 1, 1,true) << std::endl;

  
  
  // ===== linar combinations =====
  // Although fadd, fsub, fdot and fscal exist, it's often better
  // to use faxpy due to processor optimisation

  std::cout << "\n=== Linear Combinations ==="  << std::endl;
  write_field(F, std::cout << "X:=", X, 2, 1, 1,true) << std::endl;
  write_field(F, std::cout << "Y:=", Y, 2, 1, 1,true) << std::endl;

  faxpy(F,2,alpha,X,1,Y,1); // Y := alpha * X + Y

  std::cout << "\nY := alpha * X + Y"  << std::endl;
  write_field(F, std::cout << "Y:=", Y, 2, 1, 1,true) << std::endl;

  // faxpby is a more general form of faxpy

  //F.init(beta,14);
  //faxpby(F,2,alpha,X,1,beta,Y,1); // Y = alpha * X + beta * Y

  // ===== dot product (scalar product) =====
  
  std::cout << "\n=== Dot Product ==="  << std::endl;
  write_field(F, std::cout << "Z:=", Z, 2, 1, 1,true) << std::endl;
  write_field(F, std::cout << "X:=", X, 2, 1, 1,true) << std::endl;
  
  fdot(F,2,Z,1,X,1); // X := Z^T dot X

  std::cout << "\nX := Z^T dot X"  << std::endl;
  write_field(F, std::cout << "X:=", X, 2, 1, 1,true) << std::endl;

  // ===== tests =====
  bool b;
  F.init(Z[0], 0); F.init(Z[1], 120);
  
  std::cout << "\n=== Tests ==="  << std::endl;
  write_field(F, std::cout << "Z:=", Z, 2, 1, 1,true) << std::endl;
  write_field(F, std::cout << "X:=", X, 2, 1, 1,true) << std::endl;

  
  b = fiszero(F,2,Z,1); // TRUE if X is zero vector in F

  std::cout << "Z == zero_vector : " << b << std::endl;
  
  b = fequal(F,2,Z,1,X,1); // TRUE if X and Y are equivalent in F

  std::cout << "Z == X : " << b << std::endl;
  
  // ===== swap  =====

  std::cout << "\n=== Swap ==="  << std::endl;
  write_field(F, std::cout << "X:=", X, 2, 1, 1,true) << std::endl;
  write_field(F, std::cout << "Y:=", Y, 2, 1, 1,true) << std::endl;

  fswap(F,2,X,1,Y,1);   // exchanges X and Y;

  std::cout << "\nX <-> Y"  << std::endl;
  write_field(F, std::cout << "X:=", X, 2, 1, 1,true) << std::endl;
  write_field(F, std::cout << "Y:=", Y, 2, 1, 1,true) << std::endl;
  
  
} // End main
