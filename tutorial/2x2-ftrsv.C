/* Copyright (c) 2017 FFLAS-FFPACK
* Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
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


#include <fflas-ffpack/fflas-ffpack-config.h>
#include <givaro/modular-balanced.h>
#include <fflas-ffpack/fflas/fflas.h>
#include <fflas-ffpack/utils/timer.h>
#include <fflas-ffpack/utils/Matio.h>
#include <fflas-ffpack/utils/args-parser.h>

#include <iostream>
#include <array>

using namespace FFLAS;

int main(int argc, char** argv) {

    typedef Givaro::Modular<float> Ring;
    Ring F(11);
    
    Ring::Element L[4]{1,0,2,3};
    Ring::Element B[2]{4,5};
    
    size_t m(2);
    
    write_field(F, std::cout << "L:=", L, m, m, m, true) << std::endl;
    write_field(F, std::cout << "B:=", B, m, 1, 1, true) << std::endl;
    
        // In place system solve
    ftrsv (F, FflasLower,FflasNoTrans,FflasNonUnit, m, L, m, B, 1);
    
    write_field(F, std::cout << "X:=", B, m, 1, 1, true) << std::endl;
    std::cerr << "0 = L.X - B mod " << F.characteristic() << ';' << std::endl;
    
    Ring::Element U[4]{3,2,0,5}, C[2]{4,7};
    
    
    write_field(F, std::cout << "U:=", U, m, m, m, true) << std::endl;
    write_field(F, std::cout << "C:=", C, m, 1, 1, true) << std::endl;
    
        // In place system solve
    ftrsv (F, FflasUpper, FflasNoTrans ,FflasNonUnit, m, U, m, B, 1);
    
    write_field(F, std::cout << "X:=", C, m, 1, 1, true) << std::endl;
    std::cerr << "0 = U.X - C mod " << F.characteristic() << ';' << std::endl;
    
    Ring::Element D[2]{4,5};
    
    write_field(F, std::cout << "L:=", L, m, m, m, true) << std::endl;
    write_field(F, std::cout << "D:=", D, m, 1, 1, true) << std::endl;
    
        // In place system solve
    ftrsv (F, FflasLower,FflasTrans,FflasNonUnit, m, L, m, D, 1);
    
    write_field(F, std::cout << "X:=", D, m, 1, 1, true) << std::endl;
    std::cerr << "0 = Transpose(X).L - Transpose(D) mod " << F.characteristic() << ';' << std::endl;

    Ring::Element E[2]{4,7};
    
    
    write_field(F, std::cout << "U:=", U, m, m, m, true) << std::endl;
    write_field(F, std::cout << "E:=", E, m, 1, 1, true) << std::endl;
    
        // In place system solve
    ftrsv (F, FflasUpper, FflasTrans ,FflasNonUnit, m, U, m, E, 1);
    
    write_field(F, std::cout << "X:=", E, m, 1, 1, true) << std::endl;
    std::cerr << "0 = Transpose(X).U - Transpose(E) mod " << F.characteristic() << ';' << std::endl;
    
    return 0;
}

