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

using namespace FFLAS;

int main(int argc, char** argv) {

    typedef Givaro::Modular<float> Ring;
    Ring F(13);
    
    Ring::Element L[4]{1,0,2,3}, B[2]{4,5};
    
    size_t m(2);
    
    write_field(F, std::cout << "L:=", L, m, m, m, true) << std::endl;
    write_field(F, std::cout << "B:=", B, m, 1, 1, true) << std::endl;
    
        // In place system solve
    ftrsv (F, FflasLower,FflasNoTrans,FflasNonUnit, m, L, m, B, 1);
    
    write_field(F, std::cout << "X:=", B, m, 1, 1, true) << std::endl;
    std::cerr << "0 = L.X - B mod " << F.characteristic() << ';' << std::endl;
    
    return 0;
}

