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


/*        ===========================================
          |    Linar Algebra Crash Course Level 0   |
          ===========================================
          The puspose of this file is to be able to handle the element level
          of the linear algebra over finite fields used in fflas.
          In this example the finite field will be the Z over 101 Z modular ring.
          This set is created using the givaro library although modular rings are
          not the only possibilities.
          There had to be a choice between having a clean file and a clean display.
          The clean display has been adopted therefore the reader is invited
          to pay attention to what lines are used for display and what lines
          actually do the work.
          */

// =========Information=======
// This first tutorial has a lot of commentaries to enable
// the user to understand how operations on basic elements
// are performed.

// Although the elements are usually manipulated in the machine
// as floats and doubles,the display is done with integers.

// The different functions are split along four levels :
//    - Level 0 : Element operations
//    - Level 1 : Vector-Vector operations
//    - Level 2 : Matrix-Vector operations
//    - Level 3 : Matrix-Matrix operations

// This first file has many comments. As the levels increase
// the amount of comments will deminish as it will be assumed
// That the reader knows about how the programs form lower
// levels work.

// ==========Includes=========

// fflas-ffpack modules
// #include <fflas-ffpack/fflas-ffpack-config.h>
#include <fflas-ffpack/fflas/fflas.h>

// Givaro provides some finite fields such as modular rings.
// Here we include two types of modular rings.
//    1. Normal   modular rings with elements numbered from 0 to p
//    2. Balanced modular rings with elements numbered form -(p-1)/2 to (p-1)/2
// Typically p is a prime number but Givaro works with any interger n as well.
// For other ways of defining fields using Givaro please refer to
// Givaro's documentation.
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>



#include <iostream>

// All of the methods used come from the FFLAS namespace.
using namespace FFLAS;

int main(int argc, char** argv) {

    // ========== Defining finite fields  =========

    // Defining rings using Givaro.
    typedef Givaro::Modular<float> Float_Ring;
    Float_Ring F1(101);
    // F1 is a ring of floats defined modulo 101.

    typedef Givaro::ModularBalanced<double> Double_Balanced_Ring;
    Double_Balanced_Ring F2(105);
    // F2 is a balanced ring of doubles defined modulo 105;


    // ===== Elements in finite fields ======


    // Let a, b and c be three elements of F1
    // and d, e and f be three elements of F2.

    Float_Ring::Element a,b,c,h;
    Double_Balanced_Ring::Element d,e,f,g;

    // The elements are defined but not yet assigned.
    // There are two ways of assigning values to elements
    // that depend on wether these values are integers or
    // already exixting elements.
    // The first way is using init().

    F1.init(a,2);
    F1.init(b,205);
    // a and b are initialised with the elements of F1 that
    // correspond to 2 and 205.
    // Since F1 is made of floats modulo 101, these are 2.0 and 3.0.

    // This means the second parameter will be converted into an
    // actual element of the field.

    // The second way is using assign().
    // A field has two elements that are neutral with respect to
    // the addition and multiplication laws, ie : zero and one.

    F2.assign(d, F2.zero);
    F2.assign(e, F2.one);
    F2.assign(f, e);

    // Operations on elements
    // Since the modular reduction is very expensive it's important
    // to call it as rarely as possible. Therefore the following
    // operations way yield a result outside of the fields.

    F1.add(c,a,b);        // c = a + b
    F2.mul(g,e,-5);       // g = e * (-5)
    F1.mul(h,c,-3);       // h = c * (-3)

    // The opposite operations sub and div are used just the same.

    // ===== Display =====
    // The fields form Givaro provide a display method called
    // write().
    // It's important to use write() because of the fact that the
    // opertations may get out of the field. Here's why :

    std::cout << "\n===== Linear Algebra level 0 =====" << std::endl;
    std::cout << "a := " << a << " = 2      % 101" << std::endl;
    std::cout << "b := " << b << " = 205    % 101" << std::endl;
    std::cout << "c := " << c << " = a + b  % 101" << std::endl;
    std::cout << "d := " << d << " = 1      % 105" << std::endl;
    std::cout << "e := " << e << " = 0      % 105" << std::endl;
    std::cout << "f := " << f << " = e      % 105" << std::endl;
    std::cout << "g := " << g <<"  = e*(-5) % 105" << std::endl;
    std::cout << "h := " << h << " = c*(-3) % 101" << std::endl;

    std::cout << "\nThis last result is strange,"<< std::endl;
    std::cout << "h is supposed to be defined %101."<< std::endl;
    std::cout << "With the Field display :"<< std::endl;
    std::cout << "h := " << std::flush;
    F1.write(std::cout,h) << std::endl;

    std::cout << "\nThis is strange as well."<< std::endl;
    std::cout << "That's because -3 is not"<< std::endl;
    std::cout << "an element of the field."<< std::endl;
    std::cout << "Thus -3 ought to be converted first."<< std::endl;
    F1.init(h,-3);
    std::cout << "F1.init(h,-3);\nF1.mulin(h,c);"<< std::endl;
    std::cout << "With the Field display :"<< std::endl;
    std::cout << "h := " << std::flush;
    F1.write(std::cout,h) << "\n" <<std::endl;
} //end main
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
