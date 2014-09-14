/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_memory.h
 * Copyright (C) 2014 fflas-ffpack group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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

#ifndef __FFLASFFPACK_memory_H
#define __FFLASFFPACK_memory_H

#include "fflas-ffpack/utils/align-allocator.h"

namespace FFLAS{

    // template<class Field>
    // inline typename Field::Element_ptr fflas_new (const Field& F, size_t m, size_t n)
    // {
	   //  //return new typename Field::Element[m*n];
	   //  return malloc_align<typename Field::Element>(m*n, Alignment::AVX);
    // }
 
    template<class Field>
    inline typename Field::Element_ptr fflas_new (const Field& F, const size_t m, const size_t n, const Alignment align = Alignment::AVX)
    {
        //return new typename Field::Element[m*n];
        return malloc_align<typename Field::Element>(m*n, align);
    }

    template<class Element >
    inline Element* fflas_new (const size_t m, const Alignment align = Alignment::AVX)
    {
	   // return new typename Field::Element[m*n];
	   return malloc_align<Element>(m, align);
    }

    template<class Element_ptr>
    inline void fflas_delete (Element_ptr A)
    {
	    //delete[] A;
	    free(A);
    }

}
#endif // __FFLASFFPACK_memory_H
