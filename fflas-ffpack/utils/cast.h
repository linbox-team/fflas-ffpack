/*
 * This file is part of FFLAS-FFPACK
 * Copyright (C) 2011 Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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
 * License along with this library; if not, see
 * <https://www.gnu.org/licenses/>.
 * ========LICENCE========
 *.
 */

#ifndef __FFLASFFPACK_const_H
#define __FFLASFFPACK_const_H

namespace FFPACK {

    template<class T, class CT = const T>
    T fflas_const_cast (CT x)
    {
        return const_cast<T>(x);
    }

}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
