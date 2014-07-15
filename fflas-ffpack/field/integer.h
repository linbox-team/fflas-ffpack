/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* field/unparametric.h
 * Copyright (C) FFLAS group
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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
 *
 */


/*! @file field/integer.h
 * @ingroup field
 * @brief  multiprecision integer from Givaro 
 */

#include "fflas-ffpack/config.h"
#ifdef __FFLASFFPACK_HAVE_GIVARO
#ifndef __FFLASFFPACK_NOINTEGER
 
#include "givaro/givinteger.h"
namespace FFPACK {
  typedef Givaro::Integer Integer;
  typedef Givaro::Integer integer;
  #define __FFLASFFPACK_HAVE_INTEGER
};
#endif
#endif 
 
