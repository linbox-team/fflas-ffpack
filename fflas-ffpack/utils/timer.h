/* test/timer.h
 * Copyright (C) 1994-1997 Givaro Team
 *
 * Written by T. Gautier
 *
 * ------------------------------------
 * Modified by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Added _start_t member to BaseTimer, so that stop () does not clobber the
 * class' memory of its start time. This allows it to be called repeatedly to
 * get elapsed times.
 * ------------------------------------
 * Modified by Clement Pernet
 * integrated into FFLAS_FFPACK
 *
 * ------------------------------------
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
 *
 * This file implements the C++ interface to commentators (for
 * providing runtime commentary to the user)
 */

#ifndef __FFLASFFPACK_timer_H
#define __FFLASFFPACK_timer_H


#include <time.h>

#ifdef __FFLASFFPACK_USE_OPENMP
#  ifndef __GIVARO_USE_OPENMP
#    define __GIVARO_USE_OPENMP 1
#  endif
#endif

#include <givaro/givtimer.h>

#ifdef __GIVARO_USE_OPENMP
#include <givaro/givomptimer.h>
#endif

namespace FFLAS {
    typedef Givaro::Timer Timer  ;
    typedef Givaro::BaseTimer BaseTimer ;
    typedef Givaro::UserTimer UserTimer ;
    typedef Givaro::SysTimer SysTimer ;
#ifdef __GIVARO_USE_OPENMP
    typedef Givaro::OMPTimer OMPTimer ;
#endif
}

#endif // __FFLASFFPACK_timer_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
