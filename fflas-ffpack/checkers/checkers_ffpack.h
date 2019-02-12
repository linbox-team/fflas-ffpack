/* checkers/checkers.h
 * Copyright (C) 2016 Ashley Lesdalons, JG Dumas
 *
 * Written by Ashley Lesdalons <Ashley.Lesdalons@e.ujf-grenoble.fr>
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
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

#ifndef __FFLASFFPACK_checkers_ffpack_H
#define __FFLASFFPACK_checkers_ffpack_H

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "checker_empty.h"

#ifdef __FFLASFFPACK_DEBUG
#define CHECKING_MODE 1
#define ENABLE_ALL_CHECKINGS 1
#endif

#ifdef ENABLE_ALL_CHECKINGS
#define ENABLE_CHECKER_PLUQ 1
#define ENABLE_CHECKER_Det 1
#define ENABLE_CHECKER_invert 1
#define ENABLE_CHECKER_charpoly 1
#endif

#ifdef TIME_CHECKERS
#include <givaro/givtimer.h>
#define TIME_CHECKER_PLUQ
#define TIME_CHECKER_Det
#define TIME_CHECKER_INVERT
#define TIME_CHECKER_CHARPOLY
#endif


// definition of the exceptions
class FailurePLUQCheck {};
class FailureDetCheck {};
class FailureInvertCheck {};
class FailureCharpolyCheck {};

namespace FFPACK {
    template <class Field> class CheckerImplem_PLUQ;
    template <class Field> class CheckerImplem_Det;
    template <class Field> class CheckerImplem_invert;
    template <class Field, class Polynomial> class CheckerImplem_charpoly;
}


namespace FFPACK {
#ifdef ENABLE_CHECKER_PLUQ
    template <class Field> using Checker_PLUQ = CheckerImplem_PLUQ<Field>;
#else
    template <class Field> using Checker_PLUQ = FFLAS::Checker_Empty<Field>;
#endif

#ifdef ENABLE_CHECKER_Det
    template <class Field> using Checker_Det = CheckerImplem_Det<Field>;
#else
    template <class Field> using Checker_Det = FFLAS::Checker_Empty<Field>;
#endif

#ifdef ENABLE_CHECKER_invert
    template <class Field> using Checker_invert = CheckerImplem_invert<Field>;
#else
    template <class Field> using Checker_invert = FFLAS::Checker_Empty<Field>;
#endif

#ifdef ENABLE_CHECKER_charpoly
    template <class Field, class Polynomial> using Checker_charpoly = CheckerImplem_charpoly<Field,Polynomial>;
#else
    template <class Field, class Polynomial> using Checker_charpoly = FFLAS::Checker_Empty<Field>;
#endif
}

#include "fflas-ffpack/ffpack/ffpack.h"

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
