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

#ifndef __FFLASFFPACK_checkers_fflas_H
#define __FFLASFFPACK_checkers_fflas_H

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "checker_empty.h"

#ifdef __FFLASFFPACK_DEBUG
#define CHECKING_MODE 1
#define ENABLE_ALL_CHECKINGS 1
#endif

#ifdef ENABLE_ALL_CHECKINGS
#define ENABLE_CHECKER_fgemm 1
#define ENABLE_CHECKER_ftrsm 1
#endif

#ifdef TIME_CHECKERS
#include <givaro/givtimer.h>
#define TIME_CHECKER_FGEMM
#define TIME_CHECKER_FTRSM
#endif

// definition of the exceptions
class FailureFgemmCheck {};
class FailureTrsmCheck {};

namespace FFLAS {
    template <class Field> class CheckerImplem_fgemm;
    template <class Field> class CheckerImplem_ftrsm;
}

namespace FFLAS {
#ifdef ENABLE_CHECKER_fgemm
    template <class Field> using Checker_fgemm = CheckerImplem_fgemm<Field>;
#else
    template <class Field> using Checker_fgemm = FFLAS::Checker_Empty<Field>;
#endif

#ifdef ENABLE_CHECKER_ftrsm
    template <class Field> using Checker_ftrsm = CheckerImplem_ftrsm<Field>;
#else
    template <class Field> using Checker_ftrsm = FFLAS::Checker_Empty<Field>;
#endif
}

#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas/fflas_enum.h"
#include "fflas-ffpack/utils/fflas_memory.h"

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
