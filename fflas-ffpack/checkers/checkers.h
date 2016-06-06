/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* checkers/checkers.h
 * Copyright (C) 2016 Ashley Lesdalons
 *
 * Written by Ashley Lesdalons <Ashley.Lesdalons@e.ujf-grenoble.fr>
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

#ifndef __FFLASFFPACK_checkers_H
#define __FFLASFFPACK_checkers_H

#include <list>
#include <vector>
#include <iostream>
#include <algorithm>


#ifdef DEBUG
 	#define CHECKING_MODE 1
 	#define ENABLE_ALL_CHECKINGS 1
#endif

#ifdef ENABLE_ALL_CHECKINGS
	#define ENABLE_CHECKER_PLUQ 1
	#define ENABLE_CHECKER_fgemm 1
 	#define ENABLE_CHECKER_ftrsm 1
 	#define ENABLE_CHECKER_invert 1
#endif



template <class Field>
struct Checker_Empty {
	template<typename... Params> Checker_Empty(Params... parameters) {}
	template<typename... Params> bool check(Params... parameters) { return true; }
};


#ifdef ENABLE_CHECKER_PLUQ
	template <class Field> class Checker_PLUQ;
#else
	template <class Field> using Checker_PLUQ = Checker_Empty<Field>;
#endif

#ifdef ENABLE_CHECKER_fgemm
	template <class Field> class Checker_fgemm;
#else
	template <class Field> using Checker_fgemm = Checker_Empty<Field>;
#endif

#ifdef ENABLE_CHECKER_ftrsm
	template <class Field> class Checker_ftrsm;
#else
	template <class Field> using Checker_ftrsm = Checker_Empty<Field>;
#endif

#ifdef ENABLE_CHECKER_invert
	template <class Field> class Checker_invert;
#else
	template <class Field> using Checker_invert = Checker_Empty<Field>;
#endif

#endif