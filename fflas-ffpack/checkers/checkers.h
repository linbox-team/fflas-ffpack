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

#include "fflas-ffpack/fflas-ffpack.h"
#include <list>
#include <vector>
#include <iostream>
#include <algorithm>

#ifdef DEBUG
	#define ENABLE_CHECKING 1
#endif


// interface for all the checkers
template <class Field>
class Checker_Itf {
	Field F;
public:
	Checker_Itf(Field F): F(F) {}
	~Checker_Itf() = 0;
	virtual void init()  = 0;
	virtual bool check() = 0;
};


template <class Field>
class Checker_Empty : public Checker_Itf<Field> {
public:
	bool check() { return true; }
};

template <class Field> class Checker_PLUQ;

#endif