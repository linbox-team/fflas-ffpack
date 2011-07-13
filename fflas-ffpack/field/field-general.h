/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* field/field-general
 * This file is part of FFLAS-FFPACK
 * Copyright (C) 2011 Brice Boyer <bboyer@imag.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __FFLASFFPACK_field_general_H
#define __FFLASFFPACK_field_general_H

#include <ostream>

namespace FFPACK {


	template<class T>
	class UnparametricField ;

	template<class T>
	class Modular ;

	template<class T>
	class ModularBalanced ;

	template<class T>
	std::ostream & operator<<( std::ostream & o, const Modular<T> & F)
	{
		return F.write(o);
	}

	template<class T>
	std::ostream & operator<<( std::ostream & o, const ModularBalanced<T> & F)
	{
		return F.write(o);
	}

	template<class T>
	std::ostream & operator<<( std::ostream & o, const UnparametricField<T> & F)
	{
		return F.write(o);
	}

}


#endif // __FFLASFFPACK_field_general_H
