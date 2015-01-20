/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* field/field-general
 * This file is part of FFLAS-FFPACK
 * Copyright (C) 2011 Brice Boyer <bboyer@imag.fr>
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
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */


/*! @file field/modular-general.h
 * @ingroup field
 * @brief fields in FFLAS-FFPACK
 */

#ifndef __FFLASFFPACK_field_general_H
#define __FFLASFFPACK_field_general_H

#include <ostream>

namespace Givaro {

	template<class T>
	class Modular;

	template<class T>
	class ModularBalanced ;

}

namespace FFPACK {


	template<class T>
	class UnparametricField ;

	template<class T>
	class RNSInteger;

	template<class T>
	class RNSIntegerMod;


	template<class T>
	std::ostream & operator<<( std::ostream & o, const Givaro::Modular<T> & F)
	{
		return F.write(o);
	}

	template<class T>
	std::ostream & operator<<( std::ostream & o, const Givaro::ModularBalanced<T> & F)
	{
		return F.write(o);
	}

	template<class T>
	std::ostream & operator<<( std::ostream & o, const UnparametricField<T> & F)
	{
		return F.write(o);
	}

	template<class T, class CT=const T>
	T fflas_const_cast (CT x){return const_cast<T>(x);}
	
	// template<>
	// RNSInteger::Element_ptr 
	// fflas_const_cast<RNSInteger::Element_ptr, RNSInteger::ConstElement_ptr>(RNSInteger::ConstElement_ptr x){
	// 	return (RNSInteger::Element_ptr)x;
	// }
}

#endif // __FFLASFFPACK_field_general_H
