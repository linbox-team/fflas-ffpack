/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* tests/print-utils.h
 * Copyright (C) 2011, Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 * 					   Bastien Vialla <bastien.vialla@lirmm.fr>
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
 */

#ifndef __FFLASFFPACK_print_utils_H
#define __FFLASFFPACK_print_utils_H

#include <fflas-ffpack/fflas-ffpack-config.h>
#include <vector>
// #include <pair>
#include <list>
#include <set>
#include <iterator>

namespace std
{

	/*! Prints a vector on output.
	 * @param o output stream
	 * @param v vector
	 * @warning <<(ostream&,T&) exists !
	 */
	template<class T, class Alloc>
	std::ostream & operator<<(std::ostream&o, const std::vector<T, Alloc> & v)
	{
		o << '[' ;
		std::copy(v.begin(), v.end(), std::ostream_iterator<T>(o, " "));
		// if (v.size()) {
		// 	size_t i = 0  ;
		// 	for (; i < v.size()-1 ; ++i)
		// 		o << v[i] << ',' ;
		// 	o <<  v[i] ;
		// }
		return o << ']' ;
	}


	/*! Prints a pair.
	 * @param o output stream
	 * @param C a pair
	 * @warning <<(ostream&,T&) exists !
	 */
	template<class S, class T>
	std::ostream& operator<<(std::ostream& o, const std::pair<S, T> & C)
	{
		o << '(' << C.first << ", " << C.second << ')';
		return o ;
	}


	/*! Prints a list.
	 * @param o output stream
	 * @param C a pair
	 * @warning <<(ostream&,T&) exists !
	 */
	template<class T, class Alloc>
	std::ostream& operator<< (std::ostream& o, const std::list<T, Alloc> & L)
	{
		o << '{' ;
		std::copy(L.begin(), L.end(), std::ostream_iterator<T>(o, " "));
		return o << '}' ;

		// typename std::list<T>::const_iterator it = L.begin() ;
		// if (it != L.end() )
		// 	while(true) {
		// 		o << *it ;
		// 		if (++it != L.end())
		// 			o << ", " ;
		// 		else
		// 			break;
		// 	}
	}


	/*! Prints a set.
	 * @param o output stream
	 * @param C a pair
	 * @warning <<(ostream&,T&) exists !
	 */
	template<class T, class Alloc>
	std::ostream& operator<< (std::ostream& o, const std::set<T, Alloc> & S)
	{
		o << '|' ;
		std::copy(S.begin(), S.end(), std::ostream_iterator<T>(o, " "));
		return o << '|' ;
		// typename std::set<T>::const_iterator it = L.begin() ;
		// o << '|' ;
		// if (it != L.end() )
		// 	while(true) {
		// 		o << *it ;
		// 		if (++it != L.end())
		// 			o << ", " ;
		// 		else
		// 			break;
		// 	}
		// return o << '|' ;
	}


#if 0
	std::ostream &operator << (std::ostream &out, const std::vector<bool> &S)
	{
		std::vector<bool>::const_iterator i;

		for (i = S.begin (); i != S.end (); ++i) {
			out << ((*i) ? "1" : "0");
			if (i != S.end () - 1)
				out << ", ";
		}

		return out;
	}

	template<class T, template <class T> class Container>
	std::ostream& operator<< (std::ostream& o, const Container<T>& C)
	{
		for(typename Container<T>::const_iterator refs =  C.begin();
		    refs != C.end() ;
		    ++refs )
			o << (*refs) << " " ;
		return o << std::endl;
	}


#endif

}

#endif // __FFLASFFPACK_print_utils_H
