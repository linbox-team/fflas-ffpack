/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* tests/print-utils.h
 * Copyright (C) 2011, Brice Boyer <bboyer@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __FFLAFLAS_print_utils_H
#define __FFLAFLAS_print_utils_H

#include <vector>
// #include <pair>
#include <list>
#include <set>

namespace std
{

	/*! Prints a vector on output.
	 * @param o output stream
	 * @param v vector
	 * @warning <<(ostream&,T&) exists !
	 */
	template<class T>
	std::ostream & operator<<(std::ostream&o, const std::vector<T> & v)
	{
		o << '[' ;
		if (v.size()) {
			size_t i = 0  ;
			for (; i < v.size()-1 ; ++i)
				o << v[i] << ',' ;
			o <<  v[i] ;
		}
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
	template<class T>
	std::ostream& operator<< (std::ostream& o, const std::list<T> & L)
	{
		typename std::list<T>::const_iterator it = L.begin() ;
		o << '{' ;
		if (it != L.end() )
			while(true) {
				o << *it ;
				if (++it != L.end())
					o << ", " ;
				else
					break;
			}
		return o << '}' ;
	}


	/*! Prints a set.
	 * @param o output stream
	 * @param C a pair
	 * @warning <<(ostream&,T&) exists !
	 */
	template<class T>
	std::ostream& operator<< (std::ostream& o, const std::set<T> & L)
	{
		typename std::set<T>::const_iterator it = L.begin() ;
		o << '|' ;
		if (it != L.end() )
			while(true) {
				o << *it ;
				if (++it != L.end())
					o << ", " ;
				else
					break;
			}
		return o << '|' ;
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

#endif // __FFLAFLAS_print_utils_H
