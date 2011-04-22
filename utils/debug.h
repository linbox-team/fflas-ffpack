/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* utils/debug.h
 *
 * Copyright (C) 2011 Fflas-ffpack
 * Modified by BB, from LinBox
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 */

/*! @file utils/debug.h
 * @ingroup util
 * Various utilities for debugging.
 * @todo we should put vector printing elsewhere.
 */

#ifndef __FFLAFLAS_util_debug_H
#define __FFLAFLAS_util_debug_H

#include <sstream>

#include "fflasffpack-config.h"


#ifdef __FFLAFLAS_HAVE_STDINT_H
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#else
#error "you need intXX_t types"
#endif

#ifndef DEBUG
#define fflaflas_check(check)
#else
#define fflaflas_check(check) \
if (!(check)) \
throw FFPACK::Failure (__func__, __FILE__, __LINE__, #check); //BB : should work on non gnu compilers too
#endif



namespace FFPACK {
#if 0

/*!  A precondtion failed.
	 * @ingroup util
	 * The \c throw mechanism is usually used here as in
	 \code
	 if (!check)
	 throw(Failure(__func__,__LINE__,"this check just failed");
	 \endcode
	 * The parameters of the constructor help debugging.
	 */
	class Failure {//: public LinboxError BB: otherwise,  error.h:39 segfaults
		static std::ostream *_errorStream;

	public:
		/*! @internal
		 * A precondtion failed.
		 * @param function usually \c __func__, the function that threw the error
		 * @param line     usually \c __LINE__, the line where it happened
		 * @param check    a string telling what failed.
		 */
		Failure (const char *function, int line, const char *check)
		{
			if (_errorStream == (std::ostream *) 0)
				_errorStream = &std::cerr;

			(*_errorStream) << std::endl << std::endl;
			(*_errorStream) << "ERROR (" << function << ":" << line << "): ";
			(*_errorStream) << "Precondition not met:" << check << std::endl;
		}

		/*! @internal
		 * A precondtion failed.
		 * The parameter help debugging. This is not much different from the previous
		 * except we can digg faster in the file where the exception was triggered.
		 * @param function usually \c __func__, the function that threw the error
		 * @param file     usually \c __FILE__, the file where this function is
		 * @param line     usually \c __LINE__, the line where it happened
		 * @param check    a string telling what failed.
		 */
		Failure (const char* function, const char *file, int line, const char *check)
		{
			if (_errorStream == (std::ostream *) 0)
				_errorStream = &std::cerr;

			(*_errorStream) << std::endl << std::endl;
			(*_errorStream) << "ERROR (at " << function << " in " << file << ':' <<  line << "): " << std::endl;
			(*_errorStream) << "Precondition not met:" << check << std::endl;
		}

		static void setErrorStream (std::ostream &stream);

		/*! @internal overload the virtual print of LinboxError.
		 * @param o output stream
		 */
		std::ostream &print (std::ostream &o) const
		{
			if (std::ostringstream * str = dynamic_cast<std::ostringstream*>(_errorStream))
				return o << str->str() ;
			else
				throw "Fflas-Ffpack ERROR: Failure exception is not initialized correctly";
		}
	};

	/*! @internal A function is "not implemented yet(tm)".
	 * where, why ?
	 */
	class NotImplementedYet {
	protected:
		static std::ostream *_errorStream;

	public:
		/*! @internal
		 * A precondtion failed.
		 * The parameter help debugging. This is not much different from the previous
		 * except we can digg faster in the file where the exception was triggered.
		 * @param function usually \c __func__, the function that threw the error
		 * @param file     usually \c __FILE__, the file where this function is
		 * @param line     usually \c __LINE__, the line where it happened
		 * @param why      by default, lazy people don't provide an explanation.
		 */
		NotImplementedYet() {}

		NotImplementedYet(const char * function,
				  const char* file,
				  int line,
				  const char * why='\0')
		{
			if (_errorStream == (std::ostream *) 0)
				_errorStream = &std::cerr;

			(*_errorStream) << std::endl << std::endl;
			(*_errorStream) << "ERROR (at " << function << " in " << file << ':' <<  line << "): " << std::endl;
			(*_errorStream) << " This function is not implemented yet" ;
			if (why)
				(*_errorStream)	<< " (" << why << ")" <<std::endl;
			else
				(*_errorStream)	<<  "." << std::endl;

		}
	};
#endif

} // FFPACK

#endif // __FFLAFLAS_util_debug_H
