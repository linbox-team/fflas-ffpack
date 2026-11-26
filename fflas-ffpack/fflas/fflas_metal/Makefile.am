# Copyright (c) 2026 FFLAS-FFPACK group
#
#
# ========LICENCE========
# This file is part of the library FFLAS-FFPACK.
#
# FFLAS-FFPACK is free software: you can redistribute it and/or modify
# it under the terms of the  GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========
#/


pkgincludesubdir=$(pkgincludedir)/fflas/fflas_metal

pkgincludesub_HEADERS=          \
	fflas_metal.h

EXTRA_DIST=fflas_metal.m

if FFLASFFPACK_HAVE_METAL
# The one compiled file of the Metal backend (Objective-C, built by the C++ compiler's driver),
# as a static library beside the header-only library.
lib_LIBRARIES=libfflas_metal.a
libfflas_metal_a_SOURCES=
libfflas_metal_a_LIBADD=fflas_metal.$(OBJEXT)

fflas_metal.$(OBJEXT): $(srcdir)/fflas_metal.m $(srcdir)/fflas_metal.h
	$(CXX) -x objective-c -fobjc-arc -O2 -Wall -c $(srcdir)/fflas_metal.m -o $@

CLEANFILES=fflas_metal.$(OBJEXT)
endif
