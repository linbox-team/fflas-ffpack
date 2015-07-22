# Copyright (c) 2010 the LinBox group
# Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
# ========LICENCE========
# This file is part of the library LinBox.
#
# LinBox is free software: you can redistribute it and/or modify
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


pkgincludesubdir=$(pkgincludedir)/interfaces/libs
AM_CPPFLAGS=-I$(top_srcdir)
AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS += $(OPTFLAGS)  -I$(top_srcdir)/fflas-ffpack/utils/ -I$(top_srcdir)/fflas-ffpack/fflas/  -I$(top_srcdir)/fflas-ffpack/ffpack  -I$(top_srcdir)/fflas-ffpack/field $(GIVARO_CFLAGS) $(CBLAS_FLAG) $(CUDA_CFLAGS) $(PARFLAGS)
LDADD = $(CBLAS_LIBS) $(GIVARO_LIBS) $(CUDA_LIBS) $(PARFLAGS)
AM_LDFLAGS=-static



pkgincludesub_HEADERS=fflas_c.h  \
		      ffpack_c.h

lib_LTLIBRARIES=libfflas_c.la \
		libffpack_c.la

libfflas_c_la_SOURCES=fflas_lvl1.C \
		    fflas_lvl2.C \
		    fflas_lvl3.C \
		    fflas_sparse.C

libffpack_c_la_SOURCES=ffpack.C

libffpack_c_la_LDFLAGS= -version-info 1:0:0
libffpack_c_la_LDFLAGS=  \
		       -no-undefined

libfflas_c_la_LDFLAGS= -version-info 1:0:0
libfflas_c_la_LDFLAGS+= \
                      -no-undefined

EXTRA_DIST=c_libs.doxy
