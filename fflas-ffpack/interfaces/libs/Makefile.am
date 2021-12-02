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

if FFLASFFPACK_PRECOMPILED

pkgincludesubdir=$(pkgincludedir)/interfaces/libs

AM_CXXFLAGS = $(FFLASFFPACK_CXXFLAGS) $(GIVARO_CFLAGS) $(BLAS_CFLAGS) $(PARFLAGS)
AM_CPPFLAGS = -I$(top_srcdir)
#AM_LDFLAGS=-static 



pkgincludesub_HEADERS=fflas_c.h  \
		      ffpack_c.h \
		      fflas_L3_inst.h \
		      fflas_L3_inst_implem.inl \
		      fflas_L2_inst.h \
		      fflas_L2_inst_implem.inl \
		      fflas_L1_inst.h \
		      fflas_L1_inst_implem.inl \
		      ffpack_inst.h \
		      ffpack_inst_implem.inl


lib_LTLIBRARIES=libfflas.la \
		libffpack.la \
		libfflas_c.la \
		libffpack_c.la



libfflas_la_SOURCES= fflas_L1_inst.C \
		    fflas_L1_inst_implem.inl\
		    fflas_L2_inst.C \
		    fflas_L2_inst_implem.inl \
		    fflas_L3_inst.C \
		    fflas_L3_inst_implem.inl

libfflas_la_LIBADD= $(GIVARO_LIBS) $(BLAS_LIBS) $(PARLIBS)
libfflas_la_LDFLAGS= -version-info 1:0:0 -no-undefined

libffpack_la_SOURCES= ffpack_inst.C \
		      ffpack_inst_implem.inl
libffpack_la_LIBADD= libfflas.la
libffpack_la_LDFLAGS= -version-info 1:0:0 -no-undefined

libfflas_c_la_SOURCES=fflas_lvl1.C \
		    fflas_lvl2.C \
		    fflas_lvl3.C \
		    fflas_sparse.C
#libfflas_c_la_CPPFLAGS=$(AM_CPPFLAGS) -DFFLAS_COMPILED -DFFPACK_COMPILED
libfflas_c_la_LIBADD= libfflas.la
libfflas_c_la_LDFLAGS= -version-info 1:0:0 -no-undefined

libffpack_c_la_SOURCES=ffpack.C
#libffpack_c_la_CPPFLAGS=$(AM_CPPFLAGS) -DFFLAS_COMPILED -DFFPACK_COMPILED
libffpack_c_la_LIBADD= libffpack.la
libffpack_c_la_LDFLAGS= -version-info 1:0:0 -no-undefined

EXTRA_DIST=c_libs.doxy

endif
