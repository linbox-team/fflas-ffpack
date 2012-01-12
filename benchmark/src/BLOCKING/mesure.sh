#!/bin/csh -f
# Copyright (c) FFLAS-FFPACK
# Written by T. Gautier, jgdumas
# adapted from Givaro.
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


set block = $1

#foreach i (10 12 14 16 18 20 22 24 26 28 30 50 75 100 125 150)
foreach i (175 200 225)
	@ N = $i * $block
	tblockmat-$block $N 1 
end
	
