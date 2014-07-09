#!/bin/bash
# Written by Clément Pernet <clement.pernet@imag.fr>
# Copyright (c)  FFLAS-FFPACK
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

# usage : process.sh "test name"  blas-result.txt  fflas-result.txt 

prime=`basename $2 | cut -d "-" -f 4 | cut -d "." -f 1`
blas_fct=`basename $2 | cut -d "-" -f 3`
fflas_fct=`basename $3 | cut -d "-" -f 3`


echo `pwd`
echo  "<test name=\"$1\">"
echo -n "<prime> "
echo -n "$prime "
echo "</prime>"

echo -n "<function name="
echo -n "\"$blas_fct\""
echo  " blas=\"GOTO\">"
awk -F " " '{print "<run size=\42"$2"\42> "$6"</run>" }' GOTO/$2
echo "</function>"

echo -n "<function name="
echo -n "\"$fflas_fct\""
echo  " blas=\"GOTO\">"
awk -F " " '{print "<run size=\42"$2"\42> "$6"</run>" }' GOTO/$3
echo "</function>"


echo -n "<function name="
echo -n "\"$blas_fct\""
echo  " blas=\"ATLAS\">"
awk -F " " '{print "<run size=\42"$2"\42> "$6"</run>" }' ATLAS/$2
echo "</function>"

echo -n "<function name="
echo -n "\"$fflas_fct\""
echo  " blas=\"ATLAS\">"
awk -F " " '{print "<run size=\42"$2"\42> "$6"</run>" }' ATLAS/$3
echo "</function>"

echo "</test>"

