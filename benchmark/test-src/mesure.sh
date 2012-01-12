#!/bin/bash
# Written by Clément Pernet <clement.pernet@imag.fr>
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

#echo -n -e "   running $1 test  \t \t..."
printf "    running %-15s ........... " $1
while read parameter
do 
  ${BIN_PATH}/$1 $2 $parameter  2>>  ${TEST_PATH}/timing-$1-$2.txt
done < "${TEST_SRC_PATH}/parameter.in"
echo "[done]"
