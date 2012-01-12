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

primes="65521"

for p in $primes; do 
    ${TEST_SRC_PATH}/mesure.sh check-dgemm  $p
    ${TEST_SRC_PATH}/mesure.sh check-dgetrf $p
    ${TEST_SRC_PATH}/mesure.sh check-dgetri $p
    ${TEST_SRC_PATH}/mesure.sh check-dtrsm  $p
    ${TEST_SRC_PATH}/mesure.sh check-dtrtri $p
done
