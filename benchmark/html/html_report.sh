#!/bin/bash
# Copyright (c) FFLAS-FFPACK
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


CURRENT_PATH=`pwd`
HOME_PATH="${CURRENT_PATH}/.."


echo "Choose the architecture to process (`ls $HOME_PATH/Target -I CVS`):"
read ARCHI

echo "Enter the description of the architecture :"
read ARCHI_DESCR

TEST_DIR=`ls $HOME_PATH/Target/$ARCHI  -I GOTO -I ATLAS -I compilation.log -I report.xml| tail -n 1`
echo "Processing testing directory [$TEST_DIR]"

TEST_PATH="$HOME_PATH/Target/$ARCHI/$TEST_DIR"

PRIME=65521


XML_FILE="$TEST_PATH/report.xml"
HTML_FILE="$CURRENT_PATH/report-${ARCHI}.html"
cd $TEST_PATH

echo "<?xml version=\"1.0\"  encoding=\"ISO-8859-1\"?>" > ${XML_FILE}
echo "<benchmark>"  >> ${XML_FILE}
echo "<archi> $ARCHI_DESCR </archi>"  >> ${XML_FILE}
echo "<prime> $PRIME </prime>" >> ${XML_FILE}

${CURRENT_PATH}/process.sh "Matrix Multiplication"        timing-check-dgemm-${PRIME}.txt \
                                                          timing-check-fgemm-${PRIME}.txt  >> ${XML_FILE}

${CURRENT_PATH}/process.sh "Matrix Triangularization"     timing-check-dgetrf-${PRIME}.txt \
                                                          timing-check-lqup-${PRIME}.txt  >> ${XML_FILE}

${CURRENT_PATH}/process.sh "Multiple Triangular System Solving"    timing-check-dtrsm-${PRIME}.txt \
                                                                   timing-check-ftrsm-${PRIME}.txt >> ${XML_FILE}

${CURRENT_PATH}/process.sh "Matrix Inversion"             timing-check-dgetri-${PRIME}.txt \
                                                          timing-check-inverse-${PRIME}.txt  >> ${XML_FILE}

${CURRENT_PATH}/process.sh "Triangular Matrix Inversion"  timing-check-dtrtri-${PRIME}.txt \
                                                          timing-check-ftrtri-${PRIME}.txt >> ${XML_FILE}
echo "</benchmark>" >> ${XML_FILE}


xsltproc -o ${HTML_FILE} $CURRENT_PATH/html_report.xsl ${XML_FILE}
