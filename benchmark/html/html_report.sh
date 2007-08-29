#!/bin/bash

CURRENT_PATH=`pwd`
HOME_PATH="${CURRENT_PATH}/.."


echo "Choose the architecture to process (`ls $HOME_PATH/Target -I CVS`):"
read ARCHI

echo "Enter the description of the architecture :"
read ARCHI_DESCR

TEST_DIR=`ls $HOME_PATH/Target/$ARCHI -I CVS -I GOTO -I ATLAS -I compilation.log | tail -n 1`
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
                                                                   timing-check-ftrsm-${PRIME}.txt  >> ${XML_FILE}

${CURRENT_PATH}/process.sh "Matrix Inversion"             timing-check-dgetri-${PRIME}.txt \
                                                          timing-check-inverse-${PRIME}.txt  >> ${XML_FILE}

${CURRENT_PATH}/process.sh "Triangular Matrix Inversion"  timing-check-dtrtri-${PRIME}.txt \
                                                          timing-check-trinverse-${PRIME}.txt  >> ${XML_FILE}
echo "</benchmark>" >> ${XML_FILE}


xsltproc -o ${HTML_FILE} $CURRENT_PATH/html_report.xsl ${XML_FILE}