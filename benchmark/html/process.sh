#!/bin/bash

# usage : process.sh "test name"  blas-result.txt  fflas-result.txt 

prime=`basename $2 | cut -d "-" -f 4 | cut -d "." -f 1`
blas_fct=`basename $2 | cut -d "-" -f 3`
fflas_fct=`basename $3 | cut -d "-" -f 3`


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

