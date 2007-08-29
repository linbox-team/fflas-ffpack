#!/bin/bash

#echo -n -e "   running $1 test  \t \t..."
printf "    running %-15s ........... " $1
while read parameter
do 
  ${BIN_PATH}/$1 $2 $parameter  2>>  ${TEST_PATH}/timing-$1-$2.txt
done < "${TEST_SRC_PATH}/parameter.in"
echo "[done]"
