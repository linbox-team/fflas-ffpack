#!/bin/bash

primes="65521"

for p in $primes; do 
    ${TEST_SRC_PATH}/mesure.sh check-dgemm  $p
    ${TEST_SRC_PATH}/mesure.sh check-dgetrf $p
    ${TEST_SRC_PATH}/mesure.sh check-dgetri $p
    ${TEST_SRC_PATH}/mesure.sh check-dtrsm  $p
    ${TEST_SRC_PATH}/mesure.sh check-dtrtri $p
done