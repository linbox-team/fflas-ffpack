#!/bin/bash

primes="65521"

for p in $primes; do 
    ${TEST_SRC_PATH}/mesure.sh check-fgemm     $p 
    ${TEST_SRC_PATH}/mesure.sh check-lqup      $p 
    ${TEST_SRC_PATH}/mesure.sh check-inverse   $p 
    ${TEST_SRC_PATH}/mesure.sh check-ftrsm     $p 
    ${TEST_SRC_PATH}/mesure.sh check-trinverse $p 
done