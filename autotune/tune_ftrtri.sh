#!/bin/sh
echo =================================================
echo ========= FFLAS-FFPACK ftrtri Autotuning ========
echo =================================================
echo 
./ftrtri 2> ftrtri-threshold.h  | tee ftrtri-autotune.log
