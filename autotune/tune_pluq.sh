#!/bin/sh
echo =================================================
echo ========= FFLAS-FFPACK PLUQ Autotuning ==========
echo =================================================
echo 
./pluq 2> pluq-threshold.h  | tee pluq-autotune.log
