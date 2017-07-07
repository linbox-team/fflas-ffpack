#!/bin/sh
echo =================================================
echo ========= FFLAS-FFPACK fsytrf Autotuning ========
echo =================================================
echo 
./fsytrf 2> fsytrf-threshold.h  | tee fsytrf-autotune.log
