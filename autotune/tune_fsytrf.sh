#!/bin/sh
echo =================================================
echo ========= FFLAS-FFPACK fsytrf Autotuning ========
echo =================================================
echo 
(./fsytrf > fsytrf-threshold.h) 2>&1 | tee fsytrf-autotune.log
