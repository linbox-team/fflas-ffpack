#!/bin/bash
echo =================================================
echo ========= FFLAS-FFPACK fsytrf Autotuning ========
echo =================================================
echo 
(./fsytrf > fsytrf-threshold.h) 2>&1 | tee fsytrf-autotune.log
val=${PIPESTATUS[0]}; if test ${val} -ne 0 ; then exit ${val}; fi
