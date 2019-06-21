#!/bin/bash
echo =================================================
echo ========= FFLAS-FFPACK fsyrk Autotuning =========
echo =================================================
echo 
(./fsyrk > fsyrk-threshold.h) 2>&1 | tee fsyrk-autotune.log
val=${PIPESTATUS[0]}; if test ${val} -ne 0 ; then exit ${val}; fi
