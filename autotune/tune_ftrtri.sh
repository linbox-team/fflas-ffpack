#!/bin/bash
echo =================================================
echo ========= FFLAS-FFPACK ftrtri Autotuning ========
echo =================================================
echo 
(./ftrtri > ftrtri-threshold.h) 2>&1 | tee ftrtri-autotune.log
val=${PIPESTATUS[0]}; if test ${val} -ne 0 ; then exit ${val}; fi
