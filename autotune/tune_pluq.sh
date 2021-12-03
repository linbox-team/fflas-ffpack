#!/bin/bash
echo =================================================
echo ========= FFLAS-FFPACK PLUQ Autotuning ==========
echo =================================================
echo 
(./pluq > pluq-threshold.h) 2>&1 | tee pluq-autotune.log
val=${PIPESTATUS[0]}; if test ${val} -ne 0 ; then exit ${val}; fi
