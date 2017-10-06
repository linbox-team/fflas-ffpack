#!/bin/sh
echo =================================================
echo ========= FFLAS-FFPACK PLUQ Autotuning ==========
echo =================================================
echo 
(./pluq > pluq-threshold.h) 2>&1 | tee pluq-autotune.log
