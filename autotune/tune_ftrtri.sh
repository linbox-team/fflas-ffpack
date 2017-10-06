#!/bin/sh
echo =================================================
echo ========= FFLAS-FFPACK ftrtri Autotuning ========
echo =================================================
echo 
(./ftrtri > ftrtri-threshold.h) 2>&1 | tee ftrtri-autotune.log
