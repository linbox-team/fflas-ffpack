#!/bin/sh
echo =================================================
echo ========= FFLAS-FFPACK fsyrk Autotuning =========
echo =================================================
echo 
(./fsyrk > fsyrk-threshold.h) 2>&1 | tee fsyrk-autotune.log
