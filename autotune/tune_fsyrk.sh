#!/bin/sh
echo =================================================
echo ========= FFLAS-FFPACK fsyrk Autotuning =========
echo =================================================
echo 
./fsyrk 2> fsyrk-threshold.h  | tee fsyrk-autotune.log
